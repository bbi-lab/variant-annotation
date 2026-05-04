"""Add VCF-format columns (chromosome, position, ref, alt) for mapped HGVS fields.

Given mapped HGVS columns (genomic/transcript/protein), this script appends:
- <col>_chromosome: Chromosome number (1-22), X, Y, or M (mitochondria)
- <col>_start
- <col>_stop
- <col>_ref
- <col>_alt

For protein HGVS, ref/alt use one-character amino acid codes (e.g., A, R, C, etc.;
* for stop codon; empty string for deletion). For synonymous variants (ref == alt), the
amino acid is repeated.

When an HGVS column is pipe-delimited (multiple DNA candidates produced by
``reverse_translate_protein_variants`` + ``add_dna_clingen_allele_ids``), each
component is parsed independently and the output columns are also pipe-delimited,
preserving candidate cardinality so downstream steps can correlate columns.

It also adds transcript-derived boolean columns:
- touches_intronic_region: transcript HGVS contains intronic offset coordinates
- spans_intron: transcript HGVS spans both sides of an intron
"""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import csv
from itertools import islice
import logging
import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_HGVS_RESOLVER_UNAVAILABLE_LOGGED = False
_HGVS_REF_RESOLVER: Optional["_HgvsRefResolver"] = None

# Three-letter to one-letter amino acid code mapping
_AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O",  # Selenocysteine and Pyrrolysine
    "Ter": "*",  # Termination / stop codon
}

# Map accession prefixes to chromosomes
_ACCESSION_TO_CHROMOSOME = {
    "NC_000001": "1", "NC_000002": "2", "NC_000003": "3", "NC_000004": "4",
    "NC_000005": "5", "NC_000006": "6", "NC_000007": "7", "NC_000008": "8",
    "NC_000009": "9", "NC_000010": "10", "NC_000011": "11", "NC_000012": "12",
    "NC_000013": "13", "NC_000014": "14", "NC_000015": "15", "NC_000016": "16",
    "NC_000017": "17", "NC_000018": "18", "NC_000019": "19", "NC_000020": "20",
    "NC_000021": "21", "NC_000022": "22", "NC_000023": "X", "NC_000024": "Y",
    "NC_012920": "M",  # Mitochondrial genome
}


class _HgvsRefResolver:
    def __init__(self) -> None:
        import hgvs.dataproviders.uta  # type: ignore[import-untyped]
        import hgvs.normalizer  # type: ignore[import-untyped]
        import hgvs.parser  # type: ignore[import-untyped]

        uta_db_url = (os.environ.get("UTA_DB_URL") or "").strip()
        if not uta_db_url:
            raise RuntimeError("UTA_DB_URL is required to resolve missing HGVS ref alleles.")

        hdp = hgvs.dataproviders.uta.connect(uta_db_url)
        self._parser = hgvs.parser.Parser()
        self._normalizer = hgvs.normalizer.Normalizer(hdp)

    def resolve_ref(self, hgvs_value: str) -> Optional[str]:
        var = self._parser.parse_hgvs_variant(hgvs_value)
        normalized = self._normalizer.normalize(var)
        ref = getattr(normalized.posedit.edit, "ref", None)
        if isinstance(ref, str) and ref.strip() and ref.strip() != "?":
            return ref.strip()
        return None


def _get_hgvs_ref_resolver() -> Optional[_HgvsRefResolver]:
    global _HGVS_RESOLVER_UNAVAILABLE_LOGGED
    global _HGVS_REF_RESOLVER

    if _HGVS_REF_RESOLVER is not None:
        return _HGVS_REF_RESOLVER

    try:
        _HGVS_REF_RESOLVER = _HgvsRefResolver()
    except Exception as exc:
        if not _HGVS_RESOLVER_UNAVAILABLE_LOGGED:
            logger.warning(
                "Reference-aware ref allele resolution disabled: %s",
                exc,
            )
            _HGVS_RESOLVER_UNAVAILABLE_LOGGED = True
        return None

    return _HGVS_REF_RESOLVER


@lru_cache(maxsize=20000)
def _resolve_missing_ref_allele(hgvs_value: str) -> Optional[str]:
    resolver = _get_hgvs_ref_resolver()
    if resolver is None:
        return None

    try:
        return resolver.resolve_ref(hgvs_value)
    except Exception:
        return None


def _detect_separator(file_path: str) -> str:
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _is_blank(value: Optional[str]) -> bool:
    return value is None or not value.strip()


def _split_coord(coord: str) -> tuple[str, str]:
    coord = coord.strip()
    if "_" in coord:
        start, stop = coord.split("_", 1)
        return start, stop
    return coord, coord


def _is_intronic_component(component: str) -> bool:
    component = component.strip()
    if not component:
        return False
    return bool(re.match(r"^[*\-]?\d+[+\-]\d+$", component))


def _parse_intronic_component(component: str) -> Optional[tuple[str, int, str, int]]:
    m = re.match(r"^(?P<prefix>[*\-]?)(?P<base>\d+)(?P<sign>[+\-])(?P<offset>\d+)$", component.strip())
    if not m:
        return None
    return (
        m.group("prefix"),
        int(m.group("base")),
        m.group("sign"),
        int(m.group("offset")),
    )


def _spans_intron(start: Optional[str], stop: Optional[str]) -> bool:
    if start is None or stop is None:
        return False

    start_parsed = _parse_intronic_component(start)
    stop_parsed = _parse_intronic_component(stop)
    if start_parsed is None or stop_parsed is None:
        return False

    start_prefix, start_base, start_sign, _ = start_parsed
    stop_prefix, stop_base, stop_sign, _ = stop_parsed

    # Spanning requires opposite intronic sides around adjacent coding positions
    # (for example: c.76+1_77-1).
    return (
        start_prefix == stop_prefix
        and start_sign != stop_sign
        and abs(stop_base - start_base) == 1
    )


def _aa_3to1(aa_3letter: str) -> str:
    """Convert 3-letter amino acid code to 1-letter (case-insensitive)."""
    aa_3letter = aa_3letter.strip()
    if len(aa_3letter) == 1:
        # Already 1-letter
        if aa_3letter == "*":
            return "*"
        if aa_3letter == "-":
            return "-"
        return aa_3letter.upper()
    
    # Try to convert 3-letter to 1-letter
    key = aa_3letter[0].upper() + aa_3letter[1:].lower()
    if key in _AA_3TO1:
        return _AA_3TO1[key]

    # Unknown amino acid, return as-is but uppercase
    return aa_3letter.upper()


def _reverse_complement(seq: str) -> str:
    complement_table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement_table)[::-1]


def _normalize_protein_allele(allele: Optional[str], is_ref: bool, ref_aa: Optional[str]) -> Optional[str]:
    """Convert protein allele to 1-letter codes.
    
    For synonymous variants (ref == alt), repeat the amino acid.
    """
    if allele is None or allele == "":
        return "" if is_ref else None
    
    allele = allele.strip()
    if not allele:
        return ""
    
    # Handle special cases
    if allele == "dup" or allele == "fs" or allele == "inv":
        return allele
    
    # Handle range deletions like "Ala_Arg"
    if "_" in allele:
        parts = allele.split("_")
        converted_parts = [_aa_3to1(p) for p in parts]
        result = "_".join(converted_parts)
        # If this is a deletion (no ref/alt in range form), use "-"
        return result if not is_ref else result
    
    # Single amino acid
    converted_allele = _aa_3to1(allele)
    
    # For synonymous variants: if ref and alt are the same, repeat the amino acid
    # This will only apply to single-position variants where ref_aa == converted
    return converted_allele


def _parse_nucleotide_hgvs(hgvs_body: str) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    text = hgvs_body.strip()

    m = re.match(
        r"^(?P<coord>.+?)(?:(?P<ref>[A-Za-z*?]+)>(?P<alt>[A-Za-z*?]+)|del(?P<del_seq>[A-Za-z*?]*)ins(?P<delins_alt>[A-Za-z*?]+)|del(?P<del_only>[A-Za-z*?]*)|dup(?P<dup_seq>[A-Za-z*?]*)|ins(?P<ins_seq>[A-Za-z*?]+)|inv(?P<inv_seq>[A-Za-z*?]*)|(?P<eq>=))$",
        text,
    )
    if not m:
        return None, None, None, None

    coord = m.group("coord")
    start, stop = _split_coord(coord)

    if m.group("ref") is not None:
        return start, stop, m.group("ref"), m.group("alt")
    if m.group("delins_alt") is not None:
        return start, stop, (m.group("del_seq") or ""), m.group("delins_alt")
    if m.group("del_only") is not None:
        return start, stop, (m.group("del_only") or ""), ""
    if m.group("dup_seq") is not None:
        ref = m.group("dup_seq") or ""
        return start, stop, ref, (ref + ref) if ref else "dup"
    if m.group("ins_seq") is not None:
        return start, stop, "", m.group("ins_seq")
    if m.group("inv_seq") is not None:
        inv_seq = m.group("inv_seq") or ""
        return start, stop, inv_seq, "inv"
    if m.group("eq") is not None:
        return start, stop, "", ""

    return start, stop, None, None


def _parse_protein_hgvs(hgvs_body: str) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    text = hgvs_body.strip()
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1].strip()

    m_range = re.match(
        r"^(?P<aa1>[A-Za-z*]{1,3})(?P<pos1>\d+)_(?P<aa2>[A-Za-z*]{1,3})(?P<pos2>\d+)(?P<edit>.*)$",
        text,
    )
    if m_range:
        start = m_range.group("pos1")
        stop = m_range.group("pos2")
        edit = m_range.group("edit")
        aa1_1letter = _aa_3to1(m_range.group("aa1"))
        aa2_1letter = _aa_3to1(m_range.group("aa2"))
        ref_range = f"{aa1_1letter}_{aa2_1letter}"
        
        if edit.startswith("delins"):
            alt_part = edit[len("delins"):]
            alt_1letter = _normalize_protein_allele(alt_part, False, None) or alt_part
            return start, stop, ref_range, alt_1letter
        if edit.startswith("del"):
            return start, stop, ref_range, ""
        if edit.startswith("dup"):
            return start, stop, ref_range, ref_range + ref_range
        if edit == "=":
            return start, stop, ref_range, ref_range
        return start, stop, None, None

    m_single = re.match(r"^(?P<ref>[A-Za-z*]{1,3})(?P<pos>\d+)(?P<edit>.*)$", text)
    if not m_single:
        return None, None, None, None

    start = m_single.group("pos")
    stop = start
    ref_3letter = m_single.group("ref")
    ref_1letter = _aa_3to1(ref_3letter)
    edit = m_single.group("edit")

    if not edit:
        return start, stop, ref_1letter, None
    if edit == "=":
        return start, stop, ref_1letter, ref_1letter
    if edit.startswith("delins"):
        alt_part = edit[len("delins"):]
        alt_1letter = _normalize_protein_allele(alt_part, False, ref_1letter) or alt_part
        return start, stop, ref_1letter, alt_1letter
    if edit.startswith("del"):
        return start, stop, ref_1letter, ""
    if edit.startswith("dup"):
        return start, stop, ref_1letter, ref_1letter + ref_1letter
    if edit.startswith("ins"):
        alt_part = edit[len("ins"):]
        alt_1letter = _normalize_protein_allele(alt_part, False, None) or alt_part
        return start, stop, "", alt_1letter
    if edit.startswith("fs"):
        return start, stop, ref_1letter, "fs"

    # Substitution form like p.Pro656Leu -> edit is "Leu"
    if re.match(r"^[A-Za-z*]{1,3}$", edit):
        alt_3letter = edit
        alt_1letter = _aa_3to1(alt_3letter)
        # For synonymous variants, repeat the amino acid
        if ref_1letter == alt_1letter:
            return start, stop, ref_1letter, ref_1letter
        return start, stop, ref_1letter, alt_1letter

    return start, stop, ref_1letter, None


def _extract_chromosome_from_hgvs(hgvs_value: Optional[str]) -> Optional[str]:
    """Extract chromosome number from HGVS accession code.
    
    Examples:
    - "NC_000001.11:g..." -> "1"
    - "NC_000023.11:g..." -> "X"
    - "NC_012920.1:m..." -> "M"
    """
    if _is_blank(hgvs_value):
        return None
    
    text = (hgvs_value or "").strip()
    if ":" not in text:
        return None
    
    accession, _ = text.split(":", 1)
    accession = accession.strip()
    
    # Try to match known RefSeq prefixes
    for prefix, chromosome in _ACCESSION_TO_CHROMOSOME.items():
        if accession.startswith(prefix):
            return chromosome
    
    return None


def _parse_hgvs(
    hgvs_value: Optional[str],
    *,
    resolve_missing_ref_alleles: bool = False,
) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str], bool, bool, Optional[str]]:
    """Parse HGVS and return (start, stop, ref, alt, touches_intronic, spans_intron, chromosome)."""
    if _is_blank(hgvs_value):
        return None, None, None, None, False, False, None

    text = (hgvs_value or "").strip()
    if ":" not in text:
        return None, None, None, None, False, False, None

    _, body = text.split(":", 1)
    body = body.strip()

    if len(body) < 2 or body[1] != ".":
        return None, None, None, None, False, False, None

    coord_type = body[0].lower()
    posedit = body[2:]
    
    # Extract chromosome
    chromosome = _extract_chromosome_from_hgvs(hgvs_value)

    if coord_type == "p":
        start, stop, ref, alt = _parse_protein_hgvs(posedit)
        return start, stop, ref, alt, False, False, chromosome

    start, stop, ref, alt = _parse_nucleotide_hgvs(posedit)

    if (
        resolve_missing_ref_alleles
        and coord_type in {"g", "c", "n"}
        and (ref is None or ref == "")
        and ("del" in posedit or "dup" in posedit or "inv" in posedit)
    ):
        resolved_ref = _resolve_missing_ref_allele(text)
        if resolved_ref is not None:
            ref = resolved_ref

    touches_intronic_region = False
    spans_intron = False
    if coord_type in {"c", "n"} and start is not None and stop is not None:
        touches_intronic_region = _is_intronic_component(start) or _is_intronic_component(stop)
        spans_intron = _spans_intron(start, stop)

    return start, stop, ref, alt, touches_intronic_region, spans_intron, chromosome


def _set_or_append_fieldnames(fieldnames: list[str], col: str) -> None:
    if col not in fieldnames:
        fieldnames.append(col)


def _annotate_row(
    row_idx: int,
    row: dict[str, str],
    *,
    mapped_hgvs_g_col: str,
    mapped_hgvs_c_col: str,
    mapped_hgvs_p_col: str,
    touches_intronic_region_col: str,
    spans_intron_col: str,
    resolve_missing_ref_alleles: bool,
) -> tuple[int, dict[str, str]]:
    """Annotate one row with parsed position/allele/chromosome fields.

    When an HGVS column is pipe-delimited (multiple DNA candidates for a protein
    reverse translation), each segment is parsed independently and the output
    columns are likewise pipe-delimited, preserving candidate cardinality.
    """
    for base_col in (mapped_hgvs_g_col, mapped_hgvs_c_col, mapped_hgvs_p_col):
        raw_value = row.get(base_col) or ""
        segments = raw_value.split("|") if raw_value else [""]

        chromosomes, starts, stops, refs, alts = [], [], [], [], []
        for seg in segments:
            start, stop, ref, alt, _, _, chromosome = _parse_hgvs(
                seg or None,
                resolve_missing_ref_alleles=resolve_missing_ref_alleles,
            )
            chromosomes.append("" if chromosome is None else chromosome)
            starts.append("" if start is None else start)
            stops.append("" if stop is None else stop)
            refs.append("" if ref is None else ref)
            if alt == "inv" and ref:
                alt = _reverse_complement(ref)
            alts.append("" if alt is None else alt)

        row[f"{base_col}_chromosome"] = "|".join(chromosomes)
        row[f"{base_col}_start"] = "|".join(starts)
        row[f"{base_col}_stop"] = "|".join(stops)
        row[f"{base_col}_ref"] = "|".join(refs)
        row[f"{base_col}_alt"] = "|".join(alts)

    # touches_intronic_region / spans_intron: true if any c. candidate fires
    c_raw = row.get(mapped_hgvs_c_col) or ""
    c_segments = c_raw.split("|") if c_raw else [""]
    touches_any = False
    spans_any = False
    for seg in c_segments:
        _, _, _, _, touches, spans, _ = _parse_hgvs(
            seg or None,
            resolve_missing_ref_alleles=resolve_missing_ref_alleles,
        )
        if touches:
            touches_any = True
        if spans:
            spans_any = True
    row[touches_intronic_region_col] = "true" if touches_any else "false"
    row[spans_intron_col] = "true" if spans_any else "false"
    return row_idx, row


def annotate_variants(
    input_file: str,
    output_file: str,
    *,
    mapped_hgvs_g_col: str = "mapped_hgvs_g",
    mapped_hgvs_c_col: str = "mapped_hgvs_c",
    mapped_hgvs_p_col: str = "mapped_hgvs_p",
    touches_intronic_region_col: str = "touches_intronic_region",
    spans_intron_col: str = "spans_intron",
    resolve_missing_ref_alleles: bool = False,
    max_workers: int = 8,
    skip: int = 0,
    limit: Optional[int] = None,
) -> None:
    if max_workers < 1:
        raise ValueError("max_workers must be >= 1")
    if skip < 0:
        raise ValueError("skip must be >= 0")
    if limit is not None and limit < 0:
        raise ValueError("limit must be >= 0")

    in_sep = _detect_separator(input_file)
    out_sep = _detect_separator(output_file)

    with open(input_file, newline="", encoding="utf-8") as in_fh, open(
        output_file, "w", newline="", encoding="utf-8"
    ) as out_fh:
        reader = csv.DictReader(in_fh, delimiter=in_sep)
        if reader.fieldnames is None:
            raise ValueError(f"Input file {input_file!r} appears to be empty.")
        fieldnames = list(reader.fieldnames)

        for base_col in (mapped_hgvs_g_col, mapped_hgvs_c_col, mapped_hgvs_p_col):
            _set_or_append_fieldnames(fieldnames, f"{base_col}_chromosome")
            _set_or_append_fieldnames(fieldnames, f"{base_col}_start")
            _set_or_append_fieldnames(fieldnames, f"{base_col}_stop")
            _set_or_append_fieldnames(fieldnames, f"{base_col}_ref")
            _set_or_append_fieldnames(fieldnames, f"{base_col}_alt")

        _set_or_append_fieldnames(fieldnames, touches_intronic_region_col)
        _set_or_append_fieldnames(fieldnames, spans_intron_col)

        writer = csv.DictWriter(
            out_fh,
            fieldnames=fieldnames,
            delimiter=out_sep,
            extrasaction="ignore",
        )
        writer.writeheader()
        out_fh.flush()

        pending_by_index: dict[int, dict[str, str]] = {}
        next_to_write = 0
        written = 0

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            max_in_flight = max_workers * 4
            sliced_reader = islice(reader, skip, None if limit is None else skip + limit)
            row_iter = enumerate(sliced_reader)
            in_flight: dict = {}

            def _submit_until_full() -> None:
                while len(in_flight) < max_in_flight:
                    try:
                        idx, row = next(row_iter)
                    except StopIteration:
                        break
                    fut = executor.submit(
                        _annotate_row,
                        idx,
                        row,
                        mapped_hgvs_g_col=mapped_hgvs_g_col,
                        mapped_hgvs_c_col=mapped_hgvs_c_col,
                        mapped_hgvs_p_col=mapped_hgvs_p_col,
                        touches_intronic_region_col=touches_intronic_region_col,
                        spans_intron_col=spans_intron_col,
                        resolve_missing_ref_alleles=resolve_missing_ref_alleles,
                    )
                    in_flight[fut] = idx

            _submit_until_full()

            while in_flight:
                done, _ = wait(tuple(in_flight.keys()), return_when=FIRST_COMPLETED)
                for fut in done:
                    in_flight.pop(fut, None)
                    row_idx, row_out = fut.result()
                    pending_by_index[row_idx] = row_out

                while next_to_write in pending_by_index:
                    writer.writerow(pending_by_index.pop(next_to_write))
                    out_fh.flush()
                    written += 1
                    next_to_write += 1

                _submit_until_full()
    logger.info("Wrote %d rows to %s", written, output_file)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Add VCF-format columns (chromosome, position, ref, alt) for mapped HGVS variants.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="Input CSV/TSV file path.")
    p.add_argument("output", help="Output CSV/TSV file path.")
    p.add_argument("--mapped-hgvs-g", default="mapped_hgvs_g", help="Genomic HGVS column name.")
    p.add_argument("--mapped-hgvs-c", default="mapped_hgvs_c", help="Transcript HGVS column name.")
    p.add_argument("--mapped-hgvs-p", default="mapped_hgvs_p", help="Protein HGVS column name.")
    p.add_argument(
        "--touches-intronic-region-column",
        dest="touches_intronic_region_col",
        default="touches_intronic_region",
        help="Output boolean column indicating transcript HGVS touches an intronic region.",
    )
    p.add_argument(
        "--spans-intron-column",
        dest="spans_intron_col",
        default="spans_intron",
        help="Output boolean column indicating transcript HGVS spans both sides of an intron.",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity level.",
    )
    p.add_argument(
        "--resolve-missing-ref-alleles",
        dest="resolve_missing_ref_alleles",
        action="store_true",
        default=True,
        help="Use HGVS normalization to fill missing ref allele for accession-backed del/delins-like variants.",
    )
    p.add_argument(
        "--max-workers",
        type=int,
        default=8,
        help="Concurrent worker threads for row annotation.",
    )
    p.add_argument(
        "--skip",
        type=int,
        default=0,
        help="Number of data rows to skip from the start of the input.",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of data rows to process after applying --skip.",
    )
    return p


def main(argv: Optional[list[str]] = None) -> None:
    args = _build_parser().parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    annotate_variants(
        args.input,
        args.output,
        mapped_hgvs_g_col=args.mapped_hgvs_g,
        mapped_hgvs_c_col=args.mapped_hgvs_c,
        mapped_hgvs_p_col=args.mapped_hgvs_p,
        touches_intronic_region_col=args.touches_intronic_region_col,
        spans_intron_col=args.spans_intron_col,
        resolve_missing_ref_alleles=args.resolve_missing_ref_alleles,
        max_workers=args.max_workers,
        skip=args.skip,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
