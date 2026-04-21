"""Add parsed position/allele columns for mapped HGVS fields.

Given mapped HGVS columns (genomic/transcript/protein), this script appends:
- <col>_start
- <col>_stop
- <col>_ref
- <col>_alt

It also adds transcript-derived boolean columns:
- touches_intronic_region: transcript HGVS contains intronic offset coordinates
- spans_intron: transcript HGVS spans both sides of an intron
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_HGVS_RESOLVER_UNAVAILABLE_LOGGED = False
_HGVS_REF_RESOLVER: Optional["_HgvsRefResolver"] = None


class _HgvsRefResolver:
    def __init__(self) -> None:
        import hgvs.dataproviders.uta  # type: ignore[import-not-found]
        import hgvs.normalizer  # type: ignore[import-not-found]
        import hgvs.parser  # type: ignore[import-not-found]

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
    global _HGVS_REF_RESOLVER_UNAVAILABLE_LOGGED
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
        if edit.startswith("delins"):
            return start, stop, f"{m_range.group('aa1')}_{m_range.group('aa2')}", edit[len("delins") :]
        if edit.startswith("del"):
            return start, stop, f"{m_range.group('aa1')}_{m_range.group('aa2')}", ""
        if edit.startswith("dup"):
            ref = f"{m_range.group('aa1')}_{m_range.group('aa2')}"
            return start, stop, ref, "dup"
        if edit == "=":
            return start, stop, "", ""
        return start, stop, None, None

    m_single = re.match(r"^(?P<ref>[A-Za-z*]{1,3})(?P<pos>\d+)(?P<edit>.*)$", text)
    if not m_single:
        return None, None, None, None

    start = m_single.group("pos")
    stop = start
    ref = m_single.group("ref")
    edit = m_single.group("edit")

    if not edit:
        return start, stop, ref, None
    if edit == "=":
        return start, stop, ref, ref
    if edit.startswith("delins"):
        return start, stop, ref, edit[len("delins") :]
    if edit.startswith("del"):
        return start, stop, ref, ""
    if edit.startswith("dup"):
        return start, stop, ref, ref + ref
    if edit.startswith("ins"):
        return start, stop, "", edit[len("ins") :]
    if edit.startswith("fs"):
        return start, stop, ref, "fs"

    # Substitution form like p.Pro656Leu -> edit is "Leu"
    if re.match(r"^[A-Za-z*]{1,3}$", edit):
        return start, stop, ref, edit

    return start, stop, ref, None


def _parse_hgvs(
    hgvs_value: Optional[str],
    *,
    resolve_missing_ref_alleles: bool = False,
) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str], bool, bool]:
    if _is_blank(hgvs_value):
        return None, None, None, None, False, False

    text = (hgvs_value or "").strip()
    if ":" not in text:
        return None, None, None, None, False, False

    _, body = text.split(":", 1)
    body = body.strip()

    if len(body) < 2 or body[1] != ".":
        return None, None, None, None, False, False

    coord_type = body[0].lower()
    posedit = body[2:]

    if coord_type == "p":
        start, stop, ref, alt = _parse_protein_hgvs(posedit)
        return start, stop, ref, alt, False, False

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

    return start, stop, ref, alt, touches_intronic_region, spans_intron


def _set_or_append_fieldnames(fieldnames: list[str], col: str) -> None:
    if col not in fieldnames:
        fieldnames.append(col)


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
) -> None:
    in_sep = _detect_separator(input_file)
    out_sep = _detect_separator(output_file)

    with open(input_file, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=in_sep)
        if reader.fieldnames is None:
            raise ValueError(f"Input file {input_file!r} appears to be empty.")
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    for base_col in (mapped_hgvs_g_col, mapped_hgvs_c_col, mapped_hgvs_p_col):
        _set_or_append_fieldnames(fieldnames, f"{base_col}_start")
        _set_or_append_fieldnames(fieldnames, f"{base_col}_stop")
        _set_or_append_fieldnames(fieldnames, f"{base_col}_ref")
        _set_or_append_fieldnames(fieldnames, f"{base_col}_alt")

    _set_or_append_fieldnames(fieldnames, touches_intronic_region_col)
    _set_or_append_fieldnames(fieldnames, spans_intron_col)

    for row in rows:
        for base_col in (mapped_hgvs_g_col, mapped_hgvs_c_col, mapped_hgvs_p_col):
            start, stop, ref, alt, _, _ = _parse_hgvs(
                row.get(base_col),
                resolve_missing_ref_alleles=resolve_missing_ref_alleles,
            )
            row[f"{base_col}_start"] = "" if start is None else start
            row[f"{base_col}_stop"] = "" if stop is None else stop
            row[f"{base_col}_ref"] = "" if ref is None else ref
            row[f"{base_col}_alt"] = "" if alt is None else alt

        _, _, _, _, touches_intronic_region, spans_intron = _parse_hgvs(
            row.get(mapped_hgvs_c_col),
            resolve_missing_ref_alleles=resolve_missing_ref_alleles,
        )
        row[touches_intronic_region_col] = "true" if touches_intronic_region else "false"
        row[spans_intron_col] = "true" if spans_intron else "false"

    with open(output_file, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=out_sep, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Wrote %d rows to %s", len(rows), output_file)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Add parsed start/stop/ref/alt fields for mapped HGVS columns.",
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
    )


if __name__ == "__main__":
    main()
