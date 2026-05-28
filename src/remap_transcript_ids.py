"""Replace Ensembl transcript/protein accessions with RefSeq equivalents (or vice versa).

Operates on HGVS-containing columns produced by earlier pipeline steps:
  - mapped_hgvs_g, mapped_hgvs_c, mapped_hgvs_p
  - mapped_hgvs_c_transcript  (transcript accession extracted by add_vcf_identifiers)
  - mapped_hgvs_p_protein      (protein accession extracted by add_vcf_identifiers)

For each of those columns, any accession that appears in the mapping table is swapped
for its equivalent; all other values are passed through unchanged.

Pipe-delimited values (multiple candidates per cell) are handled element-by-element.

---

Mapping source
--------------
Two input modes are supported:

1. MANE Select summary file from NCBI:
     https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.*.summary.txt.gz
   This file contains RefSeq↔Ensembl equivalences for MANE Select and MANE Plus Clinical
   transcripts, which are guaranteed to be sequence-identical across namespaces.
   Pass via --mane-file.  Control the swap direction with --direction.

2. Generic two-column TSV/CSV mapping file (e.g. from another database or a custom list)
   with columns ``source_id`` and ``target_id`` (header required).
   Pass via --mapping-file.  Direction is fixed by which column is ``source_id``.

MANE file column names expected (tab-delimited, first column header starts with #):
  RefSeq_nuc, RefSeq_prot, Ensembl_nuc, Ensembl_prot
  (plus MANE_status to allow optional filtering to MANE Select only)

---

Accession version matching
--------------------------
Exact versioned match is tried first (e.g. NM_007294.4 → ENST00000357654.9).
If the versioned accession is not found, the version suffix is stripped and the
base identifier is looked up instead; the target version from the mapping file is
then used.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Mapping helpers
# ---------------------------------------------------------------------------

def _open_maybe_gzipped(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, newline="", encoding="utf-8")


def _build_mane_mapping(mane_file: str, direction: str) -> dict[str, str]:
    """Build {source_accession: target_accession} from an NCBI MANE summary file.

    Args:
        mane_file: Path to MANE.GRCh38.*.summary.txt[.gz].
        direction: ``"to-refseq"`` or ``"to-ensembl"``.

    Returns:
        Mapping dictionary (versioned accessions as keys and values).
    """
    mapping: dict[str, str] = {}
    with _open_maybe_gzipped(mane_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # NCBI uses "#NCBI_GeneID" as the first column header; csv.DictReader
        # treats the leading # as part of the name.
        for row in reader:
            refseq_nuc = (row.get("RefSeq_nuc") or "").strip()
            refseq_prot = (row.get("RefSeq_prot") or "").strip()
            ensembl_nuc = (row.get("Ensembl_nuc") or "").strip()
            ensembl_prot = (row.get("Ensembl_prot") or "").strip()

            if direction == "to-ensembl":
                if refseq_nuc and ensembl_nuc:
                    mapping[refseq_nuc] = ensembl_nuc
                if refseq_prot and ensembl_prot:
                    mapping[refseq_prot] = ensembl_prot
            else:  # to-refseq
                if ensembl_nuc and refseq_nuc:
                    mapping[ensembl_nuc] = refseq_nuc
                if ensembl_prot and refseq_prot:
                    mapping[ensembl_prot] = refseq_prot

    logger.info("Loaded %d accession mappings from MANE file (direction=%s).", len(mapping), direction)
    return mapping


def _build_custom_mapping(mapping_file: str) -> dict[str, str]:
    """Build {source_id: target_id} from a two-column TSV/CSV with headers source_id, target_id."""
    mapping: dict[str, str] = {}
    sep = "\t" if mapping_file.endswith(".tsv") else ","
    with open(mapping_file, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        required = {"source_id", "target_id"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"Custom mapping file must have 'source_id' and 'target_id' columns; "
                f"got: {reader.fieldnames}"
            )
        for row in reader:
            src = (row.get("source_id") or "").strip()
            tgt = (row.get("target_id") or "").strip()
            if src and tgt:
                mapping[src] = tgt

    logger.info("Loaded %d accession mappings from custom mapping file.", len(mapping))
    return mapping


def _build_base_mapping(mapping: dict[str, str]) -> dict[str, str]:
    """Build a version-stripped fallback mapping from a versioned mapping.

    For each entry ``ACC.N → TARGET.M``, adds ``ACC → TARGET.M`` if ``ACC``
    is not already present as a key (versioned takes priority).
    """
    base: dict[str, str] = {}
    for src, tgt in mapping.items():
        if "." in src:
            stripped = src.rsplit(".", 1)[0]
            if stripped not in base:
                base[stripped] = tgt
    return base


# ---------------------------------------------------------------------------
# Per-value remapping
# ---------------------------------------------------------------------------

# Accession prefixes that belong to a remappable namespace (RefSeq transcript/
# protein or Ensembl transcript/protein).  Accessions with these prefixes that
# are not found in the mapping table are reported as unmapped.
_REMAPPABLE_PREFIXES = ("NM_", "NR_", "NP_", "ENST", "ENSP")


def _looks_remappable(accession: str) -> bool:
    """Return True if *accession* belongs to a namespace we know how to remap."""
    return any(accession.startswith(p) for p in _REMAPPABLE_PREFIXES)


def _remap_accession(
    value: str,
    mapping: dict[str, str],
    base_mapping: dict[str, str],
    unmapped: set[str],
) -> str:
    """Return the mapped accession for *value*, or *value* unchanged.

    Accessions that look remappable but are absent from the mapping are added
    to *unmapped* for end-of-run reporting.
    """
    if not value:
        return value
    if value in mapping:
        return mapping[value]
    # Version-stripped fallback
    if "." in value:
        stripped = value.rsplit(".", 1)[0]
        if stripped in base_mapping:
            return base_mapping[stripped]
    if _looks_remappable(value):
        unmapped.add(value)
    return value


def _remap_hgvs(
    hgvs_str: str,
    mapping: dict[str, str],
    base_mapping: dict[str, str],
    unmapped: set[str],
) -> str:
    """Replace the accession prefix of a single HGVS string."""
    if not hgvs_str or ":" not in hgvs_str:
        return hgvs_str
    accession, rest = hgvs_str.split(":", 1)
    new_accession = _remap_accession(accession, mapping, base_mapping, unmapped)
    return f"{new_accession}:{rest}"


def _remap_cell(
    cell_value: str,
    mapping: dict[str, str],
    base_mapping: dict[str, str],
    is_hgvs: bool,
    unmapped: set[str],
) -> str:
    """Remap one TSV cell, handling pipe-delimited multivalues."""
    if not cell_value:
        return cell_value
    parts = cell_value.split("|")
    remapped = [
        _remap_hgvs(p, mapping, base_mapping, unmapped)
        if is_hgvs
        else _remap_accession(p, mapping, base_mapping, unmapped)
        for p in parts
    ]
    return "|".join(remapped)


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

_HGVS_COLS = ("mapped_hgvs_g", "mapped_hgvs_c", "mapped_hgvs_p")
_ACCESSION_COLS = ("mapped_hgvs_c_transcript", "mapped_hgvs_p_protein")


def remap_transcript_ids(
    input_file: str,
    output_file: str,
    mapping: dict[str, str],
    *,
    hgvs_cols: tuple[str, ...] = _HGVS_COLS,
    accession_cols: tuple[str, ...] = _ACCESSION_COLS,
) -> None:
    """Rewrite *input_file* to *output_file*, remapping transcript/protein accessions.

    Args:
        input_file: Path to input TSV/CSV.
        output_file: Path to output TSV/CSV.
        mapping: Versioned {source: target} accession mapping.
        hgvs_cols: Columns containing full HGVS strings (accession before ``:``)
        accession_cols: Columns containing bare accession values only.
    """
    base_mapping = _build_base_mapping(mapping)

    in_sep = "\t" if input_file.endswith(".tsv") else ","
    out_sep = "\t" if output_file.endswith(".tsv") else ","

    replaced_total = 0
    rows_written = 0
    unmapped: set[str] = set()

    with open(input_file, newline="", encoding="utf-8") as in_fh, \
         open(output_file, "w", newline="", encoding="utf-8") as out_fh:

        reader = csv.DictReader(in_fh, delimiter=in_sep)
        if reader.fieldnames is None:
            raise ValueError(f"Input file {input_file!r} appears to be empty.")

        writer = csv.DictWriter(out_fh, fieldnames=reader.fieldnames, delimiter=out_sep, extrasaction="ignore")
        writer.writeheader()

        present_hgvs = [c for c in hgvs_cols if c in reader.fieldnames]
        present_acc = [c for c in accession_cols if c in reader.fieldnames]

        if not present_hgvs and not present_acc:
            logger.warning(
                "None of the expected columns (%s, %s) found in input; output will be unchanged.",
                ", ".join(hgvs_cols),
                ", ".join(accession_cols),
            )

        for row in reader:
            for col in present_hgvs:
                original = row[col]
                remapped = _remap_cell(original, mapping, base_mapping, is_hgvs=True, unmapped=unmapped)
                if remapped != original:
                    replaced_total += 1
                row[col] = remapped

            for col in present_acc:
                original = row[col]
                remapped = _remap_cell(original, mapping, base_mapping, is_hgvs=False, unmapped=unmapped)
                if remapped != original:
                    replaced_total += 1
                row[col] = remapped

            writer.writerow(row)
            rows_written += 1

    logger.info("Wrote %d rows; replaced accessions in %d cells.", rows_written, replaced_total)

    if unmapped:
        logger.warning(
            "%d remappable accession(s) were not found in the mapping table and were left unchanged:\n  %s",
            len(unmapped),
            "\n  ".join(sorted(unmapped)),
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Replace Ensembl or RefSeq transcript/protein accessions in HGVS columns "
            "using a MANE Select mapping or a custom two-column mapping file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="Input TSV/CSV file.")
    p.add_argument("output", help="Output TSV/CSV file.")

    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--mane-file",
        metavar="FILE",
        help=(
            "NCBI MANE Select summary file (MANE.GRCh38.*.summary.txt[.gz]). "
            "Use --direction to control swap direction."
        ),
    )
    source.add_argument(
        "--mapping-file",
        metavar="FILE",
        help=(
            "Custom two-column TSV/CSV file with headers 'source_id' and 'target_id'. "
            "Direction is determined by which accessions appear in the source_id column."
        ),
    )

    p.add_argument(
        "--direction",
        choices=["to-refseq", "to-ensembl"],
        default="to-refseq",
        help="Swap direction when using --mane-file. Ignored for --mapping-file.",
    )
    p.add_argument(
        "--hgvs-col",
        action="append",
        dest="hgvs_cols",
        metavar="COL",
        help=(
            "HGVS column to remap (accession is the part before ':'). "
            "May be repeated. Defaults to: mapped_hgvs_g, mapped_hgvs_c, mapped_hgvs_p."
        ),
    )
    p.add_argument(
        "--accession-col",
        action="append",
        dest="accession_cols",
        metavar="COL",
        help=(
            "Plain-accession column to remap (the whole cell value is treated as an accession). "
            "May be repeated. Defaults to: mapped_hgvs_c_transcript, mapped_hgvs_p_protein."
        ),
    )
    p.add_argument(
        "--csv-field-size-limit",
        type=int,
        default=131072,
        help="Maximum field size for CSV/TSV parsing.",
    )
    p.add_argument(
        "--log-level",
        default="WARNING",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity.",
    )
    return p


def main(argv: Optional[list[str]] = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    csv.field_size_limit(args.csv_field_size_limit)

    # Build mapping
    if args.mane_file:
        mapping = _build_mane_mapping(args.mane_file, args.direction)
    else:
        mapping = _build_custom_mapping(args.mapping_file)

    hgvs_cols = tuple(args.hgvs_cols) if args.hgvs_cols else _HGVS_COLS
    accession_cols = tuple(args.accession_cols) if args.accession_cols else _ACCESSION_COLS

    remap_transcript_ids(
        args.input,
        args.output,
        mapping,
        hgvs_cols=hgvs_cols,
        accession_cols=accession_cols,
    )


if __name__ == "__main__":
    main()
