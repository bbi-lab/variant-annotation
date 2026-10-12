"""Lift genomic positions between genome assemblies.

Reads a CSV/TSV file and maps coordinates from one genome assembly to another
using a UCSC chain file.  Two modes are available and may be combined on the
same input file:

**HGVS mode** (``--input-col``): lifts genomic HGVS g. strings such as
``NC_000001.10:g.925952C>T``.  Supported variant types: substitutions,
deletions, insertions, delins, duplications, and inversions.  Input strings
may use RefSeq accession identifiers (NC_000001.10) or plain chromosome names
(chr1); the output format mirrors the input.  Pipe-delimited multi-candidate
values are lifted independently and results are pipe-aligned to the input.

**Coordinate mode** (``--chrom-col``): lifts raw chromosome + position columns
such as those produced by a BED or region file.  Positions are 1-based
inclusive.  ``--end-col`` is optional; omit it for single-position rows.

Chain files are downloaded automatically from UCSC on first use and cached in
``~/.pyliftover/``; pass ``--chain-file`` to supply a pre-downloaded file.

Requires the ``pyliftover`` package::

    pip install pyliftover
    # or: pip install variant-annotation[liftover]

Usage::

    # HGVS mode
    python -m src.liftover input.tsv output.tsv \\
        --from-assembly hg19 --to-assembly hg38 \\
        --input-col hgvs_g

    # Coordinate mode
    python -m src.liftover regions.tsv lifted.tsv \\
        --from-assembly hg19 --to-assembly hg38 \\
        --chrom-col chrom --start-col start --end-col end

    # Both modes on the same file
    python -m src.liftover variants.tsv lifted.tsv \\
        --from-assembly hg19 --to-assembly hg38 \\
        --input-col hgvs_g --output-col hgvs_g_hg38 \\
        --chrom-col chrom --start-col start --end-col end

    # CSV, custom assembly names, pre-downloaded chain file
    python -m src.liftover variants.csv lifted.csv \\
        --delimiter , \\
        --from-assembly GRCh37 --to-assembly GRCh38 \\
        --chain-file /data/hg19ToHg38.over.chain.gz \\
        --input-col mapped_hgvs_g --output-col hgvs_g_hg38
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from itertools import islice
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

load_dotenv()
logger = logging.getLogger(__name__)

PROGRESS_EVERY_ROWS = 1000

# Canonical primary-assembly RefSeq accessions per assembly + chromosome.
_REFSEQ: dict[str, dict[str, str]] = {
    "hg19": {
        "chr1": "NC_000001.10", "chr2": "NC_000002.11", "chr3": "NC_000003.11",
        "chr4": "NC_000004.11", "chr5": "NC_000005.9",  "chr6": "NC_000006.11",
        "chr7": "NC_000007.13", "chr8": "NC_000008.10", "chr9": "NC_000009.11",
        "chr10": "NC_000010.10", "chr11": "NC_000011.9", "chr12": "NC_000012.11",
        "chr13": "NC_000013.10", "chr14": "NC_000014.8", "chr15": "NC_000015.9",
        "chr16": "NC_000016.9", "chr17": "NC_000017.10", "chr18": "NC_000018.9",
        "chr19": "NC_000019.9", "chr20": "NC_000020.10", "chr21": "NC_000021.8",
        "chr22": "NC_000022.10", "chrX": "NC_000023.10", "chrY": "NC_000024.9",
        "chrM": "NC_012920.1",
    },
    "hg38": {
        "chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12",
        "chr4": "NC_000004.12", "chr5": "NC_000005.10", "chr6": "NC_000006.12",
        "chr7": "NC_000007.14", "chr8": "NC_000008.11", "chr9": "NC_000009.12",
        "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
        "chr13": "NC_000013.11", "chr14": "NC_000014.9", "chr15": "NC_000015.10",
        "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
        "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
        "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10",
        "chrM": "NC_012920.1",
    },
}

# Reverse map: RefSeq accession → chromosome name (unique across both assemblies).
_ACCESSION_TO_CHR: dict[str, str] = {
    acc: chrom
    for _asm, chrs in _REFSEQ.items()
    for chrom, acc in chrs.items()
}

_ASSEMBLY_ALIASES: dict[str, str] = {
    "hg19": "hg19", "grch37": "hg19", "grch37p13": "hg19",
    "hg38": "hg38", "grch38": "hg38", "grch38p14": "hg38",
}

# Matches reference:g.START[_END]VARIANT_SUFFIX
# Groups: (1) reference seq, (2) start pos, (3) end pos or None, (4) variant suffix
_HGVS_G_RE = re.compile(r"^([^:]+):g\.(\d+)(?:_(\d+))?(.*)", re.IGNORECASE)


def _normalize_assembly(name: str) -> str:
    key = name.lower().replace(".", "").replace("-", "")
    canonical = _ASSEMBLY_ALIASES.get(key)
    if canonical is None:
        raise ValueError(f"Unknown assembly '{name}'. Supported: hg19, GRCh37, hg38, GRCh38")
    return canonical


def _ref_to_chr(ref_seq: str) -> tuple[str, bool]:
    """Return (chr_name, used_accession). used_accession is True for NC_ inputs."""
    upper = ref_seq.upper()
    if upper.startswith("NC_"):
        chrom = _ACCESSION_TO_CHR.get(ref_seq)
        if chrom is None:
            raise ValueError(f"Unrecognized RefSeq accession '{ref_seq}'")
        return chrom, True
    if upper.startswith("CHR"):
        return ref_seq, False
    # Bare number like "1" → "chr1"
    return f"chr{ref_seq}", False


def _chr_to_ref(chrom: str, assembly: str, use_accession: bool) -> str:
    if not use_accession:
        return chrom
    acc = _REFSEQ.get(assembly, {}).get(chrom)
    if acc is None:
        raise ValueError(f"No RefSeq accession for '{chrom}' in {assembly}")
    return acc


def _normalize_chrom(value: str) -> str:
    """Return a chr-prefixed chromosome name from any common notation."""
    v = value.strip()
    if v.upper().startswith("CHR"):
        return v
    return f"chr{v}"


# ---------------------------------------------------------------------------
# HGVS mode
# ---------------------------------------------------------------------------

def _lift_hgvs_single(hgvs: str, lifter, from_asm: str, to_asm: str) -> tuple[str, str]:
    """Lift one HGVS g. string. Returns (lifted_hgvs, error)."""
    hgvs = hgvs.strip()
    if not hgvs:
        return "", ""

    m = _HGVS_G_RE.match(hgvs)
    if m is None:
        return "", f"Cannot parse as HGVS g.: {hgvs!r}"

    ref_seq, start_str, end_str, variant_suffix = m.groups()

    try:
        chrom, use_accession = _ref_to_chr(ref_seq)
    except ValueError as exc:
        return "", str(exc)

    # pyliftover expects 0-based coordinates.
    mapped_start = lifter.convert_coordinate(chrom, int(start_str) - 1)
    if not mapped_start:
        return "", f"Position {chrom}:{start_str} did not map to {to_asm}"

    new_chrom, new_start_0, new_strand = mapped_start[0][:3]
    if new_strand == "-":
        return "", f"Strand flip at {chrom}:{start_str}; manual review required"

    if end_str is not None:
        mapped_end = lifter.convert_coordinate(chrom, int(end_str) - 1)
        if not mapped_end:
            return "", f"End position {chrom}:{end_str} did not map to {to_asm}"
        end_chrom, new_end_0, end_strand = mapped_end[0][:3]
        if end_chrom != new_chrom:
            return "", "Start and end map to different chromosomes after liftover"
        if end_strand == "-":
            return "", f"Strand flip at end position {chrom}:{end_str}; manual review required"
        pos_part = f"{new_start_0 + 1}_{new_end_0 + 1}"
    else:
        pos_part = str(new_start_0 + 1)

    try:
        new_ref = _chr_to_ref(new_chrom, to_asm, use_accession)
    except ValueError as exc:
        return "", str(exc)

    return f"{new_ref}:g.{pos_part}{variant_suffix}", ""


def _lift_hgvs_field(field: str, lifter, from_asm: str, to_asm: str) -> tuple[str, str]:
    """Lift a possibly pipe-delimited HGVS g. field, preserving pipe alignment."""
    if "|" not in field:
        return _lift_hgvs_single(field, lifter, from_asm, to_asm)
    lifted_parts: list[str] = []
    error_parts: list[str] = []
    for part in field.split("|"):
        lifted, err = _lift_hgvs_single(part, lifter, from_asm, to_asm)
        lifted_parts.append(lifted)
        error_parts.append(err)
    return "|".join(lifted_parts), "|".join(error_parts)


# ---------------------------------------------------------------------------
# Coordinate mode
# ---------------------------------------------------------------------------

def _lift_coords(
    chrom: str,
    start: str,
    end: str | None,
    lifter,
    to_asm: str,
) -> tuple[str, str, str, str]:
    """Lift raw 1-based genomic coordinates.

    Returns (new_chrom, new_start, new_end, error).  new_end is empty when
    end is not provided or the input end is blank.
    """
    chrom = chrom.strip()
    start = start.strip()
    if not chrom or not start:
        return "", "", "", ""

    norm_chrom = _normalize_chrom(chrom)

    try:
        start_1 = int(start)
    except ValueError:
        return "", "", "", f"Invalid start position: {start!r}"

    mapped_start = lifter.convert_coordinate(norm_chrom, start_1 - 1)
    if not mapped_start:
        return "", "", "", f"Position {norm_chrom}:{start_1} did not map to {to_asm}"

    new_chrom, new_start_0, new_strand = mapped_start[0][:3]
    if new_strand == "-":
        logger.debug("Strand flip at %s:%d after liftover (coordinate mapped correctly)", norm_chrom, start_1)

    new_end_str = ""
    if end is not None and end.strip():
        try:
            end_1 = int(end.strip())
        except ValueError:
            return "", "", "", f"Invalid end position: {end!r}"

        mapped_end = lifter.convert_coordinate(norm_chrom, end_1 - 1)
        if not mapped_end:
            return "", "", "", f"End position {norm_chrom}:{end_1} did not map to {to_asm}"

        end_new_chrom, new_end_0, _ = mapped_end[0][:3]
        if end_new_chrom != new_chrom:
            return "", "", "", "Start and end map to different chromosomes after liftover"

        new_end_str = str(new_end_0 + 1)

    return new_chrom, str(new_start_0 + 1), new_end_str, ""


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Lift genomic positions between genome assemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input_file", help="Input CSV/TSV file")
    p.add_argument("output_file", help="Output CSV/TSV file")

    asm = p.add_argument_group("assembly")
    asm.add_argument("--from-assembly", default="hg19", metavar="ASM",
                     help="Source assembly (hg19/GRCh37 or hg38/GRCh38)")
    asm.add_argument("--to-assembly", default="hg38", metavar="ASM",
                     help="Target assembly (hg19/GRCh37 or hg38/GRCh38)")
    asm.add_argument("--chain-file", default=None, metavar="PATH",
                     help="UCSC .over.chain or .over.chain.gz file; auto-downloaded if omitted")

    hgvs = p.add_argument_group("HGVS mode (activated by --input-col)")
    hgvs.add_argument("--input-col", default=None, metavar="COL",
                      help="Column containing HGVS g. strings to lift")
    hgvs.add_argument("--output-col", default=None, metavar="COL",
                      help="Column to write lifted HGVS g. strings "
                           "(default: <input-col>_<to-assembly>)")
    hgvs.add_argument("--error-col", default="liftover_error", metavar="COL",
                      help="Column to write HGVS liftover errors")

    coord = p.add_argument_group("coordinate mode (activated by --chrom-col)")
    coord.add_argument("--chrom-col", default=None, metavar="COL",
                       help="Column containing chromosome names (chr1 or 1)")
    coord.add_argument("--start-col", default=None, metavar="COL",
                       help="Column containing 1-based start positions")
    coord.add_argument("--end-col", default=None, metavar="COL",
                       help="Column containing 1-based end positions (optional)")
    coord.add_argument("--lifted-chrom-col", default=None, metavar="COL",
                       help="Output column for lifted chromosome (default: <chrom-col>_<to-assembly>)")
    coord.add_argument("--lifted-start-col", default=None, metavar="COL",
                       help="Output column for lifted start position (default: <start-col>_<to-assembly>)")
    coord.add_argument("--lifted-end-col", default=None, metavar="COL",
                       help="Output column for lifted end position (default: <end-col>_<to-assembly>)")
    coord.add_argument("--coord-error-col", default="coord_liftover_error", metavar="COL",
                       help="Column to write coordinate liftover errors")

    io = p.add_argument_group("I/O")
    io.add_argument("--skip", type=int, default=0,
                    help="Number of data rows to skip before processing")
    io.add_argument("--limit", type=int, default=None,
                    help="Maximum number of data rows to process")
    io.add_argument("--delimiter", default="\t", metavar="DELIM",
                    help="Field delimiter (default: TAB; use ',' for CSV)")
    io.add_argument("--csv-field-size-limit", type=int, default=csv.field_size_limit(),
                    help="Maximum field size in bytes")
    io.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:  # noqa: C901
    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    csv.field_size_limit(args.csv_field_size_limit)

    hgvs_mode = args.input_col is not None
    coord_mode = args.chrom_col is not None

    if not hgvs_mode and not coord_mode:
        logger.error("Specify at least one of --input-col (HGVS mode) or --chrom-col (coordinate mode)")
        sys.exit(1)

    if coord_mode and not args.start_col:
        logger.error("--start-col is required when --chrom-col is set")
        sys.exit(1)

    if args.skip < 0:
        logger.error("--skip must be >= 0, got: %d", args.skip)
        sys.exit(1)
    if args.limit is not None and args.limit < 1:
        logger.error("--limit must be >= 1 when provided, got: %d", args.limit)
        sys.exit(1)

    try:
        from_asm = _normalize_assembly(args.from_assembly)
        to_asm = _normalize_assembly(args.to_assembly)
    except ValueError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    # Derive default output column names.
    hgvs_output_col = args.output_col or (f"{args.input_col}_{to_asm}" if hgvs_mode else None)
    lifted_chrom_col = args.lifted_chrom_col or (f"{args.chrom_col}_{to_asm}" if coord_mode else None)
    lifted_start_col = args.lifted_start_col or (f"{args.start_col}_{to_asm}" if coord_mode else None)
    lifted_end_col = (
        args.lifted_end_col or (f"{args.end_col}_{to_asm}" if coord_mode and args.end_col else None)
    )

    try:
        from pyliftover import LiftOver
    except ImportError:
        logger.error(
            "pyliftover is not installed. "
            "Install with: pip install pyliftover  (or: pip install variant-annotation[liftover])"
        )
        sys.exit(1)

    logger.info("Loading chain: %s → %s", from_asm, to_asm)
    try:
        lifter = LiftOver(args.chain_file) if args.chain_file else LiftOver(from_asm, to_asm)
    except Exception as exc:
        logger.error("Failed to load liftover chain: %s", exc)
        sys.exit(1)

    delim = "\t" if args.delimiter == "\\t" else args.delimiter

    input_path = Path(args.input_file)
    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    hgvs_lifted = 0
    hgvs_errors = 0
    coord_lifted = 0
    coord_errors = 0

    with (
        input_path.open("r", encoding="utf-8", newline="") as in_fh,
        output_path.open("w", encoding="utf-8", newline="") as out_fh,
    ):
        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears to be empty: %s", input_path)
            sys.exit(1)

        fieldnames = list(reader.fieldnames)

        # Validate that required input columns exist.
        required = []
        if hgvs_mode:
            required.append(args.input_col)
        if coord_mode:
            required += [args.chrom_col, args.start_col]
            if args.end_col:
                required.append(args.end_col)
        missing = [c for c in required if c not in fieldnames]
        if missing:
            logger.error("Required column(s) not found: %s. Available: %s",
                         ", ".join(missing), ", ".join(fieldnames))
            sys.exit(1)

        # Build output fieldnames, appending new columns in a consistent order.
        out_fieldnames = fieldnames[:]
        new_cols: list[str] = []
        if hgvs_mode:
            new_cols += [hgvs_output_col, args.error_col]
        if coord_mode:
            new_cols += [lifted_chrom_col, lifted_start_col]
            if lifted_end_col:
                new_cols.append(lifted_end_col)
            new_cols.append(args.coord_error_col)
        for col in new_cols:
            if col not in out_fieldnames:
                out_fieldnames.append(col)

        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=delim,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        window = islice(reader, args.skip, None if args.limit is None else args.skip + args.limit)
        for row in window:
            if hgvs_mode:
                field_value = row.get(args.input_col, "")
                lifted_hgvs, hgvs_err = _lift_hgvs_field(field_value, lifter, from_asm, to_asm)
                row[hgvs_output_col] = lifted_hgvs
                row[args.error_col] = hgvs_err
                if hgvs_err.replace("|", "").strip():
                    hgvs_errors += 1
                    logger.debug("Row %d HGVS: %s", total_rows + 1, hgvs_err)
                elif lifted_hgvs.replace("|", "").strip():
                    hgvs_lifted += 1

            if coord_mode:
                new_chrom, new_start, new_end, coord_err = _lift_coords(
                    row.get(args.chrom_col, ""),
                    row.get(args.start_col, ""),
                    row.get(args.end_col, "") if args.end_col else None,
                    lifter,
                    to_asm,
                )
                row[lifted_chrom_col] = new_chrom
                row[lifted_start_col] = new_start
                if lifted_end_col:
                    row[lifted_end_col] = new_end
                row[args.coord_error_col] = coord_err
                if coord_err:
                    coord_errors += 1
                    logger.debug("Row %d coords: %s", total_rows + 1, coord_err)
                elif new_chrom:
                    coord_lifted += 1

            writer.writerow(row)
            total_rows += 1
            if total_rows % PROGRESS_EVERY_ROWS == 0:
                parts = [f"{total_rows} rows"]
                if hgvs_mode:
                    parts.append(f"HGVS: {hgvs_lifted} lifted, {hgvs_errors} errors")
                if coord_mode:
                    parts.append(f"coords: {coord_lifted} lifted, {coord_errors} errors")
                logger.info("Processed %s", "; ".join(parts))

    summary_parts = [f"Processed {total_rows} rows."]
    if hgvs_mode:
        summary_parts.append(f"HGVS: {hgvs_lifted} lifted, {hgvs_errors} errors.")
    if coord_mode:
        summary_parts.append(f"Coords: {coord_lifted} lifted, {coord_errors} errors.")
    summary_parts.append(f"Output: {output_path}")
    logger.info(" ".join(summary_parts))


if __name__ == "__main__":
    main()
