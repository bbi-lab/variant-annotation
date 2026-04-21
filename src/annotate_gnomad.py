"""Annotate variants with gnomAD allele frequency metrics using a local Hail table cache.

This script annotates each input row using DNA-level ClinGen allele IDs from
``dna_clingen_allele_id`` (or a custom column). It resolves one row by trying
pipe-delimited candidate CAIDs in order and using the first gnomAD hit.

On first execution, it downloads/reads the source gnomAD Hail table and writes a
local indexed cache keyed by ``caid``. Subsequent runs reuse the local cache.

Default output columns:
  - <namespace>.<version>.minor_allele_frequency
  - <namespace>.<version>.allele_frequency
  - <namespace>.<version>.allele_count
  - <namespace>.<version>.allele_number
  - <namespace>.<version>.faf95_max
  - <namespace>.<version>.faf95_max_ancestry

Usage:
    python -m src.annotate_gnomad input.tsv output.tsv [OPTIONS]
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

GNOMAD_DATA_VERSION = os.environ.get("GNOMAD_DATA_VERSION", "v4.1")
GNOMAD_HT_URI_DEFAULT = os.environ.get(
    "GNOMAD_HT_URI",
    "gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht",
)


@dataclass
class GnomadRecord:
    caid: str
    allele_count: int
    allele_number: int
    allele_frequency: float
    minor_allele_frequency: float
    faf95_max: Optional[float]
    faf95_max_ancestry: str


def _import_hail():
    try:
        import hail as hl  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "Hail is required for gnomAD annotation. Install the gnomad extra and ensure Java is available."
        ) from exc
    return hl


def _detect_separator(file_path: str) -> str:
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _split_pipe(value: str) -> list[str]:
    if "|" not in value:
        return [value.strip()] if value.strip() else []
    return [part.strip() for part in value.split("|") if part.strip()]


def _local_ht_path(cache_dir: Path, version: str) -> Path:
    return cache_dir / f"gnomad_{version.replace('.', '_')}_indexed.ht"


def _has_path(dtype: Any, parts: list[str]) -> bool:
    hl = _import_hail()
    cur = dtype
    for part in parts:
        if not isinstance(cur, hl.tstruct) or part not in cur:
            return False
        cur = cur[part]
    return True


def _get_path(expr: Any, parts: list[str]) -> Any:
    cur = expr
    for part in parts:
        cur = cur[part]
    return cur


def _choose_expr(ht: Any, candidates: list[list[str]]) -> Optional[Any]:
    for path in candidates:
        if _has_path(ht.row.dtype, path):
            return _get_path(ht, path)
    return None


def ensure_local_gnomad_ht(
    cache_dir: Path,
    *,
    version: str,
    source_ht_uri: str,
    overwrite: bool = False,
) -> Path:
    """Create a local CAID-keyed gnomAD Hail table cache if missing."""
    hl = _import_hail()

    cache_dir.mkdir(parents=True, exist_ok=True)
    ht_path = _local_ht_path(cache_dir, version)
    if ht_path.exists() and not overwrite:
        logger.info("Using cached gnomAD Hail table: %s", ht_path)
        return ht_path

    hail_tmp = cache_dir / "hail-tmp"
    hail_tmp.mkdir(parents=True, exist_ok=True)

    logger.info("Initializing Hail and loading source table: %s", source_ht_uri)
    hl.init(tmp_dir=str(hail_tmp), quiet=True, idempotent=True)
    try:
        source_ht = hl.read_table(source_ht_uri)

        caid_expr = _choose_expr(source_ht, [["caid"]])
        if caid_expr is None:
            raise RuntimeError("Could not find 'caid' field in source gnomAD Hail table")

        ac_expr = _choose_expr(source_ht, [["joint", "freq", "all", "ac"], ["freq", "all", "ac"]])
        an_expr = _choose_expr(source_ht, [["joint", "freq", "all", "an"], ["freq", "all", "an"]])
        if ac_expr is None or an_expr is None:
            raise RuntimeError("Could not find allele count/number fields in source gnomAD Hail table")

        faf_anc_expr = _choose_expr(
            source_ht,
            [["joint", "fafmax", "faf95_max_gen_anc"], ["fafmax", "faf95_max_gen_anc"]],
        )
        faf_max_expr = _choose_expr(
            source_ht,
            [["joint", "fafmax", "faf95_max"], ["fafmax", "faf95_max"]],
        )

        prepared = source_ht.select(
            caid=hl.str(caid_expr),
            allele_count=hl.int64(ac_expr),
            allele_number=hl.int64(an_expr),
            faf95_max_ancestry=hl.if_else(
                hl.is_defined(faf_anc_expr), hl.str(faf_anc_expr), hl.str("")
            )
            if faf_anc_expr is not None
            else hl.str(""),
            faf95_max=hl.if_else(hl.is_defined(faf_max_expr), hl.float64(faf_max_expr), hl.missing(hl.tfloat64))
            if faf_max_expr is not None
            else hl.missing(hl.tfloat64),
        )

        prepared = prepared.key_by(prepared.caid)
        logger.info("Writing local gnomAD cache table: %s", ht_path)
        prepared.write(str(ht_path), overwrite=True)
    finally:
        hl.stop()

    return ht_path


def load_gnomad_records_for_caids(local_ht_path: Path, caids: set[str], cache_dir: Path) -> dict[str, GnomadRecord]:
    """Load gnomAD records for requested CAIDs from local Hail table."""
    if not caids:
        return {}

    hl = _import_hail()
    hail_tmp = cache_dir / "hail-tmp"
    hail_tmp.mkdir(parents=True, exist_ok=True)

    hl.init(tmp_dir=str(hail_tmp), quiet=True, idempotent=True)
    try:
        ht = hl.read_table(str(local_ht_path))
        caid_literal = hl.literal(caids)
        filtered = ht.filter(caid_literal.contains(ht.caid))
        rows = filtered.collect()
    finally:
        hl.stop()

    out: dict[str, GnomadRecord] = {}
    for row in rows:
        caid = str(row.caid)
        ac = int(row.allele_count)
        an = int(row.allele_number)
        if an <= 0:
            continue
        af = float(ac) / float(an)
        maf = min(af, 1.0 - af)
        faf95_max = float(row.faf95_max) if row.faf95_max is not None else None
        faf95_max_ancestry = str(row.faf95_max_ancestry or "")
        out[caid] = GnomadRecord(
            caid=caid,
            allele_count=ac,
            allele_number=an,
            allele_frequency=af,
            minor_allele_frequency=maf,
            faf95_max=faf95_max,
            faf95_max_ancestry=faf95_max_ancestry,
        )
    return out


def annotate_row(row: dict[str, str], records: dict[str, GnomadRecord], col_prefix: str, dna_col: str) -> dict[str, str]:
    out = {
        f"{col_prefix}.minor_allele_frequency": "",
        f"{col_prefix}.allele_frequency": "",
        f"{col_prefix}.allele_count": "",
        f"{col_prefix}.allele_number": "",
        f"{col_prefix}.faf95_max": "",
        f"{col_prefix}.faf95_max_ancestry": "",
    }

    caids = _split_pipe((row.get(dna_col) or "").strip())
    for caid in caids:
        rec = records.get(caid)
        if rec is None:
            continue
        out[f"{col_prefix}.minor_allele_frequency"] = str(rec.minor_allele_frequency)
        out[f"{col_prefix}.allele_frequency"] = str(rec.allele_frequency)
        out[f"{col_prefix}.allele_count"] = str(rec.allele_count)
        out[f"{col_prefix}.allele_number"] = str(rec.allele_number)
        out[f"{col_prefix}.faf95_max"] = "" if rec.faf95_max is None else str(rec.faf95_max)
        out[f"{col_prefix}.faf95_max_ancestry"] = rec.faf95_max_ancestry
        return out

    return out


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Annotate rows with gnomAD minor allele frequency using DNA-level ClinGen allele IDs "
            "and a local Hail table cache."
        )
    )
    p.add_argument("input_file", help="Input CSV/TSV file")
    p.add_argument("output_file", help="Output CSV/TSV file")
    p.add_argument("--gnomad-version", default=GNOMAD_DATA_VERSION, help="gnomAD version label for output columns")
    p.add_argument(
        "--gnomad-namespace",
        default="gnomad",
        help="Namespace for output columns (default: gnomad)",
    )
    p.add_argument(
        "--gnomad-ht-uri",
        default=GNOMAD_HT_URI_DEFAULT,
        help="Source gnomAD Hail table URI (gs://...)",
    )
    p.add_argument(
        "--cache-dir",
        default=os.environ.get("GNOMAD_CACHE_DIR", "/tmp/gnomad_cache"),
        help="Cache dir for local gnomAD Hail table",
    )
    p.add_argument(
        "--dna-clingen-allele-id-col",
        default="dna_clingen_allele_id",
        help="Column containing DNA-level ClinGen allele IDs",
    )
    p.add_argument("--delimiter", default="\t", help="Input/output delimiter (default TAB)")
    p.add_argument("--download-only", action="store_true", help="Only materialize local gnomAD cache; do not annotate")
    p.add_argument("--refresh-cache", action="store_true", help="Rebuild local gnomAD cache even if present")
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    if not args.gnomad_ht_uri:
        logger.error("--gnomad-ht-uri (or GNOMAD_HT_URI env var) is required")
        sys.exit(1)

    cache_dir = Path(args.cache_dir)
    local_ht = ensure_local_gnomad_ht(
        cache_dir,
        version=args.gnomad_version,
        source_ht_uri=args.gnomad_ht_uri,
        overwrite=args.refresh_cache,
    )

    if args.download_only:
        logger.info("Local gnomAD cache ready at %s", local_ht)
        return

    delim = "\t" if args.delimiter == "\\t" else args.delimiter
    input_path = Path(args.input_file)
    output_path = Path(args.output_file)

    with input_path.open("r", encoding="utf-8", newline="") as in_fh:
        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears empty: %s", input_path)
            sys.exit(1)
        rows = list(reader)
        fieldnames = list(reader.fieldnames)

    caids: set[str] = set()
    for row in rows:
        caids.update(_split_pipe((row.get(args.dna_clingen_allele_id_col) or "").strip()))

    logger.info("Loading gnomAD records for %d unique CAIDs", len(caids))
    records = load_gnomad_records_for_caids(local_ht, caids, cache_dir)

    prefix = f"{args.gnomad_namespace}.{args.gnomad_version}"
    ann_cols = [
        f"{prefix}.minor_allele_frequency",
        f"{prefix}.allele_frequency",
        f"{prefix}.allele_count",
        f"{prefix}.allele_number",
        f"{prefix}.faf95_max",
        f"{prefix}.faf95_max_ancestry",
    ]

    out_fieldnames = fieldnames + [c for c in ann_cols if c not in fieldnames]
    annotated = 0
    with output_path.open("w", encoding="utf-8", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_fieldnames, delimiter=delim, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            ann = annotate_row(row, records, prefix, args.dna_clingen_allele_id_col)
            row.update(ann)
            writer.writerow(row)
            if ann[f"{prefix}.minor_allele_frequency"]:
                annotated += 1

    logger.info("Done. %d/%d rows annotated -> %s", annotated, len(rows), output_path)


if __name__ == "__main__":
    main()
