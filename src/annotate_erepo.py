"""Annotate variant rows with ClinGen Expert Panel classifications from the Evidence Repository.

The ClinGen Evidence Repository (erepo) publishes expert-panel variant classifications.
This script downloads (and caches) the full TSV from::

    https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download

and then for each input row joins the classification data using up to three keys:

1. **HGVS** – each candidate in ``mapped_hgvs_c`` is looked up against every
   expression in the ``HGVS Expressions`` column of the erepo TSV (comma-separated).
   The ``HGVS Expressions`` column contains pre-normalised HGVS strings across
   multiple representations (transcript, genomic, protein) and does not include
   extraneous gene-symbol annotations, making it the most reliable join surface.
   Gene symbols in parentheses that may appear in ``mapped_hgvs_c`` values
   (e.g. ``NM_004333.5(BRAF):c.740T>C``) are stripped before comparison.

2. **ClinVar Variation ID** – matched against ``ClinVar Variation Id``.

3. **ClinGen Allele Registry ID (CAID)** – matched against ``Allele Registry Id``.

Which keys are used is controlled by ``--join-keys`` (default: all three).

Output columns (prefixed ``clingen_evidence_repository.``) carry the pipe-separated
values across each ``mapped_hgvs_c`` candidate:

    clingen_evidence_repository.ClinVar Variation Id
    clingen_evidence_repository.Allele Registry Id
    clingen_evidence_repository.Disease Mondo Id
    clingen_evidence_repository.Mode of Inheritance
    clingen_evidence_repository.Assertion
    clingen_evidence_repository.Applied Evidence Codes (Met)
    clingen_evidence_repository.Applied Evidence Codes (Not Met)
    clingen_evidence_repository.Summary of interpretation
    clingen_evidence_repository.PubMed Articles
    clingen_evidence_repository.Expert Panel
    clingen_evidence_repository.Guideline
    clingen_evidence_repository.Approval Date
    clingen_evidence_repository.Published Date
    clingen_evidence_repository.Retracted
    clingen_evidence_repository.Evidence Repo Link
    clingen_evidence_repository.Uuid

Additionally, ``clingen_evidence_repository.warnings`` records cross-key
discrepancies (e.g. variant found by HGVS but not CAID).

Usage::

    python -m src.annotate_erepo input.tsv output.tsv [OPTIONS]
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import sys
import time
from itertools import islice
from pathlib import Path
from typing import Optional

import requests

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EREPO_DOWNLOAD_URL = (
    "https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download"
)
EREPO_CACHE_FILENAME = "clingen_erepo_classifications.tsv"

# Column names exactly as they appear in the erepo TSV header.
EREPO_COL_VARIATION = "Variation"
EREPO_COL_HGVS_EXPRESSIONS = "HGVS Expressions"
EREPO_COL_CLINVAR_VARIATION_ID = "ClinVar Variation Id"
EREPO_COL_ALLELE_REGISTRY_ID = "Allele Registry Id"

# Ordered list of erepo columns that are forwarded to the output.
EREPO_OUTPUT_COLS = [
    "ClinVar Variation Id",
    "Allele Registry Id",
    "Disease Mondo Id",
    "Mode of Inheritance",
    "Assertion",
    "Applied Evidence Codes (Met)",
    "Applied Evidence Codes (Not Met)",
    "Summary of interpretation",
    "PubMed Articles",
    "Expert Panel",
    "Guideline",
    "Approval Date",
    "Published Date",
    "Retracted",
    "Evidence Repo Link",
    "Uuid",
]

OUTPUT_COL_PREFIX = "clingen_evidence_repository"
OUTPUT_WARNINGS_COL = f"{OUTPUT_COL_PREFIX}.warnings"

# Regex to strip a gene symbol in parentheses from the start of an HGVS string,
# e.g. "NM_004333.5(BRAF):c.740T>C" → "NM_004333.5:c.740T>C"
_GENE_SYMBOL_RE = re.compile(r"^([A-Za-z0-9_.]+)\([^)]+\)(:.*)$")

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Erepo TSV download / cache
# ---------------------------------------------------------------------------


def fetch_erepo_tsv(cache_dir: Path, refresh: bool = False) -> Path:
    """Download (or return cached) the ClinGen erepo classifications TSV.

    The file is cached in *cache_dir* and reused on subsequent calls unless
    *refresh* is True.

    Returns:
        Path to the cached TSV file.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / EREPO_CACHE_FILENAME

    if dest.exists() and not refresh:
        logger.info("Using cached erepo TSV: %s", dest)
        return dest

    logger.info("Downloading ClinGen erepo TSV from %s", EREPO_DOWNLOAD_URL)
    response = requests.get(EREPO_DOWNLOAD_URL, stream=True, timeout=120)
    response.raise_for_status()
    tmp = dest.with_suffix(".tmp")
    with tmp.open("wb") as fh:
        for chunk in response.iter_content(chunk_size=1 << 16):
            fh.write(chunk)
    tmp.rename(dest)
    logger.info("Saved erepo TSV to %s", dest)
    return dest


# ---------------------------------------------------------------------------
# Index building
# ---------------------------------------------------------------------------


def _strip_gene_symbol(hgvs: str) -> str:
    """Strip a parenthesised gene symbol from an input HGVS string before lookup.

    ``mapped_hgvs_c`` values may carry a gene symbol annotation, e.g.
    ``NM_004333.5(BRAF):c.740T>C``, while the erepo ``HGVS Expressions``
    column does not.  This normalises the input side so the two can be
    compared directly.
    """
    m = _GENE_SYMBOL_RE.match(hgvs.strip())
    if m:
        return m.group(1) + m.group(2)
    return hgvs.strip()


def _hgvs_candidates(raw: str) -> list[str]:
    """Return all HGVS candidates from a comma-separated expression field."""
    return [h.strip() for h in raw.split(",") if h.strip()]


ErepoRecord = dict[str, str]


class ErepoIndex:
    """In-memory index of erepo records for fast lookup by three key types."""

    def __init__(self) -> None:
        # Maps each key to a *list* of matching records (one variant can appear
        # in multiple expert-panel classifications with different assertions).
        self.by_hgvs: dict[str, list[ErepoRecord]] = {}
        self.by_clinvar_id: dict[str, list[ErepoRecord]] = {}
        self.by_caid: dict[str, list[ErepoRecord]] = {}

    @classmethod
    def build(cls, tsv_path: Path) -> "ErepoIndex":
        """Parse the erepo TSV and build the index."""
        index = cls()
        logger.info("Loading erepo TSV from %s …", tsv_path)
        row_count = 0
        with tsv_path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                row_count += 1
                record: ErepoRecord = {col: (row.get(col) or "").strip() for col in EREPO_OUTPUT_COLS}

                # Index by every HGVS expression in the HGVS Expressions column.
                # These are already clean (no gene-symbol parentheticals) and
                # cover multiple representations (transcript, genomic, protein).
                for hgvs in _hgvs_candidates(row.get(EREPO_COL_HGVS_EXPRESSIONS) or ""):
                    if ":" not in hgvs:
                        continue  # skip non-HGVS tokens
                    index.by_hgvs.setdefault(hgvs, []).append(record)

                # Index by ClinVar Variation ID
                clinvar_id = (row.get(EREPO_COL_CLINVAR_VARIATION_ID) or "").strip()
                if clinvar_id:
                    index.by_clinvar_id.setdefault(clinvar_id, []).append(record)

                # Index by Allele Registry ID (CAID)
                caid = (row.get(EREPO_COL_ALLELE_REGISTRY_ID) or "").strip()
                if caid:
                    index.by_caid.setdefault(caid, []).append(record)

        logger.info(
            "Erepo index built: %d rows; %d HGVS keys, %d ClinVar Variation ID keys, %d CAID keys",
            row_count,
            len(index.by_hgvs),
            len(index.by_clinvar_id),
            len(index.by_caid),
        )
        return index


# ---------------------------------------------------------------------------
# Per-candidate lookup
# ---------------------------------------------------------------------------


def _records_for_candidate(
    hgvs_c: str,
    clinvar_variation_id: str,
    caid: str,
    index: ErepoIndex,
    join_keys: set[str],
) -> tuple[list[ErepoRecord], str]:
    """Look up erepo records for a single candidate.

    Returns ``(records, warnings)`` where *records* is the deduplicated list of
    matching erepo records (empty if nothing found) and *warnings* is a
    human-readable string describing any cross-key discrepancies.
    """
    hgvs_c_norm = _strip_gene_symbol(hgvs_c) if hgvs_c else ""

    found_by: dict[str, list[ErepoRecord]] = {}

    if "hgvs" in join_keys and hgvs_c_norm:
        recs = index.by_hgvs.get(hgvs_c_norm)
        if recs:
            found_by["hgvs"] = recs

    if "clinvar" in join_keys and clinvar_variation_id:
        recs = index.by_clinvar_id.get(clinvar_variation_id.strip())
        if recs:
            found_by["clinvar"] = recs

    if "caid" in join_keys and caid:
        recs = index.by_caid.get(caid.strip())
        if recs:
            found_by["caid"] = recs

    if not found_by:
        return [], ""

    # Deduplicate by Uuid across all found record lists.
    seen_uuids: set[str] = set()
    merged: list[ErepoRecord] = []
    for recs in found_by.values():
        for rec in recs:
            uid = rec.get("Uuid", "")
            if uid not in seen_uuids:
                seen_uuids.add(uid)
                merged.append(rec)

    # Build warnings when different keys yield different (non-overlapping) record sets.
    warnings: list[str] = []
    if len(found_by) > 1:
        uuid_sets = {k: {r.get("Uuid", "") for r in recs} for k, recs in found_by.items()}
        key_list = list(uuid_sets)
        for i, k1 in enumerate(key_list):
            for k2 in key_list[i + 1 :]:
                if not uuid_sets[k1] & uuid_sets[k2]:
                    warnings.append(
                        f"non-overlapping results for keys {k1!r} and {k2!r}"
                    )
    elif len(join_keys) > 1:
        # Only one key yielded results; note which keys were present but missed.
        keys_present = set()
        if hgvs_c_norm:
            keys_present.add("hgvs")
        if clinvar_variation_id:
            keys_present.add("clinvar")
        if caid:
            keys_present.add("caid")
        missed = (join_keys & keys_present) - set(found_by)
        if missed:
            found_key = next(iter(found_by))
            warnings.append(
                f"found via {found_key!r} but not {', '.join(sorted(missed))!r}"
            )

    return merged, "; ".join(warnings)


# ---------------------------------------------------------------------------
# Row annotation
# ---------------------------------------------------------------------------

_EMPTY_RECORD: ErepoRecord = {col: "" for col in EREPO_OUTPUT_COLS}


def _join_records(records: list[ErepoRecord]) -> ErepoRecord:
    """Merge multiple erepo records for the same candidate into a single dict.

    When there are multiple records (e.g. the same variant classified by two
    different expert panels), values are joined with " | " within each field.
    """
    if not records:
        return _EMPTY_RECORD
    if len(records) == 1:
        return records[0]
    merged: ErepoRecord = {}
    for col in EREPO_OUTPUT_COLS:
        vals = [r.get(col, "") for r in records]
        # Deduplicate while preserving order
        seen: set[str] = set()
        deduped: list[str] = []
        for v in vals:
            if v not in seen:
                seen.add(v)
                deduped.append(v)
        merged[col] = " | ".join(deduped)
    return merged


def annotate_row(
    row: dict[str, str],
    index: ErepoIndex,
    join_keys: set[str],
    mapped_hgvs_c_col: str,
    clinvar_variation_id_col: str,
    caid_col: str,
) -> dict[str, str]:
    """Return annotation values for *row* as a flat dict keyed by output column name."""
    empty_out = {f"{OUTPUT_COL_PREFIX}.{col}": "" for col in EREPO_OUTPUT_COLS}
    empty_out[OUTPUT_WARNINGS_COL] = ""

    raw_hgvs_c = (row.get(mapped_hgvs_c_col) or "").strip()
    hgvs_candidates = [c.strip() for c in raw_hgvs_c.split("|")] if raw_hgvs_c else [""]
    n = len(hgvs_candidates)

    # Split pipe-delimited ClinVar variation IDs and CAIDs to align with hgvs_candidates.
    raw_clinvar = (row.get(clinvar_variation_id_col) or "").strip()
    clinvar_ids = [c.strip() for c in raw_clinvar.split("|")] if raw_clinvar else []
    raw_caid = (row.get(caid_col) or "").strip()
    caids = [c.strip() for c in raw_caid.split("|")] if raw_caid else []

    per_candidate_records: list[ErepoRecord] = []
    per_candidate_warnings: list[str] = []

    for i in range(n):
        hgvs_c = hgvs_candidates[i]
        cv_id = clinvar_ids[i] if i < len(clinvar_ids) else ""
        caid = caids[i] if i < len(caids) else ""
        records, warn = _records_for_candidate(hgvs_c, cv_id, caid, index, join_keys)
        per_candidate_records.append(_join_records(records))
        per_candidate_warnings.append(warn)

    out: dict[str, str] = {}
    for col in EREPO_OUTPUT_COLS:
        out[f"{OUTPUT_COL_PREFIX}.{col}"] = "|".join(r.get(col, "") for r in per_candidate_records)
    out[OUTPUT_WARNINGS_COL] = "|".join(per_candidate_warnings)
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Annotate variant rows with ClinGen Expert Panel classifications "
            "from the ClinGen Evidence Repository."
        )
    )
    p.add_argument("input_file", help="Input TSV file")
    p.add_argument("output_file", help="Output TSV file with erepo annotation columns appended")

    p.add_argument(
        "--cache-dir",
        default=os.environ.get("EREPO_CACHE_DIR", "/tmp/erepo_cache"),
        metavar="DIR",
        help=(
            "Directory for caching the downloaded erepo TSV "
            "(default: $EREPO_CACHE_DIR or /tmp/erepo_cache)."
        ),
    )
    p.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Re-download the erepo TSV even if a cached copy exists.",
    )
    p.add_argument(
        "--join-keys",
        default="hgvs,clinvar,caid",
        metavar="KEY[,KEY...]",
        help=(
            "Comma-separated list of keys to use when joining to the erepo. "
            "Allowed values: hgvs, clinvar, caid (default: all three). "
            "Example: --join-keys hgvs,caid"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-c-col",
        default="mapped_hgvs_c",
        metavar="COL",
        help="Column containing pipe-delimited transcript HGVS strings (default: mapped_hgvs_c).",
    )
    p.add_argument(
        "--clinvar-variation-id-col",
        default="clinvar.202601.variation_id",
        metavar="COL",
        help=(
            "Column containing ClinVar Variation IDs for CAID-based join "
            "(default: clinvar.202601.variation_id)."
        ),
    )
    p.add_argument(
        "--caid-col",
        default="dna_clingen_allele_id",
        metavar="COL",
        help="Column containing ClinGen Allele Registry IDs (default: dna_clingen_allele_id).",
    )
    p.add_argument("--delimiter", default="\t", help="Input/output field delimiter (default: TAB).")
    p.add_argument(
        "--skip",
        type=int,
        default=0,
        help="Number of data rows to skip before annotation.",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of data rows to annotate.",
    )
    p.add_argument(
        "--csv-field-size-limit",
        type=int,
        default=csv.field_size_limit(),
        metavar="BYTES",
        help="Maximum per-field character length for TSV parsing (default: %(default)s).",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    csv.field_size_limit(args.csv_field_size_limit)

    if args.skip < 0:
        logger.error("--skip must be >= 0, got: %d", args.skip)
        sys.exit(1)
    if args.limit is not None and args.limit < 1:
        logger.error("--limit must be >= 1 when provided, got: %d", args.limit)
        sys.exit(1)

    # Parse and validate join keys.
    allowed_keys = {"hgvs", "clinvar", "caid"}
    raw_keys = [k.strip().lower() for k in args.join_keys.split(",") if k.strip()]
    invalid = set(raw_keys) - allowed_keys
    if invalid:
        logger.error(
            "--join-keys contains invalid values: %s. Allowed: %s",
            ", ".join(sorted(invalid)),
            ", ".join(sorted(allowed_keys)),
        )
        sys.exit(1)
    if not raw_keys:
        logger.error("--join-keys must contain at least one key.")
        sys.exit(1)
    join_keys = set(raw_keys)
    logger.info("Join keys: %s", ", ".join(sorted(join_keys)))

    delim = "\t" if args.delimiter == "\\t" else args.delimiter
    cache_dir = Path(args.cache_dir)

    tsv_path = fetch_erepo_tsv(cache_dir, refresh=args.refresh_cache)
    index = ErepoIndex.build(tsv_path)

    in_path = Path(args.input_file)
    out_path = Path(args.output_file)

    ann_cols = [f"{OUTPUT_COL_PREFIX}.{col}" for col in EREPO_OUTPUT_COLS] + [OUTPUT_WARNINGS_COL]

    with in_path.open("r", encoding="utf-8", newline="") as in_fh, \
         out_path.open("w", encoding="utf-8", newline="") as out_fh:

        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears empty: %s", in_path)
            sys.exit(1)

        fieldnames = list(reader.fieldnames)
        out_fieldnames = fieldnames + [c for c in ann_cols if c not in fieldnames]

        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=delim,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        processed = 0
        annotated = 0

        selected = islice(
            reader,
            args.skip,
            None if args.limit is None else args.skip + args.limit,
        )

        for row in selected:
            ann = annotate_row(
                row,
                index,
                join_keys,
                mapped_hgvs_c_col=args.mapped_hgvs_c_col,
                clinvar_variation_id_col=args.clinvar_variation_id_col,
                caid_col=args.caid_col,
            )
            row.update(ann)
            writer.writerow(row)
            processed += 1
            # Count a row as annotated if any candidate had a non-empty Assertion.
            if ann[f"{OUTPUT_COL_PREFIX}.Assertion"].replace("|", "").strip():
                annotated += 1
            if processed % 10000 == 0:
                out_fh.flush()
                logger.info("Processed %d rows (%d annotated) …", processed, annotated)

    logger.info(
        "Done. %d rows processed, %d annotated with erepo data → %s",
        processed,
        annotated,
        out_path,
    )


if __name__ == "__main__":
    main()
