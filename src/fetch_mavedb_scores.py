# Copyright (C) 2026 Jeremy Stone
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Fetch MaveDB score sets into the pipeline's starting variant table.

Pulls scores and score set metadata for a list of score sets from the public
MaveDB REST API and writes one row per variant, ready for ``map_variants``.
The columns are those of ``cvfg_variants.0.tsv`` minus ``variant_id`` and
``mavedb_mapping_error``, which the API does not provide. MaveDB's own variant
mappings (``mavedb_mapped_hgvs_g``/``_c``/``_p``) are added only with
``--include-mavedb-mappings``. No dataset-specific handling is applied; every
score set is fetched as is.

For each score set the script calls:

1. ``GET /api/v1/score-sets/{urn}`` for the title and target gene/sequence.
2. ``GET /api/v1/score-sets/{urn}/scores`` for the per-variant scores (CSV).
   Every column other than ``accession``/``hgvs_nt``/``hgvs_splice``/``hgvs_pro``
   is packed into ``score_data`` as JSON, and ``rna_score``/``rna_score_d6``/
   ``rna_score_d20`` are pulled out of it when present.
3. With ``--include-mavedb-mappings`` only,
   ``GET /api/v1/score-sets/{urn}/variant-details`` (NDJSON) for MaveDB's mapped
   HGVS strings. Unmapped variants are absent from this export and get blank
   ``mavedb_mapped_hgvs_*`` columns. A 404 is treated as "no mapped variants".
   ``map_variants`` computes its own mappings, so these are for reference only.

Gene symbols come from a single ``GET /api/v1/score-sets/mapped-genes`` call,
falling back to the target gene's ``mappedHgncName``.

Score sets are given with ``--score-sets-file`` (a CSV/TSV with a
``score_set_urn`` column, e.g. ``data/excalibr_domainome/score_sets.tsv``)
and/or repeated ``--score-set-urn`` options. If the file has a
``preferred_transcript`` column (see ``--preferred-transcript-col``), its value
is copied to every row of that score set for ``map_variants
--preferred-transcript-col``.

Differences from a MaveDB database export (``data/cvfg/mavedb_queries.sql``):
``variant_id`` (MaveDB's internal integer ID) and ``mavedb_mapping_error`` are
omitted, since the API exposes neither, and the ``mavedb_mapped_hgvs_*``
columns are optional. Nothing downstream reads any of them.

Usage::

    python -m src.fetch_mavedb_scores output.tsv --score-sets-file score_sets.tsv [OPTIONS]
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import os
import re
import time
from pathlib import Path
from typing import Any, Iterable, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_MAVEDB_API_URL = "https://api.mavedb.org"

OUTPUT_COLS = [
    "gene_symbol",
    "score_set_title",
    "variant_urn",
    "raw_hgvs_nt",
    "raw_hgvs_pro",
    "score",
    "rna_score",
    "rna_score_d6",
    "rna_score_d20",
    "score_data",
    "reference_based",
    "target_sequence_type",
    "target_sequence",
    "preferred_transcript",
]

# Added after reference_based with --include-mavedb-mappings.
MAVEDB_MAPPING_COLS = ["mavedb_mapped_hgvs_g", "mavedb_mapped_hgvs_c", "mavedb_mapped_hgvs_p"]


def output_columns(include_mavedb_mappings: bool = False) -> list[str]:
    """Return the output columns, with the MaveDB mapping columns if requested."""
    if not include_mavedb_mappings:
        return list(OUTPUT_COLS)
    i = OUTPUT_COLS.index("reference_based") + 1
    return OUTPUT_COLS[:i] + MAVEDB_MAPPING_COLS + OUTPUT_COLS[i:]

# Columns of the /scores CSV that are not packed into score_data.
FIXED_SCORE_COLS = {"accession", "hgvs_nt", "hgvs_splice", "hgvs_pro"}

# MaveDB's CSV exports use "NA" for missing values.
_MISSING_VALUES = {"", "NA"}
_INT_RE = re.compile(r"^-?\d+$")


# ---------------------------------------------------------------------------
# Score set list
# ---------------------------------------------------------------------------


def load_score_sets(path: Path, transcript_col: str = "preferred_transcript") -> list[tuple[str, str]]:
    """Load ``(score_set_urn, preferred_transcript)`` pairs from a CSV/TSV file.

    The file must have a ``score_set_urn`` column; *transcript_col* is optional
    and yields ``""`` when absent or blank. Rows with a blank ``score_set_urn``
    are skipped, and duplicate URNs keep their first occurrence. Other columns
    (e.g. ``dataset_name``) are ignored. The delimiter is auto-detected from the
    file extension (``.tsv``/``.txt`` → tab; otherwise comma).
    """
    delim = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if reader.fieldnames is None:
            raise ValueError(f"Score sets file is empty: {path}")
        if "score_set_urn" not in reader.fieldnames:
            raise ValueError(f"Score sets file {path} is missing column: score_set_urn")

        score_sets: list[tuple[str, str]] = []
        seen: set[str] = set()
        for row in reader:
            urn = (row.get("score_set_urn") or "").strip()
            if not urn or urn in seen:
                continue
            seen.add(urn)
            score_sets.append((urn, (row.get(transcript_col) or "").strip()))
    return score_sets


# ---------------------------------------------------------------------------
# MaveDB API fetchers
# ---------------------------------------------------------------------------


def build_session(retries=8, backoff_factor=2.0):
    """Build a `requests.Session` that retries the transient 5xx errors this API returns often."""
    session = requests.Session()
    retry = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=(500, 502, 503, 504),
        allowed_methods=("GET",),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": "variant-annotation/fetch_mavedb_scores"})
    return session


def fetch_mapped_gene_symbols(
    api_url: str,
    session: requests.Session,
    timeout: int = 120,
) -> dict[str, list[str]]:
    """Return MaveDB's ``{score_set_urn: [gene_symbol, ...]}`` mapping for all score sets."""
    resp = session.get(f"{api_url}/api/v1/score-sets/mapped-genes", timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def fetch_score_set(api_url: str, urn: str, session: requests.Session, timeout: int = 120) -> dict[str, Any]:
    """Return the score set record for *urn*."""
    resp = session.get(f"{api_url}/api/v1/score-sets/{urn}", timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def fetch_scores_csv(api_url: str, urn: str, session: requests.Session, timeout: int = 120) -> str:
    """Return the score set's ``/scores`` CSV export as text."""
    resp = session.get(f"{api_url}/api/v1/score-sets/{urn}/scores", timeout=timeout)
    resp.raise_for_status()
    return resp.text


def fetch_variant_details(
    api_url: str,
    urn: str,
    session: requests.Session,
    timeout: int = 300,
) -> list[str]:
    """Return the score set's ``/variant-details`` NDJSON export as a list of lines.

    Returns an empty list when MaveDB has no mapped variants for the score set
    (HTTP 404).
    """
    resp = session.get(f"{api_url}/api/v1/score-sets/{urn}/variant-details", timeout=timeout)
    if resp.status_code == 404:
        logger.warning("%s: no mapped variants available (HTTP 404)", urn)
        return []
    resp.raise_for_status()
    return resp.text.splitlines()


# ---------------------------------------------------------------------------
# Response parsing
# ---------------------------------------------------------------------------


def score_set_metadata(score_set: dict[str, Any], mapped_gene_symbols: dict[str, list[str]]) -> dict[str, str]:
    """Extract the per-score-set output columns from a ``/score-sets/{urn}`` record.

    Only the first target gene is used; score sets with several target genes
    are logged.
    """
    urn = score_set.get("urn", "")
    target_genes = score_set.get("targetGenes") or []
    if len(target_genes) > 1:
        logger.warning("%s has %d target genes; using the first", urn, len(target_genes))
    target_gene = target_genes[0] if target_genes else {}
    target_sequence = target_gene.get("targetSequence") or {}

    symbols = mapped_gene_symbols.get(urn) or []
    gene_symbol = ",".join(symbols) or target_gene.get("mappedHgncName") or ""
    if not gene_symbol:
        logger.warning("%s: MaveDB has no mapped gene symbol; gene_symbol will be blank", urn)

    return {
        "gene_symbol": gene_symbol,
        "score_set_title": score_set.get("title") or "",
        "reference_based": "true" if target_gene.get("targetAccession") else "false",
        "target_sequence_type": target_sequence.get("sequenceType") or "",
        "target_sequence": target_sequence.get("sequence") or "",
    }


def _json_value(value: str) -> Any:
    """Convert a MaveDB CSV cell to the JSON value MaveDB stores in ``score_data``."""
    if _INT_RE.match(value):
        return int(value)
    try:
        return float(value)
    except ValueError:
        return value


def _cell(row: dict[str, Optional[str]], col: str) -> str:
    value = (row.get(col) or "").strip()
    return "" if value in _MISSING_VALUES else value


def scores_from_csv_text(csv_text: str, urn: str = "<score set>") -> list[dict[str, str]]:
    """Parse a ``/scores`` CSV export into partial output rows, one per variant.

    ``score`` and the ``rna_score*`` columns keep MaveDB's text exactly; ``NA``
    becomes blank.
    """
    reader = csv.DictReader(io.StringIO(csv_text))
    if reader.fieldnames is None or "accession" not in reader.fieldnames:
        raise ValueError(f"{urn}: /scores response has no 'accession' column")
    if "score" not in reader.fieldnames:
        raise ValueError(f"{urn}: /scores response has no 'score' column")
    score_cols = [c for c in reader.fieldnames if c not in FIXED_SCORE_COLS]

    rows: list[dict[str, str]] = []
    for row in reader:
        score_data = {c: _json_value(_cell(row, c)) for c in score_cols if _cell(row, c)}
        rows.append(
            {
                "variant_urn": _cell(row, "accession"),
                "raw_hgvs_nt": _cell(row, "hgvs_nt"),
                "raw_hgvs_pro": _cell(row, "hgvs_pro"),
                "score": _cell(row, "score"),
                "rna_score": _cell(row, "rna_score") or _cell(row, "score_rna"),
                "rna_score_d6": _cell(row, "rna_score_d6"),
                "rna_score_d20": _cell(row, "rna_score_d20"),
                "score_data": json.dumps(score_data),
            }
        )
    return rows


def mapped_hgvs_from_variant_details(lines: Iterable[str]) -> dict[str, dict[str, str]]:
    """Parse ``/variant-details`` NDJSON lines into ``{variant_urn: {mavedb_mapped_hgvs_*: ...}}``.

    Records with ``isCurrent`` false are skipped.
    """
    mapped: dict[str, dict[str, str]] = {}
    for line in lines:
        if not line.strip():
            continue
        record = json.loads(line)
        if record.get("isCurrent") is False:
            continue
        levels = record.get("mapped") or {}
        mapped[record["urn"]] = {
            "mavedb_mapped_hgvs_g": ((levels.get("genomic") or {}).get("hgvs")) or "",
            "mavedb_mapped_hgvs_c": ((levels.get("cdna") or {}).get("hgvs")) or "",
            "mavedb_mapped_hgvs_p": ((levels.get("protein") or {}).get("hgvs")) or "",
        }
    return mapped


def _variant_number(variant_urn: str) -> int:
    suffix = variant_urn.rpartition("#")[2]
    return int(suffix) if suffix.isdigit() else 0


def build_score_set_rows(
    metadata: dict[str, str],
    scores: list[dict[str, str]],
    mapped: Optional[dict[str, dict[str, str]]] = None,
    preferred_transcript: str = "",
) -> list[dict[str, str]]:
    """Combine one score set's metadata, scores, and (if *mapped* is given) MaveDB mappings into output rows."""
    columns = output_columns(include_mavedb_mappings=mapped is not None)
    rows = []
    for score_row in sorted(scores, key=lambda r: _variant_number(r["variant_urn"])):
        row = {col: "" for col in columns}
        row.update(metadata)
        row.update(score_row)
        if mapped is not None:
            row.update(mapped.get(score_row["variant_urn"], {}))
        row["preferred_transcript"] = preferred_transcript
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Fetch MaveDB score sets into the pipeline's starting variant table."
    )
    p.add_argument("output_file", help="Output CSV/TSV file")
    p.add_argument(
        "--score-sets-file",
        default=None,
        metavar="PATH",
        help=(
            "CSV/TSV file with a score_set_urn column (e.g. data/excalibr_domainome/score_sets.tsv). "
            "Other columns are ignored except --preferred-transcript-col."
        ),
    )
    p.add_argument(
        "--score-set-urn",
        action="append",
        default=[],
        metavar="URN",
        help="Score set URN to fetch. Repeatable; combined with --score-sets-file.",
    )
    p.add_argument(
        "--preferred-transcript-col",
        default="preferred_transcript",
        metavar="COLUMN",
        help=(
            "Column of --score-sets-file whose value is copied to the output's preferred_transcript "
            "column for every variant of that score set (default: preferred_transcript). "
            "Ignored if the column is absent."
        ),
    )
    p.add_argument(
        "--include-mavedb-mappings",
        action="store_true",
        help=(
            "Add mavedb_mapped_hgvs_g/_c/_p columns with MaveDB's own variant mappings "
            "(one extra, larger API request per score set). map_variants does not use them."
        ),
    )
    p.add_argument(
        "--mavedb-api-url",
        default=os.environ.get("MAVEDB_API_URL", DEFAULT_MAVEDB_API_URL),
        metavar="URL",
        help=(
            "Base URL for the MaveDB REST API. "
            "Defaults to the MAVEDB_API_URL env var, or https://api.mavedb.org."
        ),
    )
    p.add_argument(
        "--retries",
        type=int,
        default=8,
        help="HTTP retries per request on 5xx responses (default: 8)",
    )
    p.add_argument(
        "--timeout",
        type=int,
        default=300,
        metavar="SECONDS",
        help="Per-request timeout in seconds (default: 300)",
    )
    p.add_argument(
        "--sleep-between-requests",
        type=float,
        default=0.0,
        metavar="SECONDS",
        help="Seconds to sleep between score sets, to go easier on the API (default: 0)",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO)",
    )
    args = p.parse_args(argv)
    if not args.score_sets_file and not args.score_set_urn:
        p.error("give --score-sets-file and/or at least one --score-set-urn")
    return args


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    score_sets: list[tuple[str, str]] = []
    if args.score_sets_file:
        score_sets_path = Path(args.score_sets_file)
        if not score_sets_path.exists():
            logger.error("Score sets file not found: %s", score_sets_path)
            raise SystemExit(1)
        score_sets = load_score_sets(score_sets_path, args.preferred_transcript_col)
    listed = {urn for urn, _ in score_sets}
    score_sets += [(urn, "") for urn in dict.fromkeys(args.score_set_urn) if urn not in listed]
    logger.info("Fetching %d score set(s)", len(score_sets))

    api_url = args.mavedb_api_url.rstrip("/")
    session = build_session(retries=args.retries)
    mapped_gene_symbols = fetch_mapped_gene_symbols(api_url, session, timeout=args.timeout)

    rows: list[dict[str, str]] = []
    for i, (urn, preferred_transcript) in enumerate(score_sets, start=1):
        metadata = score_set_metadata(fetch_score_set(api_url, urn, session, args.timeout), mapped_gene_symbols)
        scores = scores_from_csv_text(fetch_scores_csv(api_url, urn, session, args.timeout), urn=urn)
        mapped = None
        if args.include_mavedb_mappings:
            mapped = mapped_hgvs_from_variant_details(fetch_variant_details(api_url, urn, session, args.timeout))
        rows += build_score_set_rows(metadata, scores, mapped, preferred_transcript)
        logger.info("[%d/%d] %s (%s): %d variant(s)", i, len(score_sets), urn, metadata["gene_symbol"], len(scores))
        if mapped is not None:
            logger.info("  %d variant(s) with a MaveDB mapping", sum(1 for r in scores if r["variant_urn"] in mapped))
        if args.sleep_between_requests and i < len(score_sets):
            time.sleep(args.sleep_between_requests)

    # Same order as the database export: by gene, then by score set, then by variant.
    order = {urn: i for i, (urn, _) in enumerate(score_sets)}
    rows.sort(key=lambda r: (r["gene_symbol"], order[r["variant_urn"].rpartition("#")[0]]))

    output_path = Path(args.output_file)
    out_delim = "\t" if output_path.suffix.lower() in (".tsv", ".txt") else ","
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as outf:
        writer = csv.DictWriter(
            outf, fieldnames=output_columns(args.include_mavedb_mappings), delimiter=out_delim
        )
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d row(s) to %s", len(rows), output_path)


if __name__ == "__main__":
    main()
