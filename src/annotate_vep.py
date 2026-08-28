"""Annotate variants with VEP molecular consequences.

Composition root for :mod:`variant_annotation.lib.vep`. This script owns CSV streaming, the Redis
cache, and CLI wiring; every decision about *which* consequence applies to a variant lives in the
library, so this pipeline and mavedb-api cannot drift apart.

Each input row is annotated from one or more pipe-delimited HGVS candidates. By default the columns
``mapped_hgvs_c``, ``mapped_hgvs_g``, and ``mapped_hgvs_p`` are tried in that order; the first
non-blank value in a row is used (see ``--hgvs-cols``).

For rows with multiple DNA candidates (e.g. reverse translations of a protein variant), each candidate
is annotated independently. The output columns ``vep.mutational_consequences``,
``vep.most_severe_mutational_consequence``, ``vep.consequence_source``, ``vep.access_date``, and
``vep.error`` are pipe-delimited lists aligned to the input candidate positions.

``vep.mutational_consequences`` contains a ``^``-delimited list of all consequence terms for the
matched transcript when ``source == "transcript"``, or just the single term when the top-level
``most_severe_consequence`` was used. ``vep.most_severe_mutational_consequence`` always contains the
single most-severe term, or is empty when no consequence was determined.

Transcript selection
--------------------
Ensembl VEP annotates every transcript that overlaps a variant's genomic position, not just the
transcript referenced in the input HGVS string. Its top-level ``most_severe_consequence`` is therefore
the worst outcome across **all** overlapping transcripts, which regularly differs from the consequence
on the transcript of interest.

The library resolves against the input's own transcript where one can be identified:

- **RefSeq transcripts** (``NM_``/``NR_``/``XM_``/``XR_``): the request carries ``refseq=1`` so
  ``transcript_consequences`` comes back with RefSeq ids, and the entry matching the input accession
  (versions stripped) is selected. ``source = "transcript"``.
- **Ensembl transcripts** (``ENST``): the default transcript set is used and the matching ``ENST``
  accession is looked up directly. ``source = "transcript"``.
- **Genomic HGVS** (``NC_``) and **protein HGVS** (``NP_``): no transcript is identifiable from the
  string, so VEP's cross-transcript headline is used. ``source = "most_severe"``.
- **Reference-identical input** (``c.123=``, ``p.Met4=``, or a ``delins`` restating the reference): no
  VEP call is made; the consequence is ``no_change`` with ``source = "reference_identical"``. Reads the
  transcript reference from UTA when ``UTA_DB_URL`` is set.

Column priority and transcript specificity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Transcript HGVS strings (``mapped_hgvs_c``) yield transcript-specific consequences; genomic strings
(``mapped_hgvs_g``) do not. ``mapped_hgvs_c`` should therefore come first in ``--hgvs-cols``, as the
default ordering ``mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p`` reflects.

Redis caching
-------------
When a Redis instance is reachable, each resolved consequence is cached to avoid redundant Ensembl API
calls across runs. Failed requests are never cached, so they are retried on the next run. Caching is
transparent and fails gracefully: if Redis is unreachable the script continues without it.

Relevant environment variables:

``VEP_CACHE_ENABLED``
    Set to ``0`` / ``false`` to disable caching entirely (default: enabled).
``VEP_CACHE_REDIS_URL``
    Redis connection URL (default: ``redis://redis:6379/0``; also falls back to the generic
    ``REDIS_URL`` variable).
``VEP_CACHE_PREFIX``
    Key namespace prefix (default: ``vep:v3``). The version suffix is bumped whenever the cached value
    shape or the library's ``RESOLVER_VERSION`` changes, so stale answers are ignored rather than served.
``VEP_CACHE_TTL_SECONDS``
    TTL for cached hits, in seconds (default: 86 400 — 1 day).
``VEP_CACHE_MISS_TTL_SECONDS``
    TTL for cached misses, in seconds (default: 86 400 — 1 day).
``UTA_DB_URL``
    UTA connection URL, used only for reference-identical detection. Absent means no-change rows are
    reported as having no consequence rather than as ``no_change``.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
from datetime import date
from itertools import islice
from pathlib import Path
from typing import Any, Optional, TextIO

from variant_annotation.lib.clients.ensembl import (
    ENSEMBL_API_URL_DEFAULT,
    VEP_MAX_HGVS_PER_REQUEST,
    EnsemblRestClient,
)
from variant_annotation.lib.clients.uta import UtaClient, connect_uta
from variant_annotation.lib.vep import (
    RESOLVER_VERSION,
    ConsequenceOutcome,
    ConsequenceResolution,
    ConsequenceSource,
    ReferenceSequence,
    VepConfig,
    VepInput,
    requested_transcript,
    resolve_consequences,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Redis cache configuration
# ---------------------------------------------------------------------------

VEP_CACHE_REDIS_URL_DEFAULT = "redis://redis:6379/0"
# v3: values now carry the library's resolution outcome and source enum rather than the old
# (most_severe, all_consequences, source) tuple. Entries written by earlier versions are unreadable
# under this prefix, so they are simply never consulted.
VEP_CACHE_PREFIX_DEFAULT = "vep:v3"
VEP_CACHE_TTL_SECONDS_DEFAULT = 86400  # 1 day
VEP_CACHE_MISS_TTL_SECONDS_DEFAULT = 86400  # 1 day

_VEP_MISS_SENTINEL = "__MISS__"
_VEP_REDIS_CLIENT: Any = None
_VEP_REDIS_INIT_ATTEMPTED = False
_VEP_REDIS_UNAVAILABLE_LOGGED = False


def _split_pipe_preserve_positions(value: str) -> list[str]:
    raw = value or ""
    if "|" not in raw:
        return [raw.strip()]
    return [part.strip() for part in raw.split("|")]


def _sanitize_for_tsv(msg: str, max_length: int = 300) -> str:
    """Replace TSV/pipe-delimited field-breaking characters and truncate.

    Sink-specific escaping: the library bounds an error's length but knows nothing about this file's
    delimiters, so neutralising them is this layer's job.
    """
    sanitized = msg.replace("|", ";").replace("\t", " ").replace("\r", "").replace("\n", " ").strip()
    return sanitized[:max_length] + ("..." if len(sanitized) > max_length else "")


# ---------------------------------------------------------------------------
# Redis cache helpers
# ---------------------------------------------------------------------------


def _vep_env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off", ""}


def _vep_env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        logger.warning("Ignoring non-integer %s=%r; using %d", name, raw, default)
        return default


def _vep_cache_enabled() -> bool:
    return _vep_env_bool("VEP_CACHE_ENABLED", True)


def _vep_cache_prefix() -> str:
    return (os.environ.get("VEP_CACHE_PREFIX") or VEP_CACHE_PREFIX_DEFAULT).strip()


def _vep_cache_ttl() -> int:
    return max(1, _vep_env_int("VEP_CACHE_TTL_SECONDS", VEP_CACHE_TTL_SECONDS_DEFAULT))


def _vep_cache_miss_ttl() -> int:
    return max(1, _vep_env_int("VEP_CACHE_MISS_TTL_SECONDS", VEP_CACHE_MISS_TTL_SECONDS_DEFAULT))


def _vep_cache_redis_url() -> str:
    return os.environ.get("VEP_CACHE_REDIS_URL") or os.environ.get("REDIS_URL") or VEP_CACHE_REDIS_URL_DEFAULT


def _vep_cache_key(vep_input: VepInput, ensembl_release: str) -> str:
    """Cache key for one input, namespaced by *both* versioning axes plus the resolution question.

    Keyed on the HGVS string *and* the transcript the consequence is resolved against: the same HGVS
    resolved against two transcripts is two different questions with two different answers. This script
    always derives the transcript from the HGVS string, so the transcript component is redundant
    today — it is included so that adding a transcript-hint column later cannot silently serve answers
    computed for a different transcript.

    Both versioning axes are in the key so a stored answer is only reused when it was computed under the
    same rules *and* the same upstream data (mirrors the ``(ensembl_release, RESOLVER_VERSION)`` staleness
    contract in :mod:`variant_annotation.lib.vep`):

    - ``RESOLVER_VERSION`` — a change to the resolution rule or the severity ranking.
    - ``ensembl_release`` — the Ensembl release VEP resolved against; a release boundary can change the
      answer for an unchanged input, so an answer from the prior release must not be served after a bump
      (the TTL alone would keep serving it for up to a day).
    """
    transcript = requested_transcript(vep_input) or "-"
    return f"{_vep_cache_prefix()}:r{RESOLVER_VERSION}:e{ensembl_release}:{transcript}:{vep_input.hgvs}"


def _vep_get_redis_client(*, force: bool = False):
    """Return the module-level Redis client, connecting lazily on first call.

    Returns ``None`` when caching is disabled or Redis is unreachable.
    Pass ``force=True`` to bypass the ``VEP_CACHE_ENABLED`` check.
    """
    global _VEP_REDIS_CLIENT
    global _VEP_REDIS_INIT_ATTEMPTED
    global _VEP_REDIS_UNAVAILABLE_LOGGED

    if not force and not _vep_cache_enabled():
        return None
    if _VEP_REDIS_CLIENT is not None:
        return _VEP_REDIS_CLIENT
    if _VEP_REDIS_INIT_ATTEMPTED:
        return None

    _VEP_REDIS_INIT_ATTEMPTED = True
    try:
        import redis  # type: ignore[import-not-found]

        client = redis.Redis.from_url(_vep_cache_redis_url(), decode_responses=True)
        client.ping()
        _VEP_REDIS_CLIENT = client
        logger.info("VEP Redis cache connected: %s", _vep_cache_redis_url())
        return _VEP_REDIS_CLIENT
    except Exception as exc:
        if not _VEP_REDIS_UNAVAILABLE_LOGGED:
            logger.warning("VEP Redis cache unavailable; continuing without cache: %s", exc)
            _VEP_REDIS_UNAVAILABLE_LOGGED = True
        return None


def _resolution_to_cache_value(resolution: ConsequenceResolution) -> Optional[str]:
    """Serialise a resolution for Redis, or return ``None`` when it must not be cached.

    ``ERRORED`` resolutions are never stored: they record a failed request, not an answer about the
    variant, so caching one would turn a transient outage into a persistent wrong result.
    """
    if resolution.outcome is ConsequenceOutcome.ERRORED:
        return None
    if resolution.outcome is ConsequenceOutcome.ABSENT:
        return _VEP_MISS_SENTINEL
    return json.dumps({
        "c": resolution.most_severe_consequence,
        "cs": resolution.consequence_terms,
        "s": resolution.source.value if resolution.source else None,
        "t": resolution.matched_transcript,
    })


def _cache_value_to_resolution(vep_input: VepInput, value: str) -> Optional[ConsequenceResolution]:
    """Rebuild a resolution from its Redis value, or ``None`` when the entry is unusable."""
    if value == _VEP_MISS_SENTINEL:
        return ConsequenceResolution(input=vep_input, outcome=ConsequenceOutcome.ABSENT)
    try:
        parsed = json.loads(value)
        source = parsed.get("s")
        return ConsequenceResolution(
            input=vep_input,
            outcome=ConsequenceOutcome.RESOLVED,
            consequence_terms=list(parsed.get("cs") or []),
            most_severe_consequence=parsed.get("c"),
            source=ConsequenceSource(source) if source else None,
            matched_transcript=parsed.get("t"),
        )
    except (json.JSONDecodeError, AttributeError, ValueError, TypeError):
        return None


def _vep_cache_get_many(inputs: list[VepInput], ensembl_release: str) -> dict[str, ConsequenceResolution]:
    """Batch-fetch resolutions from Redis, keyed by HGVS. Only cached inputs appear in the result."""
    client = _vep_get_redis_client()
    if client is None or not inputs:
        return {}
    try:
        values = client.mget([_vep_cache_key(i, ensembl_release) for i in inputs])
    except Exception:
        return {}

    resolutions: dict[str, ConsequenceResolution] = {}
    for vep_input, value in zip(inputs, values):
        if value is None:
            continue
        rebuilt = _cache_value_to_resolution(vep_input, value)
        if rebuilt is not None:
            resolutions[vep_input.hgvs] = rebuilt
    return resolutions


def _vep_cache_set_many(resolutions: list[ConsequenceResolution], ensembl_release: str) -> None:
    """Store resolved and confirmed-absent answers in Redis. Failures are skipped."""
    client = _vep_get_redis_client()
    if client is None or not resolutions:
        return
    try:
        pipe = client.pipeline(transaction=False)
        for resolution in resolutions:
            value = _resolution_to_cache_value(resolution)
            if value is None:
                continue
            ttl = _vep_cache_miss_ttl() if value == _VEP_MISS_SENTINEL else _vep_cache_ttl()
            pipe.set(_vep_cache_key(resolution.input, ensembl_release), value, ex=ttl)
        pipe.execute()
    except Exception as exc:
        logger.debug("VEP Redis cache write failed: %s", exc)


# ---------------------------------------------------------------------------
# Row annotation
# ---------------------------------------------------------------------------


def _get_hgvs_for_row(row: dict[str, str], hgvs_cols: list[str]) -> str:
    """Return the first non-blank HGVS value from the priority-ordered column list."""
    for col in hgvs_cols:
        val = (row.get(col) or "").strip()
        if val:
            return val
    return ""


def annotate_row(
    row: dict[str, str],
    resolutions: dict[str, ConsequenceResolution],
    *,
    col_prefix: str,
    hgvs_cols: list[str],
    access_date: str,
) -> dict[str, str]:
    """Build the annotation columns for one row from already-resolved consequences.

    Output columns are pipe-delimited and position-aligned to the row's HGVS candidates, so a blank
    candidate yields a blank in every column rather than shifting the alignment.
    """
    consequences_col = f"{col_prefix}.mutational_consequences"
    most_severe_col = f"{col_prefix}.most_severe_mutational_consequence"
    access_col = f"{col_prefix}.access_date"
    source_col = f"{col_prefix}.consequence_source"
    error_col = f"{col_prefix}.error"

    out = {
        consequences_col: "",
        most_severe_col: "",
        access_col: "",
        source_col: "",
        error_col: "",
    }

    candidates = _split_pipe_preserve_positions(_get_hgvs_for_row(row, hgvs_cols))
    if not candidates or all(c == "" for c in candidates):
        return out

    columns: dict[str, list[str]] = {
        consequences_col: [],
        most_severe_col: [],
        access_col: [],
        source_col: [],
        error_col: [],
    }

    for hgvs in candidates:
        if not hgvs:
            for values in columns.values():
                values.append("")
            continue

        resolution = resolutions.get(hgvs)
        # Every non-blank candidate was submitted for resolution, so a lookup was attempted and the
        # access date applies even when the answer is empty or unknown.
        columns[access_col].append(access_date)

        if resolution is None or resolution.outcome is ConsequenceOutcome.ABSENT:
            columns[consequences_col].append("")
            columns[most_severe_col].append("")
            columns[source_col].append("")
            columns[error_col].append("")
            continue

        if resolution.outcome is ConsequenceOutcome.ERRORED:
            columns[consequences_col].append("")
            columns[most_severe_col].append("")
            columns[source_col].append("")
            columns[error_col].append(_sanitize_for_tsv(resolution.error or "request failed"))
            continue

        columns[consequences_col].append("^".join(resolution.consequence_terms))
        columns[most_severe_col].append(resolution.most_severe_consequence or "")
        columns[source_col].append(resolution.source.value if resolution.source else "")
        columns[error_col].append("")

    return {column: "|".join(values) for column, values in columns.items()}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Annotate rows with VEP molecular consequence from HGVS candidates")
    p.add_argument("input_file", help="Input CSV/TSV file")
    p.add_argument("output_file", help="Output CSV/TSV file")
    p.add_argument(
        "--vep-namespace",
        default="vep",
        help="Namespace for output columns (default: vep)",
    )
    p.add_argument(
        "--hgvs-cols",
        default="mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p",
        metavar="COLS",
        help=(
            "Comma-separated list of columns to try for HGVS input, in priority order. "
            "The first non-blank value found is used for each row. "
            "Column order matters for transcript specificity: a transcript HGVS "
            "(e.g. NM_000049.4:c.256_257delinsAA from mapped_hgvs_c) is resolved against that "
            "specific transcript (source=transcript), while a genomic HGVS (mapped_hgvs_g) can only "
            "yield the most severe consequence across all overlapping transcripts "
            "(source=most_severe). (default: mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p)"
        ),
    )
    p.add_argument(
        "--vep-api-url",
        default=ENSEMBL_API_URL_DEFAULT,
        help="Base URL for Ensembl REST API (default from ENSEMBL_API_URL env var)",
    )
    p.add_argument(
        "--vep-batch-size",
        type=int,
        default=int(os.environ.get("VEP_BATCH_SIZE", str(VEP_MAX_HGVS_PER_REQUEST))),
        help=(
            f"Number of HGVS values per VEP request batch (default and maximum: "
            f"{VEP_MAX_HGVS_PER_REQUEST}, Ensembl's documented POST limit)"
        ),
    )
    p.add_argument(
        "--vep-workers",
        type=int,
        default=int(os.environ.get("VEP_WORKERS", "8")),
        help="Number of concurrent VEP batch requests (default: 8)",
    )
    p.add_argument(
        "--row-batch-size",
        type=int,
        default=int(os.environ.get("VEP_ROW_BATCH_SIZE", "1000")),
        help="Number of input rows per lookup/write batch (default: 1000)",
    )
    p.add_argument(
        "--vep-timeout-seconds",
        type=int,
        default=int(os.environ.get("VEP_TIMEOUT_SECONDS", "60")),
        help="HTTP timeout in seconds for VEP API calls (default: 60)",
    )
    p.add_argument(
        "--no-recoder",
        action="store_true",
        default=False,
        help=(
            "Skip the Variant Recoder fallback. Inputs VEP cannot resolve are reported as having no "
            "consequence instead of being recoded to genomic and re-queried. Use when a "
            "cross-transcript answer is worse than no answer."
        ),
    )
    p.add_argument("--delimiter", default="\t", help="Input/output delimiter (default TAB)")
    p.add_argument("--skip", type=int, default=0, help="Number of data rows to skip before annotation")
    p.add_argument("--limit", type=int, default=None, help="Maximum number of data rows to annotate")
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    p.add_argument(
        "--csv-field-size-limit",
        type=int,
        default=csv.field_size_limit(),
        metavar="BYTES",
        help="Maximum per-field character length for CSV/TSV parsing (default: %(default)s).",
    )
    p.add_argument(
        "--keep-existing",
        action="store_true",
        default=False,
        help=(
            "Preserve rows that already have a non-empty annotation value and only "
            "annotate rows with a blank consequence column. Reports how many rows "
            "were newly annotated."
        ),
    )
    return p.parse_args(argv)


def _open_reference_sequence() -> Optional[ReferenceSequence]:
    """Build the UTA-backed reference reader, or ``None`` when UTA is not configured/reachable.

    One connection for the whole run: UTA's reserved connection budget is exhausted by opening one per
    lookup. Absence is not fatal — it only means reference-identical rows cannot be labelled.
    """
    database_url = os.environ.get("UTA_DB_URL")
    if not database_url:
        logger.info("UTA_DB_URL is not set; reference-identical (no_change) rows will not be labelled.")
        return None
    try:
        return UtaClient(connect_uta(database_url))
    except Exception as exc:
        logger.warning("UTA unavailable; reference-identical rows will not be labelled: %s", exc)
        return None


def _resolve_ensembl_release(client: EnsemblRestClient) -> str:
    """Return the Ensembl release for cache keying, or ``"unknown"`` if it cannot be determined.

    The release is one of the two axes the cache key is namespaced by, so it is resolved once per run.
    Best-effort: a mirror that does not serve ``/info/software`` must not fail the run — it degrades to
    ``"unknown"``, under which cached entries simply never collide with release-keyed ones (correctness
    is preserved; only cross-run reuse is lost until the endpoint answers).
    """
    try:
        return client.software_release()
    except Exception as exc:
        logger.warning(
            "Could not determine Ensembl release for cache keying (%s); using 'unknown'. Cached answers "
            "will not be shared with release-keyed runs.",
            exc,
        )
        return "unknown"


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )
    csv.field_size_limit(args.csv_field_size_limit)

    if args.skip < 0:
        logger.error("--skip must be >= 0, got: %d", args.skip)
        sys.exit(1)
    if args.limit is not None and args.limit < 1:
        logger.error("--limit must be >= 1 when provided, got: %d", args.limit)
        sys.exit(1)
    if args.vep_batch_size < 1 or args.vep_batch_size > VEP_MAX_HGVS_PER_REQUEST:
        logger.error(
            "--vep-batch-size must be between 1 and %d, got: %d",
            VEP_MAX_HGVS_PER_REQUEST,
            args.vep_batch_size,
        )
        sys.exit(1)
    if args.vep_workers < 1:
        logger.error("--vep-workers must be >= 1, got: %d", args.vep_workers)
        sys.exit(1)
    if args.row_batch_size < 1:
        logger.error("--row-batch-size must be >= 1, got: %d", args.row_batch_size)
        sys.exit(1)
    if args.vep_timeout_seconds < 1:
        logger.error("--vep-timeout-seconds must be >= 1, got: %d", args.vep_timeout_seconds)
        sys.exit(1)

    delim = "\t" if args.delimiter == "\\t" else args.delimiter
    input_path = Path(args.input_file)
    output_path = Path(args.output_file)

    hgvs_cols = [c.strip() for c in args.hgvs_cols.split(",") if c.strip()]
    if not hgvs_cols:
        logger.error("--hgvs-cols must contain at least one column name")
        sys.exit(1)

    prefix = args.vep_namespace
    access_date = date.today().isoformat()
    ann_cols = [
        f"{prefix}.mutational_consequences",
        f"{prefix}.most_severe_mutational_consequence",
        f"{prefix}.access_date",
        f"{prefix}.consequence_source",
        f"{prefix}.error",
    ]

    ensembl = EnsemblRestClient(
        api_url=args.vep_api_url,
        timeout_seconds=args.vep_timeout_seconds,
    )
    ensembl_release = _resolve_ensembl_release(ensembl)
    reference = _open_reference_sequence()
    config = VepConfig(
        batch_size=args.vep_batch_size,
        max_workers=args.vep_workers,
        recode_misses=not args.no_recoder,
    )

    # Resolutions accumulate across row batches: the same HGVS recurs constantly in reverse-translated
    # data, and one answer per distinct string is enough for the whole run.
    resolutions: dict[str, ConsequenceResolution] = {}
    total_rows = 0
    kept_rows = 0
    newly_resolved_rows = 0
    discrepancy_rows = 0

    with (
        input_path.open("r", encoding="utf-8", newline="") as in_fh,
        output_path.open("w", encoding="utf-8", newline="") as out_fh,
    ):
        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears empty: %s", input_path)
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

        selected_reader = islice(
            reader,
            args.skip,
            None if args.limit is None else args.skip + args.limit,
        )

        def flush(batch_rows: list[dict[str, str]]) -> None:
            nonlocal total_rows, kept_rows, newly_resolved_rows, discrepancy_rows
            batch_total, batch_kept, batch_resolved, batch_discrepancies = _process_batch(
                batch_rows,
                writer,
                out_fh,
                resolutions,
                hgvs_cols=hgvs_cols,
                col_prefix=prefix,
                access_date=access_date,
                ensembl=ensembl,
                reference=reference,
                config=config,
                ensembl_release=ensembl_release,
                keep_existing=args.keep_existing,
            )
            total_rows += batch_total
            kept_rows += batch_kept
            newly_resolved_rows += batch_resolved
            discrepancy_rows += batch_discrepancies

        batch_rows: list[dict[str, str]] = []
        for row in selected_reader:
            batch_rows.append(row)
            if len(batch_rows) >= args.row_batch_size:
                flush(batch_rows)
                batch_rows = []

        if batch_rows:
            flush(batch_rows)

    logger.info(
        "Done. %d rows processed; %d kept (existing); %d rows newly annotated; %d discrepancy rows -> %s",
        total_rows,
        kept_rows,
        newly_resolved_rows,
        discrepancy_rows,
        output_path,
    )


def _process_batch(
    rows: list[dict[str, str]],
    writer: csv.DictWriter,
    out_fh: TextIO,
    resolutions: dict[str, ConsequenceResolution],
    *,
    hgvs_cols: list[str],
    col_prefix: str,
    access_date: str,
    ensembl: EnsemblRestClient,
    reference: Optional[ReferenceSequence],
    config: VepConfig,
    ensembl_release: str,
    keep_existing: bool = False,
) -> tuple[int, int, int, int]:
    """Resolve this batch's unseen HGVS candidates, then write every row.

    Returns ``(total, kept, newly_resolved, discrepancies)``. *kept* counts rows whose existing
    annotation was preserved (only non-zero when *keep_existing* is True).
    """
    most_severe_col = f"{col_prefix}.most_severe_mutational_consequence"
    error_col = f"{col_prefix}.error"

    pending: list[VepInput] = []
    seen: set[str] = set()
    for row in rows:
        if keep_existing and (row.get(most_severe_col) or "").strip():
            continue
        for hgvs in _split_pipe_preserve_positions(_get_hgvs_for_row(row, hgvs_cols)):
            if hgvs and hgvs not in resolutions and hgvs not in seen:
                seen.add(hgvs)
                pending.append(VepInput(hgvs=hgvs))

    if pending:
        # Satisfy as many lookups as possible from Redis before hitting the API.
        cached = _vep_cache_get_many(pending, ensembl_release)
        if cached:
            resolutions.update(cached)
            logger.debug("VEP Redis cache: %d/%d hits", len(cached), len(pending))
        pending = [i for i in pending if i.hgvs not in resolutions]

    if pending:
        fresh = resolve_consequences(
            pending,
            vep=ensembl,
            recoder=ensembl if config.recode_misses else None,
            reference=reference,
            config=config,
        )
        resolutions.update({r.input.hgvs: r for r in fresh})
        _vep_cache_set_many(fresh, ensembl_release)

    kept = 0
    newly_resolved = 0
    discrepancies = 0
    for row in rows:
        if keep_existing and (row.get(most_severe_col) or "").strip():
            kept += 1
            writer.writerow(row)
            continue
        ann = annotate_row(
            row,
            resolutions,
            col_prefix=col_prefix,
            hgvs_cols=hgvs_cols,
            access_date=access_date,
        )
        row.update(ann)
        writer.writerow(row)
        if ann[most_severe_col]:
            newly_resolved += 1
        elif ann[error_col]:
            discrepancies += 1

    out_fh.flush()
    return len(rows), kept, newly_resolved, discrepancies


if __name__ == "__main__":
    main()
