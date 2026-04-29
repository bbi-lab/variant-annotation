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
from itertools import islice
import logging
import os
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

GNOMAD_DATA_VERSION = os.environ.get("GNOMAD_DATA_VERSION", "v4.1")
GNOMAD_HT_URI_DEFAULT = os.environ.get(
    "GNOMAD_HT_URI",
    "gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht",
)
HAIL_GCS_CONNECTOR_JAR_DEFAULT = (
    "https://repo1.maven.org/maven2/com/google/cloud/bigdataoss/"
    "gcs-connector/hadoop3-2.2.11/gcs-connector-hadoop3-2.2.11-shaded.jar"
)
GNOMAD_ATHENA_DB_DEFAULT = os.environ.get("GNOMAD_ATHENA_DATABASE", "gnomad")
GNOMAD_ATHENA_TABLE_DEFAULT = os.environ.get("GNOMAD_ATHENA_TABLE", "")
GNOMAD_ATHENA_OUTPUT_LOCATION_DEFAULT = os.environ.get("GNOMAD_ATHENA_OUTPUT_LOCATION", "")
GNOMAD_ATHENA_WORKGROUP_DEFAULT = os.environ.get("GNOMAD_ATHENA_WORKGROUP", "")
GNOMAD_ATHENA_REGION_DEFAULT = os.environ.get("GNOMAD_ATHENA_REGION", os.environ.get("AWS_REGION", ""))
GNOMAD_ATHENA_ROW_BATCH_SIZE_DEFAULT = int(os.environ.get("GNOMAD_ATHENA_ROW_BATCH_SIZE", "1000"))


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


def _is_gs_uri(uri: str) -> bool:
    return (uri or "").strip().lower().startswith("gs://")


def _hail_init_kwargs(tmp_dir: Path, source_uri: str) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "tmp_dir": str(tmp_dir),
        "quiet": True,
        "idempotent": True,
    }
    if not _is_gs_uri(source_uri):
        return kwargs

    spark_conf = {
        "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
        "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",
    }
    gcs_auth_type = os.environ.get("HAIL_GCS_AUTH_TYPE", "UNAUTHENTICATED").strip()
    if gcs_auth_type:
        spark_conf["spark.hadoop.fs.gs.auth.type"] = gcs_auth_type
        if gcs_auth_type.upper() == "UNAUTHENTICATED":
            spark_conf["spark.hadoop.google.cloud.auth.service.account.enable"] = "false"
            spark_conf["spark.hadoop.google.cloud.auth.null.enable"] = "true"
            spark_conf["spark.hadoop.fs.gs.auth.null.enable"] = "true"

    base_jars = os.environ.get("HAIL_SPARK_JARS", "").strip()
    gcs_jar = os.environ.get("HAIL_GCS_CONNECTOR_JAR", HAIL_GCS_CONNECTOR_JAR_DEFAULT).strip()
    jars = ",".join(part for part in [base_jars, gcs_jar] if part)
    if jars:
        spark_conf["spark.jars"] = jars

    base_packages = os.environ.get("HAIL_SPARK_JARS_PACKAGES", "").strip()
    gcs_packages = os.environ.get("HAIL_GCS_CONNECTOR_PACKAGES", "").strip()
    packages = ",".join(part for part in [base_packages, gcs_packages] if part)
    if packages:
        spark_conf["spark.jars.packages"] = packages

    kwargs["spark_conf"] = spark_conf
    return kwargs


def _raise_actionable_hail_error(source_uri: str, exc: Exception) -> None:
    message = str(exc)
    if _is_gs_uri(source_uri) and "UnsupportedFileSystemException" in message and "scheme \"gs\"" in message:
        raise RuntimeError(
            "Hail cannot read gs:// paths because the Google Cloud Storage connector is unavailable "
            "in Spark/Hadoop. This script now configures a GCS connector automatically. "
            "Rebuild the annotate-gnomad image and retry. If this still fails in your environment, "
            "download the gnomAD Hail table locally and pass a local --gnomad-ht-uri path. "
            f"Original error: {message}"
        ) from exc
    if _is_gs_uri(source_uri) and "No valid credential configuration discovered" in message:
        raise RuntimeError(
            "Hail reached the Google Cloud Storage connector, but no usable credential mode was configured. "
            "For public gnomAD buckets this command defaults to unauthenticated access; if your environment "
            "overrides auth settings, set HAIL_GCS_AUTH_TYPE=UNAUTHENTICATED. "
            f"Original error: {message}"
        ) from exc
    raise exc


def _split_pipe(value: str) -> list[str]:
    if "|" not in value:
        return [value.strip()] if value.strip() else []
    return [part.strip() for part in value.split("|") if part.strip()]


def _split_pipe_preserve_positions(value: str) -> list[str]:
    raw = value or ""
    if "|" not in raw:
        return [raw.strip()]
    return [part.strip() for part in raw.split("|")]


def _batched(values: list[str], batch_size: int) -> list[list[str]]:
    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")
    return [values[i : i + batch_size] for i in range(0, len(values), batch_size)]


def _athena_table_name_from_version(version: str) -> str:
    return version.replace(".", "_")


def _load_athena_rows_for_caids(
    caids: list[str],
    *,
    database: str,
    table: str,
    output_location: str,
    workgroup: Optional[str],
    region: Optional[str],
    max_caids_per_query: int,
    poll_seconds: int,
) -> list[dict[str, Optional[str]]]:
    try:
        import boto3  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "boto3 is required for --execution-mode athena. Install dependencies with boto3 available."
        ) from exc

    if not caids:
        return []

    if not output_location:
        raise ValueError("Athena output location is required")

    client_kwargs: dict[str, Any] = {}
    if region:
        client_kwargs["region_name"] = region
    client = boto3.client("athena", **client_kwargs)

    all_rows: list[dict[str, Optional[str]]] = []
    caid_chunks = _batched(caids, max_caids_per_query)
    total_chunks = len(caid_chunks)

    for chunk_index, chunk in enumerate(caid_chunks, start=1):
        caid_values = ",".join("'" + c.replace("'", "''") + "'" for c in chunk)
        query = (
            f'SELECT "caid", '
            f'"joint.freq.all.ac", '
            f'"joint.freq.all.an", '
            f'"joint.fafmax.faf95_max_gen_anc", '
            f'"joint.fafmax.faf95_max" '
            f'FROM "{database}"."{table}" '
            f'WHERE caid IN ({caid_values})'
        )

        logger.info(
            "Athena gnomAD query %d/%d: requesting %d CAIDs",
            chunk_index,
            total_chunks,
            len(chunk),
        )

        start_kwargs: dict[str, Any] = {
            "QueryString": query,
            "QueryExecutionContext": {"Database": database},
            "ResultConfiguration": {"OutputLocation": output_location},
        }
        if workgroup:
            start_kwargs["WorkGroup"] = workgroup

        start_resp = client.start_query_execution(**start_kwargs)
        execution_id = start_resp["QueryExecutionId"]

        while True:
            execution = client.get_query_execution(QueryExecutionId=execution_id)
            status = execution["QueryExecution"]["Status"]["State"]
            if status in {"SUCCEEDED", "FAILED", "CANCELLED"}:
                break
            time.sleep(max(1, poll_seconds))

        if status != "SUCCEEDED":
            reason = execution["QueryExecution"]["Status"].get("StateChangeReason", "")
            raise RuntimeError(f"Athena query failed ({status}) for chunk {chunk_index}/{total_chunks}: {reason}")

        next_token: Optional[str] = None
        header: list[str] = []
        while True:
            result_kwargs: dict[str, Any] = {"QueryExecutionId": execution_id}
            if next_token:
                result_kwargs["NextToken"] = next_token
            result = client.get_query_results(**result_kwargs)

            result_set = result.get("ResultSet", {})
            rows = result_set.get("Rows", [])
            for row_num, row in enumerate(rows):
                values = [cell.get("VarCharValue") for cell in row.get("Data", [])]
                if not header:
                    header = [str(v or "") for v in values]
                    continue
                if len(values) < len(header):
                    values.extend([None] * (len(header) - len(values)))
                all_rows.append(dict(zip(header, values)))

            next_token = result.get("NextToken")
            if not next_token:
                break

        logger.info(
            "Athena gnomAD query %d/%d complete; cumulative rows fetched: %d",
            chunk_index,
            total_chunks,
            len(all_rows),
        )

    return all_rows


def _local_ht_path(cache_dir: Path, version: str) -> Path:
    return cache_dir / f"gnomad_{version.replace('.', '_')}_indexed.ht"


def _format_bytes(num_bytes: int) -> str:
    value = float(max(0, num_bytes))
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{value:.1f} TiB"


def _path_stats(path: Path) -> tuple[int, int]:
    if not path.exists():
        return 0, 0
    if path.is_file():
        try:
            return 1, path.stat().st_size
        except OSError:
            return 1, 0

    file_count = 0
    byte_count = 0
    for child in path.rglob("*"):
        try:
            if child.is_file():
                file_count += 1
                byte_count += child.stat().st_size
        except OSError:
            continue
    return file_count, byte_count


def _cache_progress_message(stage: str, started_at: float, cache_dir: Path, ht_path: Path) -> str:
    elapsed = max(time.monotonic() - started_at, 0.0)
    cache_files, cache_bytes = _path_stats(cache_dir)
    ht_files, ht_bytes = _path_stats(ht_path)
    success_marker = ht_path.exists() and (ht_path / "_SUCCESS").exists()
    return (
        "gnomAD cache prep progress: "
        f"stage={stage}; elapsed={elapsed / 3600.0:.2f}h; "
        f"cache_dir={cache_files} files/{_format_bytes(cache_bytes)}; "
        f"local_ht={ht_files} files/{_format_bytes(ht_bytes)}; "
        f"success_marker={'yes' if success_marker else 'no'}"
    )


class _CachePrepProgressLogger:
    def __init__(self, cache_dir: Path, ht_path: Path, interval_seconds: int) -> None:
        self._cache_dir = cache_dir
        self._ht_path = ht_path
        self._interval_seconds = interval_seconds
        self._started_at = time.monotonic()
        self._stage = "starting"
        self._stage_lock = threading.Lock()
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None

    def start(self) -> None:
        if self._interval_seconds <= 0:
            return
        self._thread = threading.Thread(target=self._run, name="gnomad-cache-progress", daemon=True)
        self._thread.start()

    def stop(self) -> None:
        if self._thread is None:
            return
        self._stop_event.set()
        self._thread.join()
        self._thread = None

    def set_stage(self, stage: str) -> None:
        with self._stage_lock:
            self._stage = stage

    def _run(self) -> None:
        while not self._stop_event.wait(self._interval_seconds):
            with self._stage_lock:
                stage = self._stage
            logger.info(_cache_progress_message(stage, self._started_at, self._cache_dir, self._ht_path))


def _has_path(dtype: Any, parts: list[Any]) -> bool:
    hl = _import_hail()
    cur = dtype
    for part in parts:
        if isinstance(part, int):
            if not isinstance(cur, hl.tarray):
                return False
            cur = cur.element_type
            continue
        if not isinstance(cur, hl.tstruct) or part not in cur:
            return False
        cur = cur[part]
    return True


def _get_path(expr: Any, parts: list[Any]) -> Any:
    cur = expr
    for part in parts:
        cur = cur[part]
    return cur


def _choose_expr(ht: Any, candidates: list[list[Any]]) -> Optional[Any]:
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
    progress_every_seconds: int = 300,
) -> Path:
    """Download and cache a local gnomAD Hail table.

    The local table is keyed by ``caid`` when the source table contains that
    field (case 1), or by a ``"chrN:pos:ref:alt"`` string (``gnomad_key``) when it
    does not (cases 2 & 3).
    """
    hl = _import_hail()

    cache_dir.mkdir(parents=True, exist_ok=True)
    ht_path = _local_ht_path(cache_dir, version)
    if ht_path.exists() and not overwrite:
        logger.info("Using cached gnomAD Hail table: %s", ht_path)
        return ht_path

    hail_tmp = cache_dir / "hail-tmp"
    hail_tmp.mkdir(parents=True, exist_ok=True)
    progress_logger = _CachePrepProgressLogger(cache_dir, ht_path, progress_every_seconds)

    logger.info("Initializing Hail and loading source table: %s", source_ht_uri)
    progress_logger.set_stage("initializing Hail")
    progress_logger.start()
    hl.init(**_hail_init_kwargs(hail_tmp, source_ht_uri))
    try:
        progress_logger.set_stage("loading source table")
        try:
            source_ht = hl.read_table(source_ht_uri)
        except Exception as exc:
            _raise_actionable_hail_error(source_ht_uri, exc)

        progress_logger.set_stage("preparing local cache projection")

        # Resolve allele counts across known gnomAD schema variants:
        # - joint v4.1 sites HT: joint.freq[0].AC / .AN
        # - browser v4.1.1 sites HT: joint.freq.all.ac / .an
        ac_expr = _choose_expr(
            source_ht,
            [
                ["joint", "freq", "all", "AC"],
                ["joint", "freq", "all", "ac"],
                ["joint", "freq", 0, "AC"],
                ["joint", "freq", 0, "ac"],
                ["freq", "all", "AC"],
                ["freq", "all", "ac"],
                ["freq", 0, "AC"],
                ["freq", 0, "ac"],
                ["joint", "AC"],
                ["joint", "ac"],
                ["AC"],
                ["ac"],
            ],
        )
        an_expr = _choose_expr(
            source_ht,
            [
                ["joint", "freq", "all", "AN"],
                ["joint", "freq", "all", "an"],
                ["joint", "freq", 0, "AN"],
                ["joint", "freq", 0, "an"],
                ["freq", "all", "AN"],
                ["freq", "all", "an"],
                ["freq", 0, "AN"],
                ["freq", 0, "an"],
                ["joint", "AN"],
                ["joint", "an"],
                ["AN"],
                ["an"],
            ],
        )

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

        common_select: dict[str, Any] = dict(
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

        caid_expr = _choose_expr(source_ht, [["caid"], ["CAID"]])
        if caid_expr is not None:
            # Case 1: source table has a caid field — build a caid-keyed local cache.
            logger.info(
                "Source gnomAD table has a 'caid' field; building caid-keyed local cache (case 1)"
            )
            prepared = source_ht.select(caid=hl.str(caid_expr), **common_select)
            prepared = prepared.key_by(prepared.caid)
        else:
            # Cases 2 & 3: source table is keyed by locus + alleles — build a
            # coordinate-keyed ("chrN:pos:ref:alt") local cache.
            logger.debug("Source gnomAD row fields when 'caid' was not found: %s", ", ".join(source_ht.row.dtype.keys()))
            logger.info(
                "Source gnomAD table has no 'caid' field; "
                "building coordinate-keyed local cache (gnomad_key) for cases 2 & 3"
            )
            gnomad_key_expr = (
                source_ht.locus.contig
                + ":"
                + hl.str(source_ht.locus.position)
                + ":"
                + source_ht.alleles[0]
                + ":"
                + source_ht.alleles[1]
            )
            prepared = source_ht.select(gnomad_key=gnomad_key_expr, **common_select)
            prepared = prepared.key_by(prepared.gnomad_key)

        logger.info("Writing local gnomAD cache table: %s", ht_path)
        progress_logger.set_stage("writing local cache table")
        prepared.write(str(ht_path), overwrite=True)
    finally:
        progress_logger.stop()
        hl.stop()

    logger.info(_cache_progress_message("completed", progress_logger._started_at, cache_dir, ht_path))

    return ht_path


def _build_caid_to_gnomad_key(
    rows: list[dict[str, str]],
    *,
    dna_col: str,
    chrom_col: str,
    pos_col: str,
    ref_col: str,
    alt_col: str,
) -> dict[str, str]:
    """Build a CAID → ``"chrN:pos:ref:alt"`` mapping from pre-computed coordinate columns.

    All four coordinate columns must hold pipe-delimited values whose order
    matches the pipe-delimited CAID candidates in *dna_col* (the cardinality
    contract maintained by ``add_vcf_identifiers``).
    """
    mapping: dict[str, str] = {}
    for row in rows:
        caids = _split_pipe(row.get(dna_col, ""))
        chroms = _split_pipe(row.get(chrom_col, ""))
        positions = _split_pipe(row.get(pos_col, ""))
        refs = _split_pipe(row.get(ref_col, ""))
        alts = _split_pipe(row.get(alt_col, ""))
        for i, caid in enumerate(caids):
            chrom = chroms[i] if i < len(chroms) else ""
            pos = positions[i] if i < len(positions) else ""
            ref = refs[i] if i < len(refs) else ""
            alt = alts[i] if i < len(alts) else ""
            if not (chrom and pos and ref and alt):
                continue
            chrom_norm = chrom if chrom.startswith("chr") else f"chr{chrom}"
            mapping[caid] = f"{chrom_norm}:{pos}:{ref}:{alt}"
    return mapping


def load_gnomad_records_for_caids(
    local_ht_path: Path,
    caids: set[str],
    cache_dir: Path,
    *,
    caid_to_gnomad_key: Optional[dict[str, str]] = None,
) -> dict[str, GnomadRecord]:
    """Load gnomAD records for requested CAIDs from local Hail table.

    Three resolution strategies are tried in order:

    **Case 1** — the local cache is keyed by ``caid`` (source table had a caid
    field).  Rows are filtered directly.

    **Case 2** — ``caid_to_gnomad_key`` mapping is provided (pre-computed
    coordinate columns from the input file).  CAIDs are translated to
    ``"chrN:pos:ref:alt"`` keys and rows are filtered by those keys.

    **Case 3** — fall back to ClinGen Allele Registry API lookups to resolve
    CAID → GRCh38 coordinates → gnomad_key.
    """
    if not caids:
        return {}

    hl = _import_hail()
    hail_tmp = cache_dir / "hail-tmp"
    hail_tmp.mkdir(parents=True, exist_ok=True)

    hl.init(**_hail_init_kwargs(hail_tmp, str(local_ht_path)))
    try:
        ht = hl.read_table(str(local_ht_path))
        # Detect the key field from the local cache schema.
        try:
            key_field = next(iter(ht.key.dtype.keys()), "gnomad_key")
        except Exception:
            key_field = "caid" if hasattr(ht, "caid") else "gnomad_key"

        if key_field == "caid":
            logger.info(
                "gnomAD lookup strategy: caid-indexed local cache (case 1) — "
                "filtering %d CAIDs directly",
                len(caids),
            )
            caid_literal = hl.literal(caids)
            key_expr = getattr(ht, key_field)
            filtered = ht.filter(caid_literal.contains(key_expr))
            rows = filtered.collect()
        else:
            # Resolve CAID → gnomad_key
            if caid_to_gnomad_key is not None:
                logger.info(
                    "gnomAD lookup strategy: pre-computed coordinate columns (case 2) — "
                    "resolved %d/%d CAIDs to gnomad_key",
                    sum(1 for c in caids if c in caid_to_gnomad_key),
                    len(caids),
                )
                resolved = {k: v for k, v in caid_to_gnomad_key.items() if k in caids}
            else:
                logger.info(
                    "gnomAD lookup strategy: ClinGen Allele Registry API lookups (case 3) — "
                    "resolving %d CAIDs to GRCh38 coordinates",
                    len(caids),
                )
                from src.lib.clingen import resolve_grch38_coordinates  # local import to keep Hail optional

                coord_cache: dict[str, Optional[tuple[str, int, str, str]]] = {}
                resolved = {}
                for caid in caids:
                    coords = resolve_grch38_coordinates(caid, coord_cache)
                    if coords is not None:
                        chrom, pos, ref, alt = coords
                        resolved[caid] = f"{chrom}:{pos}:{ref}:{alt}"
                logger.info(
                    "ClinGen resolved %d/%d CAIDs to GRCh38 coordinates",
                    len(resolved),
                    len(caids),
                )

            if not resolved:
                return {}

            gnomad_key_to_caid = {v: k for k, v in resolved.items()}
            key_literal = hl.literal(set(resolved.values()))
            key_expr = getattr(ht, key_field)
            filtered = ht.filter(key_literal.contains(key_expr))
            rows = filtered.collect()
    finally:
        hl.stop()

    out: dict[str, GnomadRecord] = {}
    for row in rows:
        if key_field == "caid":
            caid = str(row.caid)
        else:
            row_key_value = getattr(row, key_field)
            resolved_caid = gnomad_key_to_caid.get(str(row_key_value))
            if resolved_caid is None:
                continue
            caid = resolved_caid
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


def load_gnomad_records_for_caids_athena(
    caids: set[str],
    *,
    database: str,
    table: str,
    output_location: str,
    workgroup: Optional[str] = None,
    region: Optional[str] = None,
    max_caids_per_query: int = 16250,
    poll_seconds: int = 5,
) -> dict[str, GnomadRecord]:
    if not caids:
        return {}

    rows = _load_athena_rows_for_caids(
        sorted(caids),
        database=database,
        table=table,
        output_location=output_location,
        workgroup=workgroup,
        region=region,
        max_caids_per_query=max_caids_per_query,
        poll_seconds=poll_seconds,
    )

    out: dict[str, GnomadRecord] = {}
    for row in rows:
        caid = (row.get("caid") or "").strip()
        if not caid:
            continue
        try:
            ac = int(str(row.get("joint.freq.all.ac") or "0"))
            an = int(str(row.get("joint.freq.all.an") or "0"))
        except ValueError:
            continue
        if an <= 0:
            continue

        af = float(ac) / float(an)
        maf = min(af, 1.0 - af)

        faf95_max_raw = row.get("joint.fafmax.faf95_max")
        faf95_max: Optional[float]
        if faf95_max_raw is None or str(faf95_max_raw).strip() == "":
            faf95_max = None
        else:
            try:
                faf95_max = float(str(faf95_max_raw))
            except ValueError:
                faf95_max = None

        faf95_max_ancestry = str(row.get("joint.fafmax.faf95_max_gen_anc") or "")

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

    caids = _split_pipe_preserve_positions((row.get(dna_col) or "").strip())
    if not caids:
        return out

    minor_af_values: list[str] = []
    af_values: list[str] = []
    ac_values: list[str] = []
    an_values: list[str] = []
    faf95_values: list[str] = []
    faf95_anc_values: list[str] = []

    for caid in caids:
        if not caid:
            minor_af_values.append("")
            af_values.append("")
            ac_values.append("")
            an_values.append("")
            faf95_values.append("")
            faf95_anc_values.append("")
            continue
        rec = records.get(caid)
        if rec is None:
            minor_af_values.append("")
            af_values.append("")
            ac_values.append("")
            an_values.append("")
            faf95_values.append("")
            faf95_anc_values.append("")
            continue
        minor_af_values.append(str(rec.minor_allele_frequency))
        af_values.append(str(rec.allele_frequency))
        ac_values.append(str(rec.allele_count))
        an_values.append(str(rec.allele_number))
        faf95_values.append("" if rec.faf95_max is None else str(rec.faf95_max))
        faf95_anc_values.append(rec.faf95_max_ancestry)

    out[f"{col_prefix}.minor_allele_frequency"] = "|".join(minor_af_values)
    out[f"{col_prefix}.allele_frequency"] = "|".join(af_values)
    out[f"{col_prefix}.allele_count"] = "|".join(ac_values)
    out[f"{col_prefix}.allele_number"] = "|".join(an_values)
    out[f"{col_prefix}.faf95_max"] = "|".join(faf95_values)
    out[f"{col_prefix}.faf95_max_ancestry"] = "|".join(faf95_anc_values)

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
    p.add_argument(
        "--execution-mode",
        choices=["hail", "athena"],
        default="hail",
        help="Execution backend for gnomAD lookups (default: hail)",
    )
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
        "--athena-database",
        default=GNOMAD_ATHENA_DB_DEFAULT,
        help="Athena database containing gnomAD table (used in --execution-mode athena)",
    )
    p.add_argument(
        "--athena-table",
        default=GNOMAD_ATHENA_TABLE_DEFAULT,
        help="Athena table containing gnomAD data (default derived from --gnomad-version)",
    )
    p.add_argument(
        "--athena-output-location",
        default=GNOMAD_ATHENA_OUTPUT_LOCATION_DEFAULT,
        help="Athena query output S3 location (s3://...)",
    )
    p.add_argument(
        "--athena-workgroup",
        default=GNOMAD_ATHENA_WORKGROUP_DEFAULT,
        help="Optional Athena workgroup",
    )
    p.add_argument(
        "--athena-region",
        default=GNOMAD_ATHENA_REGION_DEFAULT,
        help="Optional AWS region for Athena client",
    )
    p.add_argument(
        "--athena-max-caids-per-query",
        type=int,
        default=16250,
        help="Maximum CAIDs per Athena IN query batch",
    )
    p.add_argument(
        "--athena-poll-seconds",
        type=int,
        default=5,
        help="Polling interval in seconds for Athena query execution",
    )
    p.add_argument(
        "--athena-row-batch-size",
        type=int,
        default=GNOMAD_ATHENA_ROW_BATCH_SIZE_DEFAULT,
        help="Input rows per Athena lookup/write batch (preserves input row order)",
    )
    p.add_argument(
        "--dna-clingen-allele-id-col",
        default="dna_clingen_allele_id",
        help="Column containing DNA-level ClinGen allele IDs",
    )
    # Coordinate columns (case 2): when the input file already has pre-mapped GRCh38
    # coordinates, these are used instead of ClinGen API lookups.
    p.add_argument(
        "--coord-chromosome-col",
        default="mapped_hgvs_g_chromosome",
        help="Input column with chromosome (e.g. '1' or 'chr1') for case-2 lookups "
             "(default: mapped_hgvs_g_chromosome)",
    )
    p.add_argument(
        "--coord-pos-col",
        default="mapped_hgvs_g_stop",
        help="Input column with 1-based position for case-2 lookups "
             "(default: mapped_hgvs_g_stop)",
    )
    p.add_argument(
        "--coord-ref-col",
        default="mapped_hgvs_g_ref",
        help="Input column with reference allele for case-2 lookups "
             "(default: mapped_hgvs_g_ref)",
    )
    p.add_argument(
        "--coord-alt-col",
        default="mapped_hgvs_g_alt",
        help="Input column with alternate allele for case-2 lookups "
             "(default: mapped_hgvs_g_alt)",
    )
    p.add_argument("--delimiter", default="\t", help="Input/output delimiter (default TAB)")
    p.add_argument("--skip", type=int, default=0, help="Number of data rows to skip before annotation")
    p.add_argument("--limit", type=int, default=None, help="Maximum number of data rows to annotate")
    p.add_argument("--download-only", action="store_true", help="Only materialize local gnomAD cache; do not annotate")
    p.add_argument("--refresh-cache", action="store_true", help="Rebuild local gnomAD cache even if present")
    p.add_argument(
        "--cache-progress-every-seconds",
        type=int,
        default=int(os.environ.get("GNOMAD_CACHE_PROGRESS_EVERY_SECONDS", "300")),
        help="Log cache preparation heartbeat every N seconds (0 disables heartbeat logs)",
    )
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
    if args.skip < 0:
        logger.error("--skip must be >= 0, got: %d", args.skip)
        sys.exit(1)
    if args.limit is not None and args.limit < 1:
        logger.error("--limit must be >= 1 when provided, got: %d", args.limit)
        sys.exit(1)

    local_ht: Optional[Path] = None
    cache_dir = Path(args.cache_dir)
    if args.execution_mode == "hail":
        local_ht = ensure_local_gnomad_ht(
            cache_dir,
            version=args.gnomad_version,
            source_ht_uri=args.gnomad_ht_uri,
            overwrite=args.refresh_cache,
            progress_every_seconds=args.cache_progress_every_seconds,
        )

        if args.download_only:
            logger.info("Local gnomAD cache ready at %s", local_ht)
            return
    else:
        if args.download_only:
            logger.error("--download-only is only supported in --execution-mode hail")
            sys.exit(1)
        if not args.athena_output_location:
            logger.error("--athena-output-location (or GNOMAD_ATHENA_OUTPUT_LOCATION env var) is required in athena mode")
            sys.exit(1)
        if args.athena_max_caids_per_query < 1:
            logger.error("--athena-max-caids-per-query must be >= 1, got: %d", args.athena_max_caids_per_query)
            sys.exit(1)
        if args.athena_poll_seconds < 1:
            logger.error("--athena-poll-seconds must be >= 1, got: %d", args.athena_poll_seconds)
            sys.exit(1)
        if args.athena_row_batch_size < 1:
            logger.error("--athena-row-batch-size must be >= 1, got: %d", args.athena_row_batch_size)
            sys.exit(1)
        if args.refresh_cache:
            logger.warning("--refresh-cache is ignored in --execution-mode athena")

    delim = "\t" if args.delimiter == "\\t" else args.delimiter
    input_path = Path(args.input_file)
    output_path = Path(args.output_file)

    prefix = f"{args.gnomad_namespace}.{args.gnomad_version.replace('.', '_')}"
    ann_cols = [
        f"{prefix}.minor_allele_frequency",
        f"{prefix}.allele_frequency",
        f"{prefix}.allele_count",
        f"{prefix}.allele_number",
        f"{prefix}.faf95_max",
        f"{prefix}.faf95_max_ancestry",
    ]

    if args.execution_mode == "athena":
        athena_table = args.athena_table or _athena_table_name_from_version(args.gnomad_version)
        logger.info(
            "Using Athena backend for gnomAD lookups: database=%s table=%s",
            args.athena_database,
            athena_table,
        )

        looked_up_caids: set[str] = set()
        record_cache: dict[str, GnomadRecord] = {}
        selected_rows = 0
        annotated = 0

        def _process_athena_batch(
            batch_rows: list[dict[str, str]],
            writer: csv.DictWriter,
            out_handle: Any,
        ) -> tuple[int, int]:
            if not batch_rows:
                return 0, 0

            batch_caids: set[str] = set()
            for row in batch_rows:
                batch_caids.update(_split_pipe((row.get(args.dna_clingen_allele_id_col) or "").strip()))

            missing_caids = batch_caids - looked_up_caids
            if missing_caids:
                logger.info(
                    "Athena row batch: %d rows, %d unique CAIDs (%d new lookups)",
                    len(batch_rows),
                    len(batch_caids),
                    len(missing_caids),
                )
                fetched = load_gnomad_records_for_caids_athena(
                    missing_caids,
                    database=args.athena_database,
                    table=athena_table,
                    output_location=args.athena_output_location,
                    workgroup=args.athena_workgroup or None,
                    region=args.athena_region or None,
                    max_caids_per_query=args.athena_max_caids_per_query,
                    poll_seconds=args.athena_poll_seconds,
                )
                record_cache.update(fetched)
                looked_up_caids.update(missing_caids)

            batch_annotated = 0
            for row in batch_rows:
                ann = annotate_row(row, record_cache, prefix, args.dna_clingen_allele_id_col)
                row.update(ann)
                writer.writerow(row)
                if ann[f"{prefix}.minor_allele_frequency"].replace("|", "").strip():
                    batch_annotated += 1

            out_handle.flush()
            return len(batch_rows), batch_annotated

        with input_path.open("r", encoding="utf-8", newline="") as in_fh, output_path.open(
            "w", encoding="utf-8", newline=""
        ) as out_fh:
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

            batch_rows: list[dict[str, str]] = []
            for row in selected_reader:
                batch_rows.append(row)
                if len(batch_rows) >= args.athena_row_batch_size:
                    wrote, batch_annotated = _process_athena_batch(batch_rows, writer, out_fh)
                    selected_rows += wrote
                    annotated += batch_annotated
                    batch_rows = []

            if batch_rows:
                wrote, batch_annotated = _process_athena_batch(batch_rows, writer, out_fh)
                selected_rows += wrote
                annotated += batch_annotated

        logger.info("Done. %d/%d rows annotated -> %s", annotated, selected_rows, output_path)
        return

    with input_path.open("r", encoding="utf-8", newline="") as in_fh:
        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears empty: %s", input_path)
            sys.exit(1)
        fieldnames = list(reader.fieldnames)
        rows = list(
            islice(
                reader,
                args.skip,
                None if args.limit is None else args.skip + args.limit,
            )
        )

    caids: set[str] = set()
    for row in rows:
        caids.update(_split_pipe((row.get(args.dna_clingen_allele_id_col) or "").strip()))

    logger.info("Loading gnomAD records for %d unique CAIDs", len(caids))
    if args.execution_mode == "hail":
        assert local_ht is not None
        coord_cols = [
            args.coord_chromosome_col,
            args.coord_pos_col,
            args.coord_ref_col,
            args.coord_alt_col,
        ]
        present_coord_cols = [col for col in coord_cols if col in fieldnames]
        missing_coord_cols = [col for col in coord_cols if col not in fieldnames]
        coord_cols_present = len(missing_coord_cols) == 0
        if present_coord_cols and missing_coord_cols:
            logger.warning(
                "Only some coordinate columns were found in input. Found: %s. Missing: %s. "
                "Case-2 coordinate lookup requires all four columns; falling back to case-1/3 lookup.",
                ", ".join(present_coord_cols),
                ", ".join(missing_coord_cols),
            )
        coord_mapping = (
            _build_caid_to_gnomad_key(
                rows,
                dna_col=args.dna_clingen_allele_id_col,
                chrom_col=args.coord_chromosome_col,
                pos_col=args.coord_pos_col,
                ref_col=args.coord_ref_col,
                alt_col=args.coord_alt_col,
            )
            if coord_cols_present
            else None
        )

        records = load_gnomad_records_for_caids(
            local_ht,
            caids,
            cache_dir,
            caid_to_gnomad_key=coord_mapping,
        )
    else:
        athena_table = args.athena_table or _athena_table_name_from_version(args.gnomad_version)
        logger.info(
            "Using Athena backend for gnomAD lookups: database=%s table=%s",
            args.athena_database,
            athena_table,
        )
        records = load_gnomad_records_for_caids_athena(
            caids,
            database=args.athena_database,
            table=athena_table,
            output_location=args.athena_output_location,
            workgroup=args.athena_workgroup or None,
            region=args.athena_region or None,
            max_caids_per_query=args.athena_max_caids_per_query,
            poll_seconds=args.athena_poll_seconds,
        )

    out_fieldnames = fieldnames + [c for c in ann_cols if c not in fieldnames]
    annotated = 0
    with output_path.open("w", encoding="utf-8", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_fieldnames, delimiter=delim, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            ann = annotate_row(row, records, prefix, args.dna_clingen_allele_id_col)
            row.update(ann)
            writer.writerow(row)
            if ann[f"{prefix}.minor_allele_frequency"].replace("|", "").strip():
                annotated += 1

    logger.info("Done. %d/%d rows annotated -> %s", annotated, len(rows), output_path)


if __name__ == "__main__":
    main()
