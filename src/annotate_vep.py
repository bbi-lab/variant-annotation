"""Annotate variants with VEP mutational consequences.

This script annotates each input row with VEP mutational consequences derived
from one or more pipe-delimited HGVS candidates.  By default the columns
``mapped_hgvs_c``, ``mapped_hgvs_g``, and ``mapped_hgvs_p`` are tried in that
order; the first non-blank value in a row is used (see ``--hgvs-cols``).

For rows with multiple DNA candidates (e.g. reverse translations of a protein
variant), each candidate is annotated independently.  The output columns
``vep.mutational_consequences``, ``vep.most_severe_mutational_consequence``,
``vep.consequence_source``, and ``vep.access_date`` are pipe-delimited lists
aligned to the input candidate positions.

``vep.mutational_consequences`` contains a ``^``-delimited list of all consequence
terms for the matched transcript when ``source == "transcript"``, or just the
single most-severe term when the top-level ``most_severe_consequence`` is used.
``vep.most_severe_mutational_consequence`` always contains the single most-severe
consequence term, or is empty when VEP returned no result.

Transcript selection
--------------------
Ensembl VEP annotates every transcript that overlaps the genomic position of a
variant, not just the transcript referenced in the input HGVS string.  The
top-level ``most_severe_consequence`` field therefore reflects the worst outcome
across **all** overlapping transcripts, which may differ from the consequence on
the transcript of interest.

This script addresses that by selecting the consequence from the specific input
transcript when possible:

- **RefSeq transcripts** (``NM_``, ``NR_``): the request is sent with
  ``refseq=1``, which makes VEP populate ``transcript_consequences`` with RefSeq
  transcript IDs.  The entry whose ``transcript_id`` matches the accession in the
  input HGVS string (version-stripped) is selected.  If no match is found, the
  script falls back to ``most_severe_consequence``.
- **Ensembl transcripts** (``ENST``): the default VEP transcript set is used;
  the matching ``ENST`` accession is looked up in ``transcript_consequences``
  directly.
- **Genomic HGVS** (``NC_``, ``chr:g.``) and protein HGVS (``NP_``): no
  specific transcript can be identified, so ``most_severe_consequence`` is used.

The output column ``vep.consequence_source`` records which path was taken:
``"transcript"`` for a matched transcript entry, ``"most_severe"`` for the
global fallback.

Column priority and transcript specificity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Because transcript HGVS strings (``mapped_hgvs_c``) yield more specific
consequences than genomic HGVS strings (``mapped_hgvs_g``), ``mapped_hgvs_c``
should appear first in ``--hgvs-cols``.  The default ordering
``mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p`` reflects this preference.

Redis caching
-------------
When a Redis instance is reachable, each HGVS → (consequence, source) mapping is
cached to avoid redundant Ensembl API calls across runs.  Hits are stored for 1
day by default; misses (no consequence returned by VEP) are stored for 1 day.
Caching is transparent and fails gracefully: if Redis is unreachable the script
continues without it.

Relevant environment variables:

``VEP_CACHE_ENABLED``
    Set to ``0`` / ``false`` to disable caching entirely (default: enabled).
``VEP_CACHE_REDIS_URL``
    Redis connection URL (default: ``redis://redis:6379/0``; also falls back to
    the generic ``REDIS_URL`` variable).
``VEP_CACHE_PREFIX``
    Key namespace prefix (default: ``vep:v1``).  Bump the version suffix to
    invalidate all cached entries after a significant Ensembl / VEP release.
``VEP_CACHE_TTL_SECONDS``
    TTL for cached hits, in seconds (default: 86 400 — 1 day).
``VEP_CACHE_MISS_TTL_SECONDS``
    TTL for cached misses, in seconds (default: 86 400 — 1 day).
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from itertools import islice
import logging
import os
from pathlib import Path
import re
import sys
from typing import Any, Optional, TextIO
from urllib.parse import urlsplit

import json
import requests  # type: ignore[import-untyped]

logger = logging.getLogger(__name__)

ENSEMBL_API_URL_DEFAULT = os.environ.get("ENSEMBL_API_URL", "https://rest.ensembl.org")

# ---------------------------------------------------------------------------
# Redis cache configuration
# ---------------------------------------------------------------------------

VEP_CACHE_REDIS_URL_DEFAULT = "redis://redis:6379/0"
VEP_CACHE_PREFIX_DEFAULT = "vep:v1"
VEP_CACHE_TTL_SECONDS_DEFAULT = 86400 # 1 day
VEP_CACHE_MISS_TTL_SECONDS_DEFAULT = 86400 # 1 day

_VEP_MISS_SENTINEL = "__MISS__"
_VEP_REDIS_CLIENT: Any = None
_VEP_REDIS_INIT_ATTEMPTED = False
_VEP_REDIS_UNAVAILABLE_LOGGED = False

_UTA_CONN: Any = None
_UTA_CONN_INIT_ATTEMPTED = False
_UTA_CONN_UNAVAILABLE_LOGGED = False

_C_DELINS_RE = re.compile(
    r"^(?P<tx>[^:]+):c\.(?P<start>\d+)_(?P<stop>\d+)delins(?P<alt>[ACGTNacgtn]+)$"
)

# Ordered from most to least severe (mirrors MaveDB worker logic).
VEP_CONSEQUENCES = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
]


def _split_pipe_preserve_positions(value: str) -> list[str]:
    raw = value or ""
    if "|" not in raw:
        return [raw.strip()]
    return [part.strip() for part in raw.split("|")]


def _sanitize_for_tsv(msg: str, max_length: int = 300) -> str:
    """Replace TSV/pipe-delimited field-breaking characters and truncate."""
    sanitized = msg.replace("|", ";").replace("\t", " ").replace("\r", "").replace("\n", " ").strip()
    return sanitized[:max_length] + ("..." if len(sanitized) > max_length else "")


def _build_pg_connection_kwargs(database_url: str) -> dict[str, Any]:
    parsed = urlsplit(database_url)
    path_parts = [part for part in parsed.path.split("/") if part]
    if not path_parts:
        raise RuntimeError(f"UTA_DB_URL is missing a database name: {database_url}")

    kwargs: dict[str, Any] = {"dbname": path_parts[0]}
    if parsed.hostname:
        kwargs["host"] = parsed.hostname
    if parsed.port:
        kwargs["port"] = parsed.port
    if parsed.username:
        kwargs["user"] = parsed.username
    if parsed.password:
        kwargs["password"] = parsed.password
    return kwargs


def _get_uta_connection() -> Optional[Any]:
    """Return a reusable psycopg2 connection to UTA, or None when unavailable."""
    global _UTA_CONN
    global _UTA_CONN_INIT_ATTEMPTED
    global _UTA_CONN_UNAVAILABLE_LOGGED

    if _UTA_CONN is not None:
        return _UTA_CONN
    if _UTA_CONN_INIT_ATTEMPTED:
        return None

    _UTA_CONN_INIT_ATTEMPTED = True
    uta_db_url = (os.environ.get("UTA_DB_URL") or "").strip()
    if not uta_db_url:
        if not _UTA_CONN_UNAVAILABLE_LOGGED:
            logger.warning("UTA_DB_URL not set; unchanged-variant no_change fallback disabled.")
            _UTA_CONN_UNAVAILABLE_LOGGED = True
        return None

    try:
        import psycopg2  # type: ignore[import-untyped]

        _UTA_CONN = psycopg2.connect(**_build_pg_connection_kwargs(uta_db_url))
        return _UTA_CONN
    except Exception as exc:
        if not _UTA_CONN_UNAVAILABLE_LOGGED:
            logger.warning("Cannot connect to UTA; unchanged-variant no_change fallback disabled: %s", exc)
            _UTA_CONN_UNAVAILABLE_LOGGED = True
        return None


def _lookup_transcript_ref_for_c_interval(transcript: str, c_start: int, c_stop: int) -> Optional[str]:
    """Return transcript CDS ref sequence for c. interval [c_start, c_stop], or None."""
    if c_start < 1 or c_stop < c_start:
        return None
    conn = _get_uta_connection()
    if conn is None:
        return None

    length = c_stop - c_start + 1
    try:
        with conn.cursor() as cursor:
            cursor.execute(
                """
                SELECT SUBSTRING(s.seq FROM t.cds_start_i + %s FOR %s)
                FROM uta_20241220.transcript t
                JOIN uta_20241220.seq_anno sa ON sa.ac = t.ac
                JOIN uta_20241220.seq s ON s.seq_id = sa.seq_id
                WHERE t.ac = %s
                LIMIT 1
                """,
                (c_start, length, transcript),
            )
            row = cursor.fetchone()
    except Exception as exc:
        logger.debug("UTA ref lookup failed for %s:c.%d_%d: %s", transcript, c_start, c_stop, exc)
        return None

    if row and row[0]:
        ref = str(row[0]).upper()
        if len(ref) == length:
            return ref
    return None


def _is_unchanged_transcript_delins(hgvs: str) -> bool:
    """Return True when transcript c. delins encodes no sequence change (ref == alt)."""
    m = _C_DELINS_RE.match((hgvs or "").strip())
    if not m:
        return False
    tx = m.group("tx")
    c_start = int(m.group("start"))
    c_stop = int(m.group("stop"))
    alt = m.group("alt").upper()
    ref = _lookup_transcript_ref_for_c_interval(tx, c_start, c_stop)
    return ref is not None and ref == alt


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
        logger.warning("Invalid integer for %s=%r; using default=%d", name, raw, default)
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
    return (
        os.environ.get("VEP_CACHE_REDIS_URL")
        or os.environ.get("REDIS_URL")
        or VEP_CACHE_REDIS_URL_DEFAULT
    )


def _vep_cache_key(hgvs: str) -> str:
    return f"{_vep_cache_prefix()}:{hgvs}"


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


def _vep_cache_get_many(hgvs_list: list[str]) -> dict[str, tuple[Optional[str], Optional[list[str]], str]]:
    """Batch-fetch VEP consequences from Redis.

    Returns a dict containing only the HGVS strings that were found in the
    cache.  Each value is a ``(most_severe, all_consequences, source)`` tuple.
    ``most_severe`` is ``None`` for a stored miss.  ``all_consequences`` is the
    full list of consequence terms from the matched transcript entry, or ``None``
    when only the most-severe consequence was available.  Cache entries written
    by earlier versions of this script (which lack the ``cs`` field) are
    silently discarded so they will be re-queried and overwritten in the new
    format.
    """
    client = _vep_get_redis_client()
    if client is None or not hgvs_list:
        return {}
    keys = [_vep_cache_key(h) for h in hgvs_list]
    try:
        values = client.mget(keys)
    except Exception:
        return {}
    result: dict[str, tuple[Optional[str], Optional[list[str]], str]] = {}
    for hgvs, value in zip(hgvs_list, values):
        if value is None:
            continue
        if value == _VEP_MISS_SENTINEL:
            result[hgvs] = (None, None, "most_severe")
        else:
            try:
                parsed = json.loads(value)
                if "cs" not in parsed:
                    # Old-format entry (no all_consequences field) — discard so it
                    # is re-queried and overwritten in the new format.
                    continue
                cs = parsed.get("cs") or None
                result[hgvs] = (parsed.get("c"), cs, parsed.get("s", "most_severe"))
            except (json.JSONDecodeError, AttributeError):
                # Unrecognised format — discard.
                continue
    return result


def _vep_cache_set_many(results: dict[str, tuple[Optional[str], Optional[list[str]], str]]) -> None:
    """Store a batch of HGVS → (most_severe, all_consequences, source) mappings in Redis.

    ``None`` consequences (VEP returned nothing) are stored under a sentinel
    value so that repeated no-hit queries are also short-circuited.

    Entries with ``source == "api_error"`` are skipped entirely — they result
    from a failed HTTP request rather than a definitive VEP response, and must
    not be cached so that subsequent runs can retry them.
    """
    client = _vep_get_redis_client()
    if client is None or not results:
        return
    try:
        pipe = client.pipeline(transaction=False)
        for hgvs, (consequence, all_consequences, source) in results.items():
            if source == "api_error":
                continue
            key = _vep_cache_key(hgvs)
            if consequence is not None:
                pipe.set(key, json.dumps({"c": consequence, "cs": all_consequences, "s": source}), ex=_vep_cache_ttl())
            else:
                pipe.set(key, _VEP_MISS_SENTINEL, ex=_vep_cache_miss_ttl())
        pipe.execute()
    except Exception as exc:
        logger.debug("VEP Redis cache write failed: %s", exc)


# ---------------------------------------------------------------------------


def run_variant_recoder(
    hgvs_strings: list[str],
    *,
    api_url: str,
    timeout_seconds: int,
) -> tuple[Optional[dict[str, list[str]]], str]:
    """Return ``(mapping, error)`` for a batch of HGVS strings.

    *mapping* maps each input HGVS to a list of recoded genomic HGVS strings.
    On request failure *mapping* is ``None`` and *error* contains a sanitized
    description of the failure; on success *error* is an empty string.
    """
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    out: dict[str, list[str]] = {}
    try:
        response = requests.post(
            f"{api_url}/variant_recoder/human",
            headers=headers,
            json={"ids": hgvs_strings},
            timeout=timeout_seconds,
        )
    except requests.RequestException as exc:
        logger.error("Variant Recoder request failed: %s", exc)
        return None, f"Recoder network error: {_sanitize_for_tsv(str(exc))}"

    if response.status_code != 200:
        logger.error("Variant Recoder request failed: %s %s", response.status_code, response.text)
        return None, f"Recoder HTTP {response.status_code}: {_sanitize_for_tsv(response.text)}"

    data = response.json()
    for entry in data:
        for _, variant_data in entry.items():
            # The API returns either a single result dict or a list of result dicts.
            records = variant_data if isinstance(variant_data, list) else [variant_data]
            for record in records:
                input_hgvs = record.get("input")
                if not input_hgvs:
                    continue
                genomic_hgvs_list: list[str] = []
                for genomic_hgvs in record.get("hgvsg") or []:
                    if genomic_hgvs.startswith("NC_"):
                        genomic_hgvs_list.append(genomic_hgvs)
                if genomic_hgvs_list:
                    out.setdefault(input_hgvs, []).extend(genomic_hgvs_list)
    return out, ""


def _is_refseq_transcript_hgvs(hgvs: str) -> bool:
    """Return True if the HGVS string refers to a RefSeq transcript (NM_ or NR_)."""
    return hgvs.startswith(("NM_", "NR_"))


def _extract_transcript_accession(hgvs: str) -> Optional[str]:
    """Return the versionless transcript accession from a transcript HGVS string.

    ``'NM_000049.4:c.256_257delinsAA'`` → ``'NM_000049'``
    ``'ENST00000366667.8:c.803C>T'``    → ``'ENST00000366667'``
    ``'9:g.22125504G>C'``               → ``None`` (genomic, no transcript)
    """
    if ":" not in hgvs:
        return None
    prefix = hgvs.split(":")[0]
    if "." in prefix:
        prefix = prefix.rsplit(".", 1)[0]
    if prefix.startswith(("NM_", "NR_", "ENST", "LRG_")):
        return prefix
    return None


def _pick_most_severe_from_terms(terms: list[str]) -> Optional[str]:
    """Return the most severe consequence term from a list, using VEP_CONSEQUENCES order."""
    for severity in VEP_CONSEQUENCES:
        if severity in terms:
            return severity
    return terms[0] if terms else None


def _vep_lookup_batch(
    hgvs_strings: list[str],
    *,
    api_url: str,
    timeout_seconds: int,
    refseq: bool = False,
) -> tuple[dict[str, tuple[Optional[str], Optional[list[str]], str]], Optional[str]]:
    """Look up a batch of HGVS strings against the VEP API.

    Returns ``(results, error_message)``.  On a batch-level request failure
    *results* is an empty dict and *error_message* is a sanitized description
    of the failure.  On success *error_message* is ``None``.
    """
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    out: dict[str, tuple[Optional[str], Optional[list[str]], str]] = {}
    body: dict[str, Any] = {"hgvs_notations": hgvs_strings}
    if refseq:
        body["refseq"] = 1
    try:
        response = requests.post(
            f"{api_url}/vep/human/hgvs",
            headers=headers,
            json=body,
            timeout=timeout_seconds,
        )
    except requests.RequestException as exc:
        logger.error("VEP request failed: %s", exc)
        return out, f"VEP network error: {_sanitize_for_tsv(str(exc))}"

    if response.status_code != 200:
        logger.error("VEP request failed: %s %s", response.status_code, response.text)
        return out, f"VEP HTTP {response.status_code}: {_sanitize_for_tsv(response.text)}"

    for entry in response.json():
        hgvs = entry.get("input")
        if not hgvs:
            continue
        transcript_accession = _extract_transcript_accession(hgvs)
        transcript_consequences = entry.get("transcript_consequences") or []
        matched = False
        if transcript_accession and transcript_consequences:
            for tc in transcript_consequences:
                tc_id = (tc.get("transcript_id") or "").split(".")[0]
                if tc_id == transcript_accession:
                    terms = tc.get("consequence_terms") or []
                    out[hgvs] = (_pick_most_severe_from_terms(terms), list(terms) if terms else None, "transcript")
                    matched = True
                    break
        if not matched:
            out[hgvs] = (entry.get("most_severe_consequence"), None, "most_severe")
    return out, None


def _run_batches_concurrent(
    hgvs_list: list[str],
    *,
    api_url: str,
    timeout_seconds: int,
    batch_size: int,
    max_workers: int,
) -> tuple[dict[str, tuple[Optional[str], Optional[list[str]], str]], dict[str, str]]:
    """Submit all VEP batches concurrently and merge results.

    Returns ``(results, failed_errors)`` where *failed_errors* maps each HGVS
    string whose batch failed to a sanitized error description.

    RefSeq transcript HGVS strings (NM_/NR_) are sent with ``refseq=1`` so that
    ``transcript_consequences`` entries use RefSeq IDs and can be matched back to
    the input accession.  Ensembl transcript and genomic HGVS strings are sent
    without it.
    """
    refseq_hgvs = [h for h in hgvs_list if _is_refseq_transcript_hgvs(h)]
    other_hgvs = [h for h in hgvs_list if not _is_refseq_transcript_hgvs(h)]
    batches: list[tuple[list[str], bool]] = []
    for i in range(0, len(refseq_hgvs), batch_size):
        batches.append((refseq_hgvs[i : i + batch_size], True))
    for i in range(0, len(other_hgvs), batch_size):
        batches.append((other_hgvs[i : i + batch_size], False))
    if not batches:
        return {}, {}
    result: dict[str, tuple[Optional[str], Optional[list[str]], str]] = {}
    failed_errors: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=min(max_workers, len(batches))) as pool:
        futures = {
            pool.submit(
                _vep_lookup_batch, batch, api_url=api_url, timeout_seconds=timeout_seconds, refseq=refseq
            ): batch
            for batch, refseq in batches
        }
        for fut in as_completed(futures):
            batch_results, batch_error = fut.result()
            result.update(batch_results)
            if batch_error is not None:
                for hgvs in futures[fut]:
                    if hgvs not in result:
                        failed_errors[hgvs] = batch_error
    return result, failed_errors


def get_functional_consequence(
    hgvs_strings: list[str],
    *,
    api_url: str,
    timeout_seconds: int,
    batch_size: int,
    max_workers: int = 1,
) -> dict[str, tuple[Optional[str], Optional[list[str]], str]]:
    """Return HGVS -> (most_severe, all_consequences, source) with recoder fallback.

    ``most_severe`` is the single most-severe VEP consequence term, or ``None``
    when VEP returned no result.  ``all_consequences`` is the full list of
    consequence terms taken from the matched transcript entry when
    ``source == "transcript"``; it is ``None`` when the top-level
    ``most_severe_consequence`` was used or when the recoder fallback was needed.
    ``source`` is ``"transcript"``, ``"most_severe"``, or ``"api_error"`` /
    ``"api_error:<sanitized error message>"`` when a request failed.

    Mirrors MaveDB worker behavior:
    1) VEP lookup for input HGVS (batches run concurrently when max_workers > 1)
    2) For unresolved entries, Variant Recoder to genomic HGVS
    3) VEP lookup for recoded genomic HGVS (also concurrent)
    4) Choose the most severe consequence by fixed priority order
    """
    result: dict[str, tuple[Optional[str], Optional[list[str]], str]] = {}
    if not hgvs_strings:
        return result

    # Preserve caller order while de-duplicating.
    unique_hgvs = list(dict.fromkeys(hgvs_strings))

    batch1_results, batch1_errors = _run_batches_concurrent(
        unique_hgvs,
        api_url=api_url,
        timeout_seconds=timeout_seconds,
        batch_size=batch_size,
        max_workers=max_workers,
    )
    result.update(batch1_results)

    missing_hgvs = [h for h in unique_hgvs if h not in result]
    if not missing_hgvs:
        return result

    recoded, recoder_error = run_variant_recoder(missing_hgvs, api_url=api_url, timeout_seconds=timeout_seconds)
    if recoded is None:
        # The Recoder request itself failed (network error, timeout, non-200 response).
        # Mark all missing HGVS as api_error so they are not cached and will be retried
        # on the next run.
        for missing in missing_hgvs:
            error = batch1_errors.get(missing) or recoder_error
            result[missing] = (None, None, f"api_error:{error}" if error else "api_error")
        return result

    for missing in missing_hgvs:
        if missing not in recoded:
            # Recoder succeeded but found no genomic equivalent — genuine miss.
            result[missing] = (None, None, "most_severe")

    all_recoded_hgvs: list[str] = []
    for input_hgvs in missing_hgvs:
        all_recoded_hgvs.extend(recoded.get(input_hgvs, []))

    recoded_results: dict[str, tuple[Optional[str], Optional[list[str]], str]] = {}
    batch2_errors: dict[str, str] = {}
    if all_recoded_hgvs:
        recoded_results, batch2_errors = _run_batches_concurrent(
            all_recoded_hgvs,
            api_url=api_url,
            timeout_seconds=timeout_seconds,
            batch_size=batch_size,
            max_workers=max_workers,
        )

    for input_hgvs, recoded_hgvs_list in recoded.items():
        found = {rh: recoded_results[rh] for rh in recoded_hgvs_list if rh in recoded_results}
        if not found and recoded_hgvs_list:
            # The second VEP batch returned nothing for all recoded genomic HGVS —
            # most likely a batch-level API error.  Don't cache.
            err: Optional[str] = next((batch2_errors[rh] for rh in recoded_hgvs_list if rh in batch2_errors), None)
            if not err:
                err = batch1_errors.get(input_hgvs)
            result[input_hgvs] = (None, None, f"api_error:{err}" if err else "api_error")
        else:
            consequences = [cons for cons, _cs, _src in found.values()]
            chosen = _pick_most_severe_from_terms([c for c in consequences if c])
            result[input_hgvs] = (chosen, None, "most_severe")

    return result


def _get_hgvs_for_row(row: dict[str, str], hgvs_cols: list[str]) -> str:
    """Return the first non-blank HGVS value from the priority-ordered column list."""
    for col in hgvs_cols:
        val = (row.get(col) or "").strip()
        if val:
            return val
    return ""


def annotate_row(
    row: dict[str, str],
    consequence_cache: dict[str, tuple[Optional[str], Optional[list[str]], str]],
    *,
    col_prefix: str,
    hgvs_cols: list[str],
    access_date: str,
) -> dict[str, str]:
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

    consequences_values: list[str] = []
    most_severe_values: list[str] = []
    access_values: list[str] = []
    source_values: list[str] = []
    error_values: list[str] = []

    for hgvs in candidates:
        if not hgvs:
            consequences_values.append("")
            most_severe_values.append("")
            access_values.append("")
            source_values.append("")
            error_values.append("")
            continue
        cached = consequence_cache.get(hgvs)
        if cached is None:
            access_values.append(access_date)
            if _is_unchanged_transcript_delins(hgvs):
                consequences_values.append("no_change")
                most_severe_values.append("no_change")
                source_values.append("no_change")
                error_values.append("")
            else:
                consequences_values.append("")
                most_severe_values.append("")
                source_values.append("")
                error_values.append("")
            continue
        most_severe, all_consequences, source = cached
        if source.startswith("api_error"):
            access_values.append(access_date)
            if _is_unchanged_transcript_delins(hgvs):
                consequences_values.append("no_change")
                most_severe_values.append("no_change")
                source_values.append("no_change")
                error_values.append("")
            else:
                consequences_values.append("")
                most_severe_values.append("")
                source_values.append("")
                error_values.append(source)
            continue
        most_severe_str = most_severe or ""
        if not most_severe_str and _is_unchanged_transcript_delins(hgvs):
            access_values.append(access_date)
            consequences_values.append("no_change")
            most_severe_values.append("no_change")
            source_values.append("no_change")
            error_values.append("")
            continue
        if source == "transcript" and all_consequences:
            sorted_terms = sorted(
                all_consequences,
                key=lambda t: VEP_CONSEQUENCES.index(t) if t in VEP_CONSEQUENCES else len(VEP_CONSEQUENCES),
            )
            cs_str = "^".join(sorted_terms)
        else:
            cs_str = most_severe_str
        consequences_values.append(cs_str)
        most_severe_values.append(most_severe_str)
        access_values.append(access_date)
        source_values.append(source if most_severe_str or cs_str else "")
        error_values.append("")

    out[consequences_col] = "|".join(consequences_values)
    out[most_severe_col] = "|".join(most_severe_values)
    out[access_col] = "|".join(access_values)
    out[source_col] = "|".join(source_values)
    out[error_col] = "|".join(error_values)
    return out


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Annotate rows with VEP mutational consequence using genomic HGVS candidates"
        )
    )
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
            "Column order matters for transcript specificity: when a transcript HGVS "
            "(e.g. NM_000049.4:c.256_257delinsAA from mapped_hgvs_c) is used, the "
            "script queries VEP with refseq=1 and selects the consequence for the "
            "specific input transcript (source=transcript). When a genomic HGVS "
            "(mapped_hgvs_g) is used instead, VEP returns the most severe consequence "
            "across all overlapping transcripts (source=most_severe). "
            "(default: mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p)"
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
        default=int(os.environ.get("VEP_BATCH_SIZE", "500")),
        help="Number of HGVS values per VEP/recoder request batch (default: 500)",
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


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    csv.field_size_limit(args.csv_field_size_limit)

    if args.skip < 0:
        logger.error("--skip must be >= 0, got: %d", args.skip)
        sys.exit(1)
    if args.limit is not None and args.limit < 1:
        logger.error("--limit must be >= 1 when provided, got: %d", args.limit)
        sys.exit(1)
    if args.vep_batch_size < 1:
        logger.error("--vep-batch-size must be >= 1, got: %d", args.vep_batch_size)
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

    consequence_cache: dict[str, tuple[Optional[str], Optional[list[str]], str]] = {}
    total_rows = 0
    kept_rows = 0
    newly_resolved_rows = 0
    discrepancy_rows = 0

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
            if len(batch_rows) < args.row_batch_size:
                continue

            batch_total, batch_kept, batch_resolved, batch_discrepancies = _process_batch(
                batch_rows,
                writer,
                out_fh,
                consequence_cache,
                hgvs_cols=hgvs_cols,
                col_prefix=prefix,
                access_date=access_date,
                api_url=args.vep_api_url,
                timeout_seconds=args.vep_timeout_seconds,
                vep_batch_size=args.vep_batch_size,
                vep_workers=args.vep_workers,
                keep_existing=args.keep_existing,
            )
            total_rows += batch_total
            kept_rows += batch_kept
            newly_resolved_rows += batch_resolved
            discrepancy_rows += batch_discrepancies
            batch_rows = []

        if batch_rows:
            batch_total, batch_kept, batch_resolved, batch_discrepancies = _process_batch(
                batch_rows,
                writer,
                out_fh,
                consequence_cache,
                hgvs_cols=hgvs_cols,
                col_prefix=prefix,
                access_date=access_date,
                api_url=args.vep_api_url,
                timeout_seconds=args.vep_timeout_seconds,
                vep_batch_size=args.vep_batch_size,
                vep_workers=args.vep_workers,
                keep_existing=args.keep_existing,
            )
            total_rows += batch_total
            kept_rows += batch_kept
            newly_resolved_rows += batch_resolved
            discrepancy_rows += batch_discrepancies

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
    consequence_cache: dict[str, tuple[Optional[str], Optional[list[str]], str]],
    *,
    hgvs_cols: list[str],
    col_prefix: str,
    access_date: str,
    api_url: str,
    timeout_seconds: int,
    vep_batch_size: int,
    vep_workers: int = 1,
    keep_existing: bool = False,
) -> tuple[int, int, int, int]:
    """Process a batch of rows.

    Returns ``(total, kept, newly_resolved, discrepancies)``.
    *kept* is the number of rows whose existing annotation was preserved
    (only non-zero when *keep_existing* is True).
    """
    most_severe_col = f"{col_prefix}.most_severe_mutational_consequence"
    error_col = f"{col_prefix}.error"

    batch_hgvs: list[str] = []
    for row in rows:
        if keep_existing and (row.get(most_severe_col) or "").strip():
            continue
        for hgvs in _split_pipe_preserve_positions(_get_hgvs_for_row(row, hgvs_cols)):
            if hgvs and hgvs not in consequence_cache:
                batch_hgvs.append(hgvs)

    # Satisfy as many lookups as possible from Redis before hitting the API.
    if batch_hgvs:
        unique_batch = list(dict.fromkeys(batch_hgvs))
        redis_hits = _vep_cache_get_many(unique_batch)
        if redis_hits:
            consequence_cache.update(redis_hits)
            logger.debug("VEP Redis cache: %d/%d hits", len(redis_hits), len(unique_batch))
        batch_hgvs = [h for h in batch_hgvs if h not in consequence_cache]

    if batch_hgvs:
        looked_up = get_functional_consequence(
            batch_hgvs,
            api_url=api_url,
            timeout_seconds=timeout_seconds,
            batch_size=vep_batch_size,
            max_workers=vep_workers,
        )
        consequence_cache.update(looked_up)
        _vep_cache_set_many(looked_up)

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
            consequence_cache,
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
