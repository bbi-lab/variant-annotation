"""Map variants from a CSV/TSV input file to human-genome reference HGVS strings.

Provenance
----------
This module is intentionally modeled on MaveDB variant-processing behavior and
its workflow decomposition. If this file contains copied or adapted logic from
MaveDB, treat it as AGPL-coupled when assigning a project license or when
extracting code into a separate repository.

This script is a standalone approximation of two MaveDB pipeline jobs:

1. The VRS mapping job (``map_variants_for_score_set`` in ``mavedb-api``): aligns a
   target sequence to GRCh38, maps each variant to a VRS allele using ``dcd_mapping``,
   and records an assay-level HGVS string.
2. The HGVS population job (``populate_mapped_hgvs`` in ``mavedb-api``): queries the
   ClinGen Allele Registry with the assay-level HGVS string and extracts genomic (g.),
   transcript (c.), and protein (p.) HGVS strings, plus the ClinGen allele identifier.

Variant rows fall into three categories, detected automatically:

1. **Reference-based nucleotide** – ``raw_hgvs_nt`` contains a c./n. HGVS string with
   a transcript accession prefix (e.g. ``NM_000277.3:c.1218G>A`` or
   ``ENST00000316054.9:c.1142G>A``). To mirror MaveDB behavior, the script runs this
   through the dcd_mapping accession-based normalization path to derive an assay-level
   genomic HGVS on GRCh38, then queries the ClinGen Allele Registry to populate
   ``mapped_hgvs_g``, ``mapped_hgvs_c``, and ``mapped_hgvs_p``.

2. **Sequence-based nucleotide** – ``raw_hgvs_nt`` is present but lacks a transcript
   prefix (e.g. ``c.1218G>A``). The ``target_sequence`` column is required. Rows that
   share the same ``--group-by`` column value are aligned together to GRCh38 via
   ``dcd_mapping`` (BLAT + VRS), and the resulting genomic HGVS is used to query
   ClinGen for all three strings.

3. **Protein-level only** – ``raw_hgvs_nt`` is absent or a sentinel value, but
   ``raw_hgvs_pro`` is present (e.g. ``p.Ala406Thr``). The ``target_sequence`` column
   is required; sequences may be nucleotide or amino acid. The ``dcd_mapping`` pipeline
   is invoked at the protein level and ClinGen is queried to recover
   ``mapped_hgvs_p``.

Cases 2 and 3 require ``dcd_mapping`` and its data dependencies (BLAT binary, SeqRepo
snapshot, UTA database, Gene Normalizer database). The reference genome 2-bit file
required by BLAT is downloaded automatically to the dcd_mapping local store on first
use.

Usage::

    python -m src.map_variants input.tsv output.tsv [OPTIONS]

Example::

    python -m src.map_variants variants.tsv annotated.tsv \\
        --group-by gene_symbol \\
        --raw-hgvs-nt raw_hgvs_nt \\
        --raw-hgvs-pro raw_hgvs_pro

    # Process only rows 1000–1999 (0-based, skip 1000, limit to 1000):
    python -m src.map_variants variants.tsv annotated.tsv \\
        --skip 1000 --limit 1000

Library usage::

    from src.map_variants import map_variants

    map_variants(
        input_file="variants.tsv",
        output_file="annotated.tsv",
        group_by_col="gene_symbol",
        skip=1000,
        limit=1000,
    )
"""

import asyncio
import re
import csv
import logging
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Optional

import click
from dotenv import load_dotenv
import requests
from src.lib.clingen import query_clingen_by_hgvs


load_dotenv()

logger = logging.getLogger(__name__)

ENSEMBL_REST_URL = "https://rest.ensembl.org"
PROGRESS_EVERY_ROWS = 1000


# ---------------------------------------------------------------------------
# File / format utilities
# ---------------------------------------------------------------------------


def _detect_separator(file_path: str) -> str:
    """Return the field delimiter appropriate for *file_path* based on its suffix."""
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _normalize_merge_value(value: Optional[str]) -> str:
    """Normalize row values used for merge-key lookups."""
    return (value or "").strip()


def _build_merge_key(row: dict, key_columns: tuple[str, ...]) -> tuple[str, ...]:
    """Build a stable key tuple from *row* using *key_columns*."""
    return tuple(_normalize_merge_value(row.get(col)) for col in key_columns)


def _format_exc(exc: Exception) -> str:
    """Return a readable exception string with nested message context.

    Some dependency exceptions (for example ``DataLookupError``) may carry useful
    detail in attributes like ``msg`` while ``str(exc)`` is empty.
    """

    def _one(e: BaseException) -> str:
        details: list[str] = []

        text = str(e).strip()
        if text:
            details.append(text)

        for attr in ("msg", "message", "reason", "detail"):
            value = getattr(e, attr, None)
            if isinstance(value, str):
                value = value.strip()
            if value and str(value) not in details:
                details.append(str(value))

        if not details and getattr(e, "args", None):
            details.extend(str(a) for a in e.args if a)

        if details:
            return f"{e.__class__.__name__}: {' | '.join(details)}"
        return e.__class__.__name__

    parts = [_one(exc)]
    if exc.__cause__ is not None:
        parts.append(f"cause={_one(exc.__cause__)}")
    if exc.__context__ is not None and exc.__context__ is not exc.__cause__:
        parts.append(f"context={_one(exc.__context__)}")
    return " ; ".join(parts)


# Sentinel values that should be treated as "no data" for HGVS columns.
# ---------------------------------------------------------------------------
# Protein HGVS 1-letter → 3-letter normalization
# ---------------------------------------------------------------------------

_AA1_TO_AA3: dict[str, str] = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
}

# Matches a 1-letter-coded HGVS p. expression:
#   p.{REF1}{POS}{ALT_TOKEN}{SUFFIX}
# Multi-character keywords (del, ins, dup, etc.) are tried before single-char
# so they are not split in half. The ref group is a single uppercase letter
# immediately followed by one or more digits; if the first letter is followed
# by lowercase letters (e.g. "Ala") the regex will not match, leaving already-
# normalized strings untouched.
_HGVS_P_1LETTER_RE = re.compile(
    r"^(?P<prefix>.*?p\.)(?P<ref>[A-Z])(?P<pos>\d+)"
    r"(?P<alt>del|ins|dup|ext|fs|[A-Z*=\-])"
    r"(?P<suffix>.*)$"
)


def normalize_protein_hgvs(hgvs: str) -> str:
    """Convert a protein HGVS string from 1-letter to 3-letter amino acid codes.

    Handles the following conversions:

    * Substitution: ``p.A300T`` → ``p.Ala300Thr``
    * Synonymous:   ``p.A300=`` → ``p.Ala300Ala``  (ref AA repeated)
    * Stop gain:    ``p.A300*`` → ``p.Ala300Ter``
    * Deletion:     ``p.A300-`` or ``p.A300del`` → ``p.Ala300del``
    * Frameshift:   ``p.A300fs`` → ``p.Ala300fs``  (suffix kept verbatim)
    * Insertion / duplication: keyword token kept verbatim; ref AA converted.

    Strings already using 3-letter codes (e.g. ``p.Ala300Thr``) are returned
    unchanged.  Any string that does not match the expected pattern is also
    returned unchanged.
    """
    if not hgvs:
        return hgvs
    m = _HGVS_P_1LETTER_RE.match(hgvs)
    if m is None:
        return hgvs
    ref1 = m.group("ref")
    ref3 = _AA1_TO_AA3.get(ref1)
    if ref3 is None:
        return hgvs
    alt = m.group("alt")
    if alt == "=":
        alt_out = ref3          # synonymous – repeat the reference AA
    elif alt == "*":
        alt_out = "Ter"
    elif alt == "-":
        alt_out = "del"
    elif len(alt) == 1 and alt.isupper():
        alt_out = _AA1_TO_AA3.get(alt, alt)
    else:
        alt_out = alt           # del, ins, dup, fs, ext – keep as-is
    return f"{m.group('prefix')}{ref3}{m.group('pos')}{alt_out}{m.group('suffix')}"


# Sentinel values that should be treated as "no data" for HGVS columns.
_HGVS_BLANK_SENTINELS = frozenset({"", "NA", "N/A", "None", "none", "_wt", "_sy", "="})


def _is_blank(value: Optional[str]) -> bool:
    """Return True if *value* is None, empty, or a MaveDB/HGVS empty sentinel."""
    return value is None or value.strip() in _HGVS_BLANK_SENTINELS


def _has_transcript_reference(hgvs: str) -> bool:
    """Return True if *hgvs* contains an accession prefix before a colon.

    Distinguishes ``NM_000277.3:c.1218G>A`` (True) from ``c.1218G>A`` (False).
    """
    return ":" in hgvs.strip()


def _detect_case(raw_nt: Optional[str], raw_pro: Optional[str]) -> Optional[int]:
    """Return the handling category (1, 2, 3) for a variant row, or None.

    1 – reference-based c./n. HGVS with transcript accession prefix
    2 – bare c./g. HGVS (no reference); target_sequence required
    3 – protein-only variant; target_sequence required
    None – no usable variant information
    """
    raw_nt_text = raw_nt or ""
    raw_pro_text = raw_pro or ""
    if not _is_blank(raw_nt_text):
        return 1 if _has_transcript_reference(raw_nt_text.strip()) else 2
    if not _is_blank(raw_pro_text):
        return 3
    return None


def _load_targets_file(
    targets_file: str,
    target_name_col: str,
) -> dict[str, dict[str, str]]:
    """Load a targets file and return a dict mapping target_name → row dict.

    The targets file must contain a *target_name_col* column.  All other columns
    are carried along so they can be merged onto input rows before processing.
    """
    sep = _detect_separator(targets_file)
    targets: dict[str, dict[str, str]] = {}
    with open(targets_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        if reader.fieldnames is None:
            raise ValueError(f"Targets file {targets_file!r} is empty or has no header.")
        if target_name_col not in reader.fieldnames:
            raise ValueError(
                f"Targets file {targets_file!r} is missing the required "
                f"target_name column {target_name_col!r}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )
        for row in reader:
            name = (row.get(target_name_col) or "").strip()
            if name:
                targets[name] = dict(row)
    logger.info("Loaded %d targets from %s", len(targets), targets_file)
    return targets


def _load_existing_results(
    existing_files: tuple[str, ...],
    mapped_hgvs_g_col: str,
    mapped_hgvs_c_col: str,
    mapped_hgvs_p_col: str,
    mapping_error_col: str,
    clingen_allele_id_col: str,
    key_columns: tuple[str, ...],
) -> dict[
    tuple[str, ...],
    tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]],
]:
    """Load previously computed mapping results from one or more annotated files.

    Returns a dict mapping merge-key tuples to
    ``(hgvs_c, hgvs_g, hgvs_p, error, clingen_allele_id)``.
    Rows with all-empty mapped columns and no error are ignored.
    """
    if not key_columns:
        return {}

    merged: dict[
        tuple[str, ...],
        tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]],
    ] = {}

    for file_path in existing_files:
        sep = _detect_separator(file_path)
        loaded_from_file = 0
        with open(file_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            if reader.fieldnames is None:
                logger.warning("Merge source %s is empty; skipping.", file_path)
                continue

            required = [
                *key_columns,
                mapped_hgvs_g_col,
                mapped_hgvs_c_col,
                mapped_hgvs_p_col,
                mapping_error_col,
            ]
            missing = [c for c in required if c not in reader.fieldnames]
            if missing:
                logger.warning(
                    "Merge source %s missing required columns %s; skipping.",
                    file_path,
                    ", ".join(missing),
                )
                continue

            for row in reader:
                hgvs_c = _normalize_merge_value(row.get(mapped_hgvs_c_col)) or None
                hgvs_g = _normalize_merge_value(row.get(mapped_hgvs_g_col)) or None
                hgvs_p = _normalize_merge_value(row.get(mapped_hgvs_p_col)) or None
                error = _normalize_merge_value(row.get(mapping_error_col)) or None
                clingen_allele_id = _normalize_merge_value(row.get(clingen_allele_id_col)) or None
                if hgvs_c is None and hgvs_g is None and hgvs_p is None and error is None:
                    continue

                key = _build_merge_key(row, key_columns)
                if key in merged:
                    # Keep first-seen entry to ensure deterministic precedence by file order.
                    continue
                merged[key] = (hgvs_c, hgvs_g, hgvs_p, error, clingen_allele_id)
                loaded_from_file += 1

        logger.info("Loaded %d reusable mapped rows from %s", loaded_from_file, file_path)

    return merged


# ---------------------------------------------------------------------------
# ClinGen Allele Registry helpers
# ---------------------------------------------------------------------------


def _query_clingen_by_hgvs(hgvs_string: str, max_retries: int = 3) -> Optional[dict]:
    """Query the ClinGen Allele Registry with an HGVS string.

    Returns the full JSON response dict, or None on failure.
    """
    return query_clingen_by_hgvs(hgvs_string, max_retries=max_retries, log_404=True)


async def _query_clingen_by_hgvs_batch(
    hgvs_strings: list[str],
    max_concurrency: int = 5,
) -> dict[str, Optional[dict]]:
    """Query the ClinGen Allele Registry concurrently for a batch of HGVS strings.
    
    Returns a dict mapping each hgvs_string to its result (or None if failed).
    Limits concurrent requests to max_concurrency using asyncio.Semaphore.
    """
    semaphore = asyncio.Semaphore(max_concurrency)
    
    async def _query_with_limit(hgvs: str) -> tuple[str, Optional[dict]]:
        async with semaphore:
            # Run the sync query in a thread pool to avoid blocking
            loop = asyncio.get_event_loop()
            result = await loop.run_in_executor(None, _query_clingen_by_hgvs, hgvs)
            return hgvs, result
    
    # Run all queries concurrently
    tasks = [_query_with_limit(hgvs) for hgvs in hgvs_strings]
    results = await asyncio.gather(*tasks)
    
    return {hgvs: data for hgvs, data in results}


def _clingen_allele_type(data: dict) -> str:
    """Return 'CA', 'PA', or 'unknown' for the allele type encoded in *data*."""
    at_id: str = data.get("@id", "") or ""
    fragment = at_id.rstrip("/").rsplit("/", 1)[-1]
    if fragment.startswith("CA"):
        return "CA"
    if fragment.startswith("PA"):
        return "PA"
    return "unknown"


def _extract_clingen_allele_id(data: dict) -> Optional[str]:
    """Extract ClinGen allele identifier (for example ``CA123456``) from response.

    Returns None for placeholder values such as ``_:PA...`` or ``_:CA...`` that
    ClinGen emits when no real allele record exists.
    """
    def _is_real_id(value: str) -> bool:
        # Reject blank-node-style identifiers like "_:PA..." or "_:CA..."
        return bool(value) and not value.startswith("_:")

    at_id: str = data.get("@id", "") or ""
    if at_id:
        fragment = at_id.rstrip("/").rsplit("/", 1)[-1]
        return fragment if _is_real_id(fragment) else None

    fallback = data.get("id")
    if isinstance(fallback, str):
        stripped = fallback.strip()
        if _is_real_id(stripped):
            return stripped
    return None


def _extract_hgvs_ca(
    data: dict, transcript_accession: Optional[str]
) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Extract (hgvs_g, hgvs_c, hgvs_p) from a ClinGen CA (chromosomal allele) response.

    When *transcript_accession* is provided, the matching transcript allele is
    selected directly. Otherwise the MANE transcript is used as a fallback.
    """
    hgvs_g: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None

    # Genomic – prefer GRCh38 NC_ accession.
    for allele in data.get("genomicAlleles", []):
        if allele.get("referenceGenome") == "GRCh38":
            for h in allele.get("hgvs", []):
                if h.startswith("NC_"):
                    hgvs_g = h
                    break
        if hgvs_g:
            break

    # Transcript → c. and p.
    for allele in data.get("transcriptAlleles", []):
        if transcript_accession:
            for h in allele.get("hgvs", []):
                if h.split(":")[0] == transcript_accession:
                    hgvs_c = h
                    break
            if hgvs_c:
                pe = allele.get("proteinEffect")
                if pe:
                    hgvs_p = pe.get("hgvs")
                break
        else:
            # No specific transcript requested – use MANE if available.
            mane = allele.get("MANE")
            if mane:
                hgvs_c = mane.get("nucleotide", {}).get("RefSeq", {}).get("hgvs")
                hgvs_p = mane.get("protein", {}).get("RefSeq", {}).get("hgvs")
                break

    return hgvs_g, hgvs_c, hgvs_p


def _extract_hgvs_pa(data: dict) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Extract (None, None, hgvs_p) from a ClinGen PA (protein allele) response."""
    for allele in data.get("aminoAcidAlleles", []):
        for h in allele.get("hgvs", []):
            if h:
                return None, None, h
    return None, None, None


def _extract_hgvs_from_clingen(
    data: dict, transcript_accession: Optional[str]
) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Extract (hgvs_g, hgvs_c, hgvs_p) from a ClinGen Allele Registry response."""
    if _clingen_allele_type(data) == "PA":
        return _extract_hgvs_pa(data)
    return _extract_hgvs_ca(data, transcript_accession)


# ---------------------------------------------------------------------------
# Ensembl → RefSeq transcript lookup
# ---------------------------------------------------------------------------


def _enst_to_refseq(enst_id: str) -> Optional[str]:
    """Return a RefSeq NM_ accession for *enst_id* using the Ensembl REST API.

    The version suffix (e.g. ``.9``) is stripped before the lookup. Returns None if
    no RefSeq accession is found or the request fails.
    """
    enst_base = enst_id.split(".")[0]
    try:
        url = (
            f"{ENSEMBL_REST_URL}/xrefs/id/{enst_base}"
            "?content-type=application/json&external_db=RefSeq_mRNA"
        )
        resp = requests.get(url, timeout=20)
        if resp.status_code == 200:
            for xref in resp.json():
                primary = xref.get("primary_id", "")
                if primary.startswith("NM_"):
                    return primary
    except Exception as exc:
        logger.warning("Ensembl REST lookup failed for %s: %s", enst_id, exc)
    return None


def _normalize_transcript_accession(accession: str) -> str:
    """If *accession* is an Ensembl ENST id, attempt to resolve a RefSeq NM_ id."""
    if accession.startswith("ENST"):
        refseq = _enst_to_refseq(accession)
        if refseq:
            logger.info("Mapped Ensembl %s → RefSeq %s", accession, refseq)
            return refseq
        logger.info("No RefSeq equivalent found for %s; keeping Ensembl ID.", accession)
    return accession


# ---------------------------------------------------------------------------
# Case 1 – reference-based c. HGVS
# ---------------------------------------------------------------------------


def _process_case1(
    raw_hgvs_nt: str,
    dcd: Optional[dict],
) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]:
    """Handle a reference-based c./n. HGVS variant (case 1).

    Returns ``(mapped_hgvs_c, mapped_hgvs_g, mapped_hgvs_p, error, clingen_allele_id)``.

    This path intentionally mirrors MaveDB's dcd_mapping-based handling for
    accession-referenced variants, even when no sequence alignment is required.
    """
    raw = raw_hgvs_nt.strip()
    colon_pos = raw.find(":")
    if colon_pos <= 0:
        return raw, None, None, f"Expected 'accession:variant' format; got: {raw!r}", None

    original_accession = raw[:colon_pos]
    if dcd is not None:
        fetch_clingen_genomic_hgvs = dcd["fetch_clingen_genomic_hgvs"]
        assay_level_hgvs = fetch_clingen_genomic_hgvs(raw)
    else:
        assay_level_hgvs = raw

    if assay_level_hgvs is None:
        return raw, None, None, f"ClinGen returned no data for {raw!r}", None

    data = _query_clingen_by_hgvs(assay_level_hgvs)
    if data is None:
        return raw, None, None, f"ClinGen returned no data for {assay_level_hgvs!r}", None

    hgvs_g, hgvs_c_from_clingen, hgvs_p = _extract_hgvs_from_clingen(data, original_accession)
    clingen_allele_id = _extract_clingen_allele_id(data)
    final_hgvs_c = hgvs_c_from_clingen or raw
    if hgvs_g is None and assay_level_hgvs.startswith("NC_"):
        hgvs_g = assay_level_hgvs
    return final_hgvs_c, hgvs_g, hgvs_p, None, clingen_allele_id


def _process_case1_batch(
    raw_hgvs_nt_values: list[str],
    dcd: Optional[dict],
    max_concurrency: int,
) -> list[
    tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]
]:
    """Process a batch of class-1 HGVS values, querying ClinGen concurrently.

    Returns one ``(mapped_hgvs_c, mapped_hgvs_g, mapped_hgvs_p, error,
    clingen_allele_id)`` tuple per input string in the same order.
    """
    fetch_clingen_genomic_hgvs = dcd["fetch_clingen_genomic_hgvs"] if dcd is not None else None
    prepared: list[dict] = []
    assays_to_query: list[str] = []
    assays_seen: set[str] = set()

    for raw_hgvs_nt in raw_hgvs_nt_values:
        raw = raw_hgvs_nt.strip()
        colon_pos = raw.find(":")
        if colon_pos <= 0:
            prepared.append(
                {
                    "raw": raw,
                    "original_accession": None,
                    "assay_level_hgvs": None,
                    "error": f"Expected 'accession:variant' format; got: {raw!r}",
                }
            )
            continue

        assay_level_hgvs = raw
        if fetch_clingen_genomic_hgvs is not None:
            assay_level_hgvs = fetch_clingen_genomic_hgvs(raw)

        if assay_level_hgvs is None:
            prepared.append(
                {
                    "raw": raw,
                    "original_accession": raw[:colon_pos],
                    "assay_level_hgvs": None,
                    "error": f"ClinGen returned no data for {raw!r}",
                }
            )
            continue

        prepared.append(
            {
                "raw": raw,
                "original_accession": raw[:colon_pos],
                "assay_level_hgvs": assay_level_hgvs,
                "error": None,
            }
        )
        if assay_level_hgvs not in assays_seen:
            assays_seen.add(assay_level_hgvs)
            assays_to_query.append(assay_level_hgvs)

    clingen_results: dict[str, Optional[dict]] = {}
    if assays_to_query:
        loop = asyncio.new_event_loop()
        try:
            asyncio.set_event_loop(loop)
            clingen_results = loop.run_until_complete(
                _query_clingen_by_hgvs_batch(assays_to_query, max_concurrency=max_concurrency)
            )
        finally:
            asyncio.set_event_loop(None)
            loop.close()

    results: list[
        tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]
    ] = []
    for item in prepared:
        error = item["error"]
        if error is not None:
            raw = item["raw"]
            results.append((raw, None, None, error, None))
            continue

        raw = item["raw"]
        original_accession = item["original_accession"]
        assay_level_hgvs = item["assay_level_hgvs"]
        data = clingen_results.get(assay_level_hgvs)
        if data is None:
            results.append((raw, None, None, f"ClinGen returned no data for {assay_level_hgvs!r}", None))
            continue

        hgvs_g, hgvs_c_from_clingen, hgvs_p = _extract_hgvs_from_clingen(data, original_accession)
        clingen_allele_id = _extract_clingen_allele_id(data)
        final_hgvs_c = hgvs_c_from_clingen or raw
        if hgvs_g is None and assay_level_hgvs.startswith("NC_"):
            hgvs_g = assay_level_hgvs
        results.append((final_hgvs_c, hgvs_g, hgvs_p, None, clingen_allele_id))

    return results


# ---------------------------------------------------------------------------
# Cases 2 and 3 – sequence-based mapping via dcd_mapping
# ---------------------------------------------------------------------------


# Standard GRCh38 and GRCh37 chromosome-name → RefSeq NC_ accession mappings.
# SeqRepo snapshots sometimes lack these aliases, causing dcd_mapping to fail
# when trying to convert ``GRCh38:chrN`` to an NC_ identifier during VRS mapping.
_GRCH38_CHR_NC: dict[str, str] = {
    "GRCh38:chr1": "NC_000001.11",
    "GRCh38:chr2": "NC_000002.12",
    "GRCh38:chr3": "NC_000003.12",
    "GRCh38:chr4": "NC_000004.12",
    "GRCh38:chr5": "NC_000005.10",
    "GRCh38:chr6": "NC_000006.12",
    "GRCh38:chr7": "NC_000007.14",
    "GRCh38:chr8": "NC_000008.11",
    "GRCh38:chr9": "NC_000009.12",
    "GRCh38:chr10": "NC_000010.11",
    "GRCh38:chr11": "NC_000011.10",
    "GRCh38:chr12": "NC_000012.12",
    "GRCh38:chr13": "NC_000013.11",
    "GRCh38:chr14": "NC_000014.9",
    "GRCh38:chr15": "NC_000015.10",
    "GRCh38:chr16": "NC_000016.10",
    "GRCh38:chr17": "NC_000017.11",
    "GRCh38:chr18": "NC_000018.10",
    "GRCh38:chr19": "NC_000019.10",
    "GRCh38:chr20": "NC_000020.11",
    "GRCh38:chr21": "NC_000021.9",
    "GRCh38:chr22": "NC_000022.11",
    "GRCh38:chrX": "NC_000023.11",
    "GRCh38:chrY": "NC_000024.10",
    "GRCh38:chrMT": "NC_012920.1",
    "GRCh37:chr1": "NC_000001.10",
    "GRCh37:chr2": "NC_000002.11",
    "GRCh37:chr3": "NC_000003.11",
    "GRCh37:chr4": "NC_000004.11",
    "GRCh37:chr5": "NC_000005.9",
    "GRCh37:chr6": "NC_000006.11",
    "GRCh37:chr7": "NC_000007.13",
    "GRCh37:chr8": "NC_000008.10",
    "GRCh37:chr9": "NC_000009.11",
    "GRCh37:chr10": "NC_000010.10",
    "GRCh37:chr11": "NC_000011.9",
    "GRCh37:chr12": "NC_000012.11",
    "GRCh37:chr13": "NC_000013.10",
    "GRCh37:chr14": "NC_000014.8",
    "GRCh37:chr15": "NC_000015.9",
    "GRCh37:chr16": "NC_000016.9",
    "GRCh37:chr17": "NC_000017.10",
    "GRCh37:chr18": "NC_000018.9",
    "GRCh37:chr19": "NC_000019.9",
    "GRCh37:chr20": "NC_000020.10",
    "GRCh37:chr21": "NC_000021.8",
    "GRCh37:chr22": "NC_000022.10",
    "GRCh37:chrX": "NC_000023.10",
    "GRCh37:chrY": "NC_000024.9",
    "GRCh37:chrMT": "NC_012920.1",
}

_seqrepo_chr_patch_applied = False


def _patch_seqrepo_chr_lookup() -> None:
    """Monkey-patch ``SeqRepoAccess.translate_identifier`` to fall back to a
    hardcoded chr→NC_ table when the SeqRepo volume lacks chromosome alias records.

    Applied at most once per process. Safe to call multiple times.
    """
    global _seqrepo_chr_patch_applied
    if _seqrepo_chr_patch_applied:
        return

    from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess  # noqa: PLC0415

    _original_translate = SeqRepoAccess.translate_identifier

    def _patched_translate(self, identifier: str, *args, **kwargs):
        result, err = _original_translate(self, identifier, *args, **kwargs)
        if not result and identifier in _GRCH38_CHR_NC:
            nc = _GRCH38_CHR_NC[identifier]
            logger.debug(
                "SeqRepo chr alias fallback: %s → %s", identifier, nc
            )
            return [nc], None
        return result, err

    SeqRepoAccess.translate_identifier = _patched_translate
    _seqrepo_chr_patch_applied = True
    logger.debug("Applied SeqRepo chromosome-alias fallback patch.")


def _try_import_dcd_mapping() -> dict:
    """Lazily import dcd_mapping and return a dict of required symbols.

    Raises ``ImportError`` with a user-friendly message if dcd_mapping or
    cool_seq_tool is not installed.
    """
    try:
        # dcd_mapping2 imports legacy constants from cool_seq_tool.app that were
        # removed in cool-seq-tool 0.16+. Populate them as None so imports succeed
        # while still using CoolSeqTool defaults.
        import cool_seq_tool.app as cst_app  # noqa: PLC0415

        for legacy_name in (
            "LRG_REFSEQGENE_PATH",
            "MANE_SUMMARY_PATH",
            "TRANSCRIPT_MAPPINGS_PATH",
        ):
            if not hasattr(cst_app, legacy_name):
                setattr(cst_app, legacy_name, None)

        from cool_seq_tool.schemas import AnnotationLayer  # noqa: PLC0415
        from dcd_mapping.align import build_alignment_result  # noqa: PLC0415
        from dcd_mapping.annotate import annotate  # noqa: PLC0415
        from dcd_mapping.schemas import (  # noqa: PLC0415
            ScoreRow,
            ScoresetMetadata,
            TargetGene,
            TargetSequenceType,
            TargetType,
            VrsVersion,
        )
        from dcd_mapping.transcripts import TxSelectError, select_transcripts  # noqa: PLC0415
        from dcd_mapping.vrs_map import fetch_clingen_genomic_hgvs, vrs_map  # noqa: PLC0415
    except ImportError as exc:
        raise ImportError(
            "dcd_mapping (and cool_seq_tool) are required for sequence-based mapping "
            "(cases 2 and 3). Install dcd_mapping as a Python package, e.g.:\n"
            "  pip install 'dcd-mapping @ git+https://github.com/VariantEffect/dcd_mapping2.git@mavedb-main'\n"
            f"Original ImportError: {exc}"
        ) from exc

    _patch_seqrepo_chr_lookup()

    return {
        "AnnotationLayer": AnnotationLayer,
        "build_alignment_result": build_alignment_result,
        "annotate": annotate,
        "ScoreRow": ScoreRow,
        "ScoresetMetadata": ScoresetMetadata,
        "TargetGene": TargetGene,
        "TargetSequenceType": TargetSequenceType,
        "TargetType": TargetType,
        "VrsVersion": VrsVersion,
        "TxSelectError": TxSelectError,
        "select_transcripts": select_transcripts,
        "fetch_clingen_genomic_hgvs": fetch_clingen_genomic_hgvs,
        "vrs_map": vrs_map,
    }


def _is_dna_sequence(sequence: str) -> bool:
    """Return True if *sequence* appears to be nucleotide rather than protein."""
    nucleotide_chars = frozenset("ATGCNUatgcnuRYSWKMBDHVryswkmbdhv \t\n")
    return all(c in nucleotide_chars for c in sequence)


def _hgvs_from_annotation(annotation) -> Optional[str]:
    """Extract the HGVS string from the ``post_mapped`` VRS object of a dcd_mapping
    ``ScoreAnnotationWithLayer``.

    Handles GA4GH VRS Allele (single expression) and Haplotype/CisPhasedBlock with
    exactly one member. Returns None for multi-member haplotypes (not yet supported
    downstream).
    """
    pm = annotation.post_mapped
    if pm is None:
        return None
    # GA4GH VRS Allele
    if hasattr(pm, "expressions") and pm.expressions:
        return pm.expressions[0].value
    # GA4GH VRS Haplotype / CisPhasedBlock
    if hasattr(pm, "members") and pm.members:
        if len(pm.members) != 1:
            return None  # multi-variant haplotype; not supported
        member = pm.members[0]
        if hasattr(member, "expressions") and member.expressions:
            return member.expressions[0].value
    return None


async def _run_dcd_mapping_pipeline(
    group_name: str,
    target_sequence: str,
    row_entries: list[tuple[int, str, str, int]],
    dcd: dict,
    allow_row_fallback: bool = True,
) -> tuple[list[tuple[int, Optional[str], str]], Optional[str]]:
    """Run the full dcd_mapping pipeline for one group of rows sharing a target sequence.

    Builds a synthetic ``ScoresetMetadata``, runs BLAT alignment, selects transcripts,
    performs VRS mapping, and annotates the results. Each input row is represented by a
    ``ScoreRow`` whose ``accession`` field is set to the original row index (as a string)
    so that results can be matched back unambiguously.

    Args:
        group_name: A string key identifying this group (used as the target gene name and
            for logging).
        target_sequence: The nucleotide or amino acid target sequence for the group.
        row_entries: List of ``(orig_idx, hgvs_nt, hgvs_pro, case)`` tuples. *hgvs_nt*
            and *hgvs_pro* are the raw HGVS strings from the input file. *case* is 2
            (nucleotide) or 3 (protein).
        dcd: Dict of imported dcd_mapping symbols from ``_try_import_dcd_mapping``.

    Returns:
        A 2-tuple of:
        - A list of ``(orig_idx, hgvs_assay_level, error)`` per row. *hgvs_assay_level*
          is the post-mapped HGVS string; *error* is an empty string on success.
        - The NM_ transcript accession selected by dcd_mapping, or None if unavailable.
    """
    ScoresetMetadata = dcd["ScoresetMetadata"]
    TargetGene = dcd["TargetGene"]
    TargetSequenceType = dcd["TargetSequenceType"]
    TargetType = dcd["TargetType"]
    VrsVersion = dcd["VrsVersion"]
    ScoreRow = dcd["ScoreRow"]
    TxSelectError = dcd["TxSelectError"]
    build_alignment_result = dcd["build_alignment_result"]
    select_transcripts = dcd["select_transcripts"]
    vrs_map = dcd["vrs_map"]
    annotate = dcd["annotate"]

    def _fail_all(error: str):
        return [(orig_idx, None, error) for orig_idx, _, _, _ in row_entries], None

    seq_type = (
        TargetSequenceType.DNA if _is_dna_sequence(target_sequence) else TargetSequenceType.PROTEIN
    )
    metadata = ScoresetMetadata(
        urn=f"local:{group_name}",
        target_genes={
            group_name: TargetGene(
                target_gene_name=group_name,
                target_gene_category=TargetType.PROTEIN_CODING,
                target_sequence=target_sequence,
                target_sequence_type=seq_type,
            )
        },
    )

    # Build ScoreRow objects; use the original row index as the accession for matching.
    score_rows = []
    for orig_idx, hgvs_nt, hgvs_pro, case in row_entries:
        if case == 2:
            nt = hgvs_nt if hgvs_nt else "_wt"
            pro = "_wt"
        else:  # case 3
            nt = "_wt"
            pro = hgvs_pro if hgvs_pro else "_wt"
        score_rows.append(ScoreRow(hgvs_nt=nt, hgvs_pro=pro, score=None, accession=str(orig_idx)))

    records = {group_name: score_rows}

    logger.info("Aligning target sequence for group %r (%d rows).", group_name, len(row_entries))
    try:
        alignment_results = build_alignment_result(metadata, silent=True)
    except Exception as exc:
        logger.error("BLAT alignment failed for group %r: %s", group_name, _format_exc(exc))
        logger.debug("BLAT alignment failure details:", exc_info=True)
        return _fail_all(f"BLAT alignment failed: {_format_exc(exc)}")

    logger.info("Selecting transcripts for group %r.", group_name)
    try:
        transcripts = await select_transcripts(metadata, records, alignment_results)
    except Exception as exc:
        logger.error(
            "Transcript selection failed for group %r: %s",
            group_name,
            _format_exc(exc),
        )
        logger.debug("Transcript selection full-group failure details:", exc_info=True)

        # Retry transcript selection with a single representative row. Some
        # dcd_mapping paths can fail on large record sets even when one-row
        # selection succeeds and can be reused for the full group.
        try:
            representative_records = {group_name: score_rows[:1]}
            transcripts = await select_transcripts(metadata, representative_records, alignment_results)
            logger.info(
                "Transcript selection retry (single representative row) succeeded for group %r.",
                group_name,
            )
        except Exception as retry_exc:
            logger.error(
                "Transcript selection retry failed for group %r: %s",
                group_name,
                _format_exc(retry_exc),
            )
            logger.debug("Transcript selection retry failure details:", exc_info=True)
            if allow_row_fallback and len(row_entries) > 1:
                logger.warning(
                    "Falling back to per-row mapping for group %r after transcript selection failure.",
                    group_name,
                )
                merged_results: list[tuple[int, Optional[str], str]] = []
                selected_transcript_nm: Optional[str] = None
                for orig_idx, hgvs_nt, hgvs_pro, case in row_entries:
                    one_row_results, one_tx_nm = await _run_dcd_mapping_pipeline(
                        group_name=f"{group_name}#row{orig_idx}",
                        target_sequence=target_sequence,
                        row_entries=[(orig_idx, hgvs_nt, hgvs_pro, case)],
                        dcd=dcd,
                        allow_row_fallback=False,
                    )
                    merged_results.extend(one_row_results)
                    if selected_transcript_nm is None and one_tx_nm is not None:
                        selected_transcript_nm = one_tx_nm
                return merged_results, selected_transcript_nm

            return _fail_all(f"Transcript selection failed: {_format_exc(retry_exc)}")

    transcript = transcripts.get(group_name)
    transcript_nm: Optional[str] = None
    if transcript is not None and not isinstance(transcript, TxSelectError):
        transcript_nm = getattr(transcript, "nm", None)

    align_result = alignment_results.get(group_name)
    metadata_tg = metadata.target_genes[group_name]

    logger.info("Running VRS mapping for group %r.", group_name)
    try:
        vrs_results = vrs_map(
            metadata=metadata_tg,
            align_result=align_result,
            records=score_rows,
            transcript=transcript,
            silent=True,
        )
    except Exception as exc:
        logger.error("VRS mapping failed for group %r: %s", group_name, _format_exc(exc))
        logger.debug("VRS mapping failure details:", exc_info=True)
        return _fail_all(f"VRS mapping failed: {_format_exc(exc)}")

    try:
        annotated = annotate(vrs_results, transcript, metadata_tg, metadata.urn, VrsVersion.V_2)
    except Exception as exc:
        logger.error("Annotation failed for group %r: %s", group_name, _format_exc(exc))
        logger.debug("Annotation failure details:", exc_info=True)
        return _fail_all(f"Annotation failed: {_format_exc(exc)}")

    # Index annotated results by accession (original row index as string).
    ann_by_accession: dict[str, list] = defaultdict(list)
    for ann in annotated:
        ann_by_accession[ann.mavedb_id].append(ann)

    per_row: list[tuple[int, Optional[str], str]] = []
    for orig_idx, _, _, _ in row_entries:
        ann_list = ann_by_accession.get(str(orig_idx), [])
        if not ann_list:
            per_row.append((orig_idx, None, "No mapping result returned by dcd_mapping"))
            continue
        ann = ann_list[0]
        if ann.error_message:
            per_row.append((orig_idx, None, ann.error_message))
            continue
        hgvs_assay = _hgvs_from_annotation(ann)
        per_row.append((orig_idx, hgvs_assay, ""))

    return per_row, transcript_nm


async def _process_all_groups(
    groups: dict[str, list],
    dcd: dict,
) -> dict[int, tuple[Optional[str], Optional[str], str]]:
    """Run the dcd_mapping pipeline for every group and collect per-row results.

    Args:
        groups: Dict mapping group key → list of
            ``(orig_idx, hgvs_nt, hgvs_pro, case, target_sequence)`` tuples.
        dcd: Dict of imported dcd_mapping symbols.

    Returns:
        Dict mapping ``orig_idx`` → ``(hgvs_assay_level, transcript_nm, error)``.
    """
    all_results: dict[int, tuple[Optional[str], Optional[str], str]] = {}
    for group_name, group_rows in groups.items():
        target_seq = group_rows[0][4]
        target_seqs = {r[4] for r in group_rows}
        if len(target_seqs) > 1:
            logger.warning(
                "Group %r contains %d different target_sequence values; using the first one.",
                group_name,
                len(target_seqs),
            )

        row_entries = [(r[0], r[1], r[2], r[3]) for r in group_rows]
        per_row, transcript_nm = await _run_dcd_mapping_pipeline(
            group_name, target_seq, row_entries, dcd
        )
        for orig_idx, hgvs_assay, error in per_row:
            all_results[orig_idx] = (hgvs_assay, transcript_nm, error)

    return all_results


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def map_variants(
    input_file: str,
    output_file: str,
    raw_hgvs_nt_col: str = "raw_hgvs_nt",
    raw_hgvs_pro_col: str = "raw_hgvs_pro",
    target_sequence_col: str = "target_sequence",
    target_type_col: Optional[str] = "target_type",
    group_by_col: str = "target_sequence",
    targets_file: Optional[str] = None,
    target_name_col: str = "target_name",
    mapped_hgvs_g_col: str = "mapped_hgvs_g",
    mapped_hgvs_c_col: str = "mapped_hgvs_c",
    mapped_hgvs_p_col: str = "mapped_hgvs_p",
    mapping_error_col: str = "mapping_error",
    clingen_allele_id_col: str = "clingen_allele_id",
    skip: int = 0,
    limit: Optional[int] = None,
    drop_columns: tuple[str, ...] = (),
    max_clingen_concurrency: int = 5,
    dcd_chunk_on_137: bool = True,
    dcd_chunk_size_on_137: int = 500,
    dcd_max_retry_attempts: int = 3,
    preserve_order: str = "groups",
    merge_existing_files: tuple[str, ...] = (),
    merge_match_columns: tuple[str, ...] = (),
) -> None:
    """Map variants in *input_file* to human-genome reference HGVS strings.

    Reads a CSV or TSV file (detected by file extension), processes each row according
    to the case-detection logic described in the module docstring, and writes an output
    file in the same format with additional mapped columns appended.

    Args:
        input_file: Path to the input CSV or TSV file.
        output_file: Path to write the annotated output file.
        raw_hgvs_nt_col: Column containing raw nucleotide HGVS strings.
        raw_hgvs_pro_col: Column containing raw protein HGVS strings.
        target_sequence_col: Column containing the target sequence (required for cases
            2 and 3).
        target_type_col: Optional column that distinguishes 'reference' from
            'target-sequence' rows. Informational; currently not used by the
            case-detection logic.
        group_by_col: Column used to group sequence-based rows for batch alignment.
            Rows sharing the same value are aligned together in a single BLAT run.
            Defaults to ``target_sequence``.
        mapped_hgvs_g_col: Output column name for genomic (g.) HGVS strings.
        mapped_hgvs_c_col: Output column name for transcript (c.) HGVS strings.
        mapped_hgvs_p_col: Output column name for protein (p.) HGVS strings.
        mapping_error_col: Output column name for error messages.
        clingen_allele_id_col: Output column name for ClinGen Allele Registry IDs.
        skip: Number of data rows to skip from the beginning of the file (after the
            header). Useful for resuming an interrupted run.
        limit: Maximum number of data rows to process. When None, all remaining rows
            after *skip* are processed.
        drop_columns: Column names to omit from the output file. Useful for removing
            large columns such as ``target_sequence`` that are not needed downstream.
        dcd_chunk_on_137: If True, automatically retry groups that fail with BLAT
            error code 137 using chunked processing. Disabled by default for
            static chunking.
        dcd_chunk_size_on_137: Initial chunk size applied on first BLAT 137 retry
            (default: 500).
        dcd_max_retry_attempts: Maximum number of retry attempts with progressively
            smaller chunks (default: 3).
        merge_existing_files: One or more existing annotated CSV/TSV files whose
            mapped results should be reused. Matching rows are emitted directly and
            skipped from fresh processing.
        merge_match_columns: Additional columns that must match for merge reuse,
            in addition to ``raw_hgvs_nt_col`` and ``raw_hgvs_pro_col``.
        preserve_order: Order guarantee for output rows. Options are:
            'no': Write immediately as results arrive; groups may appear
                out-of-order. Fastest mode.
            'index': Buffer by row index to guarantee exact input order. Rows emit
                incrementally when the next contiguous index is available.
            'groups' (default): Preserve input order while processing groups; buffers
                by row index and emits incrementally.
    """
    in_sep = _detect_separator(input_file)
    out_sep = _detect_separator(output_file)
    started = time.monotonic()

    # Load optional targets file and build name → target-row lookup.
    targets_lookup: dict[str, dict[str, str]] = {}
    if targets_file is not None:
        targets_lookup = _load_targets_file(targets_file, target_name_col)

    def _write_result_row(
        writer: csv.DictWriter,
        out_handle,
        row: dict,
        hgvs_c: Optional[str],
        hgvs_g: Optional[str],
        hgvs_p: Optional[str],
        error: Optional[str],
        clingen_allele_id: Optional[str],
    ) -> None:
        out_row = dict(row)
        out_row[mapped_hgvs_c_col] = hgvs_c or ""
        out_row[mapped_hgvs_g_col] = hgvs_g or ""
        out_row[mapped_hgvs_p_col] = hgvs_p or ""
        out_row[mapping_error_col] = error or ""
        out_row[clingen_allele_id_col] = clingen_allele_id or ""
        writer.writerow(out_row)
        # Keep output visible/usable during long runs.
        out_handle.flush()

    with open(input_file, newline="") as in_fh, open(output_file, "w", newline="") as out_fh:
        reader = csv.DictReader(in_fh, delimiter=in_sep)
        in_fieldnames: list[str] = list(reader.fieldnames or [])

        if not in_fieldnames:
            logger.warning("Input file %s is empty or missing a header; nothing to do.", input_file)
            return

        # Determine extra columns introduced by the targets file (not already in input).
        targets_extra_cols: list[str] = []
        if targets_lookup:
            # Use the first target row to discover available columns.
            sample_target = next(iter(targets_lookup.values()))
            targets_extra_cols = [
                c for c in sample_target if c not in in_fieldnames and c != target_name_col
            ]

        # Output columns: preserve original order and append new columns if absent.
        out_fieldnames = list(in_fieldnames)
        for col in targets_extra_cols:
            if col not in out_fieldnames:
                out_fieldnames.append(col)
        for col in [
            mapped_hgvs_g_col,
            mapped_hgvs_c_col,
            mapped_hgvs_p_col,
            mapping_error_col,
            clingen_allele_id_col,
        ]:
            if col not in out_fieldnames:
                out_fieldnames.append(col)
        if drop_columns:
            drop_set = set(drop_columns)
            out_fieldnames = [c for c in out_fieldnames if c not in drop_set]

        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=out_sep,
            extrasaction="ignore",
        )
        writer.writeheader()
        out_fh.flush()

        # Sequence-based rows (cases 2 and 3).
        # In preserve_order="groups" mode, contiguous blocks are processed as soon
        # as the block ends so downstream consumers can see output incrementally.
        # Otherwise rows are accumulated by group_key for end-of-scan processing.
        groups: dict[str, list] = defaultdict(list)

        n_total = 0
        n_written = 0
        n_errors = 0
        n_merged = 0
        next_to_write = 0
        dcd_for_case1: Optional[dict] = None
        dcd_for_case1_failed = False

        dedup_merge_match_columns = tuple(dict.fromkeys(merge_match_columns))
        merge_key_columns = (raw_hgvs_nt_col, raw_hgvs_pro_col) + tuple(
            c for c in dedup_merge_match_columns if c not in (raw_hgvs_nt_col, raw_hgvs_pro_col)
        )
        existing_results_by_key: dict[
            tuple[str, ...],
            tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]],
        ] = {}
        if merge_existing_files:
            missing_in_input = [c for c in merge_key_columns if c not in in_fieldnames]
            if missing_in_input:
                logger.warning(
                    "Merge requested but input is missing merge key columns (%s); merge disabled.",
                    ", ".join(missing_in_input),
                )
            else:
                logger.info(
                    "Merge key columns: %s",
                    ", ".join(merge_key_columns),
                )
                existing_results_by_key = _load_existing_results(
                    existing_files=merge_existing_files,
                    mapped_hgvs_g_col=mapped_hgvs_g_col,
                    mapped_hgvs_c_col=mapped_hgvs_c_col,
                    mapped_hgvs_p_col=mapped_hgvs_p_col,
                    mapping_error_col=mapping_error_col,
                    clingen_allele_id_col=clingen_allele_id_col,
                    key_columns=merge_key_columns,
                )
                logger.info(
                    "Loaded %d reusable variant results from %d merge file(s).",
                    len(existing_results_by_key),
                    len(merge_existing_files),
                )

        # Used by preserve-order modes to hold rows/results until they can be written.
        rows_by_idx: dict[int, dict] = {}
        pending_results: dict[
            int,
            tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]],
        ] = {}
        pending_case1_rows: list[tuple[int, dict, str]] = []
        current_group_key: Optional[str] = None
        current_group_rows: list[tuple[int, dict, str, str, int, str]] = []
        processed_group_count = 0
        sequence_loop: Optional[asyncio.AbstractEventLoop] = None
        dcd_for_groups: Optional[dict] = None

        def _record_result(
            idx: int,
            row: dict,
            hgvs_c: Optional[str],
            hgvs_g: Optional[str],
            hgvs_p: Optional[str],
            error: Optional[str],
            clingen_allele_id: Optional[str] = None,
            group_id: Optional[str] = None,
        ) -> None:
            nonlocal n_written, n_errors, next_to_write
            
            if preserve_order == "no":
                _write_result_row(
                    writer,
                    out_fh,
                    row,
                    hgvs_c,
                    hgvs_g,
                    hgvs_p,
                    error,
                    clingen_allele_id,
                )
                n_written += 1
                if error:
                    n_errors += 1
                return
            
            # Both "index" and "groups" modes use the same buffering logic:
            # buffer by index and emit in order. The only difference is semantic—
            # "groups" is an optimization hint that groups are contiguous.
            rows_by_idx[idx] = row
            pending_results[idx] = (hgvs_c, hgvs_g, hgvs_p, error, clingen_allele_id)
            while next_to_write in pending_results:
                next_row = rows_by_idx.pop(next_to_write)
                (
                    next_hgvs_c,
                    next_hgvs_g,
                    next_hgvs_p,
                    next_error,
                    next_clingen_allele_id,
                ) = pending_results.pop(next_to_write)
                _write_result_row(
                    writer,
                    out_fh,
                    next_row,
                    next_hgvs_c,
                    next_hgvs_g,
                    next_hgvs_p,
                    next_error,
                    next_clingen_allele_id,
                )
                n_written += 1
                if next_error:
                    n_errors += 1
                next_to_write += 1

        def _flush_case1_rows(case1_rows: list[tuple[int, dict, str]]) -> None:
            nonlocal dcd_for_case1, dcd_for_case1_failed

            if not case1_rows:
                return

            if dcd_for_case1 is None and not dcd_for_case1_failed:
                try:
                    dcd_for_case1 = _try_import_dcd_mapping()
                except ImportError as exc:
                    dcd_for_case1_failed = True
                    logger.warning(
                        "dcd_mapping unavailable for class-1 normalization; falling back to direct ClinGen query: %s",
                        exc,
                    )

            raw_hgvs_nt_values = [entry[2] for entry in case1_rows]
            if max_clingen_concurrency > 1 and len(raw_hgvs_nt_values) > 1:
                batch_results = _process_case1_batch(
                    raw_hgvs_nt_values,
                    dcd_for_case1,
                    max_concurrency=max_clingen_concurrency,
                )
            else:
                batch_results = [_process_case1(raw_nt, dcd_for_case1) for raw_nt in raw_hgvs_nt_values]

            for (idx, row, _), (hgvs_c, hgvs_g, hgvs_p, err, clingen_allele_id) in zip(case1_rows, batch_results):
                _record_result(
                    idx,
                    row,
                    hgvs_c,
                    hgvs_g,
                    hgvs_p,
                    err,
                    clingen_allele_id,
                    group_id=f"__case1_{idx}__",
                )

            case1_rows.clear()

        def _process_sequence_group(
            group_name: str,
            group_rows: list[tuple[int, dict, str, str, int, str]],
        ) -> None:
            nonlocal processed_group_count, sequence_loop, dcd_for_groups

            if not group_rows:
                return

            if dcd_for_groups is None:
                dcd_for_groups = _try_import_dcd_mapping()
            if sequence_loop is None:
                sequence_loop = asyncio.new_event_loop()

            processed_group_count += 1
            target_seq = group_rows[0][5]
            target_seqs = {r[5] for r in group_rows}
            if len(target_seqs) > 1:
                logger.warning(
                    "Group %r contains %d different target_sequence values; using the first one.",
                    group_name,
                    len(target_seqs),
                )

            group_rows_by_idx = {r[0]: r[1] for r in group_rows}
            row_entries = [(r[0], r[2], r[3], r[4]) for r in group_rows]

            per_row = None
            transcript_nm = None
            last_error = None

            try:
                asyncio.set_event_loop(sequence_loop)
                for attempt in range(1, dcd_max_retry_attempts + 1):
                    try:
                        per_row, transcript_nm = sequence_loop.run_until_complete(
                            _run_dcd_mapping_pipeline(group_name, target_seq, row_entries, dcd_for_groups)
                        )
                        break
                    except Exception as exc:
                        exc_str = str(exc)
                        last_error = exc
                        if dcd_chunk_on_137 and "error code 137" in exc_str and attempt < dcd_max_retry_attempts:
                            chunk_sz = max(1, dcd_chunk_size_on_137 // (2 ** (attempt - 1)))
                            logger.warning(
                                "Group %r failed with BLAT error 137 (attempt %d/%d); retrying with chunk_size=%d.",
                                group_name,
                                attempt,
                                dcd_max_retry_attempts,
                                chunk_sz,
                            )
                            per_row = []
                            for start_idx in range(0, len(row_entries), chunk_sz):
                                chunk_entries = row_entries[start_idx : start_idx + chunk_sz]
                                chunk_name = f"{group_name}#retry{attempt}_chunk{start_idx // chunk_sz + 1}"
                                try:
                                    chunk_per_row, chunk_tx = sequence_loop.run_until_complete(
                                        _run_dcd_mapping_pipeline(chunk_name, target_seq, chunk_entries, dcd_for_groups)
                                    )
                                    per_row.extend(chunk_per_row)
                                    if transcript_nm is None:
                                        transcript_nm = chunk_tx
                                except Exception as chunk_exc:
                                    logger.error(
                                        "Chunk %s failed: %s",
                                        chunk_name,
                                        _format_exc(chunk_exc),
                                    )
                                    for orig_idx, _, _, _ in chunk_entries:
                                        per_row.append((orig_idx, None, f"BLAT chunk failed: {_format_exc(chunk_exc)}"))
                            break
                        raise

                if per_row is None:
                    logger.error(
                        "Group %r failed after %d attempt(s): %s",
                        group_name,
                        dcd_max_retry_attempts,
                        _format_exc(last_error) if last_error else "unknown error",
                    )
                    fallback_error = _format_exc(last_error) if last_error is not None else "unknown error"
                    per_row = [(orig_idx, None, fallback_error) for orig_idx, _, _, _ in row_entries]

                hgvs_queries = []
                for orig_idx, hgvs_assay, error in per_row:
                    row = group_rows_by_idx[orig_idx]

                    if error:
                        _record_result(orig_idx, row, None, None, None, error, group_id=group_name)
                        continue
                    if not hgvs_assay:
                        _record_result(orig_idx, row, None, None, None, "No HGVS string produced by mapper", group_id=group_name)
                        continue

                    hgvs_queries.append((orig_idx, row, hgvs_assay))

                if hgvs_queries:
                    hgvs_strings = [h[2] for h in hgvs_queries]
                    logger.debug(
                        "Batch querying ClinGen for %d variants in group %r (transcript: %s)",
                        len(hgvs_strings),
                        group_name,
                        transcript_nm,
                    )
                    clingen_batch_start = time.monotonic()
                    clingen_results = sequence_loop.run_until_complete(
                        _query_clingen_by_hgvs_batch(hgvs_strings, max_concurrency=max_clingen_concurrency)
                    )
                    clingen_batch_elapsed = time.monotonic() - clingen_batch_start
                    logger.debug(
                        "ClinGen batch query completed for %d variants in %.2f seconds",
                        len(hgvs_strings),
                        clingen_batch_elapsed,
                    )

                    for orig_idx, row, hgvs_assay in hgvs_queries:
                        data = clingen_results.get(hgvs_assay)
                        if data is None:
                            _assay_is_protein = (
                                hgvs_assay.startswith("p.") or ":p." in hgvs_assay
                            )
                            _record_result(
                                orig_idx,
                                row,
                                None,
                                None if _assay_is_protein else hgvs_assay,
                                hgvs_assay if _assay_is_protein else None,
                                f"ClinGen returned no data for {hgvs_assay!r}",
                                group_id=group_name,
                            )
                            continue

                        hgvs_g, hgvs_c, hgvs_p = _extract_hgvs_from_clingen(data, transcript_nm)
                        if _clingen_allele_type(data) == "PA" and hgvs_p is None:
                            hgvs_p = hgvs_assay

                        _record_result(
                            orig_idx,
                            row,
                            hgvs_c,
                            hgvs_g,
                            hgvs_p,
                            None,
                            _extract_clingen_allele_id(data),
                            group_id=group_name,
                        )

                        if n_written % PROGRESS_EVERY_ROWS == 0:
                            elapsed = max(time.monotonic() - started, 1e-9)
                            logger.info(
                                "Output progress: %d/%d rows written (%d errors, %.1f rows/s).",
                                n_written,
                                n_total,
                                n_errors,
                                n_written / elapsed,
                            )
            finally:
                asyncio.set_event_loop(None)

            elapsed = max(time.monotonic() - started, 1e-9)
            logger.info(
                "Group %d complete (%d rows in group). %d/%d rows written so far (%.1f rows/s).",
                processed_group_count,
                len(group_rows),
                n_written,
                n_total,
                n_written / elapsed,
            )

        def _flush_current_group() -> None:
            nonlocal current_group_key, current_group_rows

            if not current_group_rows or current_group_key is None:
                return

            _process_sequence_group(current_group_key, current_group_rows)
            current_group_key = None
            current_group_rows = []

        for src_idx, row in enumerate(reader):
            if src_idx < skip:
                continue
            if limit is not None and n_total >= limit:
                break

            idx = n_total
            n_total += 1

            # Merge target-file columns onto the row (input columns take precedence).
            if targets_lookup:
                target_name = (row.get(target_name_col) or "").strip()
                target_row = targets_lookup.get(target_name)
                if target_row is None:
                    logger.warning(
                        "Row %d: target_name %r not found in targets file; "
                        "target_sequence and other target columns will be empty.",
                        src_idx,
                        target_name,
                    )
                else:
                    for col, val in target_row.items():
                        if col == target_name_col:
                            continue
                        if not row.get(col):
                            row[col] = val

            raw_nt = row.get(raw_hgvs_nt_col) or ""
            raw_pro = row.get(raw_hgvs_pro_col) or ""
            if raw_pro and not _is_blank(raw_pro):
                normalized = normalize_protein_hgvs(raw_pro)
                if normalized != raw_pro:
                    raw_pro = normalized
                    row[raw_hgvs_pro_col] = normalized

            if existing_results_by_key:
                merge_key = _build_merge_key(row, merge_key_columns)
                merged_result = existing_results_by_key.get(merge_key)
                if merged_result is not None:
                    if preserve_order == "groups":
                        _flush_current_group()
                    hgvs_c, hgvs_g, hgvs_p, err, clingen_allele_id = merged_result
                    _record_result(
                        idx,
                        row,
                        hgvs_c,
                        hgvs_g,
                        hgvs_p,
                        err,
                        clingen_allele_id,
                        group_id="__merged__",
                    )
                    n_merged += 1
                    continue

            case = _detect_case(raw_nt, raw_pro)

            # Process each contiguous class-1 block as soon as we encounter the
            # first non-class-1 row that follows it.
            if case != 1 and pending_case1_rows:
                _flush_case1_rows(pending_case1_rows)

            if preserve_order == "groups" and case not in (2, 3):
                _flush_current_group()

            if case is None:
                _record_result(idx, row, None, None, None, "No usable HGVS variant data in this row", group_id=None)
            elif case == 1:
                if idx < 10 or idx % PROGRESS_EVERY_ROWS == 0:
                    logger.debug("Row %d: case 1 (reference-based c./n. HGVS)", idx)
                pending_case1_rows.append((idx, row, raw_nt))
                if len(pending_case1_rows) >= max(2, max_clingen_concurrency):
                    _flush_case1_rows(pending_case1_rows)
            else:
                # Cases 2 and 3 require a target_sequence.
                target_seq = (row.get(target_sequence_col) or "").strip()
                if not target_seq:
                    _record_result(
                        idx,
                        row,
                        None,
                        None,
                        None,
                        (
                            f"target_sequence is required for case {case} but column "
                            f"{target_sequence_col!r} is missing or blank"
                        ),
                        group_id=None,
                    )
                else:
                    group_key = (row.get(group_by_col) or target_seq).strip() or target_seq
                    if idx < 10 or idx % PROGRESS_EVERY_ROWS == 0:
                        logger.debug("Row %d: case %d, group %r", idx, case, group_key)
                    if preserve_order == "groups":
                        if current_group_key is not None and group_key != current_group_key:
                            _flush_current_group()
                        current_group_key = group_key
                        current_group_rows.append((idx, row, raw_nt.strip(), raw_pro.strip(), case, target_seq))
                    else:
                        groups[group_key].append((idx, row, raw_nt.strip(), raw_pro.strip(), case, target_seq))

            if n_total % PROGRESS_EVERY_ROWS == 0:
                elapsed = max(time.monotonic() - started, 1e-9)
                logger.info(
                    "Scan progress: %d rows read, %d rows written, %d groups queued (%.1f rows/s).",
                    n_total,
                    n_written,
                    len(groups),
                    n_total / elapsed,
                )

        if n_total == 0:
            logger.warning("No rows to process after applying skip=%d, limit=%s.", skip, limit)
            return

        logger.info(
            "Processing rows %d-%d (selected=%d). Immediate output rows written: %d.",
            skip,
            skip + n_total - 1,
            n_total,
            n_written,
        )

        # Flush any class-1 rows queued at end-of-input.
        if pending_case1_rows:
            _flush_case1_rows(pending_case1_rows)

        if preserve_order == "groups":
            _flush_current_group()

        # Run dcd_mapping pipeline for all sequence-based groups (cases 2 and 3).
        if preserve_order != "groups" and groups:
            dcd = _try_import_dcd_mapping()
            n_groups = len(groups)
            loop = asyncio.new_event_loop()
            try:
                asyncio.set_event_loop(loop)
                for group_num, (group_name, group_rows) in enumerate(groups.items(), start=1):
                    target_seq = group_rows[0][5]
                    target_seqs = {r[5] for r in group_rows}
                    if len(target_seqs) > 1:
                        logger.warning(
                            "Group %r contains %d different target_sequence values; using the first one.",
                            group_name,
                            len(target_seqs),
                        )

                    group_rows_by_idx = {r[0]: r[1] for r in group_rows}
                    row_entries = [(r[0], r[2], r[3], r[4]) for r in group_rows]

                    # Try to process the full group; retry with chunking on BLAT error 137
                    per_row = None
                    transcript_nm = None
                    last_error = None

                    for attempt in range(1, dcd_max_retry_attempts + 1):
                        try:
                            per_row, transcript_nm = loop.run_until_complete(
                                _run_dcd_mapping_pipeline(group_name, target_seq, row_entries, dcd)
                            )
                            break  # Success; exit retry loop
                        except Exception as exc:
                            exc_str = str(exc)
                            last_error = exc
                            if dcd_chunk_on_137 and "error code 137" in exc_str and attempt < dcd_max_retry_attempts:
                                # Retry with chunking: half the remaining size each retry
                                chunk_sz = max(1, dcd_chunk_size_on_137 // (2 ** (attempt - 1)))
                                logger.warning(
                                    "Group %r failed with BLAT error 137 (attempt %d/%d); retrying with chunk_size=%d.",
                                    group_name,
                                    attempt,
                                    dcd_max_retry_attempts,
                                    chunk_sz,
                                )
                                # Process row_entries as chunks
                                per_row = []
                                for start_idx in range(0, len(row_entries), chunk_sz):
                                    chunk_entries = row_entries[start_idx : start_idx + chunk_sz]
                                    chunk_name = f"{group_name}#retry{attempt}_chunk{start_idx // chunk_sz + 1}"
                                    try:
                                        chunk_per_row, chunk_tx = loop.run_until_complete(
                                            _run_dcd_mapping_pipeline(chunk_name, target_seq, chunk_entries, dcd)
                                        )
                                        per_row.extend(chunk_per_row)
                                        if transcript_nm is None:
                                            transcript_nm = chunk_tx
                                    except Exception as chunk_exc:
                                        logger.error(
                                            "Chunk %s failed: %s",
                                            chunk_name,
                                            _format_exc(chunk_exc),
                                        )
                                        # Fail individual rows in this chunk
                                        for orig_idx, _, _, _ in chunk_entries:
                                            per_row.append((orig_idx, None, f"BLAT chunk failed: {_format_exc(chunk_exc)}"))
                                break  # Chunked attempt completed; exit retry loop
                            else:
                                # No retry, or non-137 error
                                raise

                    # If all retries failed, record errors for all rows
                    if per_row is None:
                        logger.error(
                            "Group %r failed after %d attempt(s): %s",
                            group_name,
                            dcd_max_retry_attempts,
                            _format_exc(last_error) if last_error else "unknown error",
                        )
                        fallback_error = _format_exc(last_error) if last_error is not None else "unknown error"
                        per_row = [(orig_idx, None, fallback_error) for orig_idx, _, _, _ in row_entries]

                    # Collect all HGVS strings that need ClinGen queries
                    hgvs_queries = []  # List of (orig_idx, row, hgvs_assay) tuples
                    for orig_idx, hgvs_assay, error in per_row:
                        row = group_rows_by_idx[orig_idx]

                        if error:
                            _record_result(orig_idx, row, None, None, None, error, group_id=group_name)
                            continue
                        if not hgvs_assay:
                            _record_result(orig_idx, row, None, None, None, "No HGVS string produced by mapper", group_id=group_name)
                            continue

                        hgvs_queries.append((orig_idx, row, hgvs_assay))

                    # Batch query ClinGen for all HGVS strings concurrently
                    if hgvs_queries:
                        hgvs_strings = [h[2] for h in hgvs_queries]
                        logger.debug(
                            "Batch querying ClinGen for %d variants in group %r (transcript: %s)",
                            len(hgvs_strings),
                            group_name,
                            transcript_nm,
                        )
                        clingen_batch_start = time.monotonic()
                        clingen_results = loop.run_until_complete(
                            _query_clingen_by_hgvs_batch(hgvs_strings, max_concurrency=max_clingen_concurrency)
                        )
                        clingen_batch_elapsed = time.monotonic() - clingen_batch_start
                        logger.debug(
                            "ClinGen batch query completed for %d variants in %.2f seconds",
                            len(hgvs_strings),
                            clingen_batch_elapsed,
                        )

                        # Process results in the same order as queries
                        for orig_idx, row, hgvs_assay in hgvs_queries:
                            data = clingen_results.get(hgvs_assay)
                            if data is None:
                                # Preserve the assay-level HGVS even when ClinGen has no record.
                                # Protein assays (p. strings) go to hgvs_p; genomic/transcript
                                # assays go to hgvs_g so the mapping work is not discarded.
                                _assay_is_protein = (
                                    hgvs_assay.startswith("p.") or ":p." in hgvs_assay
                                )
                                _record_result(
                                    orig_idx,
                                    row,
                                    None,
                                    None if _assay_is_protein else hgvs_assay,
                                    hgvs_assay if _assay_is_protein else None,
                                    f"ClinGen returned no data for {hgvs_assay!r}",
                                    group_id=group_name,
                                )
                                continue

                            hgvs_g, hgvs_c, hgvs_p = _extract_hgvs_from_clingen(data, transcript_nm)

                            # For protein-layer variants where ClinGen returns a PA allele, the
                            # assay-level p. string is the best available protein HGVS if ClinGen
                            # did not populate hgvs_p from aminoAcidAlleles.
                            if _clingen_allele_type(data) == "PA" and hgvs_p is None:
                                hgvs_p = hgvs_assay

                            _record_result(
                                orig_idx,
                                row,
                                hgvs_c,
                                hgvs_g,
                                hgvs_p,
                                None,
                                _extract_clingen_allele_id(data),
                                group_id=group_name,
                            )

                            if n_written % PROGRESS_EVERY_ROWS == 0:
                                elapsed = max(time.monotonic() - started, 1e-9)
                                logger.info(
                                    "Output progress: %d/%d rows written (%d errors, %.1f rows/s).",
                                    n_written,
                                    n_total,
                                    n_errors,
                                    n_written / elapsed,
                                )

                    elapsed = max(time.monotonic() - started, 1e-9)
                    logger.info(
                        "Group %d/%d complete (%d rows in group). %d/%d rows written so far (%.1f rows/s).",
                        group_num,
                        n_groups,
                        len(group_rows),
                        n_written,
                        n_total,
                        n_written / elapsed,
                    )
            finally:
                asyncio.set_event_loop(None)
                loop.close()

        if sequence_loop is not None:
            sequence_loop.close()

        # Safety net: emit any buffered rows that could not be flushed during normal flow.
        if preserve_order in {"index", "groups"} and pending_results:
            logger.warning(
                "Preserve-order flush: %d buffered rows remained at end of run.",
                len(pending_results),
            )
            for idx in sorted(pending_results):
                row = rows_by_idx[idx]
                (
                    pending_hgvs_c,
                    pending_hgvs_g,
                    pending_hgvs_p,
                    pending_error,
                    pending_clingen_allele_id,
                ) = pending_results[idx]
                _write_result_row(
                    writer,
                    out_fh,
                    row,
                    pending_hgvs_c,
                    pending_hgvs_g,
                    pending_hgvs_p,
                    pending_error,
                    pending_clingen_allele_id,
                )
                n_written += 1
                if pending_error:
                    n_errors += 1
            pending_results.clear()
            rows_by_idx.clear()

        # Ensure all buffered output is visible to downstream tail/readers.
        out_fh.flush()
        try:
            os.fsync(out_fh.fileno())
        except OSError:
            # Some filesystems may not support fsync; regular flush already happened.
            pass

    logger.info(
        "Done. %d/%d rows had errors; reused %d rows from merge sources. Output written to %s.",
        n_errors,
        n_total,
        n_merged,
        output_file,
    )
    print(
        f"Done. {n_errors}/{n_total} rows had errors; reused {n_merged} rows from merge sources. "
        f"Output written to {output_file}",
        file=sys.stderr,
    )
    sys.stderr.flush()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.command()
@click.argument("input_file")
@click.argument("output_file")
@click.option(
    "--raw-hgvs-nt",
    "raw_hgvs_nt_col",
    default="raw_hgvs_nt",
    show_default=True,
    help="Column name for raw nucleotide HGVS strings.",
)
@click.option(
    "--raw-hgvs-pro",
    "raw_hgvs_pro_col",
    default="raw_hgvs_pro",
    show_default=True,
    help="Column name for raw protein HGVS strings.",
)
@click.option(
    "--target-sequence",
    "target_sequence_col",
    default="target_sequence",
    show_default=True,
    help="Column name for target sequences (required for cases 2 and 3).",
)
@click.option(
    "--target-type",
    "target_type_col",
    default="target_type",
    show_default=True,
    help="Column distinguishing 'reference' from 'target-sequence' rows (informational).",
)
@click.option(
    "--group-by",
    "group_by_col",
    default="target_sequence",
    show_default=True,
    help="Column used to batch sequence-based rows for a single BLAT alignment run.",
)
@click.option(
    "--mapped-hgvs-g",
    "mapped_hgvs_g_col",
    default="mapped_hgvs_g",
    show_default=True,
    help="Output column name for genomic (g.) HGVS strings.",
)
@click.option(
    "--mapped-hgvs-c",
    "mapped_hgvs_c_col",
    default="mapped_hgvs_c",
    show_default=True,
    help="Output column name for transcript (c.) HGVS strings.",
)
@click.option(
    "--mapped-hgvs-p",
    "mapped_hgvs_p_col",
    default="mapped_hgvs_p",
    show_default=True,
    help="Output column name for protein (p.) HGVS strings.",
)
@click.option(
    "--mapping-error",
    "mapping_error_col",
    default="mapping_error",
    show_default=True,
    help="Output column name for error messages.",
)
@click.option(
    "--clingen-allele-id",
    "clingen_allele_id_col",
    default="clingen_allele_id",
    show_default=True,
    help="Output column name for ClinGen Allele Registry allele IDs.",
)
@click.option(
    "--skip",
    default=0,
    show_default=True,
    type=int,
    help="Number of data rows to skip from the start of the file (useful for resuming).",
)
@click.option(
    "--limit",
    default=None,
    type=int,
    help="Maximum number of data rows to process. Processes all remaining rows when omitted.",
)
@click.option(
    "--drop-columns",
    "drop_columns",
    multiple=True,
    metavar="COLUMN",
    help="Column to omit from the output. May be repeated to drop multiple columns.",
)
@click.option(
    "--max-clingen-concurrency",
    "max_clingen_concurrency",
    default=5,
    show_default=True,
    type=int,
    help="Maximum number of concurrent ClinGen Allele Registry API requests.",
)
@click.option(
    "--dcd-chunk-on-137/--no-dcd-chunk-on-137",
    "dcd_chunk_on_137",
    default=True,
    show_default=True,
    help="Automatically retry groups with BLAT error 137 using chunked processing.",
)
@click.option(
    "--dcd-chunk-size-on-137",
    "dcd_chunk_size_on_137",
    default=500,
    show_default=True,
    type=int,
    help="Initial chunk size for BLAT 137 retry.",
)
@click.option(
    "--dcd-max-retry-attempts",
    "dcd_max_retry_attempts",
    default=3,
    show_default=True,
    type=int,
    help="Maximum number of retry attempts with progressively smaller chunks.",
)
@click.option(
    "--preserve-order",
    type=click.Choice(["no", "index", "groups"], case_sensitive=False),
    default="groups",
    show_default=True,
    help=(
        "Order guarantee mode: "
        "'no' writes immediately (potentially out of order); "
        "'index' buffers by row index to guarantee exact input order; "
        "'groups' preserves input order while processing groups (recommended default)."
    ),
)
@click.option(
    "--merge-existing",
    "merge_existing_files",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=str),
    metavar="FILE",
    help=(
        "Existing annotated CSV/TSV file to reuse mapped results from. "
        "May be repeated; first match by merge key wins."
    ),
)
@click.option(
    "--merge-match-col",
    "merge_match_columns",
    multiple=True,
    metavar="COLUMN",
    help=(
        "Additional column that must match for merge reuse. "
        "Rows are merged only when raw_hgvs_nt, raw_hgvs_pro, and all specified columns match."
    ),
)
@click.option(
    "--targets-file",
    "targets_file",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=str),
    help=(
        "Optional CSV/TSV file containing target sequences and other target-level columns. "
        "Joined to the input file on the column specified by --target-name."
    ),
)
@click.option(
    "--target-name",
    "target_name_col",
    default="target_name",
    show_default=True,
    help="Column name used to join the input file to --targets-file.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Enable DEBUG-level logging.",
)
def main(
    input_file: str,
    output_file: str,
    raw_hgvs_nt_col: str,
    raw_hgvs_pro_col: str,
    target_sequence_col: str,
    target_type_col: str,
    group_by_col: str,
    mapped_hgvs_g_col: str,
    mapped_hgvs_c_col: str,
    mapped_hgvs_p_col: str,
    mapping_error_col: str,
    clingen_allele_id_col: str,
    skip: int,
    limit: Optional[int],
    drop_columns: tuple[str, ...],
    max_clingen_concurrency: int,
    dcd_chunk_on_137: bool,
    dcd_chunk_size_on_137: int,
    dcd_max_retry_attempts: int,
    preserve_order: str,
    merge_existing_files: tuple[str, ...],
    merge_match_columns: tuple[str, ...],
    targets_file: Optional[str],
    target_name_col: str,
    verbose: bool,
) -> None:
    """Map variants in INPUT_FILE to human-genome HGVS strings and write OUTPUT_FILE.

    Adds mapped_hgvs_g, mapped_hgvs_c, mapped_hgvs_p, mapping_error, and
    clingen_allele_id columns.
    Three variant categories are handled automatically:

    \b
    1. Reference-based c. HGVS with transcript accession → ClinGen query directly.
    2. Bare c. HGVS + target_sequence → dcd_mapping (BLAT, SeqRepo, UTA) + ClinGen.
    3. Protein-only HGVS + target_sequence → dcd_mapping (protein level) + ClinGen.
    """
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stderr,
    )
    map_variants(
        input_file=input_file,
        output_file=output_file,
        raw_hgvs_nt_col=raw_hgvs_nt_col,
        raw_hgvs_pro_col=raw_hgvs_pro_col,
        target_sequence_col=target_sequence_col,
        target_type_col=target_type_col,
        group_by_col=group_by_col,
        mapped_hgvs_g_col=mapped_hgvs_g_col,
        mapped_hgvs_c_col=mapped_hgvs_c_col,
        mapped_hgvs_p_col=mapped_hgvs_p_col,
        mapping_error_col=mapping_error_col,
        clingen_allele_id_col=clingen_allele_id_col,
        skip=skip,
        limit=limit,
        drop_columns=drop_columns,
        max_clingen_concurrency=max_clingen_concurrency,
        dcd_chunk_on_137=dcd_chunk_on_137,
        dcd_chunk_size_on_137=dcd_chunk_size_on_137,
        dcd_max_retry_attempts=dcd_max_retry_attempts,
        preserve_order=preserve_order,
        merge_existing_files=merge_existing_files,
        merge_match_columns=merge_match_columns,
        targets_file=targets_file,
        target_name_col=target_name_col,
    )


if __name__ == "__main__":
    main()
