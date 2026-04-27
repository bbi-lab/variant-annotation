"""ClinGen Allele Registry client helpers."""

from __future__ import annotations

import logging
import time
import urllib.parse
from typing import Optional

import requests

CLINGEN_API_URL = "https://reg.genome.network/allele"
CLINGEN_MAX_RETRIES = 3
CLINGEN_RETRY_DELAY = 2.0

logger = logging.getLogger(__name__)


def resolve_clinvar_allele_id(
    clingen_id: str,
    cache: dict[str, str],
    max_retries: int = CLINGEN_MAX_RETRIES,
    retry_delay: float = CLINGEN_RETRY_DELAY,
) -> str:
    """Return ClinVar allele ID for a ClinGen ID, or empty string if unavailable.

    Results are cached in *cache* for the duration of the process.
    """
    if clingen_id in cache:
        return cache[clingen_id]

    url = f"{CLINGEN_API_URL}/{urllib.parse.quote(clingen_id, safe='')}"
    last_exc: Optional[Exception] = None

    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 404:
                cache[clingen_id] = ""
                return ""
            response.raise_for_status()
            api_data = response.json()
            allele_id = (
                api_data.get("externalRecords", {})
                .get("ClinVarAlleles", [{}])[0]
                .get("alleleId")
            )
            result = str(allele_id) if allele_id is not None else ""
            cache[clingen_id] = result
            return result
        except (requests.HTTPError, requests.ConnectionError, requests.Timeout) as exc:
            last_exc = exc
            if isinstance(exc, requests.HTTPError) and exc.response is not None:
                status = exc.response.status_code
                # Don't retry client errors other than rate-limiting
                if status != 429 and 400 <= status < 500:
                    logger.warning(
                        "ClinGen API returned %s for %s; skipping", status, clingen_id
                    )
                    cache[clingen_id] = ""
                    return ""
            if attempt < max_retries - 1:
                logger.warning(
                    "ClinGen API error for %s (attempt %d/%d): %s; retrying in %.1fs",
                    clingen_id,
                    attempt + 1,
                    max_retries,
                    exc,
                    retry_delay,
                )
                time.sleep(retry_delay)

    logger.error(
        "ClinGen API failed for %s after %d attempts: %s",
        clingen_id,
        max_retries,
        last_exc,
    )
    cache[clingen_id] = ""
    return ""


def query_clingen_by_hgvs(
    hgvs_string: str,
    max_retries: int = CLINGEN_MAX_RETRIES,
    retry_delay: float = CLINGEN_RETRY_DELAY,
    *,
    log_404: bool = False,
) -> Optional[dict]:
    """Query ClinGen Allele Registry by HGVS string.

    Returns parsed JSON response on success, or None when lookup fails.
    """
    for attempt in range(max_retries):
        try:
            resp = requests.get(
                CLINGEN_API_URL,
                params={"hgvs": hgvs_string},
                timeout=30,
                headers={"Accept": "application/json"},
            )
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 404:
                if log_404:
                    logger.warning("ClinGen 404 for %s", hgvs_string)
                return None
            if resp.status_code == 429:
                wait = retry_delay * (2**attempt)
                logger.warning(
                    "ClinGen rate-limited for %s; waiting %.1f s", hgvs_string, wait
                )
                time.sleep(wait)
                continue
            logger.warning("ClinGen returned HTTP %d for %s", resp.status_code, hgvs_string)
            return None
        except requests.exceptions.RequestException as exc:
            logger.warning(
                "ClinGen request failed for %s (attempt %d/%d): %s",
                hgvs_string,
                attempt + 1,
                max_retries,
                exc,
            )
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
    return None
