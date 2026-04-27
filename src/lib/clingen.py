"""ClinGen Allele Registry client helpers with optional Redis caching."""

from __future__ import annotations

import json
import logging
import os
import time
import urllib.parse
from typing import Any, Optional

import requests

CLINGEN_API_URL = os.environ.get("CLINGEN_API_URL", "https://reg.genome.network/allele")
CLINGEN_MAX_RETRIES = 3
CLINGEN_RETRY_DELAY = 2.0

CLINGEN_CACHE_REDIS_URL_DEFAULT = "redis://redis:6379/0"
CLINGEN_CACHE_PREFIX_DEFAULT = "clingen:v1"
CLINGEN_CACHE_TTL_SECONDS_DEFAULT = 86400
CLINGEN_CACHE_MISS_TTL_SECONDS_DEFAULT = 86400

_CACHE_MISS_SENTINEL = "__MISS__"
_REDIS_CLIENT: Any = None
_REDIS_INIT_ATTEMPTED = False
_REDIS_UNAVAILABLE_LOGGED = False

logger = logging.getLogger(__name__)


def _env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off", ""}


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        parsed = int(raw)
    except ValueError:
        logger.warning("Invalid integer for %s=%r; using default=%d", name, raw, default)
        return default
    return parsed


def _cache_enabled() -> bool:
    return _env_bool("CLINGEN_CACHE_ENABLED", True)


def _cache_prefix() -> str:
    return (os.environ.get("CLINGEN_CACHE_PREFIX") or CLINGEN_CACHE_PREFIX_DEFAULT).strip()


def _cache_ttl_seconds() -> int:
    return max(1, _env_int("CLINGEN_CACHE_TTL_SECONDS", CLINGEN_CACHE_TTL_SECONDS_DEFAULT))


def _cache_miss_ttl_seconds() -> int:
    return max(
        1,
        _env_int("CLINGEN_CACHE_MISS_TTL_SECONDS", CLINGEN_CACHE_MISS_TTL_SECONDS_DEFAULT),
    )


def _cache_redis_url() -> str:
    return (
        os.environ.get("CLINGEN_CACHE_REDIS_URL")
        or os.environ.get("REDIS_URL")
        or CLINGEN_CACHE_REDIS_URL_DEFAULT
    )


def _allele_cache_key(allele_id: str) -> str:
    return f"{_cache_prefix()}:allele:{allele_id}"


def _hgvs_map_key(hgvs: str) -> str:
    return f"{_cache_prefix()}:hgvs:{hgvs}"


def _get_redis_client(*, force: bool = False):
    global _REDIS_CLIENT
    global _REDIS_INIT_ATTEMPTED
    global _REDIS_UNAVAILABLE_LOGGED

    if not force and not _cache_enabled():
        return None

    if _REDIS_CLIENT is not None:
        return _REDIS_CLIENT
    if _REDIS_INIT_ATTEMPTED:
        return None

    _REDIS_INIT_ATTEMPTED = True
    try:
        import redis  # type: ignore[import-not-found]

        client = redis.Redis.from_url(_cache_redis_url(), decode_responses=True)
        client.ping()
        _REDIS_CLIENT = client
        return _REDIS_CLIENT
    except Exception as exc:
        if not _REDIS_UNAVAILABLE_LOGGED:
            logger.warning("ClinGen Redis cache unavailable; continuing without cache: %s", exc)
            _REDIS_UNAVAILABLE_LOGGED = True
        return None


def _cache_get(key: str) -> tuple[bool, Optional[str]]:
    client = _get_redis_client()
    if client is None:
        return False, None
    try:
        value = client.get(key)
    except Exception:
        return False, None
    return (value is not None), value


def _cache_set(key: str, value: str, *, miss: bool = False) -> None:
    client = _get_redis_client()
    if client is None:
        return
    ttl = _cache_miss_ttl_seconds() if miss else _cache_ttl_seconds()
    try:
        client.set(key, value, ex=ttl)
    except Exception:
        return


def _extract_clingen_allele_id(data: dict) -> Optional[str]:
    def _normalize(value: str) -> str:
        text = (value or "").strip()
        if not text or text.startswith("_:"):
            return ""
        return text

    at_id: str = data.get("@id", "") or ""
    if at_id:
        fragment = at_id.rstrip("/").rsplit("/", 1)[-1]
        normalized = _normalize(fragment)
        return normalized or None

    fallback = data.get("id")
    if isinstance(fallback, str):
        normalized = _normalize(fallback)
        return normalized or None
    return None


def _extract_clinvar_allele_id(data: dict) -> str:
    allele_id = (
        data.get("externalRecords", {})
        .get("ClinVarAlleles", [{}])[0]
        .get("alleleId")
    )
    return str(allele_id) if allele_id is not None else ""


def _extract_clinvar_variation_id(data: dict) -> str:
    variation_id = (
        data.get("externalRecords", {})
        .get("ClinVarVariations", [{}])[0]
        .get("variationId")
    )
    return str(variation_id) if variation_id is not None else ""


def resolve_clinvar_ids(
    clingen_id: str,
    cache: dict[str, tuple[str, str]],
    max_retries: int = CLINGEN_MAX_RETRIES,
    retry_delay: float = CLINGEN_RETRY_DELAY,
) -> tuple[str, str]:
    """Return (ClinVar variation ID, ClinVar allele ID) for a ClinGen ID.

    Returns ``("", "")`` when unavailable. Results are cached in *cache* for
    the duration of the process. ClinGen allele responses are also cached in
    Redis when enabled.
    """
    if clingen_id in cache:
        return cache[clingen_id]

    key = _allele_cache_key(clingen_id)
    found, raw = _cache_get(key)
    if found and raw is not None:
        if raw == _CACHE_MISS_SENTINEL:
            cache[clingen_id] = ("", "")
            return cache[clingen_id]
        try:
            cached_data = json.loads(raw)
            if isinstance(cached_data, dict):
                resolved = (
                    _extract_clinvar_variation_id(cached_data),
                    _extract_clinvar_allele_id(cached_data),
                )
                cache[clingen_id] = resolved
                return resolved
        except Exception:
            pass

    data = _fetch_allele_response_by_id(
        clingen_id,
        max_retries=max_retries,
        retry_delay=retry_delay,
    )
    if data is None:
        cache[clingen_id] = ("", "")
        return cache[clingen_id]

    result = (
        _extract_clinvar_variation_id(data),
        _extract_clinvar_allele_id(data),
    )
    cache[clingen_id] = result
    return result


def _fetch_allele_response_by_id(
    clingen_id: str,
    *,
    max_retries: int,
    retry_delay: float,
) -> Optional[dict]:
    url = f"{CLINGEN_API_URL}/{urllib.parse.quote(clingen_id, safe='')}"

    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 404:
                _cache_set(_allele_cache_key(clingen_id), _CACHE_MISS_SENTINEL, miss=True)
                return None
            response.raise_for_status()
            data = response.json()
            _cache_set(_allele_cache_key(clingen_id), json.dumps(data))
            return data
        except (requests.HTTPError, requests.ConnectionError, requests.Timeout) as exc:
            if isinstance(exc, requests.HTTPError) and exc.response is not None:
                status = exc.response.status_code
                if status != 429 and 400 <= status < 500:
                    logger.warning(
                        "ClinGen API returned %s for %s; skipping", status, clingen_id
                    )
                    _cache_set(_allele_cache_key(clingen_id), _CACHE_MISS_SENTINEL, miss=True)
                    return None
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

    return None


def resolve_clinvar_allele_id(
    clingen_id: str,
    cache: dict[str, str],
    max_retries: int = CLINGEN_MAX_RETRIES,
    retry_delay: float = CLINGEN_RETRY_DELAY,
) -> str:
    """Return ClinVar allele ID for a ClinGen ID, or empty string if unavailable.

    Results are cached in *cache* for the duration of the process. ClinGen allele
    responses are also cached in Redis when enabled.
    """
    if clingen_id in cache:
        return cache[clingen_id]

    tuple_cache: dict[str, tuple[str, str]] = {}
    result = resolve_clinvar_ids(
        clingen_id,
        tuple_cache,
        max_retries=max_retries,
        retry_delay=retry_delay,
    )
    cache[clingen_id] = result[1]
    return result[1]


def query_clingen_by_hgvs(
    hgvs_string: str,
    max_retries: int = CLINGEN_MAX_RETRIES,
    retry_delay: float = CLINGEN_RETRY_DELAY,
    *,
    log_404: bool = False,
) -> Optional[dict]:
    """Query ClinGen Allele Registry by HGVS string.

    HGVS lookups are cached in Redis as HGVS->allele_id mappings and allele
    responses are cached by allele ID.
    """
    hgvs = (hgvs_string or "").strip()
    if not hgvs:
        return None

    map_key = _hgvs_map_key(hgvs)
    found_map, cached_allele_id = _cache_get(map_key)
    if found_map and cached_allele_id is not None:
        if cached_allele_id == _CACHE_MISS_SENTINEL:
            return None
        allele_key = _allele_cache_key(cached_allele_id)
        found_allele, cached_response = _cache_get(allele_key)
        if found_allele and cached_response is not None:
            if cached_response == _CACHE_MISS_SENTINEL:
                return None
            try:
                payload = json.loads(cached_response)
                if isinstance(payload, dict):
                    return payload
            except Exception:
                pass

        # Mapping exists but allele response was evicted/corrupt; refetch by allele id.
        return _fetch_allele_response_by_id(
            cached_allele_id,
            max_retries=max_retries,
            retry_delay=retry_delay,
        )

    for attempt in range(max_retries):
        try:
            resp = requests.get(
                CLINGEN_API_URL,
                params={"hgvs": hgvs},
                timeout=30,
                headers={"Accept": "application/json"},
            )
            if resp.status_code == 200:
                data = resp.json()
                allele_id = _extract_clingen_allele_id(data)
                if allele_id:
                    _cache_set(_allele_cache_key(allele_id), json.dumps(data))
                    _cache_set(map_key, allele_id)
                return data
            if resp.status_code == 404:
                if log_404:
                    logger.warning("ClinGen 404 for %s", hgvs)
                _cache_set(map_key, _CACHE_MISS_SENTINEL, miss=True)
                return None
            if resp.status_code == 429:
                wait = retry_delay * (2**attempt)
                logger.warning("ClinGen rate-limited for %s; waiting %.1f s", hgvs, wait)
                time.sleep(wait)
                continue
            logger.warning("ClinGen returned HTTP %d for %s", resp.status_code, hgvs)
            return None
        except requests.exceptions.RequestException as exc:
            logger.warning(
                "ClinGen request failed for %s (attempt %d/%d): %s",
                hgvs,
                attempt + 1,
                max_retries,
                exc,
            )
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
    return None


def clear_clingen_cache(prefix: Optional[str] = None) -> int:
    """Delete all ClinGen cache keys for the given prefix. Returns deleted count."""
    client = _get_redis_client(force=True)
    if client is None:
        raise RuntimeError("Redis client unavailable; cannot clear ClinGen cache")

    key_prefix = (prefix or _cache_prefix()).strip()
    pattern = f"{key_prefix}:*"
    cursor = 0
    deleted = 0

    while True:
        cursor, keys = client.scan(cursor=cursor, match=pattern, count=500)
        if keys:
            deleted += int(client.delete(*keys))
        if cursor == 0:
            break

    return deleted
