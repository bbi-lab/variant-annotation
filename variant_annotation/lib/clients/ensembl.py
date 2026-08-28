"""Ensembl REST-backed VepLookup / VariantRecoder implementations."""

from __future__ import annotations

import logging
import os
import time
from typing import Any, Mapping, Sequence

import requests  # type: ignore[import-untyped]

logger = logging.getLogger(__name__)

ENSEMBL_API_URL_DEFAULT = os.environ.get("ENSEMBL_API_URL", "https://rest.ensembl.org")

#: Ensembl's documented ceiling for POST /vep/human/hgvs.
VEP_MAX_HGVS_PER_REQUEST = 200

_JSON_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


class EnsemblRestClient:
    """Satisfies VepLookup and VariantRecoder against the Ensembl REST API.

    One client per job, reused across batches: it holds a ``requests.Session`` so TLS handshakes and
    connection setup are paid once rather than per batch.

    Retries are applied to the transport failures worth retrying — timeouts, connection resets, and the
    5xx/429 responses a rate-limited public service emits under load. A 4xx is *not* retried: VEP
    returning 400 for a protein HGVS it cannot parse is a settled answer about the input, and retrying
    it burns quota to receive the same rejection. Callers distinguish the two through the port contract:
    a raise means "unknown, retry"; the resolution layer converts an exhausted-retry raise into an
    ``ERRORED`` outcome.
    """

    def __init__(
        self,
        *,
        api_url: str = ENSEMBL_API_URL_DEFAULT,
        timeout_seconds: int = 60,
        recoder_timeout_seconds: int = 600,
        max_attempts: int = 4,
        backoff_seconds: float = 1.0,
        session: Any = None,
    ) -> None:
        self._api_url = api_url.rstrip("/")
        self._timeout = timeout_seconds
        # Variant Recoder is markedly slower than VEP for large batches and 504s are routine, so it gets
        # its own, far more generous budget rather than sharing VEP's.
        self._recoder_timeout = recoder_timeout_seconds
        self._max_attempts = max_attempts
        self._backoff = backoff_seconds
        self._session = session if session is not None else requests.Session()

    # --- VepLookup ---

    def annotate_hgvs(self, hgvs: Sequence[str], *, refseq: bool) -> Sequence[Mapping[str, Any]]:
        """POST a batch of HGVS notations to ``/vep/human/hgvs``.

        ``refseq`` adds Ensembl's ``refseq=1``, which switches ``transcript_consequences`` to the RefSeq
        transcript set so ``NM_``/``NR_`` accessions can be matched back.
        """
        notations = list(hgvs)
        if len(notations) > VEP_MAX_HGVS_PER_REQUEST:
            raise ValueError(
                f"VEP accepts at most {VEP_MAX_HGVS_PER_REQUEST} HGVS notations per request, got "
                f"{len(notations)}. Batching is the caller's responsibility."
            )

        body: dict[str, Any] = {"hgvs_notations": notations}
        if refseq:
            body["refseq"] = 1

        payload = self._post_json("/vep/human/hgvs", body, timeout=self._timeout)
        return payload if isinstance(payload, list) else []

    # --- VariantRecoder ---

    def to_genomic(self, hgvs: Sequence[str]) -> Mapping[str, Sequence[str]]:
        """POST to ``/variant_recoder/human`` and collect each input's genomic (``NC_``) equivalents.

        Only ``NC_`` expressions are kept: the fallback exists to obtain a genomic form VEP can resolve
        against a reference assembly, and Recoder also returns transcript- and protein-level spellings
        that would just re-pose the question VEP already declined.
        """
        payload = self._post_json("/variant_recoder/human", {"ids": list(hgvs)}, timeout=self._recoder_timeout)
        if not isinstance(payload, list):
            return {}

        recoded: dict[str, list[str]] = {}
        for record in payload:
            if not isinstance(record, Mapping):
                continue
            # Recoder nests results under arbitrary allele keys, each echoing the submitted string in
            # "input"; a top-level "input" key is metadata, not an allele. A key's value is either a
            # single result object or a list of them, so both shapes are flattened.
            for key, value in record.items():
                if key == "input":
                    continue
                alleles = value if isinstance(value, list) else [value]
                for allele in alleles:
                    if not isinstance(allele, Mapping):
                        continue
                    submitted = allele.get("input")
                    if not submitted:
                        continue
                    genomic = [
                        str(expression)
                        for expression in (allele.get("hgvsg") or ())
                        if str(expression).startswith("NC_")
                    ]
                    if genomic:
                        recoded.setdefault(str(submitted), []).extend(genomic)
        return recoded

    # --- release metadata ---

    def software_release(self) -> str:
        """Return the Ensembl release the REST API is serving, e.g. ``"116"`` (``/info/software``).

        An Ensembl release is coordinated — software, transcript set, and consequence vocabulary all
        bump together under one number — so it version-keys every consequence resolved against it.
        """
        payload = self._get_json("/info/software", timeout=30)
        return str(payload["release"])

    # --- transport ---

    def _post_json(self, path: str, body: Mapping[str, Any], *, timeout: int) -> Any:
        return self._request_json("POST", path, timeout=timeout, json=body)

    def _get_json(self, path: str, *, timeout: int) -> Any:
        return self._request_json("GET", path, timeout=timeout)

    def _request_json(self, method: str, path: str, *, timeout: int, **kwargs: Any) -> Any:
        """Issue a request with bounded retries, raising on final failure.

        Raising rather than returning a sentinel is the port contract: the resolution layer needs
        "unknown" to be distinguishable from "answered with nothing", and only an exception carries that.
        """
        url = f"{self._api_url}{path}"
        last_error: Exception | None = None

        for attempt in range(1, self._max_attempts + 1):
            try:
                response = self._session.request(method, url, headers=_JSON_HEADERS, timeout=timeout, **kwargs)
            except requests.RequestException as exc:
                last_error = exc
            else:
                if response.status_code < 400:
                    return response.json()
                # 4xx other than 429 is a statement about the request itself; retrying cannot change it.
                if response.status_code != 429 and response.status_code < 500:
                    raise requests.HTTPError(
                        f"{method} {url} returned {response.status_code}: {response.text[:500]}",
                        response=response,
                    )
                last_error = requests.HTTPError(f"{method} {url} returned {response.status_code}", response=response)
                # Honour Retry-After when the service tells us how long to wait.
                retry_after = response.headers.get("Retry-After")
                if retry_after and retry_after.strip().isdigit():
                    time.sleep(min(float(retry_after), 60.0))
                    continue

            if attempt < self._max_attempts:
                time.sleep(self._backoff * (2 ** (attempt - 1)))

        assert last_error is not None  # the loop either returns, raises, or records an error
        raise last_error
