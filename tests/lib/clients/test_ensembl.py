"""Unit tests for the Ensembl REST client — request shaping and response parsing, no network."""

import pytest
import requests

from variant_annotation.lib.clients.ensembl import VEP_MAX_HGVS_PER_REQUEST, EnsemblRestClient

pytestmark = pytest.mark.unit


class FakeResponse:
    def __init__(self, payload=None, *, status_code=200, text="", headers=None):
        self._payload = payload
        self.status_code = status_code
        self.text = text
        self.headers = headers or {}

    def json(self):
        return self._payload


class FakeSession:
    """Records requests and replays a scripted sequence of responses (or exceptions)."""

    def __init__(self, *responses):
        self._responses = list(responses)
        self.calls = []

    def request(self, method, url, **kwargs):
        self.calls.append((method, url, kwargs))
        response = self._responses.pop(0) if self._responses else FakeResponse([])
        if isinstance(response, Exception):
            raise response
        return response


def client(*responses, **kwargs):
    return EnsemblRestClient(session=FakeSession(*responses), backoff_seconds=0, **kwargs)


# --- annotate_hgvs ---


def test_refseq_flag_is_sent_when_requested():
    ensembl = client(FakeResponse([]))
    ensembl.annotate_hgvs(["NM_1.1:c.1A>G"], refseq=True)

    _, _, kwargs = ensembl._session.calls[0]
    assert kwargs["json"] == {"hgvs_notations": ["NM_1.1:c.1A>G"], "refseq": 1}


def test_refseq_flag_is_omitted_when_not_requested():
    ensembl = client(FakeResponse([]))
    ensembl.annotate_hgvs(["ENST1.1:c.1A>G"], refseq=False)

    _, _, kwargs = ensembl._session.calls[0]
    assert kwargs["json"] == {"hgvs_notations": ["ENST1.1:c.1A>G"]}
    assert "refseq" not in kwargs["json"]


def test_a_batch_over_ensembls_limit_is_refused_rather_than_silently_truncated():
    ensembl = client()
    with pytest.raises(ValueError, match="at most 200"):
        ensembl.annotate_hgvs([f"NM_1.1:c.{n}A>G" for n in range(VEP_MAX_HGVS_PER_REQUEST + 1)], refseq=True)


def test_a_non_list_payload_yields_no_entries():
    ensembl = client(FakeResponse({"error": "something"}))
    assert ensembl.annotate_hgvs(["NM_1.1:c.1A>G"], refseq=False) == []


# --- retry policy ---


def test_a_transport_error_is_retried_then_raised():
    ensembl = client(
        requests.ConnectionError("reset"),
        requests.ConnectionError("reset"),
        requests.ConnectionError("reset"),
        requests.ConnectionError("reset"),
        max_attempts=4,
    )
    with pytest.raises(requests.ConnectionError):
        ensembl.annotate_hgvs(["NM_1.1:c.1A>G"], refseq=False)

    assert len(ensembl._session.calls) == 4


def test_a_transient_failure_followed_by_success_returns_the_success():
    ensembl = client(requests.ConnectionError("reset"), FakeResponse([{"input": "NM_1.1:c.1A>G"}]))
    entries = ensembl.annotate_hgvs(["NM_1.1:c.1A>G"], refseq=False)

    assert entries == [{"input": "NM_1.1:c.1A>G"}]
    assert len(ensembl._session.calls) == 2


def test_a_server_error_is_retried():
    ensembl = client(FakeResponse(status_code=503), FakeResponse([]), max_attempts=3)
    ensembl.annotate_hgvs(["NM_1.1:c.1A>G"], refseq=False)
    assert len(ensembl._session.calls) == 2


def test_a_client_error_is_not_retried():
    """VEP rejecting a protein HGVS is a settled answer; retrying burns quota for the same rejection."""
    ensembl = client(FakeResponse(status_code=400, text="cannot parse"), max_attempts=4)

    with pytest.raises(requests.HTTPError, match="400"):
        ensembl.annotate_hgvs(["NP_1.1:p.Gly15Asn"], refseq=False)

    assert len(ensembl._session.calls) == 1


def test_rate_limiting_is_retried():
    ensembl = client(FakeResponse(status_code=429, headers={"Retry-After": "0"}), FakeResponse([]))
    ensembl.annotate_hgvs(["NM_1.1:c.1A>G"], refseq=False)
    assert len(ensembl._session.calls) == 2


# --- to_genomic ---


def test_only_genomic_expressions_are_kept():
    """Transcript and protein spellings would just re-pose the question VEP already declined."""
    ensembl = client(
        FakeResponse(
            [
                {
                    "input": "NP_1.1:p.Gly15Asn",
                    "A": {
                        "input": "NP_1.1:p.Gly15Asn",
                        "hgvsg": ["NC_000001.11:g.5G>A", "NM_1.1:c.43G>A", "ENSP1:p.Gly15Asn"],
                    },
                }
            ]
        )
    )

    assert ensembl.to_genomic(["NP_1.1:p.Gly15Asn"]) == {"NP_1.1:p.Gly15Asn": ["NC_000001.11:g.5G>A"]}


def test_multiple_allele_keys_accumulate_for_one_input():
    ensembl = client(
        FakeResponse(
            [
                {
                    "A": {"input": "NP_1.1:p.Gly15Asn", "hgvsg": ["NC_1:g.5G>A"]},
                    "C": {"input": "NP_1.1:p.Gly15Asn", "hgvsg": ["NC_1:g.6G>C"]},
                }
            ]
        )
    )

    assert ensembl.to_genomic(["NP_1.1:p.Gly15Asn"]) == {"NP_1.1:p.Gly15Asn": ["NC_1:g.5G>A", "NC_1:g.6G>C"]}


def test_a_list_valued_allele_key_is_flattened():
    """Recoder returns either a single result object or a list of them under an allele key."""
    ensembl = client(
        FakeResponse(
            [
                {
                    "A": [
                        {"input": "NP_1.1:p.Gly15Asn", "hgvsg": ["NC_1:g.5G>A"]},
                        {"input": "NP_1.1:p.Gly15Asn", "hgvsg": ["NC_1:g.7G>A"]},
                    ]
                }
            ]
        )
    )

    assert ensembl.to_genomic(["NP_1.1:p.Gly15Asn"]) == {"NP_1.1:p.Gly15Asn": ["NC_1:g.5G>A", "NC_1:g.7G>A"]}


def test_an_input_with_no_genomic_equivalent_is_absent_from_the_mapping():
    ensembl = client(FakeResponse([{"input": "NP_1.1:p.Met4=", "A": {"input": "NP_1.1:p.Met4="}}]))
    assert ensembl.to_genomic(["NP_1.1:p.Met4="]) == {}


def test_the_top_level_input_key_is_not_treated_as_an_allele():
    ensembl = client(FakeResponse([{"input": "NP_1.1:p.Gly15Asn"}]))
    assert ensembl.to_genomic(["NP_1.1:p.Gly15Asn"]) == {}


def test_the_recoder_gets_its_own_longer_timeout():
    """Recoder is markedly slower than VEP for large batches and 504s are routine."""
    ensembl = client(FakeResponse([]), timeout_seconds=60, recoder_timeout_seconds=600)
    ensembl.to_genomic(["NP_1.1:p.Gly15Asn"])

    _, _, kwargs = ensembl._session.calls[0]
    assert kwargs["timeout"] == 600


# --- software_release ---


def test_software_release_is_returned_as_a_string():
    ensembl = client(FakeResponse({"release": 116}))
    assert ensembl.software_release() == "116"
