"""Unit tests for the VEP resolution orchestration, driven through stub ports.

No HTTP: the ports exist so the batching, fallback, and failure-attribution behaviour can be exercised
against controlled responses.
"""

import pytest

from variant_annotation.lib.vep import (
    NO_CHANGE_TERM,
    ConsequenceOutcome,
    ConsequenceSource,
    VepConfig,
    VepInput,
    resolve_consequences,
)

pytestmark = pytest.mark.unit


class StubVep:
    """A VepLookup returning canned entries, recording every request it received."""

    def __init__(self, entries=None, *, fail_on=(), fail_all=False):
        # entries: hgvs -> the response entry body (without "input", which is added here)
        self._entries = entries or {}
        self._fail_on = set(fail_on)
        self._fail_all = fail_all
        self.requests = []

    def annotate_hgvs(self, hgvs, *, refseq):
        self.requests.append((tuple(hgvs), refseq))
        if self._fail_all or (self._fail_on & set(hgvs)):
            raise RuntimeError("VEP is down")
        return [{"input": h, **self._entries[h]} for h in hgvs if h in self._entries]


class StubRecoder:
    def __init__(self, mapping=None, *, fail=False):
        self._mapping = mapping or {}
        self._fail = fail
        self.calls = []

    def to_genomic(self, hgvs):
        self.calls.append(tuple(hgvs))
        if self._fail:
            raise RuntimeError("Recoder is down")
        return {h: self._mapping[h] for h in hgvs if h in self._mapping}


class StubReference:
    def __init__(self, sequences=None):
        # (transcript, start, stop) -> reference bases
        self._sequences = sequences or {}

    def coding_interval_reference(self, transcript, start, stop):
        return self._sequences.get((transcript, start, stop))


def transcript_entry(transcript, terms, *, headline=None):
    entry = {"transcript_consequences": [{"transcript_id": transcript, "consequence_terms": terms}]}
    if headline:
        entry["most_severe_consequence"] = headline
    return entry


def test_no_inputs_makes_no_requests():
    vep = StubVep()
    assert resolve_consequences([], vep=vep) == []
    assert vep.requests == []


def test_resolutions_come_back_in_input_order():
    inputs = [
        VepInput("NM_1.1:c.3A>G"),
        VepInput("NM_1.1:c.1A>G"),
        VepInput("NM_1.1:c.2A>G"),
    ]
    vep = StubVep(
        {
            "NM_1.1:c.1A>G": transcript_entry("NM_1.1", ["missense_variant"]),
            "NM_1.1:c.2A>G": transcript_entry("NM_1.1", ["synonymous_variant"]),
            "NM_1.1:c.3A>G": transcript_entry("NM_1.1", ["stop_gained"]),
        }
    )

    results = resolve_consequences(inputs, vep=vep)

    assert [r.input.hgvs for r in results] == [i.hgvs for i in inputs]
    assert [r.most_severe_consequence for r in results] == [
        "stop_gained",
        "missense_variant",
        "synonymous_variant",
    ]


def test_refseq_and_ensembl_inputs_are_requested_separately():
    inputs = [VepInput("NM_1.1:c.1A>G"), VepInput("ENST1.1:c.1A>G")]
    vep = StubVep(
        {
            "NM_1.1:c.1A>G": transcript_entry("NM_1.1", ["missense_variant"]),
            "ENST1.1:c.1A>G": transcript_entry("ENST1.1", ["synonymous_variant"]),
        }
    )

    resolve_consequences(inputs, vep=vep)

    assert sorted(vep.requests) == [
        (("ENST1.1:c.1A>G",), False),
        (("NM_1.1:c.1A>G",), True),
    ]


def test_a_genomic_input_with_a_transcript_hint_is_sent_with_refseq():
    inputs = [VepInput("NC_1:g.100A>G", transcript="NM_1.1")]
    vep = StubVep({"NC_1:g.100A>G": transcript_entry("NM_1.1", ["5_prime_UTR_variant"])})

    results = resolve_consequences(inputs, vep=vep)

    assert vep.requests == [(("NC_1:g.100A>G",), True)]
    assert results[0].source is ConsequenceSource.TRANSCRIPT
    assert results[0].most_severe_consequence == "5_prime_UTR_variant"


def test_duplicate_hgvs_are_sent_once_but_answered_individually():
    inputs = [VepInput("NM_1.1:c.1A>G"), VepInput("NM_1.1:c.1A>G")]
    vep = StubVep({"NM_1.1:c.1A>G": transcript_entry("NM_1.1", ["missense_variant"])})

    results = resolve_consequences(inputs, vep=vep)

    assert vep.requests == [(("NM_1.1:c.1A>G",), True)]
    assert len(results) == 2
    assert all(r.most_severe_consequence == "missense_variant" for r in results)


def test_the_same_hgvs_with_different_transcripts_gets_different_answers():
    """Two inputs sharing an HGVS string but naming different transcripts are different questions."""
    inputs = [
        VepInput("NC_1:g.100A>G", transcript="NM_1.1"),
        VepInput("NC_1:g.100A>G", transcript="NM_2.1"),
    ]
    vep = StubVep(
        {
            "NC_1:g.100A>G": {
                "most_severe_consequence": "stop_gained",
                "transcript_consequences": [
                    {"transcript_id": "NM_1.1", "consequence_terms": ["synonymous_variant"]},
                    {"transcript_id": "NM_2.1", "consequence_terms": ["missense_variant"]},
                ],
            }
        }
    )

    results = resolve_consequences(inputs, vep=vep)

    assert [r.most_severe_consequence for r in results] == ["synonymous_variant", "missense_variant"]


def test_an_input_vep_omits_is_absent_when_there_is_no_recoder():
    vep = StubVep({})

    results = resolve_consequences([VepInput("NM_1.1:c.1A>G")], vep=vep)

    assert results[0].outcome is ConsequenceOutcome.ABSENT
    assert results[0].error is None


def test_a_failed_batch_reports_its_inputs_errored_not_absent():
    vep = StubVep(fail_all=True)

    results = resolve_consequences([VepInput("NM_1.1:c.1A>G")], vep=vep)

    assert results[0].outcome is ConsequenceOutcome.ERRORED
    assert "VEP is down" in results[0].error


def test_one_failed_batch_does_not_discard_another_batch_s_results():
    """A transport blip must not throw away work already paid for against a rate-limited service."""
    inputs = [VepInput("NM_1.1:c.1A>G"), VepInput("ENST1.1:c.1A>G")]
    vep = StubVep(
        {"ENST1.1:c.1A>G": transcript_entry("ENST1.1", ["missense_variant"])},
        fail_on=["NM_1.1:c.1A>G"],
    )

    results = resolve_consequences(inputs, vep=vep)

    assert results[0].outcome is ConsequenceOutcome.ERRORED
    assert results[1].outcome is ConsequenceOutcome.RESOLVED


def test_an_errored_input_is_never_sent_to_the_recoder():
    """The answer is unknown, so a fallback would be answering a question that was never asked."""
    recoder = StubRecoder({"NM_1.1:c.1A>G": ["NC_1:g.5A>G"]})
    vep = StubVep(fail_all=True)

    results = resolve_consequences([VepInput("NM_1.1:c.1A>G")], vep=vep, recoder=recoder)

    assert recoder.calls == []
    assert results[0].outcome is ConsequenceOutcome.ERRORED


def test_a_protein_input_is_recovered_through_the_recoder():
    recoder = StubRecoder({"NP_1.1:p.Gly15Asn": ["NC_1:g.5G>A"]})
    vep = StubVep({"NC_1:g.5G>A": {"most_severe_consequence": "missense_variant"}})

    results = resolve_consequences([VepInput("NP_1.1:p.Gly15Asn")], vep=vep, recoder=recoder)

    assert results[0].outcome is ConsequenceOutcome.RESOLVED
    assert results[0].most_severe_consequence == "missense_variant"
    assert results[0].source is ConsequenceSource.MOST_SEVERE


def test_multiple_recoded_forms_combine_to_the_most_severe():
    recoder = StubRecoder({"NP_1.1:p.Gly15Asn": ["NC_1:g.5G>A", "NC_1:g.6G>A"]})
    vep = StubVep(
        {
            "NC_1:g.5G>A": {"most_severe_consequence": "missense_variant"},
            "NC_1:g.6G>A": {"most_severe_consequence": "stop_gained"},
        }
    )

    results = resolve_consequences([VepInput("NP_1.1:p.Gly15Asn")], vep=vep, recoder=recoder)

    assert results[0].most_severe_consequence == "stop_gained"


def test_the_recoder_is_skipped_when_disabled():
    recoder = StubRecoder({"NP_1.1:p.Gly15Asn": ["NC_1:g.5G>A"]})
    vep = StubVep({})

    results = resolve_consequences(
        [VepInput("NP_1.1:p.Gly15Asn")],
        vep=vep,
        recoder=recoder,
        config=VepConfig(recode_misses=False),
    )

    assert recoder.calls == []
    assert results[0].outcome is ConsequenceOutcome.ABSENT


def test_a_failed_recoder_reports_errored_not_absent():
    recoder = StubRecoder(fail=True)
    vep = StubVep({})

    results = resolve_consequences([VepInput("NP_1.1:p.Gly15Asn")], vep=vep, recoder=recoder)

    assert results[0].outcome is ConsequenceOutcome.ERRORED
    assert "Recoder is down" in results[0].error


def test_an_input_the_recoder_has_no_genomic_form_for_stays_absent():
    recoder = StubRecoder({})
    vep = StubVep({})

    results = resolve_consequences([VepInput("NP_1.1:p.Met4=")], vep=vep, recoder=recoder)

    assert results[0].outcome is ConsequenceOutcome.ABSENT


def test_a_vep_entry_with_no_consequence_is_still_recoded():
    """VEP knowing the variant but classifying nothing is a miss, not a settled negative."""
    recoder = StubRecoder({"NM_1.1:c.1A>G": ["NC_1:g.5A>G"]})
    vep = StubVep(
        {
            "NM_1.1:c.1A>G": {},
            "NC_1:g.5A>G": {"most_severe_consequence": "intron_variant"},
        }
    )

    results = resolve_consequences([VepInput("NM_1.1:c.1A>G")], vep=vep, recoder=recoder)

    assert recoder.calls == [("NM_1.1:c.1A>G",)]
    assert results[0].most_severe_consequence == "intron_variant"


def test_the_recoder_call_is_batched_by_batch_size():
    """The Recoder POST is chunked like VEP, so a large miss set cannot exceed Ensembl's per-POST limit."""
    recoder = StubRecoder(
        {
            "NP_1.1:p.Gly15Asn": ["NC_1:g.5G>A"],
            "NP_1.1:p.Gly16Asn": ["NC_1:g.8G>A"],
        }
    )
    vep = StubVep(
        {
            "NC_1:g.5G>A": {"most_severe_consequence": "missense_variant"},
            "NC_1:g.8G>A": {"most_severe_consequence": "missense_variant"},
        }
    )

    results = resolve_consequences(
        [VepInput("NP_1.1:p.Gly15Asn"), VepInput("NP_1.1:p.Gly16Asn")],
        vep=vep,
        recoder=recoder,
        config=VepConfig(batch_size=1),
    )

    # Two misses at batch size one -> two separate Recoder calls, each carrying a single HGVS.
    assert len(recoder.calls) == 2
    assert all(len(call) == 1 for call in recoder.calls)
    assert all(r.outcome is ConsequenceOutcome.RESOLVED for r in results)


def test_one_failed_recoder_batch_does_not_sink_another():
    """A failed Recoder chunk marks only its own inputs ERRORED; a sibling chunk still resolves."""

    class PartialRecoder:
        def __init__(self):
            self.calls = []

        def to_genomic(self, hgvs):
            self.calls.append(tuple(hgvs))
            if any("Bad" in h for h in hgvs):
                raise RuntimeError("Recoder is down")
            return {"NP_1.1:p.Gly15Asn": ["NC_1:g.5G>A"]}

    recoder = PartialRecoder()
    vep = StubVep({"NC_1:g.5G>A": {"most_severe_consequence": "missense_variant"}})

    results = resolve_consequences(
        [VepInput("NP_1.1:p.Gly15Asn"), VepInput("NP_1.1:p.BadXyz")],
        vep=vep,
        recoder=recoder,
        config=VepConfig(batch_size=1),
    )

    by_hgvs = {r.input.hgvs: r for r in results}
    assert by_hgvs["NP_1.1:p.Gly15Asn"].outcome is ConsequenceOutcome.RESOLVED
    assert by_hgvs["NP_1.1:p.BadXyz"].outcome is ConsequenceOutcome.ERRORED
    assert "Recoder is down" in by_hgvs["NP_1.1:p.BadXyz"].error


# --- reference-identical (#760) ---


def test_an_explicit_no_change_expression_is_labelled_rather_than_left_null():
    vep = StubVep({})

    results = resolve_consequences([VepInput("NP_000537.3:p.Met4=")], vep=vep, reference=StubReference())

    assert results[0].outcome is ConsequenceOutcome.RESOLVED
    assert results[0].most_severe_consequence == NO_CHANGE_TERM
    assert results[0].source is ConsequenceSource.REFERENCE_IDENTICAL


def test_a_delins_that_restates_the_reference_is_labelled_no_change():
    vep = StubVep({})
    reference = StubReference({("NM_1.1", 12, 14): "GCT"})

    results = resolve_consequences([VepInput("NM_1.1:c.12_14delinsGCT")], vep=vep, reference=reference)

    assert results[0].most_severe_consequence == NO_CHANGE_TERM
    assert results[0].source is ConsequenceSource.REFERENCE_IDENTICAL


def test_a_delins_that_changes_the_reference_stays_absent():
    vep = StubVep({})
    reference = StubReference({("NM_1.1", 12, 14): "GCT"})

    results = resolve_consequences([VepInput("NM_1.1:c.12_14delinsAAA")], vep=vep, reference=reference)

    assert results[0].outcome is ConsequenceOutcome.ABSENT


def test_the_reference_is_not_consulted_for_inputs_vep_resolved():
    """The no-change check must cost nothing for the overwhelming majority of inputs."""

    class ExplodingReference:
        def coding_interval_reference(self, transcript, start, stop):
            raise AssertionError("reference must not be consulted for a resolved input")

    vep = StubVep({"NM_1.1:c.12_14delinsGCT": transcript_entry("NM_1.1", ["synonymous_variant"])})

    results = resolve_consequences([VepInput("NM_1.1:c.12_14delinsGCT")], vep=vep, reference=ExplodingReference())

    assert results[0].most_severe_consequence == "synonymous_variant"


def test_an_errored_input_is_not_relabelled_as_no_change():
    """An unknown answer must stay unknown; labelling it no_change would fabricate a fact."""
    vep = StubVep(fail_all=True)
    reference = StubReference({("NM_1.1", 12, 14): "GCT"})

    results = resolve_consequences([VepInput("NM_1.1:c.12_14delinsGCT")], vep=vep, reference=reference)

    assert results[0].outcome is ConsequenceOutcome.ERRORED


def test_concurrent_batches_produce_the_same_answers_as_serial_ones():
    inputs = [VepInput(f"NM_1.1:c.{n}A>G") for n in range(1, 8)]
    entries = {f"NM_1.1:c.{n}A>G": transcript_entry("NM_1.1", ["missense_variant"]) for n in range(1, 8)}

    serial = resolve_consequences(inputs, vep=StubVep(entries), config=VepConfig(batch_size=2))
    concurrent = resolve_consequences(inputs, vep=StubVep(entries), config=VepConfig(batch_size=2, max_workers=4))

    assert [r.most_severe_consequence for r in serial] == [r.most_severe_consequence for r in concurrent]
    assert [r.input.hgvs for r in concurrent] == [i.hgvs for i in inputs]
