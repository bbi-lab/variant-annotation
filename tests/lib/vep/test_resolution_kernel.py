"""Unit tests for the VEP consequence-resolution kernel.

The kernel is pure, so every test here drives it with literal VEP response payloads — the shape the
real API returns, trimmed to the fields resolution reads.
"""

import pytest

from variant_annotation.lib.vep import (
    NO_CHANGE_TERM,
    ConsequenceOutcome,
    ConsequenceSource,
    VepInput,
    combine_recoded,
    needs_refseq_transcripts,
    partition_batches,
    requested_transcript,
    resolve_entry,
)

pytestmark = pytest.mark.unit


def vep_entry(submitted, *, most_severe=None, transcripts=()):
    """A VEP response entry: the echoed input, the cross-transcript headline, and per-transcript terms.

    ``transcripts`` is a sequence of ``(transcript_id, [consequence_terms])``.
    """
    entry = {"input": submitted}
    if most_severe is not None:
        entry["most_severe_consequence"] = most_severe
    if transcripts:
        entry["transcript_consequences"] = [
            {"transcript_id": tx, "consequence_terms": list(terms)} for tx, terms in transcripts
        ]
    return entry


# --- requested_transcript ---


def test_coding_hgvs_infers_its_own_transcript():
    assert requested_transcript(VepInput("NM_000049.4:c.256A>G")) == "NM_000049"


def test_ensembl_coding_hgvs_infers_its_own_transcript():
    assert requested_transcript(VepInput("ENST00000366667.8:c.803C>T")) == "ENST00000366667"


def test_genomic_hgvs_infers_nothing():
    """An NC_ accession is a chromosome, so there is no transcript to resolve against."""
    assert requested_transcript(VepInput("NC_000010.11:g.43086481G>T")) is None


def test_explicit_transcript_overrides_inference():
    """The genomic-allele fix: pass the paired coding transcript and the match becomes possible."""
    resolved = requested_transcript(VepInput("NC_000010.11:g.43086481G>T", transcript="NM_020975.6"))
    assert resolved == "NM_020975"


def test_explicit_transcript_wins_over_a_coding_accession():
    resolved = requested_transcript(VepInput("NM_000049.4:c.256A>G", transcript="ENST00000366667.8"))
    assert resolved == "ENST00000366667"


def test_protein_hgvs_infers_nothing():
    """NP_ identifies the product, but VEP keys transcript_consequences on transcript ids."""
    assert requested_transcript(VepInput("NP_071890.2:p.Gly15Asn")) is None


# --- needs_refseq_transcripts / partition_batches ---


def test_refseq_coding_input_requests_the_refseq_transcript_set():
    assert needs_refseq_transcripts(VepInput("NM_000049.4:c.256A>G")) is True


def test_ensembl_coding_input_does_not_request_refseq():
    assert needs_refseq_transcripts(VepInput("ENST00000366667.8:c.803C>T")) is False


def test_genomic_input_with_a_refseq_transcript_hint_requests_refseq():
    """The bug a prefix-based implementation has.

    Keying the refseq flag on the HGVS prefix sends this to the Ensembl transcript set, so the NM_ hint
    can never match anything in transcript_consequences and the answer silently degrades to most_severe.
    """
    assert needs_refseq_transcripts(VepInput("NC_000010.11:g.43086481G>T", transcript="NM_020975.6")) is True


def test_bare_genomic_input_does_not_request_refseq():
    assert needs_refseq_transcripts(VepInput("NC_000010.11:g.43086481G>T")) is False


def test_batches_separate_the_two_transcript_sets():
    """refseq is a per-request flag, so mixing kinds would query one group against the wrong set."""
    inputs = [
        VepInput("NM_000049.4:c.1A>G"),
        VepInput("ENST00000366667.8:c.2A>G"),
        VepInput("NR_003051.3:c.3A>G"),
    ]
    batches = partition_batches(inputs, batch_size=10)

    assert [refseq for _, refseq in batches] == [True, False]
    assert [i.hgvs for i in batches[0][0]] == ["NM_000049.4:c.1A>G", "NR_003051.3:c.3A>G"]
    assert [i.hgvs for i in batches[1][0]] == ["ENST00000366667.8:c.2A>G"]


def test_batches_respect_the_size_cap_within_a_group():
    inputs = [VepInput(f"NM_000049.4:c.{n}A>G") for n in range(1, 6)]
    batches = partition_batches(inputs, batch_size=2)

    assert [len(batch) for batch, _ in batches] == [2, 2, 1]
    assert all(refseq for _, refseq in batches)


def test_batches_preserve_input_order_within_a_group():
    inputs = [VepInput(f"NM_000049.4:c.{n}A>G") for n in range(1, 4)]
    ((batch, _),) = partition_batches(inputs, batch_size=10)
    assert [i.hgvs for i in batch] == [i.hgvs for i in inputs]


def test_empty_input_produces_no_batches():
    assert partition_batches([], batch_size=10) == []


def test_zero_batch_size_is_rejected():
    with pytest.raises(ValueError, match="batch_size must be at least 1"):
        partition_batches([VepInput("NM_000049.4:c.1A>G")], batch_size=0)


# --- resolve_entry: the transcript match ---


def test_matching_transcript_wins_over_the_cross_transcript_headline():
    """The core correctness claim of the whole module.

    A synonymous change on the assayed transcript must not be reported as splice_acceptor_variant
    because some other overlapping isoform reads that position differently.
    """
    entry = vep_entry(
        "NM_000049.4:c.256A>G",
        most_severe="splice_acceptor_variant",
        transcripts=[
            ("NM_000049.4", ["synonymous_variant"]),
            ("NM_999999.1", ["splice_acceptor_variant"]),
        ],
    )

    resolution = resolve_entry(VepInput("NM_000049.4:c.256A>G"), entry)

    assert resolution.outcome is ConsequenceOutcome.RESOLVED
    assert resolution.most_severe_consequence == "synonymous_variant"
    assert resolution.source is ConsequenceSource.TRANSCRIPT
    assert resolution.matched_transcript == "NM_000049.4"


def test_transcript_match_ignores_accession_versions():
    entry = vep_entry(
        "NM_000049.4:c.256A>G",
        most_severe="stop_gained",
        transcripts=[("NM_000049.7", ["missense_variant"])],
    )

    resolution = resolve_entry(VepInput("NM_000049.4:c.256A>G"), entry)

    assert resolution.most_severe_consequence == "missense_variant"
    assert resolution.source is ConsequenceSource.TRANSCRIPT
    # The versioned id VEP actually used is recorded, not the version we asked with.
    assert resolution.matched_transcript == "NM_000049.7"


def test_all_terms_from_the_matched_transcript_are_kept_severity_ordered():
    entry = vep_entry(
        "NM_000049.4:c.256A>G",
        most_severe="stop_gained",
        transcripts=[("NM_000049.4", ["synonymous_variant", "splice_region_variant"])],
    )

    resolution = resolve_entry(VepInput("NM_000049.4:c.256A>G"), entry)

    assert resolution.consequence_terms == ["splice_region_variant", "synonymous_variant"]
    assert resolution.most_severe_consequence == "splice_region_variant"


def test_a_genomic_input_with_a_transcript_hint_resolves_against_that_transcript():
    """What lets genomic-level alleles escape the cross-transcript headline."""
    entry = vep_entry(
        "NC_000010.11:g.43086481G>T",
        most_severe="transcript_ablation",
        transcripts=[
            ("NM_020975.6", ["5_prime_UTR_variant"]),
            ("NM_888888.1", ["transcript_ablation"]),
        ],
    )

    resolution = resolve_entry(VepInput("NC_000010.11:g.43086481G>T", transcript="NM_020975.6"), entry)

    assert resolution.most_severe_consequence == "5_prime_UTR_variant"
    assert resolution.source is ConsequenceSource.TRANSCRIPT


def test_a_bare_genomic_input_falls_back_to_the_headline():
    entry = vep_entry(
        "NC_000010.11:g.43086481G>T",
        most_severe="intron_variant",
        transcripts=[("NM_020975.6", ["5_prime_UTR_variant"])],
    )

    resolution = resolve_entry(VepInput("NC_000010.11:g.43086481G>T"), entry)

    assert resolution.most_severe_consequence == "intron_variant"
    assert resolution.source is ConsequenceSource.MOST_SEVERE
    assert resolution.matched_transcript is None
    assert resolution.consequence_terms == ["intron_variant"]


def test_no_matching_transcript_falls_back_to_the_headline():
    entry = vep_entry(
        "NM_000049.4:c.256A>G",
        most_severe="missense_variant",
        transcripts=[("NM_111111.1", ["intron_variant"])],
    )

    resolution = resolve_entry(VepInput("NM_000049.4:c.256A>G"), entry)

    assert resolution.most_severe_consequence == "missense_variant"
    assert resolution.source is ConsequenceSource.MOST_SEVERE


def test_a_matched_transcript_with_no_terms_falls_back_rather_than_reporting_absent():
    """The match proves the transcript is in VEP's set; an empty term list is a gap in that entry."""
    entry = vep_entry(
        "NM_000049.4:c.256A>G",
        most_severe="missense_variant",
        transcripts=[("NM_000049.4", [])],
    )

    resolution = resolve_entry(VepInput("NM_000049.4:c.256A>G"), entry)

    assert resolution.outcome is ConsequenceOutcome.RESOLVED
    assert resolution.most_severe_consequence == "missense_variant"
    assert resolution.source is ConsequenceSource.MOST_SEVERE


def test_an_entry_with_nothing_at_all_is_absent():
    resolution = resolve_entry(VepInput("NM_000049.4:c.256A>G"), vep_entry("NM_000049.4:c.256A>G"))

    assert resolution.outcome is ConsequenceOutcome.ABSENT
    assert resolution.most_severe_consequence is None
    assert resolution.source is None
    assert resolution.error is None


def test_unrankable_upstream_terms_still_resolve():
    """A term Ensembl adds after this library's ranking was refreshed must not become a missing answer."""
    entry = vep_entry(
        "NM_000049.4:c.256A>G",
        transcripts=[("NM_000049.4", ["some_future_ensembl_term"])],
    )

    resolution = resolve_entry(VepInput("NM_000049.4:c.256A>G"), entry)

    assert resolution.outcome is ConsequenceOutcome.RESOLVED
    assert resolution.most_severe_consequence == "some_future_ensembl_term"


def test_malformed_transcript_entries_are_skipped_not_fatal():
    entry = {
        "input": "NM_000049.4:c.256A>G",
        "most_severe_consequence": "missense_variant",
        "transcript_consequences": [
            "not-a-mapping",
            {"consequence_terms": ["stop_gained"]},  # no transcript_id
            {"transcript_id": "NM_000049.4", "consequence_terms": ["synonymous_variant"]},
        ],
    }

    resolution = resolve_entry(VepInput("NM_000049.4:c.256A>G"), entry)

    assert resolution.most_severe_consequence == "synonymous_variant"
    assert resolution.source is ConsequenceSource.TRANSCRIPT


# --- combine_recoded ---


def _resolved(hgvs, term):
    return resolve_entry(VepInput(hgvs), vep_entry(hgvs, most_severe=term))


def _absent(hgvs):
    return resolve_entry(VepInput(hgvs), vep_entry(hgvs))


def test_recoded_forms_combine_to_the_most_severe():
    combined = combine_recoded(
        VepInput("NP_071890.2:p.Gly15Asn"),
        [_resolved("NC_1:g.1A>G", "missense_variant"), _resolved("NC_1:g.2A>G", "stop_gained")],
    )

    assert combined.outcome is ConsequenceOutcome.RESOLVED
    assert combined.most_severe_consequence == "stop_gained"


def test_a_combined_answer_is_never_claimed_as_transcript_specific():
    """A recoded genomic query carries no transcript, so no transcript claim can be made."""
    combined = combine_recoded(
        VepInput("NP_071890.2:p.Gly15Asn", transcript="NM_030952.2"),
        [_resolved("NC_1:g.1A>G", "missense_variant")],
    )

    assert combined.source is ConsequenceSource.MOST_SEVERE
    assert combined.matched_transcript is None


def test_all_absent_recoded_forms_combine_to_absent():
    combined = combine_recoded(VepInput("NP_1:p.Met1Val"), [_absent("NC_1:g.1A>G")])
    assert combined.outcome is ConsequenceOutcome.ABSENT


def test_a_failed_recoded_form_makes_the_answer_errored_not_absent():
    """A partial failure is unknown, not a confirmed negative — storing it as absent would be a lie."""
    from variant_annotation.lib.vep.types import ConsequenceResolution

    errored = ConsequenceResolution(input=VepInput("NC_1:g.2A>G"), outcome=ConsequenceOutcome.ERRORED, error="boom")

    combined = combine_recoded(VepInput("NP_1:p.Met1Val"), [_absent("NC_1:g.1A>G"), errored])

    assert combined.outcome is ConsequenceOutcome.ERRORED
    assert combined.error == "boom"


def test_a_resolved_form_outweighs_a_failed_sibling():
    from variant_annotation.lib.vep.types import ConsequenceResolution

    errored = ConsequenceResolution(input=VepInput("NC_1:g.2A>G"), outcome=ConsequenceOutcome.ERRORED, error="boom")

    combined = combine_recoded(VepInput("NP_1:p.Met1Val"), [_resolved("NC_1:g.1A>G", "missense_variant"), errored])

    assert combined.outcome is ConsequenceOutcome.RESOLVED
    assert combined.most_severe_consequence == "missense_variant"


def test_combining_nothing_is_absent():
    assert combine_recoded(VepInput("NP_1:p.Met1Val"), []).outcome is ConsequenceOutcome.ABSENT


# --- the regression this module exists to prevent ---


def test_a_protein_missense_is_not_reported_as_transcript_ablation():
    """Regression for the observed production defect.

    108 single-residue protein missense variants (NP_071890.2:p.Gly15*) were stored as
    transcript_ablation — the most severe term in the ranking — because a protein HGVS was pushed through
    Variant Recoder to genomic and then read off VEP's cross-transcript headline. With the protein
    variant's coding transcript supplied, resolution reads the assayed transcript instead.
    """
    entry = vep_entry(
        "NM_030952.2:c.43G>A",
        most_severe="transcript_ablation",
        transcripts=[
            ("NM_030952.2", ["missense_variant"]),
            ("NM_777777.1", ["transcript_ablation"]),
        ],
    )

    resolution = resolve_entry(VepInput("NM_030952.2:c.43G>A"), entry)

    assert resolution.most_severe_consequence == "missense_variant"
    assert resolution.source is ConsequenceSource.TRANSCRIPT


def test_no_change_term_is_not_an_so_term():
    """#760: 'no change' is a distinct claim from synonymous_variant, and must not borrow its label."""
    from variant_annotation.lib import sequence_ontology as so

    assert not so.is_known_term(NO_CHANGE_TERM)
