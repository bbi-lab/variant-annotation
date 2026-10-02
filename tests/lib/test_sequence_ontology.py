"""Unit tests for the Sequence Ontology consequence vocabulary."""

import pytest

from variant_annotation.lib import sequence_ontology as so

pytestmark = pytest.mark.unit


def test_ranking_matches_ensembls_published_term_count():
    assert len(so.ENSEMBL_CONSEQUENCE_RANKING) == 41


def test_terms_are_unique():
    assert len(set(so.CONSEQUENCE_TERMS)) == len(so.CONSEQUENCE_TERMS)


def test_accessions_are_unique_and_well_formed():
    accessions = [accession for _, accession in so.ENSEMBL_CONSEQUENCE_RANKING]
    assert len(set(accessions)) == len(accessions)
    assert all(accession.startswith("SO:") and accession[3:].isdigit() for accession in accessions)


def test_ranks_are_dense_and_one_based():
    assert sorted(so.CONSEQUENCE_RANKS.values()) == list(range(1, len(so.CONSEQUENCE_TERMS) + 1))
    assert so.CONSEQUENCE_RANKS["transcript_ablation"] == 1
    assert so.CONSEQUENCE_RANKS["sequence_variant"] == len(so.CONSEQUENCE_TERMS)


@pytest.mark.parametrize(
    "term",
    [
        # The three fine-grained splice terms Ensembl added in release 105. Their absence from a
        # hand-maintained list makes real VEP output unrankable.
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
        # Canonical spellings — the *_variant suffix is not optional.
        "start_retained_variant",
        "stop_retained_variant",
        "NMD_transcript_variant",
        "coding_transcript_variant",
    ],
)
def test_terms_absent_from_the_legacy_hand_rolled_list_are_present(term):
    assert so.is_known_term(term)


@pytest.mark.parametrize(
    "term",
    [
        "start_retained",  # non-canonical: missing the _variant suffix
        "stop_retained",
        "disruptive_inframe_insertion",  # SnpEff/ANNOVAR vocabulary, never emitted by VEP
        "disruptive_inframe_deletion",
        "nc_transcript_variant",  # retired VEP alias
        "non_coding_exon_variant",
        "gene_variant",  # SO parent terms, not calculated consequences
        "exon_variant",
        "variant_of_unknown_significance",  # not a consequence at all
    ],
)
def test_terms_ensembl_does_not_emit_are_absent(term):
    assert not so.is_known_term(term)


def test_feature_elongation_outranks_missense():
    """Ensembl ranks feature_elongation/truncation 9th and 10th — above missense at 13th.

    The legacy hand-rolled list placed them ~30 positions too low, so a variant reported as both
    feature_elongation and a milder term resolved to the milder one.
    """
    assert so.rank("feature_elongation") < so.rank("missense_variant")
    assert so.rank("feature_truncation") < so.rank("missense_variant")


def test_sequence_variant_is_the_least_severe_term():
    """sequence_variant is Ensembl's catch-all floor; the legacy list ranked it above intron_variant."""
    assert so.rank("sequence_variant") > so.rank("intron_variant")
    assert so.rank("sequence_variant") > so.rank("intergenic_variant")
    assert so.rank("sequence_variant") == max(so.CONSEQUENCE_RANKS.values())


def test_splice_region_outranks_synonymous():
    assert so.rank("splice_region_variant") < so.rank("synonymous_variant")


def test_unknown_terms_rank_last_without_raising():
    assert so.rank("some_future_ensembl_term") == so.UNRANKED
    assert so.rank("some_future_ensembl_term") > so.rank("sequence_variant")


def test_most_severe_picks_by_rank_not_input_order():
    assert so.most_severe(["synonymous_variant", "stop_gained", "intron_variant"]) == "stop_gained"


def test_most_severe_of_empty_is_none():
    assert so.most_severe([]) is None


def test_most_severe_keeps_an_unknown_term_rather_than_returning_none():
    """An unrecognised term is still a real answer from VEP — dropping it would fabricate a miss."""
    assert so.most_severe(["a_brand_new_term"]) == "a_brand_new_term"


def test_most_severe_prefers_a_known_term_over_an_unknown_one():
    assert so.most_severe(["a_brand_new_term", "intron_variant"]) == "intron_variant"


def test_sort_by_severity_orders_and_deduplicates():
    ordered = so.sort_by_severity(["synonymous_variant", "splice_region_variant", "synonymous_variant", "stop_gained"])
    assert ordered == ["stop_gained", "splice_region_variant", "synonymous_variant"]


def test_sort_by_severity_is_stable_across_unknown_terms():
    assert so.sort_by_severity(["zeta_term", "alpha_term"]) == ["zeta_term", "alpha_term"]
