"""Unit tests for reference-identical HGVS recognition."""

import pytest

from variant_annotation.lib.hgvs.no_change import parse_coding_delins, states_no_change

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    "hgvs",
    [
        "NC_000003.12:g.10141916C=",
        "NM_000049.4:c.123=",
        "NM_000049.4:c.123_125=",
        # Residues with no codon redundancy: these cannot be recoded to an alternate synonymous
        # codon at all, so recognising the notation is the only way to classify them.
        "NP_000537.3:p.Met4=",
        "NP_000537.3:p.Trp45=",
        "NR_003051.3:n.55=",
    ],
)
def test_explicit_no_change_notation_is_recognised(hgvs):
    assert states_no_change(hgvs)


@pytest.mark.parametrize(
    "hgvs",
    [
        "NM_000049.4:c.256A>G",
        "NP_000537.3:p.Arg175His",
        "NC_000003.12:g.10141916C>T",
        "NM_000049.4:c.256_257delinsAA",
        "",
        "not-hgvs",
        "NM_000049.4",
    ],
)
def test_ordinary_variants_are_not_recognised_as_no_change(hgvs):
    assert not states_no_change(hgvs)


def test_surrounding_whitespace_is_tolerated():
    assert states_no_change("  NM_000049.4:c.123=  ")


def test_coding_delins_is_parsed_into_its_parts():
    assert parse_coding_delins("NM_000049.4:c.256_257delinsAA") == ("NM_000049.4", 256, 257, "AA")


def test_lowercase_replacement_bases_are_preserved_verbatim():
    """Case folding is the comparison's job, not the parser's."""
    assert parse_coding_delins("NM_000049.4:c.12_14delinsgct") == ("NM_000049.4", 12, 14, "gct")


@pytest.mark.parametrize(
    "hgvs",
    [
        "NM_000049.4:c.256A>G",  # not a delins
        "NM_000049.4:c.256delinsAA",  # no interval
        "NC_000003.12:g.100_102delinsAAA",  # genomic, not coding
        "NM_000049.4:c.256_257delinsXX",  # not nucleotides
        "",
    ],
)
def test_non_coding_interval_delins_returns_none(hgvs):
    assert parse_coding_delins(hgvs) is None
