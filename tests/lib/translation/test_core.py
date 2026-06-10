import pytest

from variant_annotation.lib.translation.types import (
    TranslationConfig,
    TranslationResult,
    TranslationError,
    VariantInput,
    WtCodonMode,
)

pytestmark = pytest.mark.unit


class _StubTranscripts:
    def transcript_for_protein(self, acc):
        return None

    def codon_at(self, tx, pos):
        return None


class _StubCoordinates:
    def c_to_p(self, c):
        raise NotImplementedError

    def g_to_c(self, g, tx):
        raise NotImplementedError

    def c_to_g(self, c):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# _classify_kind
# ---------------------------------------------------------------------------


def test_classify_kind_protein():
    from variant_annotation.lib.translation._core import _classify_kind, _Kind

    assert _classify_kind("NP_000001.1:p.Ala1Val") is _Kind.PROTEIN


def test_classify_kind_coding():
    from variant_annotation.lib.translation._core import _classify_kind, _Kind

    assert _classify_kind("NM_000001.1:c.368C>T") is _Kind.CODING


def test_classify_kind_genomic():
    from variant_annotation.lib.translation._core import _classify_kind, _Kind

    assert _classify_kind("NC_000001.11:g.12345C>T") is _Kind.GENOMIC


def test_classify_kind_unknown():
    from variant_annotation.lib.translation._core import _classify_kind, _Kind

    assert _classify_kind("not_an_hgvs_string") is _Kind.UNKNOWN


# ---------------------------------------------------------------------------
# _parse_protein_aa_change
# ---------------------------------------------------------------------------


def test_parse_protein_aa_change_missense():
    from variant_annotation.lib.translation._core import _parse_protein_aa_change

    assert _parse_protein_aa_change("NP_000001.1:p.Ala1Val") == ("Ala", 1, "Val")


def test_parse_protein_aa_change_synonymous_shorthand():
    from variant_annotation.lib.translation._core import _parse_protein_aa_change

    assert _parse_protein_aa_change("NP_000001.1:p.Ala5=") == ("Ala", 5, "Ala")


def test_parse_protein_aa_change_tolerates_prediction_parens():
    # c_to_p renders inferred consequences parenthesized; the WT codon path must still parse them.
    from variant_annotation.lib.translation._core import _parse_protein_aa_change

    assert _parse_protein_aa_change("NP_000050.3:p.(Phe2335=)") == ("Phe", 2335, "Phe")
    assert _parse_protein_aa_change("NP_000001.1:p.(Ala1Val)") == ("Ala", 1, "Val")


def test_parse_protein_aa_change_invalid():
    from variant_annotation.lib.translation._core import _parse_protein_aa_change

    assert _parse_protein_aa_change("not_hgvs") is None
    assert _parse_protein_aa_change("") is None


# ---------------------------------------------------------------------------
# _build_wt_c_hgvs
# ---------------------------------------------------------------------------


def test_build_wt_c_hgvs():
    from variant_annotation.lib.translation._core import _build_wt_c_hgvs

    assert _build_wt_c_hgvs("NM_000001.1", 1, "ATG") == "NM_000001.1:c.1_3delinsATG"
    assert _build_wt_c_hgvs("NM_000001.1", 2, "GAA") == "NM_000001.1:c.4_6delinsGAA"


# ---------------------------------------------------------------------------
# construct_equivalent_variants
# ---------------------------------------------------------------------------


def test_construct_equivalent_variants_empty_input():
    from variant_annotation.lib.translation._core import construct_equivalent_variants

    results, errors = construct_equivalent_variants([], transcripts=_StubTranscripts(), coordinates=_StubCoordinates())
    assert results == [] and errors == []


def test_construct_equivalent_variants_unknown_hgvs():
    from variant_annotation.lib.translation._core import construct_equivalent_variants

    results, errors = construct_equivalent_variants(
        [VariantInput(hgvs="not_valid")], transcripts=_StubTranscripts(), coordinates=_StubCoordinates()
    )
    assert results == [] and len(errors) == 1 and "Unrecognised" in errors[0].error


def test_construct_equivalent_variants_errors_when_subprocess_returns_no_rows(monkeypatch):
    """An input with a valid consequence must surface as an error — never vanish —
    when the subprocess exits cleanly but emits no output rows."""
    from variant_annotation.lib.translation import _core

    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    monkeypatch.setattr(_core, "_run_reverse_translate_batch", lambda cli, cons, *, config: ([], []))

    results, errors = _core.construct_equivalent_variants(
        [VariantInput(hgvs="NP_000001.1:p.Ala1Val", transcript="NM_000001.1")],
        transcripts=_StubTranscripts(),
        coordinates=_StubCoordinates(),
    )
    assert results == []
    assert len(errors) == 1
    assert "returned 0 rows for 1 inputs" in errors[0].error


# ---------------------------------------------------------------------------
# construct_one
# ---------------------------------------------------------------------------


def test_construct_one_returns_result_for_single_input(monkeypatch):
    from variant_annotation.lib.translation import _core

    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    monkeypatch.setattr(
        _core,
        "_run_reverse_translate_batch",
        lambda cli, cons, *, config: ([_core._BatchOutputRow(hgvs_c=["NM_000001.1:c.1A>G"], hgvs_g=[])], []),
    )

    result = _core.construct_one(
        VariantInput(hgvs="NP_000001.1:p.Ala1Val", transcript="NM_000001.1"),
        transcripts=_StubTranscripts(),
        coordinates=_StubCoordinates(),
    )
    assert isinstance(result, TranslationResult)
    assert result.hgvs_c_candidates == ["NM_000001.1:c.1A>G"]


def test_wt_codon_member_added_for_parenthesized_synonymous(monkeypatch):
    # Regression: c_to_p emits p.(Phe2335=) (parenthesized), and WtCodonMode.ALL must still
    # add the reference-identical WT codon delins -- previously the parens made it vanish.
    from variant_annotation.lib.translation import _core

    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    monkeypatch.setattr(
        _core,
        "_run_reverse_translate_batch",
        lambda cli, cons, *, config: ([_core._BatchOutputRow(hgvs_c=["NM_000059.4:c.7005T>C"], hgvs_g=[])], []),
    )

    class _Transcripts:
        def transcript_for_protein(self, acc):
            return None

        def codon_at(self, tx, pos):
            return "TTT"

    class _Coordinates:
        def c_to_p(self, c):
            raise NotImplementedError

        def g_to_c(self, g, tx):
            raise NotImplementedError

        def c_to_g(self, c):
            return "NC_000013.11:g.32346894delinsTTT"

    result = _core.construct_one(
        VariantInput(hgvs="NP_000050.3:p.(Phe2335=)", transcript="NM_000059.4"),
        transcripts=_Transcripts(),
        coordinates=_Coordinates(),
        config=TranslationConfig(wt_codon_mode=WtCodonMode.ALL, include_indels=True),
    )
    assert isinstance(result, TranslationResult)
    assert "NM_000059.4:c.7003_7005delinsTTT" in result.hgvs_c_candidates
    assert "NC_000013.11:g.32346894delinsTTT" in result.hgvs_g_candidates


def test_construct_one_returns_error_for_unresolvable_input():
    from variant_annotation.lib.translation._core import construct_one

    result = construct_one(
        VariantInput(hgvs="not_valid"), transcripts=_StubTranscripts(), coordinates=_StubCoordinates()
    )
    assert isinstance(result, TranslationError)
    assert "Unrecognised" in result.error
