import pytest

from variant_annotation.lib.translation.types import (
    ProjectionPair,
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
        lambda cli, cons, *, config: (
            [_core._BatchOutputRow(projection_pairs=[ProjectionPair(hgvs_c="NM_000001.1:c.1A>G", hgvs_g=None)])],
            [],
        ),
    )

    result = _core.construct_one(
        VariantInput(hgvs="NP_000001.1:p.Ala1Val", transcript="NM_000001.1"),
        transcripts=_StubTranscripts(),
        coordinates=_StubCoordinates(),
    )
    assert isinstance(result, TranslationResult)
    assert [c.hgvs_c for c in result.projection_pairs] == ["NM_000001.1:c.1A>G"]


def test_wt_codon_member_added_for_parenthesized_synonymous(monkeypatch):
    # Regression: c_to_p emits p.(Phe2335=) (parenthesized), and WtCodonMode.ALL must still
    # add the reference-identical WT codon delins -- previously the parens made it vanish.
    from variant_annotation.lib.translation import _core

    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    monkeypatch.setattr(
        _core,
        "_run_reverse_translate_batch",
        lambda cli, cons, *, config: (
            [_core._BatchOutputRow(projection_pairs=[ProjectionPair(hgvs_c="NM_000059.4:c.7005T>C", hgvs_g=None)])],
            [],
        ),
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
    # The WT codon rides in as a whole projection pair: its coding and genomic
    # expressions live on the same object (never on two independently-appended lists).
    wt = next(c for c in result.projection_pairs if c.hgvs_c == "NM_000059.4:c.7003_7005delinsTTT")
    assert wt.hgvs_g == "NC_000013.11:g.32346894delinsTTT"
    assert wt.variant_type == "delins"


def test_construct_one_returns_error_for_unresolvable_input():
    from variant_annotation.lib.translation._core import construct_one

    result = construct_one(
        VariantInput(hgvs="not_valid"), transcripts=_StubTranscripts(), coordinates=_StubCoordinates()
    )
    assert isinstance(result, TranslationError)
    assert "Unrecognised" in result.error


# ---------------------------------------------------------------------------
# _parse_projection_pairs — row-wise, position-aligned parse (A2)
# ---------------------------------------------------------------------------


def test_parse_projection_pairs_degenerate_consequence_yields_n_pairs():
    # A degenerate protein consequence fans out into N coding candidates, each
    # paired with its genomic projection and variant type, in CLI column order.
    from variant_annotation.lib.translation._core import _parse_projection_pairs

    row = {
        "variant_type": "snv|snv|snv",
        "hgvs_c": "NM_1:c.1A>G|NM_1:c.1A>C|NM_1:c.1A>T",
        "hgvs_g": "NC_1:g.1A>G|NC_1:g.1A>C|NC_1:g.1A>T",
    }
    candidates = _parse_projection_pairs(row)

    assert [c.hgvs_c for c in candidates] == ["NM_1:c.1A>G", "NM_1:c.1A>C", "NM_1:c.1A>T"]
    assert [c.hgvs_g for c in candidates] == ["NC_1:g.1A>G", "NC_1:g.1A>C", "NC_1:g.1A>T"]
    assert all(c.variant_type == "snv" for c in candidates)


def test_parse_projection_pairs_projection_failed_is_none_not_desync():
    # The bug this refactor fixes: a candidate whose genomic projection failed
    # leaves an empty g cell. Under the old filter-empties split the g list went
    # short and every later candidate paired with the wrong g. Row-wise parsing
    # keeps the coding candidate with hgvs_g=None and preserves alignment of the
    # candidates on either side.
    from variant_annotation.lib.translation._core import _parse_projection_pairs

    row = {
        "variant_type": "snv|snv|snv",
        "hgvs_c": "NM_1:c.1A>G|NM_1:c.1A>C|NM_1:c.1A>T",
        "hgvs_g": "NC_1:g.1A>G||NC_1:g.1A>T",  # middle candidate's projection failed
    }
    candidates = _parse_projection_pairs(row)

    assert len(candidates) == 3
    assert candidates[0].hgvs_c == "NM_1:c.1A>G" and candidates[0].hgvs_g == "NC_1:g.1A>G"
    assert candidates[1].hgvs_c == "NM_1:c.1A>C" and candidates[1].hgvs_g is None
    assert candidates[2].hgvs_c == "NM_1:c.1A>T" and candidates[2].hgvs_g == "NC_1:g.1A>T"


def test_parse_projection_pairs_no_candidates_is_empty():
    from variant_annotation.lib.translation._core import _parse_projection_pairs

    assert _parse_projection_pairs({"variant_type": "", "hgvs_c": "", "hgvs_g": ""}) == []


# ---------------------------------------------------------------------------
# _apply_wt_codon — WT-genomic-fails preserves alignment by construction (A3)
# ---------------------------------------------------------------------------


def test_apply_wt_codon_genomic_failure_yields_paired_none_not_desync():
    # Regression for the latent desync: when the WT codon's c->g projection
    # fails, the coding candidate must still be appended with hgvs_g=None as a
    # single projection pair -- never a coding entry with no matching genomic
    # entry (which the old two-list append produced).
    from variant_annotation.lib.translation._core import _apply_wt_codon
    from variant_annotation.lib.consequence import ProteinConsequence

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
            raise RuntimeError("intronic projection unavailable")

    existing = [ProjectionPair(hgvs_c="NM_000059.4:c.7005T>C", hgvs_g="NC_000013.11:g.32346894T>C")]
    out = _apply_wt_codon(
        existing,
        ProteinConsequence(hgvs_p="NP_000050.3:p.(Phe2335=)", transcript="NM_000059.4"),
        config=TranslationConfig(wt_codon_mode=WtCodonMode.ALL, include_indels=True),
        transcripts=_Transcripts(),
        coordinates=_Coordinates(),
    )

    assert len(out) == 2
    wt = out[-1]
    assert wt.hgvs_c == "NM_000059.4:c.7003_7005delinsTTT"
    assert wt.hgvs_g is None  # projection failed, but the pair object is intact
    assert wt.variant_type == "delins"


def test_construct_equivalent_variants_end_to_end_returns_projection_pairs(monkeypatch):
    # End-to-end (subprocess stubbed): the equivalence class arrives as
    # ProjectionPair objects, projection-failed members carrying hgvs_g=None.
    from variant_annotation.lib.translation import _core

    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    monkeypatch.setattr(
        _core,
        "_run_reverse_translate_batch",
        lambda cli, cons, *, config: (
            [
                _core._BatchOutputRow(
                    projection_pairs=[
                        ProjectionPair(hgvs_c="NM_1:c.3G>A", hgvs_g="NC_1:g.3G>A", variant_type="snv"),
                        ProjectionPair(hgvs_c="NM_1:c.3G>T", hgvs_g=None, variant_type="snv"),
                    ]
                )
            ],
            [],
        ),
    )

    results, errors = _core.construct_equivalent_variants(
        [VariantInput(hgvs="NP_1:p.Trp1Cys", transcript="NM_1")],
        transcripts=_StubTranscripts(),
        coordinates=_StubCoordinates(),
    )
    assert errors == []
    assert len(results) == 1
    candidates = results[0].projection_pairs
    assert [(c.hgvs_c, c.hgvs_g) for c in candidates] == [
        ("NM_1:c.3G>A", "NC_1:g.3G>A"),
        ("NM_1:c.3G>T", None),
    ]
    # Protein input -> protein apex is authoritative, not re-emitted as a derived member.
    assert results[0].hgvs_p is None
