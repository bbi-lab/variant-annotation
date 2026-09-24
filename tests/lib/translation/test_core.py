import pytest

from variant_annotation.lib.consequence import ProteinConsequence
from variant_annotation.lib.translation import _core
from variant_annotation.lib.translation._core import (
    _Kind,
    _apply_wt_codon,
    _build_wt_c_hgvs,
    _classify_kind,
    _parse_projection_pairs,
    _parse_protein_aa_change,
    _untranslatable_edit_reason,
    construct_equivalent_variants,
    construct_one,
)
from variant_annotation.lib.translation.types import (
    ProjectionPair,
    TranslationConfig,
    TranslationError,
    TranslationErrorReason,
    TranslationResult,
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
    assert _classify_kind("NP_000001.1:p.Ala1Val") is _Kind.PROTEIN


def test_classify_kind_coding():
    assert _classify_kind("NM_000001.1:c.368C>T") is _Kind.CODING


def test_classify_kind_genomic():
    assert _classify_kind("NC_000001.11:g.12345C>T") is _Kind.GENOMIC


def test_classify_kind_unknown():
    assert _classify_kind("not_an_hgvs_string") is _Kind.UNKNOWN


# ---------------------------------------------------------------------------
# _parse_protein_aa_change
# ---------------------------------------------------------------------------


def test_parse_protein_aa_change_missense():
    assert _parse_protein_aa_change("NP_000001.1:p.Ala1Val") == ("Ala", 1, "Val")


def test_parse_protein_aa_change_synonymous_shorthand():
    assert _parse_protein_aa_change("NP_000001.1:p.Ala5=") == ("Ala", 5, "Ala")


def test_parse_protein_aa_change_tolerates_prediction_parens():
    # c_to_p renders inferred consequences parenthesized; the WT codon path must still parse them.
    assert _parse_protein_aa_change("NP_000050.3:p.(Phe2335=)") == ("Phe", 2335, "Phe")
    assert _parse_protein_aa_change("NP_000001.1:p.(Ala1Val)") == ("Ala", 1, "Val")


def test_parse_protein_aa_change_invalid():
    assert _parse_protein_aa_change("not_hgvs") is None
    assert _parse_protein_aa_change("") is None


# ---------------------------------------------------------------------------
# _untranslatable_edit_reason — edit-type translatability boundary
# ---------------------------------------------------------------------------

# Representative forms the reverse-translate tool CAN express, and forms it cannot. The gate defers to
# the tool's own ``parse_hgvs_protein_change`` (see ``_untranslatable_edit_reason``), so these assert
# the classification the tool actually produces for each edit type — a deletion is included precisely
# because it is translatable, contrary to a "substitutions only" reading.
_TRANSLATABLE_EDITS = [
    "NP_000001.1:p.Ala1Val",  # missense
    "NP_000001.1:p.Arg175His",  # missense
    "NP_000001.1:p.Ala5=",  # synonymous
    "NP_000050.3:p.(Phe2335=)",  # synonymous, parenthesized prediction
    "NP_000001.1:p.Arg97Ter",  # nonsense (3-letter stop alt)
    "NP_000001.1:p.Arg97*",  # nonsense (1-letter stop alt)
    "p.A334D",  # 1-letter missense, no accession
    "NP_000001.1:p.Arg175del",  # single-residue deletion — expressly translatable
    "NP_000001.1:p.Arg175-",  # single-residue deletion, dash notation
]
_NON_TRANSLATABLE_EDITS = [
    "NP_000001.1:p.Gln39fs",  # frameshift
    "NP_000001.1:p.Arg97ThrfsTer5",  # frameshift, long form
    "NP_000001.1:p.Met1_Cys2insTrp",  # insertion
    "NP_000001.1:p.Arg175_Gly176del",  # multi-residue deletion
    "NP_000001.1:p.Arg175_Gly176delinsSer",  # multi-residue delins
    "NP_000001.1:p.Arg175dup",  # duplication
    "NP_000001.1:p.Ter810Argext*3",  # extension
    "NP_000001.1:p.Ter810Arg",  # stop-loss (3-letter ref stop)
    "NP_000001.1:p.*810Arg",  # stop-loss (1-letter ref stop)
    "NP_000001.1:p.Met1?",  # unknown effect
]


@pytest.mark.parametrize("hgvs_p", _TRANSLATABLE_EDITS)
def test_untranslatable_edit_reason_none_for_translatable_edits(hgvs_p):
    assert _untranslatable_edit_reason(hgvs_p) is None


@pytest.mark.parametrize("hgvs_p", _NON_TRANSLATABLE_EDITS)
def test_untranslatable_edit_reason_message_for_non_translatable_edits(hgvs_p):
    assert _untranslatable_edit_reason(hgvs_p) is not None


# ---------------------------------------------------------------------------
# _build_wt_c_hgvs
# ---------------------------------------------------------------------------


def test_build_wt_c_hgvs():
    assert _build_wt_c_hgvs("NM_000001.1", 1, "ATG") == "NM_000001.1:c.1_3delinsATG"
    assert _build_wt_c_hgvs("NM_000001.1", 2, "GAA") == "NM_000001.1:c.4_6delinsGAA"


# ---------------------------------------------------------------------------
# construct_equivalent_variants
# ---------------------------------------------------------------------------


def test_construct_equivalent_variants_empty_input():
    results, errors = construct_equivalent_variants([], transcripts=_StubTranscripts(), coordinates=_StubCoordinates())
    assert results == [] and errors == []


def test_construct_equivalent_variants_unknown_hgvs():
    results, errors = construct_equivalent_variants(
        [VariantInput(hgvs="not_valid")], transcripts=_StubTranscripts(), coordinates=_StubCoordinates()
    )
    assert results == [] and len(errors) == 1 and "Unrecognised" in errors[0].error
    # A consequence that cannot even be collapsed is a genuine failure, not a benign skip.
    assert errors[0].reason is TranslationErrorReason.FAILED


def test_construct_equivalent_variants_errors_when_subprocess_returns_no_rows(monkeypatch):
    """An input with a valid consequence must surface as an error — never vanish —
    when the subprocess exits cleanly but emits no output rows."""
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
    # A translatable consequence that yields no candidate is a genuine failure, not a benign skip.
    assert errors[0].reason is TranslationErrorReason.FAILED


def test_construct_equivalent_variants_non_translatable_edit_is_skipped_before_subprocess(monkeypatch):
    """A consequence the tool cannot express (here a frameshift) is screened out up front as a
    benign NOT_TRANSLATABLE skip — the reverse-translate subprocess is never invoked for it."""
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")

    def _fail_if_called(*args, **kwargs):
        raise AssertionError("subprocess must not run for a non-translatable input")

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _fail_if_called)

    results, errors = _core.construct_equivalent_variants(
        [VariantInput(hgvs="NP_000001.1:p.Gln39fs", transcript="NM_000001.1")],
        transcripts=_StubTranscripts(),
        coordinates=_StubCoordinates(),
    )
    assert results == []
    assert len(errors) == 1
    assert errors[0].reason is TranslationErrorReason.NOT_TRANSLATABLE
    assert errors[0].input.hgvs == "NP_000001.1:p.Gln39fs"


def test_construct_equivalent_variants_empty_result_for_translatable_edit_is_failed(monkeypatch):
    """An empty tool result for a *translatable* edit (a substitution) is a genuine FAILED, not a
    benign skip — the edit had a synonymous class, so producing nothing is a real failure."""
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    monkeypatch.setattr(
        _core,
        "_run_reverse_translate_batch",
        lambda cli, cons, *, config: ([_core._BatchOutputRow(projection_pairs=[])], []),
    )

    results, errors = _core.construct_equivalent_variants(
        [VariantInput(hgvs="NP_000001.1:p.Ala1Val", transcript="NM_000001.1")],
        transcripts=_StubTranscripts(),
        coordinates=_StubCoordinates(),
    )
    assert results == []
    assert len(errors) == 1
    assert errors[0].reason is TranslationErrorReason.FAILED


def test_construct_equivalent_variants_single_residue_deletion_the_tool_can_translate_succeeds(monkeypatch):
    """Regression: the tool DOES reverse-translate a single-residue deletion (into its codon
    deletion), so a deletion must never be pre-judged non-translatable and skipped. When the tool
    returns a candidate, it is a success — the classifier never gets to veto it."""
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    monkeypatch.setattr(
        _core,
        "_run_reverse_translate_batch",
        lambda cli, cons, *, config: (
            [_core._BatchOutputRow(projection_pairs=[ProjectionPair(hgvs_c="NM_000001.1:c.523_525del", hgvs_g=None)])],
            [],
        ),
    )

    result = _core.construct_one(
        VariantInput(hgvs="NP_000001.1:p.Arg175del", transcript="NM_000001.1"),
        transcripts=_StubTranscripts(),
        coordinates=_StubCoordinates(),
    )
    assert isinstance(result, TranslationResult)
    assert [c.hgvs_c for c in result.projection_pairs] == ["NM_000001.1:c.523_525del"]


def test_construct_equivalent_variants_mixes_translatable_and_not(monkeypatch):
    """A batch with one translatable and one non-translatable input feeds ONLY the translatable one
    to the tool; the non-translatable one is screened out as a NOT_TRANSLATABLE skip beside the
    result. Crucially, the tool receives exactly the translatable consequence — the screen never
    lets a non-translatable input consume a subprocess slot, and never withholds a translatable one."""
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")

    def _one_translatable(cli, cons, *, config):
        assert [c.hgvs_p for c in cons] == ["NP_000001.1:p.Ala1Val"], "only the translatable input reaches the tool"
        return (
            [_core._BatchOutputRow(projection_pairs=[ProjectionPair(hgvs_c="NM_000001.1:c.1A>G", hgvs_g=None)])],
            [],
        )

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _one_translatable)

    results, errors = _core.construct_equivalent_variants(
        [
            VariantInput(hgvs="NP_000001.1:p.Ala1Val", transcript="NM_000001.1"),
            VariantInput(hgvs="NP_000001.1:p.Met1_Cys2insTrp", transcript="NM_000001.1"),
        ],
        transcripts=_StubTranscripts(),
        coordinates=_StubCoordinates(),
    )
    assert len(results) == 1 and results[0].input.hgvs == "NP_000001.1:p.Ala1Val"
    assert len(errors) == 1
    assert errors[0].input.hgvs == "NP_000001.1:p.Met1_Cys2insTrp"
    assert errors[0].reason is TranslationErrorReason.NOT_TRANSLATABLE


# ---------------------------------------------------------------------------
# construct_one
# ---------------------------------------------------------------------------


def test_construct_one_returns_result_for_single_input(monkeypatch):
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
    assert _parse_projection_pairs({"variant_type": "", "hgvs_c": "", "hgvs_g": ""}) == []


# ---------------------------------------------------------------------------
# _apply_wt_codon — WT-genomic-fails preserves alignment by construction (A3)
# ---------------------------------------------------------------------------


def test_apply_wt_codon_genomic_failure_yields_paired_none_not_desync():
    # Regression for the latent desync: when the WT codon's c->g projection
    # fails, the coding candidate must still be appended with hgvs_g=None as a
    # single projection pair -- never a coding entry with no matching genomic
    # entry (which the old two-list append produced).
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


# ---------------------------------------------------------------------------
# transient UTA failures
# ---------------------------------------------------------------------------

_UTA_DROP = (
    'connection to server at "uta.biocommons.org" (35.95.195.14), port 5432 failed: '
    "server closed the connection unexpectedly"
)
_NO_WAIT = TranslationConfig(upstream_retry_backoff_seconds=0)


def _row(hgvs_c):
    return _core._BatchOutputRow(projection_pairs=[ProjectionPair(hgvs_c=hgvs_c, hgvs_g=None)])


def _empty_row():
    return _core._BatchOutputRow(projection_pairs=[])


def _two_inputs():
    return [
        VariantInput(hgvs="NP_000001.1:p.Ala1Val", transcript="NM_000001.1"),
        VariantInput(hgvs="NP_000001.1:p.Cys2Trp", transcript="NM_000001.1"),
    ]


@pytest.mark.parametrize(
    "message",
    [
        _UTA_DROP,
        "could not connect to server: Connection refused",
        'connection to server at "uta" port 5432 failed: timeout expired',
        "SSL SYSCALL error: EOF detected",
        "Failed reverse translation (terminating connection due to administrator command)",
    ],
)
def test_is_transient_upstream_error_recognises_connection_failures(message):
    assert _core._is_transient_upstream_error(message)


@pytest.mark.parametrize(
    "message",
    [
        'connection to server at "uta" port 5432 failed: FATAL:  password authentication failed for user "x"',
        'FATAL:  database "uta" does not exist',
        "Reference amino acid mismatch: HGVS p. requests A at position 1",
    ],
)
def test_is_transient_upstream_error_rejects_non_transient_failures(message):
    assert not _core._is_transient_upstream_error(message)


def test_subprocess_exit_on_uta_drop_is_retried_until_it_succeeds(monkeypatch):
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    calls = []

    def _flaky(cli, cons, *, config):
        calls.append(len(cons))
        if len(calls) == 1:
            raise RuntimeError(f"reverse-translate-variants failed: {_UTA_DROP}")
        return [_row("NM_000001.1:c.2C>T"), _row("NM_000001.1:c.5G>T")], []

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _flaky)

    results, errors = construct_equivalent_variants(
        _two_inputs(), transcripts=_StubTranscripts(), coordinates=_StubCoordinates(), config=_NO_WAIT
    )
    assert calls == [2, 2]
    assert errors == []
    assert [r.projection_pairs[0].hgvs_c for r in results] == ["NM_000001.1:c.2C>T", "NM_000001.1:c.5G>T"]


def test_subprocess_exit_on_uta_drop_exhausting_retries_is_upstream_unavailable(monkeypatch):
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    calls = []

    def _down(cli, cons, *, config):
        calls.append(len(cons))
        raise RuntimeError(f"reverse-translate-variants failed: {_UTA_DROP}")

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _down)

    results, errors = construct_equivalent_variants(
        _two_inputs(), transcripts=_StubTranscripts(), coordinates=_StubCoordinates(), config=_NO_WAIT
    )
    assert len(calls) == _NO_WAIT.upstream_max_attempts
    assert results == []
    assert [e.reason for e in errors] == [TranslationErrorReason.UPSTREAM_UNAVAILABLE] * 2


def test_non_transient_subprocess_exit_is_failed_without_retry(monkeypatch):
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    calls = []

    def _broken(cli, cons, *, config):
        calls.append(len(cons))
        raise RuntimeError("reverse-translate-variants failed: Unsupported --uniprot-target value")

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _broken)

    _, errors = construct_equivalent_variants(
        _two_inputs(), transcripts=_StubTranscripts(), coordinates=_StubCoordinates(), config=_NO_WAIT
    )
    assert calls == [2]
    assert [e.reason for e in errors] == [TranslationErrorReason.FAILED] * 2


def test_rows_dropped_mid_run_are_retried_alone_and_spliced_back(monkeypatch):
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    batches = []

    def _drop_second_row_once(cli, cons, *, config):
        batches.append([c.hgvs_p for c in cons])
        if len(batches) == 1:
            dropped = _core._BatchErrorRow(
                transcript="NM_000001.1", hgvs_p=cons[1].hgvs_p, error=f"Failed reverse translation ({_UTA_DROP})"
            )
            return [_row("NM_000001.1:c.2C>T"), _empty_row()], [dropped]
        return [_row("NM_000001.1:c.5G>T")], []

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _drop_second_row_once)

    results, errors = construct_equivalent_variants(
        _two_inputs(), transcripts=_StubTranscripts(), coordinates=_StubCoordinates(), config=_NO_WAIT
    )
    assert batches == [["NP_000001.1:p.Ala1Val", "NP_000001.1:p.Cys2Trp"], ["NP_000001.1:p.Cys2Trp"]]
    assert errors == []
    assert [(r.input.hgvs, r.projection_pairs[0].hgvs_c) for r in results] == [
        ("NP_000001.1:p.Ala1Val", "NM_000001.1:c.2C>T"),
        ("NP_000001.1:p.Cys2Trp", "NM_000001.1:c.5G>T"),
    ]


def test_rows_dropped_on_every_attempt_are_upstream_unavailable(monkeypatch):
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")

    def _always_drop_second(cli, cons, *, config):
        dropped = _core._BatchErrorRow(
            transcript="NM_000001.1", hgvs_p="NP_000001.1:p.Cys2Trp", error=f"Failed reverse translation ({_UTA_DROP})"
        )
        rows = [_empty_row() if c.hgvs_p == "NP_000001.1:p.Cys2Trp" else _row("NM_000001.1:c.2C>T") for c in cons]
        return rows, [dropped]

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _always_drop_second)

    results, errors = construct_equivalent_variants(
        _two_inputs(), transcripts=_StubTranscripts(), coordinates=_StubCoordinates(), config=_NO_WAIT
    )
    assert [r.input.hgvs for r in results] == ["NP_000001.1:p.Ala1Val"]
    assert len(errors) == 1
    assert errors[0].input.hgvs == "NP_000001.1:p.Cys2Trp"
    assert errors[0].reason is TranslationErrorReason.UPSTREAM_UNAVAILABLE


def test_failed_row_rerun_keeps_rows_from_the_earlier_run(monkeypatch):
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    calls = []

    def _rerun_crashes(cli, cons, *, config):
        calls.append(len(cons))
        if len(calls) == 1:
            dropped = _core._BatchErrorRow(
                transcript="NM_000001.1", hgvs_p=cons[1].hgvs_p, error=f"Failed reverse translation ({_UTA_DROP})"
            )
            return [_row("NM_000001.1:c.2C>T"), _empty_row()], [dropped]
        raise RuntimeError("reverse-translate-variants failed: unexpected crash")

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _rerun_crashes)

    results, errors = construct_equivalent_variants(
        _two_inputs(), transcripts=_StubTranscripts(), coordinates=_StubCoordinates(), config=_NO_WAIT
    )
    assert calls == [2, 1]
    assert [r.input.hgvs for r in results] == ["NP_000001.1:p.Ala1Val"]
    assert [e.reason for e in errors] == [TranslationErrorReason.UPSTREAM_UNAVAILABLE]


def test_empty_row_with_non_transient_error_is_not_retried(monkeypatch):
    monkeypatch.setattr(_core, "_find_reverse_translate_cli", lambda: "/bin/true")
    calls = []

    def _ref_mismatch(cli, cons, *, config):
        calls.append(len(cons))
        mismatch = _core._BatchErrorRow(
            transcript="NM_000001.1", hgvs_p=cons[0].hgvs_p, error="Reference amino acid mismatch"
        )
        return [_empty_row()], [mismatch]

    monkeypatch.setattr(_core, "_run_reverse_translate_batch", _ref_mismatch)

    _, errors = construct_equivalent_variants(
        _two_inputs()[:1], transcripts=_StubTranscripts(), coordinates=_StubCoordinates(), config=_NO_WAIT
    )
    assert calls == [1]
    assert errors[0].reason is TranslationErrorReason.FAILED
