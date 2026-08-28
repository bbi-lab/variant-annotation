"""Tests for src/annotate_vep.py — the CSV/CLI composition root.

Consequence *resolution* is the library's job and is tested in tests/lib/vep/. What is left here is
what this script actually owns: HGVS column priority, pipe alignment across candidates, the Redis cache
adapter, and end-to-end CLI wiring.
"""

from __future__ import annotations

import csv
from datetime import date

import pytest

import annotate_vep as mod
from variant_annotation.lib.vep import (
    NO_CHANGE_TERM,
    ConsequenceOutcome,
    ConsequenceResolution,
    ConsequenceSource,
    VepInput,
)


def resolved(hgvs, most_severe, *, terms=None, source=ConsequenceSource.MOST_SEVERE, transcript=None):
    return ConsequenceResolution(
        input=VepInput(hgvs),
        outcome=ConsequenceOutcome.RESOLVED,
        consequence_terms=list(terms) if terms is not None else [most_severe],
        most_severe_consequence=most_severe,
        source=source,
        matched_transcript=transcript,
    )


def absent(hgvs):
    return ConsequenceResolution(input=VepInput(hgvs), outcome=ConsequenceOutcome.ABSENT)


def errored(hgvs, message="VEP request failed: boom"):
    return ConsequenceResolution(input=VepInput(hgvs), outcome=ConsequenceOutcome.ERRORED, error=message)


def annotate(
    row, resolutions, *, hgvs_cols=("mapped_hgvs_c", "mapped_hgvs_g", "mapped_hgvs_p"), access_date="2026-04-30"
):
    return mod.annotate_row(
        row,
        {r.input.hgvs: r for r in resolutions},
        col_prefix="vep",
        hgvs_cols=list(hgvs_cols),
        access_date=access_date,
    )


# --- HGVS column priority ---


@pytest.mark.unit
def test_annotate_row_uses_first_non_blank_hgvs_col():
    """mapped_hgvs_c takes priority over mapped_hgvs_g when non-blank."""
    row = {"mapped_hgvs_c": "NM_000000.1:c.1A>T", "mapped_hgvs_g": "NC_000001.11:g.1A>T"}
    out = annotate(
        row,
        [
            resolved("NM_000000.1:c.1A>T", "synonymous_variant"),
            resolved("NC_000001.11:g.1A>T", "missense_variant"),
        ],
    )
    assert out["vep.most_severe_mutational_consequence"] == "synonymous_variant"


@pytest.mark.unit
def test_annotate_row_falls_back_when_first_col_blank():
    row = {"mapped_hgvs_c": "", "mapped_hgvs_g": "NC_000001.11:g.1A>T"}
    out = annotate(row, [resolved("NC_000001.11:g.1A>T", "missense_variant")])
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant"


@pytest.mark.unit
def test_annotate_row_with_no_hgvs_at_all_writes_blank_columns():
    out = annotate({"mapped_hgvs_c": "", "mapped_hgvs_g": ""}, [])
    assert out == {
        "vep.mutational_consequences": "",
        "vep.most_severe_mutational_consequence": "",
        "vep.access_date": "",
        "vep.consequence_source": "",
        "vep.error": "",
    }


# --- pipe alignment ---


@pytest.mark.unit
def test_annotate_row_pipe_aligns_matching_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    out = annotate(
        row,
        [
            resolved("NC_000001.11:g.1A>T", "missense_variant"),
            resolved("NC_000001.11:g.2C>G", "missense_variant"),
        ],
        hgvs_cols=["mapped_hgvs_g"],
    )
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant|missense_variant"


@pytest.mark.unit
def test_annotate_row_pipe_aligns_discordant_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    out = annotate(
        row,
        [
            resolved("NC_000001.11:g.1A>T", "missense_variant"),
            resolved("NC_000001.11:g.2C>G", "synonymous_variant"),
        ],
        hgvs_cols=["mapped_hgvs_g"],
    )
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant|synonymous_variant"


@pytest.mark.unit
def test_annotate_row_preserves_empty_slot_for_unresolved_candidate():
    """Alignment is positional, so a candidate with no answer must still occupy its slot."""
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    out = annotate(
        row,
        [resolved("NC_000001.11:g.1A>T", "missense_variant"), absent("NC_000001.11:g.2C>G")],
        hgvs_cols=["mapped_hgvs_g"],
    )
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant|"


@pytest.mark.unit
def test_access_date_is_pipe_aligned_for_multiple_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    out = annotate(
        row,
        [
            resolved("NC_000001.11:g.1A>T", "missense_variant"),
            resolved("NC_000001.11:g.2C>G", "synonymous_variant"),
        ],
        hgvs_cols=["mapped_hgvs_g"],
        access_date="2026-05-29",
    )
    assert out["vep.access_date"] == "2026-05-29|2026-05-29"


@pytest.mark.unit
def test_access_date_preserves_empty_candidate_slots():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T||NC_000001.11:g.3G>A"}
    out = annotate(
        row,
        [
            resolved("NC_000001.11:g.1A>T", "missense_variant"),
            resolved("NC_000001.11:g.3G>A", "intron_variant"),
        ],
        hgvs_cols=["mapped_hgvs_g"],
        access_date="2026-05-29",
    )
    assert out["vep.access_date"] == "2026-05-29||2026-05-29"


# --- how each outcome is written ---


@pytest.mark.unit
def test_all_transcript_terms_are_caret_joined_severity_ordered():
    row = {"mapped_hgvs_c": "NM_000049.4:c.256A>G"}
    out = annotate(
        row,
        [
            resolved(
                "NM_000049.4:c.256A>G",
                "splice_region_variant",
                terms=["splice_region_variant", "synonymous_variant"],
                source=ConsequenceSource.TRANSCRIPT,
                transcript="NM_000049.4",
            )
        ],
    )
    assert out["vep.mutational_consequences"] == "splice_region_variant^synonymous_variant"
    assert out["vep.most_severe_mutational_consequence"] == "splice_region_variant"
    assert out["vep.consequence_source"] == "transcript"


@pytest.mark.unit
def test_an_errored_candidate_writes_the_error_and_no_consequence():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T"}
    out = annotate(row, [errored("NC_000001.11:g.1A>T")], hgvs_cols=["mapped_hgvs_g"])

    assert out["vep.most_severe_mutational_consequence"] == ""
    assert out["vep.mutational_consequences"] == ""
    assert "boom" in out["vep.error"]


@pytest.mark.unit
def test_an_absent_candidate_writes_no_error():
    """A confirmed empty is not a failure — putting it in the error column would misreport it."""
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T"}
    out = annotate(row, [absent("NC_000001.11:g.1A>T")], hgvs_cols=["mapped_hgvs_g"])

    assert out["vep.most_severe_mutational_consequence"] == ""
    assert out["vep.error"] == ""
    assert out["vep.access_date"] == "2026-04-30"


@pytest.mark.unit
def test_a_reference_identical_candidate_is_labelled_no_change():
    row = {"mapped_hgvs_c": "NM_000049.4:c.123="}
    out = annotate(
        row,
        [
            resolved(
                "NM_000049.4:c.123=",
                NO_CHANGE_TERM,
                source=ConsequenceSource.REFERENCE_IDENTICAL,
            )
        ],
    )
    assert out["vep.most_severe_mutational_consequence"] == NO_CHANGE_TERM
    assert out["vep.consequence_source"] == "reference_identical"


@pytest.mark.unit
def test_error_text_cannot_break_the_delimited_columns():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T"}
    out = annotate(
        row,
        [errored("NC_000001.11:g.1A>T", "bad|thing\twith\nnewlines")],
        hgvs_cols=["mapped_hgvs_g"],
    )
    assert "|" not in out["vep.error"]
    assert "\t" not in out["vep.error"]
    assert "\n" not in out["vep.error"]


# --- Redis cache adapter ---


@pytest.mark.unit
def test_cache_round_trips_a_resolved_answer():
    original = resolved(
        "NM_000049.4:c.256A>G",
        "synonymous_variant",
        terms=["splice_region_variant", "synonymous_variant"],
        source=ConsequenceSource.TRANSCRIPT,
        transcript="NM_000049.4",
    )
    value = mod._resolution_to_cache_value(original)
    rebuilt = mod._cache_value_to_resolution(original.input, value)

    assert rebuilt.outcome is ConsequenceOutcome.RESOLVED
    assert rebuilt.consequence_terms == original.consequence_terms
    assert rebuilt.most_severe_consequence == original.most_severe_consequence
    assert rebuilt.source is ConsequenceSource.TRANSCRIPT
    assert rebuilt.matched_transcript == "NM_000049.4"


@pytest.mark.unit
def test_cache_round_trips_a_confirmed_absent_answer():
    original = absent("NC_000001.11:g.1A>T")
    value = mod._resolution_to_cache_value(original)
    assert mod._cache_value_to_resolution(original.input, value).outcome is ConsequenceOutcome.ABSENT


@pytest.mark.unit
def test_a_failed_request_is_never_cached():
    """Caching an outage would turn a transient failure into a persistent wrong answer."""
    assert mod._resolution_to_cache_value(errored("NC_000001.11:g.1A>T")) is None


@pytest.mark.unit
def test_unreadable_cache_entries_are_discarded_rather_than_trusted():
    assert mod._cache_value_to_resolution(VepInput("NM_1.1:c.1A>G"), "not-json") is None
    assert mod._cache_value_to_resolution(VepInput("NM_1.1:c.1A>G"), '{"s": "bogus_source"}') is None


@pytest.mark.unit
def test_cache_key_carries_the_resolver_version_and_transcript():
    """A resolution-rule change must invalidate stored answers without operator action."""
    from variant_annotation.lib.vep import RESOLVER_VERSION

    key = mod._vep_cache_key(VepInput("NM_000049.4:c.256A>G"), "116")

    assert f"r{RESOLVER_VERSION}" in key
    assert "NM_000049" in key


@pytest.mark.unit
def test_cache_key_distinguishes_the_same_hgvs_resolved_against_different_transcripts():
    first = mod._vep_cache_key(VepInput("NC_1:g.100A>G", transcript="NM_1.1"), "116")
    second = mod._vep_cache_key(VepInput("NC_1:g.100A>G", transcript="NM_2.1"), "116")
    assert first != second


@pytest.mark.unit
def test_cache_key_distinguishes_ensembl_releases():
    """A release boundary can change the answer for an unchanged input, so an answer resolved under one
    release must not be served under another — the key carries the release, not just the resolver rule."""
    vep_input = VepInput("NM_000049.4:c.256A>G")
    assert mod._vep_cache_key(vep_input, "116") != mod._vep_cache_key(vep_input, "117")


# --- CLI end to end ---


def fake_resolver(answers, *, record=None):
    """Stand in for library resolution, answering from a hgvs -> consequence map."""

    def _resolve(inputs, *, vep, recoder=None, reference=None, config=None):
        if record is not None:
            record.extend(i.hgvs for i in inputs)
        return [resolved(i.hgvs, answers[i.hgvs]) if i.hgvs in answers else absent(i.hgvs) for i in inputs]

    return _resolve


@pytest.mark.integration
def test_main_applies_skip_and_limit(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["variant_urn", "mapped_hgvs_g"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({"variant_urn": "v1", "mapped_hgvs_g": "NC_000001.11:g.1A>T"})
        writer.writerow({"variant_urn": "v2", "mapped_hgvs_g": "NC_000001.11:g.2C>G"})
        writer.writerow({"variant_urn": "v3", "mapped_hgvs_g": "NC_000001.11:g.3G>A"})

    monkeypatch.setenv("VEP_CACHE_ENABLED", "0")
    monkeypatch.setattr(
        mod,
        "resolve_consequences",
        fake_resolver(
            {
                "NC_000001.11:g.1A>T": "synonymous_variant",
                "NC_000001.11:g.2C>G": "missense_variant",
                "NC_000001.11:g.3G>A": "synonymous_variant",
            }
        ),
    )

    mod.main([str(in_path), str(out_path), "--skip", "1", "--limit", "1", "--row-batch-size", "1"])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert [r["variant_urn"] for r in rows] == ["v2"]
    assert rows[0]["vep.most_severe_mutational_consequence"] == "missense_variant"
    assert rows[0]["vep.error"] == ""
    assert rows[0]["vep.access_date"] == date.today().isoformat()


@pytest.mark.integration
def test_main_keep_existing_preserves_annotated_fills_blank(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    today = date.today().isoformat()

    fieldnames = [
        "variant_urn",
        "mapped_hgvs_g",
        "vep.most_severe_mutational_consequence",
        "vep.access_date",
        "vep.error",
    ]
    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "variant_urn": "v1",
                "mapped_hgvs_g": "NC_000001.11:g.1A>T",
                "vep.most_severe_mutational_consequence": "synonymous_variant",
                "vep.access_date": "2025-01-01",
                "vep.error": "",
            }
        )
        writer.writerow(
            {
                "variant_urn": "v2",
                "mapped_hgvs_g": "NC_000001.11:g.2C>G",
                "vep.most_severe_mutational_consequence": "",
                "vep.access_date": "",
                "vep.error": "",
            }
        )

    looked_up: list[str] = []
    monkeypatch.setenv("VEP_CACHE_ENABLED", "0")
    monkeypatch.setattr(
        mod,
        "resolve_consequences",
        fake_resolver({"NC_000001.11:g.2C>G": "missense_variant"}, record=looked_up),
    )

    mod.main([str(in_path), str(out_path), "--keep-existing"])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert looked_up == ["NC_000001.11:g.2C>G"]

    assert rows[0]["variant_urn"] == "v1"
    assert rows[0]["vep.most_severe_mutational_consequence"] == "synonymous_variant"
    assert rows[0]["vep.access_date"] == "2025-01-01"

    assert rows[1]["variant_urn"] == "v2"
    assert rows[1]["vep.most_severe_mutational_consequence"] == "missense_variant"
    assert rows[1]["vep.access_date"] == today


@pytest.mark.integration
def test_main_keep_existing_makes_no_api_calls_when_all_annotated(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    fieldnames = [
        "variant_urn",
        "mapped_hgvs_g",
        "vep.most_severe_mutational_consequence",
        "vep.access_date",
        "vep.error",
    ]
    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "variant_urn": "v1",
                "mapped_hgvs_g": "NC_000001.11:g.1A>T",
                "vep.most_severe_mutational_consequence": "synonymous_variant",
                "vep.access_date": "2025-01-01",
                "vep.error": "",
            }
        )

    calls: list[int] = []

    def exploding_resolver(inputs, **kwargs):
        calls.append(len(inputs))
        return []

    monkeypatch.setenv("VEP_CACHE_ENABLED", "0")
    monkeypatch.setattr(mod, "resolve_consequences", exploding_resolver)

    mod.main([str(in_path), str(out_path), "--keep-existing"])

    assert calls == [], "resolution should not run when every row is already annotated"

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert rows[0]["vep.most_severe_mutational_consequence"] == "synonymous_variant"
    assert rows[0]["vep.access_date"] == "2025-01-01"


@pytest.mark.integration
def test_main_rejects_a_batch_size_above_ensembls_limit(tmp_path):
    in_path = tmp_path / "in.tsv"
    in_path.write_text("variant_urn\tmapped_hgvs_g\nv1\tNC_000001.11:g.1A>T\n", encoding="utf-8")

    with pytest.raises(SystemExit):
        mod.main([str(in_path), str(tmp_path / "out.tsv"), "--vep-batch-size", "500"])
