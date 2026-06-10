"""Tests for src/annotate_vep.py."""

from __future__ import annotations

import csv
from datetime import date

import pytest

import annotate_vep as mod


@pytest.mark.unit
def test_annotate_row_uses_first_non_blank_hgvs_col():
    """mapped_hgvs_c takes priority over mapped_hgvs_g when non-blank."""
    row = {
        "mapped_hgvs_c": "NM_000000.1:c.1A>T",
        "mapped_hgvs_g": "NC_000001.11:g.1A>T",
    }
    consequence_cache = {
        "NM_000000.1:c.1A>T": ("synonymous_variant", None, "most_severe"),
        "NC_000001.11:g.1A>T": ("missense_variant", None, "most_severe"),
    }
    out = mod.annotate_row(
        row,
        consequence_cache,
        col_prefix="vep",
        hgvs_cols=["mapped_hgvs_c", "mapped_hgvs_g", "mapped_hgvs_p"],
        access_date="2026-04-30",
    )
    assert out["vep.most_severe_mutational_consequence"] == "synonymous_variant"


@pytest.mark.unit
def test_annotate_row_falls_back_when_first_col_blank():
    """Falls back to mapped_hgvs_g when mapped_hgvs_c is blank."""
    row = {
        "mapped_hgvs_c": "",
        "mapped_hgvs_g": "NC_000001.11:g.1A>T",
    }
    consequence_cache = {"NC_000001.11:g.1A>T": ("missense_variant", None, "most_severe")}
    out = mod.annotate_row(
        row,
        consequence_cache,
        col_prefix="vep",
        hgvs_cols=["mapped_hgvs_c", "mapped_hgvs_g", "mapped_hgvs_p"],
        access_date="2026-04-30",
    )
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant"


@pytest.mark.unit
def test_annotate_row_pipe_aligns_matching_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    consequence_cache = {
        "NC_000001.11:g.1A>T": ("missense_variant", None, "most_severe"),
        "NC_000001.11:g.2C>G": ("missense_variant", None, "most_severe"),
    }
    out = mod.annotate_row(
        row,
        consequence_cache,
        col_prefix="vep",
        hgvs_cols=["mapped_hgvs_g"],
        access_date="2026-04-30",
    )

    assert out["vep.mutational_consequences"] == "missense_variant|missense_variant"
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant|missense_variant"
    assert out["vep.access_date"] == "2026-04-30|2026-04-30"
    assert out["vep.error"] == "|"


@pytest.mark.unit
def test_annotate_row_preserves_empty_slot_for_unresolved_candidate():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    consequence_cache = {
        "NC_000001.11:g.1A>T": ("synonymous_variant", None, "most_severe"),
        # second candidate unresolved by VEP
        "NC_000001.11:g.2C>G": None,
    }
    out = mod.annotate_row(
        row,
        consequence_cache,
        col_prefix="vep",
        hgvs_cols=["mapped_hgvs_g"],
        access_date="2026-04-30",
    )

    assert out["vep.most_severe_mutational_consequence"] == "synonymous_variant|"
    assert out["vep.consequence_source"] == "most_severe|"
    assert out["vep.error"] == "|"


@pytest.mark.unit
def test_annotate_row_pipe_aligns_discordant_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G|NC_000001.11:g.3G>A"}
    consequence_cache = {
        "NC_000001.11:g.1A>T": ("missense_variant", None, "most_severe"),
        "NC_000001.11:g.2C>G": ("synonymous_variant", None, "most_severe"),
        "NC_000001.11:g.3G>A": None,
    }
    out = mod.annotate_row(
        row,
        consequence_cache,
        col_prefix="vep",
        hgvs_cols=["mapped_hgvs_g"],
        access_date="2026-04-30",
    )

    assert out["vep.mutational_consequences"] == "missense_variant|synonymous_variant|"
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant|synonymous_variant|"


@pytest.mark.unit
def test_get_functional_consequence_uses_recoder_fallback(monkeypatch):
    # VEP resolves the genomic HGVS directly but not the transcript HGVS, which
    # must fall back to Variant Recoder -> genomic -> VEP. (Transcript HGVS are
    # sent in their own refseq=1 batch, so VEP batches are keyed per-HGVS here
    # rather than by exact batch composition.)
    vep_calls: list[tuple[tuple[str, ...], bool]] = []
    recoder_calls: list[list[str]] = []

    resolved = {
        "NC_000001.11:g.1A>T": ("intron_variant", None, "most_severe"),
        "NC_000001.11:g.5A>T": ("synonymous_variant", None, "most_severe"),
        "NC_000001.11:g.6A>T": ("missense_variant", None, "most_severe"),
    }

    def fake_vep_lookup(hgvs_strings, *, api_url, timeout_seconds, refseq=False):
        vep_calls.append((tuple(hgvs_strings), refseq))
        return {h: resolved[h] for h in hgvs_strings if h in resolved}, None

    def fake_recoder(hgvs_strings, *, api_url, timeout_seconds):
        recoder_calls.append(list(hgvs_strings))
        return {"NM_000000.1:c.1A>T": ["NC_000001.11:g.5A>T", "NC_000001.11:g.6A>T"]}, None

    monkeypatch.setattr(mod, "_vep_lookup_batch", fake_vep_lookup)
    monkeypatch.setattr(mod, "run_variant_recoder", fake_recoder)

    out = mod.get_functional_consequence(
        ["NC_000001.11:g.1A>T", "NM_000000.1:c.1A>T"],
        api_url="https://rest.ensembl.org",
        timeout_seconds=10,
        batch_size=200,
    )

    # The transcript HGVS that VEP could not resolve is the one recoded.
    assert recoder_calls == [["NM_000000.1:c.1A>T"]]
    # The recoded genomic HGVS were looked up, and the transcript batch used refseq=1.
    assert (("NC_000001.11:g.5A>T", "NC_000001.11:g.6A>T"), False) in vep_calls
    assert (("NM_000000.1:c.1A>T",), True) in vep_calls

    assert out["NC_000001.11:g.1A>T"] == ("intron_variant", None, "most_severe")
    # most severe of synonymous vs missense should be missense
    assert out["NM_000000.1:c.1A>T"] == ("missense_variant", None, "most_severe")


@pytest.mark.integration
def test_main_applies_skip_and_limit(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "mapped_hgvs_g"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v1", "mapped_hgvs_g": "NC_000001.11:g.1A>T"})
        writer.writerow({"variant_urn": "v2", "mapped_hgvs_g": "NC_000001.11:g.2C>G"})
        writer.writerow({"variant_urn": "v3", "mapped_hgvs_g": "NC_000001.11:g.3G>A"})

    def fake_get_functional_consequence(hgvs_strings, *, api_url, timeout_seconds, batch_size, max_workers=1):
        return {
            hgvs: (
                ("missense_variant" if hgvs.endswith("2C>G") else "synonymous_variant"),
                None,
                "most_severe",
            )
            for hgvs in dict.fromkeys(hgvs_strings)
        }

    monkeypatch.setattr(mod, "get_functional_consequence", fake_get_functional_consequence)

    mod.main([
        str(in_path),
        str(out_path),
        "--skip",
        "1",
        "--limit",
        "1",
        "--row-batch-size",
        "1",
    ])

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
        # Row that already has an annotation — should be preserved unchanged.
        writer.writerow({
            "variant_urn": "v1",
            "mapped_hgvs_g": "NC_000001.11:g.1A>T",
            "vep.most_severe_mutational_consequence": "synonymous_variant",
            "vep.access_date": "2025-01-01",
            "vep.error": "",
        })
        # Row with blank annotation — should be newly annotated.
        writer.writerow({
            "variant_urn": "v2",
            "mapped_hgvs_g": "NC_000001.11:g.2C>G",
            "vep.most_severe_mutational_consequence": "",
            "vep.access_date": "",
            "vep.error": "",
        })

    looked_up: list[str] = []

    def fake_get_functional_consequence(hgvs_strings, *, api_url, timeout_seconds, batch_size, max_workers=1):
        looked_up.extend(hgvs_strings)
        return {h: ("missense_variant", None, "most_severe") for h in hgvs_strings}

    monkeypatch.setattr(mod, "get_functional_consequence", fake_get_functional_consequence)

    mod.main([str(in_path), str(out_path), "--keep-existing"])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    # Only v2's HGVS should have been looked up.
    assert looked_up == ["NC_000001.11:g.2C>G"]

    # v1 is unchanged.
    assert rows[0]["variant_urn"] == "v1"
    assert rows[0]["vep.most_severe_mutational_consequence"] == "synonymous_variant"
    assert rows[0]["vep.access_date"] == "2025-01-01"

    # v2 got a new annotation.
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
        writer.writerow({
            "variant_urn": "v1",
            "mapped_hgvs_g": "NC_000001.11:g.1A>T",
            "vep.most_severe_mutational_consequence": "synonymous_variant",
            "vep.access_date": "2025-01-01",
            "vep.error": "",
        })

    call_count: list[int] = []

    def fake_get_functional_consequence(hgvs_strings, *, api_url, timeout_seconds, batch_size, max_workers=1):
        call_count.append(1)
        return {}

    monkeypatch.setattr(mod, "get_functional_consequence", fake_get_functional_consequence)

    mod.main([str(in_path), str(out_path), "--keep-existing"])

    assert call_count == [], "get_functional_consequence should not be called when all rows already annotated"

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert rows[0]["vep.most_severe_mutational_consequence"] == "synonymous_variant"
    assert rows[0]["vep.access_date"] == "2025-01-01"


@pytest.mark.unit
def test_access_date_is_pipe_aligned_for_multiple_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    cache = {
        "NC_000001.11:g.1A>T": ("missense_variant", None, "most_severe"),
        "NC_000001.11:g.2C>G": ("synonymous_variant", None, "most_severe"),
    }

    out = mod.annotate_row(
        row,
        cache,
        col_prefix="vep",
        hgvs_cols=["mapped_hgvs_g"],
        access_date="2026-05-29",
    )

    assert out["vep.access_date"] == "2026-05-29|2026-05-29"


@pytest.mark.unit
def test_access_date_preserves_empty_candidate_slots():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T||NC_000001.11:g.3G>A"}
    cache = {
        "NC_000001.11:g.1A>T": ("missense_variant", None, "most_severe"),
        "NC_000001.11:g.3G>A": ("intron_variant", None, "most_severe"),
    }

    out = mod.annotate_row(
        row,
        cache,
        col_prefix="vep",
        hgvs_cols=["mapped_hgvs_g"],
        access_date="2026-05-29",
    )

    assert out["vep.access_date"] == "2026-05-29||2026-05-29"
