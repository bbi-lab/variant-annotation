"""Unit tests for src/annotate_vep.py."""

from __future__ import annotations

import csv
from datetime import date

import src.annotate_vep as mod


def test_annotate_row_emits_single_consequence_for_matching_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    consequence_cache = {
        "NC_000001.11:g.1A>T": "missense_variant",
        "NC_000001.11:g.2C>G": "missense_variant",
    }
    out = mod.annotate_row(
        row,
        consequence_cache,
        col_prefix="vep",
        mapped_hgvs_g_col="mapped_hgvs_g",
        access_date="2026-04-30",
    )

    assert out["vep.mutational_consequence"] == "missense_variant"
    assert out["vep.access_date"] == "2026-04-30"
    assert out["vep.error"] == ""


def test_annotate_row_assumes_missing_equals_other_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    consequence_cache = {
        "NC_000001.11:g.1A>T": "synonymous_variant",
        # second candidate unresolved by VEP
        "NC_000001.11:g.2C>G": None,
    }
    out = mod.annotate_row(
        row,
        consequence_cache,
        col_prefix="vep",
        mapped_hgvs_g_col="mapped_hgvs_g",
        access_date="2026-04-30",
    )

    assert out["vep.mutational_consequence"] == "synonymous_variant"
    assert out["vep.error"] == ""


def test_annotate_row_reports_discrepancy_with_pipe_aligned_values():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G|NC_000001.11:g.3G>A"}
    consequence_cache = {
        "NC_000001.11:g.1A>T": "missense_variant",
        "NC_000001.11:g.2C>G": "synonymous_variant",
        "NC_000001.11:g.3G>A": None,
    }
    out = mod.annotate_row(
        row,
        consequence_cache,
        col_prefix="vep",
        mapped_hgvs_g_col="mapped_hgvs_g",
        access_date="2026-04-30",
    )

    assert out["vep.mutational_consequence"] == ""
    assert "missense_variant|synonymous_variant|" in out["vep.error"]


def test_get_functional_consequence_uses_recoder_fallback(monkeypatch):
    calls = []

    def fake_vep_lookup(hgvs_strings, *, api_url, timeout_seconds):
        calls.append(tuple(hgvs_strings))
        if hgvs_strings == ["NC_000001.11:g.1A>T", "NM_000000.1:c.1A>T"]:
            return {"NC_000001.11:g.1A>T": "intron_variant"}
        if hgvs_strings == ["NC_000001.11:g.5A>T", "NC_000001.11:g.6A>T"]:
            return {
                "NC_000001.11:g.5A>T": "synonymous_variant",
                "NC_000001.11:g.6A>T": "missense_variant",
            }
        return {}

    def fake_recoder(hgvs_strings, *, api_url, timeout_seconds):
        assert hgvs_strings == ["NM_000000.1:c.1A>T"]
        return {"NM_000000.1:c.1A>T": ["NC_000001.11:g.5A>T", "NC_000001.11:g.6A>T"]}

    monkeypatch.setattr(mod, "_vep_lookup_batch", fake_vep_lookup)
    monkeypatch.setattr(mod, "run_variant_recoder", fake_recoder)

    out = mod.get_functional_consequence(
        ["NC_000001.11:g.1A>T", "NM_000000.1:c.1A>T"],
        api_url="https://rest.ensembl.org",
        timeout_seconds=10,
        batch_size=200,
    )

    assert calls == [
        ("NC_000001.11:g.1A>T", "NM_000000.1:c.1A>T"),
        ("NC_000001.11:g.5A>T", "NC_000001.11:g.6A>T"),
    ]
    assert out["NC_000001.11:g.1A>T"] == "intron_variant"
    # most severe of synonymous vs missense should be missense
    assert out["NM_000000.1:c.1A>T"] == "missense_variant"


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

    def fake_get_functional_consequence(hgvs_strings, *, api_url, timeout_seconds, batch_size):
        return {
            hgvs: "missense_variant" if hgvs.endswith("2C>G") else "synonymous_variant"
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
    assert rows[0]["vep.mutational_consequence"] == "missense_variant"
    assert rows[0]["vep.error"] == ""
    assert rows[0]["vep.access_date"] == date.today().isoformat()
