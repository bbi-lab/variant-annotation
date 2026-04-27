"""Unit tests for src/annotate_spliceai.py."""

from __future__ import annotations

import csv

import src.annotate_spliceai as mod

from src.annotate_spliceai import (
    annotate_row_with_scores,
    parse_spliceai_info,
    split_pipe_preserve_positions,
)


def test_parse_spliceai_info_extracts_all_fields_and_max_delta():
    info = (
        "NS=1;SpliceAI="
        "A|GENE1|0.10|0.20|0.30|0.40|5|6|7|8,"
        "A|GENE2|0.15|0.25|0.35|0.05|9|10|11|12"
    )
    parsed = parse_spliceai_info(info, alt="A")

    assert parsed["spliceai.ds_ag"] == "0.1500"
    assert parsed["spliceai.ds_al"] == "0.2500"
    assert parsed["spliceai.ds_dg"] == "0.3500"
    assert parsed["spliceai.ds_dl"] == "0.4000"
    assert parsed["spliceai.dp_ag"] == "9.0000"
    assert parsed["spliceai.dp_al"] == "10.0000"
    assert parsed["spliceai.dp_dg"] == "11.0000"
    assert parsed["spliceai.dp_dl"] == "12.0000"
    assert parsed["spliceai.max_delta_score"] == "0.4000"


def test_split_pipe_preserve_positions_keeps_empty_slots():
    assert split_pipe_preserve_positions("a||c") == ["a", "", "c"]


def test_annotate_row_with_scores_emits_pipe_aligned_columns():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T||NC_000001.11:g.2C>G"}
    score_map = {
        "NC_000001.11:g.1A>T": {
            "spliceai.ds_ag": "0.1000",
            "spliceai.ds_al": "0.0000",
            "spliceai.ds_dg": "0.0000",
            "spliceai.ds_dl": "0.2000",
            "spliceai.dp_ag": "3.0000",
            "spliceai.dp_al": "0.0000",
            "spliceai.dp_dg": "0.0000",
            "spliceai.dp_dl": "4.0000",
            "spliceai.max_delta_score": "0.2000",
        },
        "NC_000001.11:g.2C>G": {
            "spliceai.ds_ag": "0.4000",
            "spliceai.ds_al": "0.1000",
            "spliceai.ds_dg": "0.2000",
            "spliceai.ds_dl": "0.3000",
            "spliceai.dp_ag": "8.0000",
            "spliceai.dp_al": "9.0000",
            "spliceai.dp_dg": "10.0000",
            "spliceai.dp_dl": "11.0000",
            "spliceai.max_delta_score": "0.4000",
        },
    }

    out = annotate_row_with_scores(row, "mapped_hgvs_g", score_map)


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
            writer.writerow({"variant_urn": "v4", "mapped_hgvs_g": "NC_000001.11:g.4T>C"})

        monkeypatch.setattr(mod, "_cache_and_index_precomputed_files", lambda args: [tmp_path / "scores.vcf.gz"])
        monkeypatch.setattr(mod, "get_nc_to_chrom_map", lambda annotation: {"NC_000001.11": "1"})
        monkeypatch.setattr(mod, "parse_hgvs_g_to_vcf", lambda hgvs_g, nc_to_chrom, fasta: None)
        monkeypatch.setattr(
            mod,
            "lookup_precomputed_scores",
            lambda hgvs_to_exact_vcf, nc_to_chrom, vcf_paths: {
                hgvs_g: {
                    "spliceai.ds_ag": hgvs_g[-1],
                    "spliceai.ds_al": "",
                    "spliceai.ds_dg": "",
                    "spliceai.ds_dl": "",
                    "spliceai.dp_ag": "",
                    "spliceai.dp_al": "",
                    "spliceai.dp_dg": "",
                    "spliceai.dp_dl": "",
                    "spliceai.max_delta_score": hgvs_g[-1],
                }
                for hgvs_g in hgvs_to_exact_vcf
            },
        )

        mod.main([
            str(in_path),
            str(out_path),
            "--precomputed-vcf",
            str(tmp_path / "source.vcf.gz"),
            "--skip",
            "1",
            "--limit",
            "2",
        ])

        with out_path.open("r", encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))

        assert [r["variant_urn"] for r in rows] == ["v2", "v3"]
        assert [r["spliceai.max_delta_score"] for r in rows] == ["G", "A"]

    assert out["spliceai.ds_ag"] == "0.1000||0.4000"
    assert out["spliceai.ds_dl"] == "0.2000||0.3000"
    assert out["spliceai.dp_ag"] == "3.0000||8.0000"
    assert out["spliceai.max_delta_score"] == "0.2000||0.4000"
