"""Tests for src/filter_dna_variants.py."""

from __future__ import annotations
import pytest

import csv
from pathlib import Path

import filter_dna_variants as mod



pytestmark = pytest.mark.integration


def _write_input_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _read_output_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_filter_dna_variants_keeps_snv_and_wt_and_compacts_aligned_fields(tmp_path):
    in_path = tmp_path / "input.tsv"
    out_path = tmp_path / "output.tsv"

    _write_input_tsv(
        in_path,
        [
            {
                "variant_urn": "v1",
                "mapped_hgvs_c": "NM_000001.1:c.76A>T|NM_000001.1:c.90_91del|NM_000001.1:c.100=",
                "mapped_hgvs_g": "NC_000001.11:g.76A>T|NC_000001.11:g.90_91del|NC_000001.11:g.100=",
                "mapped_hgvs_p": "NP_000001.1:p.Gly1Val",
                "dna_clingen_allele_id": "CA1|CA2|CA3",
                "touches_intronic_region": "false|true|false",
                "spans_intron": "false|true|false",
                "vep.access_date": "2026-05-29|2026-05-29|2026-05-29",
            }
        ],
    )

    stats = mod.filter_dna_variants(in_path, out_path)
    rows = _read_output_tsv(out_path)

    assert stats.rows_processed == 1
    assert stats.candidates_removed == 1
    assert stats.rows_without_dna_and_protein == 0
    assert len(rows) == 1
    assert rows[0]["mapped_hgvs_c"] == "NM_000001.1:c.76A>T|NM_000001.1:c.100="
    assert rows[0]["mapped_hgvs_g"] == "NC_000001.11:g.76A>T|NC_000001.11:g.100="
    assert rows[0]["dna_clingen_allele_id"] == "CA1|CA3"
    assert rows[0]["touches_intronic_region"] == "false|false"
    assert rows[0]["spans_intron"] == "false|false"
    assert rows[0]["vep.access_date"] == "2026-05-29|2026-05-29"
    assert rows[0]["mapped_hgvs_p"] == "NP_000001.1:p.Gly1Val"


def test_filter_dna_variants_falls_back_to_g_when_c_is_blank_and_warns_on_empty_rows(tmp_path):
    in_path = tmp_path / "input.tsv"
    out_path = tmp_path / "output.tsv"

    _write_input_tsv(
        in_path,
        [
            {
                "variant_urn": "v1",
                "mapped_hgvs_c": "",
                "mapped_hgvs_g": "NC_000001.11:g.76A>T|NC_000001.11:g.90_91del",
                "mapped_hgvs_p": "NP_000001.1:p.Gly1Val",
                "dna_clingen_allele_id": "CA1|CA2",
                "vep.access_date": "2026-05-29|2026-05-29",
            },
            {
                "variant_urn": "v2",
                "mapped_hgvs_c": "NM_000001.1:c.10_11del|NM_000001.1:c.20_21dup",
                "mapped_hgvs_g": "NC_000001.11:g.10_11del|NC_000001.11:g.20_21dup",
                "mapped_hgvs_p": "",
                "dna_clingen_allele_id": "CA3|CA4",
                "vep.access_date": "2026-05-29|2026-05-29",
            },
        ],
    )

    stats = mod.filter_dna_variants(in_path, out_path)
    rows = _read_output_tsv(out_path)

    assert stats.rows_processed == 2
    assert stats.candidates_removed == 3
    assert stats.rows_without_dna_and_protein == 1
    assert rows[0]["mapped_hgvs_c"] == ""
    assert rows[0]["mapped_hgvs_g"] == "NC_000001.11:g.76A>T"
    assert rows[0]["dna_clingen_allele_id"] == "CA1"
    assert rows[0]["vep.access_date"] == "2026-05-29"
    assert rows[0]["mapped_hgvs_p"] == "NP_000001.1:p.Gly1Val"
    assert rows[1]["mapped_hgvs_c"] == ""
    assert rows[1]["mapped_hgvs_g"] == ""
    assert rows[1]["dna_clingen_allele_id"] == ""
    assert rows[1]["mapped_hgvs_p"] == ""


def test_filter_dna_variants_can_keep_only_snvs(tmp_path):
    in_path = tmp_path / "input.tsv"
    out_path = tmp_path / "output.tsv"

    _write_input_tsv(
        in_path,
        [
            {
                "variant_urn": "v1",
                "mapped_hgvs_c": "NM_000001.1:c.76A>T|NM_000001.1:c.100=",
                "mapped_hgvs_g": "NC_000001.11:g.76A>T|NC_000001.11:g.100=",
                "mapped_hgvs_p": "NP_000001.1:p.Gly1Val",
                "dna_clingen_allele_id": "CA1|CA2",
                "vep.access_date": "2026-05-29|2026-05-29",
            }
        ],

    )

    stats = mod.filter_dna_variants(in_path, out_path, keep_classes=["snv"])
    rows = _read_output_tsv(out_path)

    assert stats.rows_processed == 1
    assert stats.candidates_removed == 1
    assert rows[0]["mapped_hgvs_c"] == "NM_000001.1:c.76A>T"
    assert rows[0]["mapped_hgvs_g"] == "NC_000001.11:g.76A>T"
    assert rows[0]["dna_clingen_allele_id"] == "CA1"
    assert rows[0]["vep.access_date"] == "2026-05-29"
