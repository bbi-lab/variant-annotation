import csv

import pytest

from src import reverse_translate_protein_variants as rtpv
from variant_annotation.lib.translation.types import (
    TranslationError,
    VariantInput,
    TranslationConfig,
    WtCodonMode,
)



class _FakeConnection:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def _write_tsv(path, rows):
    fieldnames = [
        "variant_urn",
        "mapped_hgvs_p",
        "mapped_hgvs_c",
        "mapped_hgvs_g",
        "mapping_error",
        "assayed_variant_level",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _patch_clients(monkeypatch):
    monkeypatch.setenv("UTA_DB_URL", "postgresql://user:pass@localhost:5432/uta")
    monkeypatch.setattr("src.reverse_translate_protein_variants.connect_uta", lambda url: _FakeConnection())
    monkeypatch.setattr("src.reverse_translate_protein_variants.UtaClient", lambda conn: object())
    monkeypatch.setattr(
        "src.reverse_translate_protein_variants.HgvsMapper",
        type("M", (), {"from_url": staticmethod(lambda url, assembly: object())})(),
    )


# ---------------------------------------------------------------------------
# CLI integration tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_reverse_translate_preserves_row_order_and_updates_targets(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {
            "variant_urn": "v0",
            "mapped_hgvs_p": "NP_000001.1:p.Ala1Val",
            "mapped_hgvs_c": "",
            "mapped_hgvs_g": "",
            "mapping_error": "",
            "assayed_variant_level": "",
        },
        {
            "variant_urn": "v1",
            "mapped_hgvs_p": "NP_000001.1:p.Ala2Val",
            "mapped_hgvs_c": "NM_000001.1:c.2A>G",
            "mapped_hgvs_g": "NC_000001.11:g.2A>G",
            "mapping_error": "",
            "assayed_variant_level": "dna",
        },
        {
            "variant_urn": "v2",
            "mapped_hgvs_p": "NP_000001.1:p.Ala3Val",
            "mapped_hgvs_c": "",
            "mapped_hgvs_g": "",
            "mapping_error": "",
            "assayed_variant_level": "",
        },
    ]
    _write_tsv(input_path, rows)
    _patch_clients(monkeypatch)

    def fake_process_rows(rows, transcripts, coordinates, config=None, columns=None, **kwargs):
        c_col = columns.mapped_hgvs_c if columns else "mapped_hgvs_c"
        p_col = columns.mapped_hgvs_p if columns else "mapped_hgvs_p"
        g_col = columns.mapped_hgvs_g if columns else "mapped_hgvs_g"
        for row in rows:
            if (row.get(p_col) or "").strip() and not (row.get(c_col) or "").strip():
                row[c_col] = "NM_000001.1:c.123A>G"
                row[g_col] = "NC_000001.11:g.123A>G"
        return rows

    monkeypatch.setattr("src.reverse_translate_protein_variants.process_rows", fake_process_rows)

    rtpv.reverse_translate_protein_variants(str(input_path), str(output_path))

    out_rows = _read_tsv(output_path)
    assert [r["variant_urn"] for r in out_rows] == ["v0", "v1", "v2"]
    assert out_rows[0]["mapped_hgvs_c"] == "NM_000001.1:c.123A>G"
    assert out_rows[1]["mapped_hgvs_c"] == "NM_000001.1:c.2A>G"
    assert out_rows[2]["mapped_hgvs_c"] == "NM_000001.1:c.123A>G"


@pytest.mark.integration
def test_reverse_translate_adds_error_when_translation_fails(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {
            "variant_urn": "v0",
            "mapped_hgvs_p": "NP_000001.1:p.Ala1Val",
            "mapped_hgvs_c": "",
            "mapped_hgvs_g": "",
            "mapping_error": "existing",
            "assayed_variant_level": "",
        }
    ]
    _write_tsv(input_path, rows)
    _patch_clients(monkeypatch)

    def fake_process_rows(rows, transcripts, coordinates, config=None, columns=None, **kwargs):
        for row in rows:
            columns.write_error(
                row, TranslationError(input=VariantInput(hgvs=""), error="Unable to resolve transcript accession")
            )
        return rows

    monkeypatch.setattr("src.reverse_translate_protein_variants.process_rows", fake_process_rows)

    rtpv.reverse_translate_protein_variants(str(input_path), str(output_path))

    out_rows = _read_tsv(output_path)
    assert "existing" in out_rows[0]["mapping_error"]
    assert "Unable to resolve transcript accession" in out_rows[0]["reverse_translation_error"]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_wt_codon_mode_requires_include_indels(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(input_path, [])
    monkeypatch.setenv("UTA_DB_URL", "postgresql://user:pass@localhost:5432/uta")

    with pytest.raises(ValueError, match="include_indels=True"):
        rtpv.reverse_translate_protein_variants(
            str(input_path),
            str(output_path),
            TranslationConfig(wt_codon_mode=WtCodonMode.UNAMBIGUOUS, include_indels=False),
        )


# ---------------------------------------------------------------------------
# Parser tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_build_parser_defaults_and_flags():
    args = rtpv._build_parser().parse_args(["input.tsv", "output.tsv"])
    assert args.input_file == "input.tsv"
    assert args.assembly == "GRCh38"
    assert args.strict_ref_aa is True
    assert args.resolve_missing_ref_alleles is True
    assert args.allow_length_changing_stop_candidates is True
    assert args.wt_codon_mode == "none"


@pytest.mark.unit
def test_build_parser_wt_codon_mode_choices():
    for choice in ("none", "unambiguous", "all"):
        args = rtpv._build_parser().parse_args(["input.tsv", "output.tsv", "--wt-codon-mode", choice])
        assert args.wt_codon_mode == choice
