# Copyright (C) 2026 Jeremy Stone
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Tests for src/fetch_mavedb_scores.py."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from src.fetch_mavedb_scores import (
    MAVEDB_MAPPING_COLS,
    OUTPUT_COLS,
    build_score_set_rows,
    fetch_variant_details,
    load_score_sets,
    main,
    mapped_hgvs_from_variant_details,
    output_columns,
    score_set_metadata,
    scores_from_csv_text,
)

URN = "urn:mavedb:00000895-a-1"

SCORE_SET = {
    "urn": URN,
    "title": "MECP2 Methyl-CpG DNA binding domain domainome 1.0",
    "targetGenes": [
        {
            "name": "MECP2 Methyl-CpG DNA binding domain",
            "targetSequence": {"sequenceType": "protein", "sequence": "RGPMYDD"},
        }
    ],
}

SCORES_CSV = (
    "accession,hgvs_nt,hgvs_splice,hgvs_pro,score,sigma,raw_score\n"
    f"{URN}#2,NA,NA,p.Arg1Ala,-0.164216716233428,0.0465965450755631,4\n"
    f"{URN}#10,NA,NA,p.Gly2Ala,NA,NA,NA\n"
    f"{URN}#1,NA,NA,p.Arg1Ter,-0.930191246530987,0.120676697665939,-0.007047666\n"
)


def _detail(variant_urn: str, protein_hgvs: str, **extra) -> str:
    record = {
        "urn": variant_urn,
        "isCurrent": True,
        "mapped": {"genomic": None, "cdna": None, "protein": {"hgvs": protein_hgvs}},
    }
    record.update(extra)
    return json.dumps(record)


# ---------------------------------------------------------------------------
# load_score_sets
# ---------------------------------------------------------------------------


def test_load_score_sets_with_preferred_transcript(tmp_path: Path):
    path = tmp_path / "score_sets.tsv"
    path.write_text(
        "dataset_name\tscore_set_urn\tpreferred_transcript\n"
        f"MECP2\t{URN}\tNM_001110792.2\n"
        "CRX\turn:mavedb:00000736-a-1\t\n"
        "blank\t\tNM_1\n"
        f"dup\t{URN}\tNM_2\n"
    )
    assert load_score_sets(path) == [(URN, "NM_001110792.2"), ("urn:mavedb:00000736-a-1", "")]


def test_load_score_sets_without_transcript_column_csv(tmp_path: Path):
    path = tmp_path / "score_sets.csv"
    path.write_text(f"score_set_urn,requested_calibration_urn\n{URN},\n")
    assert load_score_sets(path) == [(URN, "")]


def test_load_score_sets_missing_urn_column(tmp_path: Path):
    path = tmp_path / "score_sets.tsv"
    path.write_text("dataset_name\nX\n")
    with pytest.raises(ValueError, match="score_set_urn"):
        load_score_sets(path)


# ---------------------------------------------------------------------------
# score_set_metadata
# ---------------------------------------------------------------------------


def test_score_set_metadata_uses_mapped_gene_symbol():
    meta = score_set_metadata(SCORE_SET, {URN: ["MECP2"]})
    assert meta == {
        "gene_symbol": "MECP2",
        "score_set_title": "MECP2 Methyl-CpG DNA binding domain domainome 1.0",
        "reference_based": "false",
        "target_sequence_type": "protein",
        "target_sequence": "RGPMYDD",
    }


def test_score_set_metadata_falls_back_to_mapped_hgnc_name_not_target_name():
    score_set = json.loads(json.dumps(SCORE_SET))
    score_set["targetGenes"][0]["mappedHgncName"] = "MECP2"
    assert score_set_metadata(score_set, {})["gene_symbol"] == "MECP2"
    # Without either, the long target name is not used as a gene symbol.
    assert score_set_metadata(SCORE_SET, {})["gene_symbol"] == ""


def test_score_set_metadata_reference_based():
    score_set = json.loads(json.dumps(SCORE_SET))
    score_set["targetGenes"][0]["targetAccession"] = {"accession": "NM_004992.4"}
    assert score_set_metadata(score_set, {})["reference_based"] == "true"


# ---------------------------------------------------------------------------
# scores_from_csv_text
# ---------------------------------------------------------------------------


def test_scores_from_csv_text_keeps_score_text_and_packs_score_data():
    rows = scores_from_csv_text(SCORES_CSV)
    assert rows[0]["variant_urn"] == f"{URN}#2"
    assert rows[0]["raw_hgvs_nt"] == ""
    assert rows[0]["raw_hgvs_pro"] == "p.Arg1Ala"
    assert rows[0]["score"] == "-0.164216716233428"
    assert json.loads(rows[0]["score_data"]) == {
        "score": -0.164216716233428,
        "sigma": 0.0465965450755631,
        "raw_score": 4,
    }
    assert rows[1]["score"] == ""
    assert json.loads(rows[1]["score_data"]) == {}


def test_scores_from_csv_text_rna_score_aliases():
    text = "accession,hgvs_nt,hgvs_pro,score,score_rna,rna_score_d6\nu#1,c.1A>G,NA,1.5,0.2,0.3\n"
    row = scores_from_csv_text(text)[0]
    assert (row["rna_score"], row["rna_score_d6"], row["rna_score_d20"]) == ("0.2", "0.3", "")


def test_scores_from_csv_text_requires_score_column():
    with pytest.raises(ValueError, match="score"):
        scores_from_csv_text("accession,hgvs_pro\nu#1,p.=\n")


# ---------------------------------------------------------------------------
# variant details
# ---------------------------------------------------------------------------


def test_mapped_hgvs_from_variant_details():
    lines = [
        _detail(f"{URN}#1", "NP_004983.1:p.Arg91Ter"),
        "",
        _detail(f"{URN}#2", "NP_004983.1:p.Arg91Ala", isCurrent=False),
        json.dumps(
            {
                "urn": f"{URN}#3",
                "mapped": {
                    "genomic": {"hgvs": "NC_000023.11:g.1A>G"},
                    "cdna": {"hgvs": "NM_004992.4:c.1A>G"},
                    "protein": None,
                },
            }
        ),
    ]
    assert mapped_hgvs_from_variant_details(lines) == {
        f"{URN}#1": {
            "mavedb_mapped_hgvs_g": "",
            "mavedb_mapped_hgvs_c": "",
            "mavedb_mapped_hgvs_p": "NP_004983.1:p.Arg91Ter",
        },
        f"{URN}#3": {
            "mavedb_mapped_hgvs_g": "NC_000023.11:g.1A>G",
            "mavedb_mapped_hgvs_c": "NM_004992.4:c.1A>G",
            "mavedb_mapped_hgvs_p": "",
        },
    }


def test_fetch_variant_details_404_returns_empty():
    session = MagicMock()
    session.get.return_value = MagicMock(status_code=404)
    assert fetch_variant_details("https://api.example", URN, session) == []
    session.get.return_value.raise_for_status.assert_not_called()


# ---------------------------------------------------------------------------
# build_score_set_rows
# ---------------------------------------------------------------------------


def test_output_columns():
    assert "variant_id" not in OUTPUT_COLS
    assert "mavedb_mapping_error" not in OUTPUT_COLS
    assert output_columns() == OUTPUT_COLS
    with_mappings = output_columns(include_mavedb_mappings=True)
    i = with_mappings.index("reference_based")
    assert with_mappings[i + 1 : i + 4] == MAVEDB_MAPPING_COLS
    assert [c for c in with_mappings if c not in MAVEDB_MAPPING_COLS] == OUTPUT_COLS


def test_build_score_set_rows_sorts_by_variant_number_and_fills_columns():
    meta = score_set_metadata(SCORE_SET, {URN: ["MECP2"]})
    rows = build_score_set_rows(meta, scores_from_csv_text(SCORES_CSV), preferred_transcript="NM_001110792.2")

    assert [r["variant_urn"] for r in rows] == [f"{URN}#1", f"{URN}#2", f"{URN}#10"]
    assert all(list(r) == OUTPUT_COLS for r in rows)
    assert {r["preferred_transcript"] for r in rows} == {"NM_001110792.2"}
    assert {r["gene_symbol"] for r in rows} == {"MECP2"}


def test_build_score_set_rows_with_mavedb_mappings():
    meta = score_set_metadata(SCORE_SET, {URN: ["MECP2"]})
    mapped = mapped_hgvs_from_variant_details([_detail(f"{URN}#1", "NP_004983.1:p.Arg91Ter")])
    rows = build_score_set_rows(meta, scores_from_csv_text(SCORES_CSV), mapped)

    assert all(list(r) == output_columns(include_mavedb_mappings=True) for r in rows)
    assert rows[0]["mavedb_mapped_hgvs_p"] == "NP_004983.1:p.Arg91Ter"
    assert rows[1]["mavedb_mapped_hgvs_p"] == ""


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def _response(status_code=200, *, json_body=None, text=""):
    resp = MagicMock(status_code=status_code, text=text)
    resp.json.return_value = json_body
    return resp


def _run_main(tmp_path: Path, *extra_args: str):
    score_sets_file = tmp_path / "score_sets.tsv"
    score_sets_file.write_text(f"score_set_urn\tpreferred_transcript\n{URN}\tNM_001110792.2\n")
    output = tmp_path / "out" / "variants.0.tsv"

    def fake_get(url, timeout):
        if url.endswith("/score-sets/mapped-genes"):
            return _response(json_body={URN: ["MECP2"]})
        if url.endswith("/scores"):
            return _response(text=SCORES_CSV)
        if url.endswith("/variant-details"):
            return _response(text=_detail(f"{URN}#1", "NP_004983.1:p.Arg91Ter") + "\n")
        return _response(json_body=SCORE_SET)

    session = MagicMock()
    session.get.side_effect = fake_get
    with patch("src.fetch_mavedb_scores.build_session", return_value=session):
        main(
            [str(output), "--score-sets-file", str(score_sets_file), "--mavedb-api-url", "https://api.example"]
            + list(extra_args)
        )

    with open(output, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)
    requested_urls = [c.args[0] for c in session.get.call_args_list]
    return fieldnames, rows, requested_urls


def test_main_writes_tsv_without_mavedb_mappings(tmp_path: Path):
    fieldnames, rows, requested_urls = _run_main(tmp_path)

    assert fieldnames == OUTPUT_COLS
    assert [r["variant_urn"] for r in rows] == [f"{URN}#1", f"{URN}#2", f"{URN}#10"]
    assert rows[0]["preferred_transcript"] == "NM_001110792.2"
    assert not any(url.endswith("/variant-details") for url in requested_urls)


def test_main_writes_tsv_with_mavedb_mappings(tmp_path: Path):
    fieldnames, rows, requested_urls = _run_main(tmp_path, "--include-mavedb-mappings")

    assert fieldnames == output_columns(include_mavedb_mappings=True)
    assert rows[0]["mavedb_mapped_hgvs_p"] == "NP_004983.1:p.Arg91Ter"
    assert any(url.endswith("/variant-details") for url in requested_urls)


def test_main_requires_score_sets(tmp_path: Path):
    with pytest.raises(SystemExit):
        main([str(tmp_path / "out.tsv")])
