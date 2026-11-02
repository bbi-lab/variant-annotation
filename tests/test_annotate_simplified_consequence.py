"""Unit tests for src/annotate_simplified_consequence.py."""

from __future__ import annotations

import csv
import gzip

import src.annotate_simplified_consequence as mod
from src.annotate_simplified_consequence import (
    CONFLICT_CONSEQUENCE,
    build_mapping_dict,
    compute_condensed_consequence,
    compute_splice_var_amino,
    compute_splice_variant,
    get_simplified_consequence,
    simplify_mutational_consequences,
)

_MAP_CSV = """VEP output term,SO summary term,Importance
splice_acceptor_variant,splice_site_variant,400
missense_variant,missense_variant,290
synonymous_variant,synonymous_variant,200
protein_altering_variant,,280
"""


def _write_map_file(tmp_path):
    map_path = tmp_path / "consequence_map.csv.gz"
    with gzip.open(map_path, "wt", encoding="utf-8") as fh:
        fh.write(_MAP_CSV)
    return map_path


# ---------------------------------------------------------------------------
# build_mapping_dict
# ---------------------------------------------------------------------------


def test_build_mapping_dict_maps_known_terms(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    assert mapping_dict["missense_variant"] == "missense_variant"
    assert mapping_dict["splice_acceptor_variant"] == "splice_site_variant"


def test_build_mapping_dict_adds_extra_terms(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    assert mapping_dict["UTR_variant"] == "UTR_variant"
    assert mapping_dict["splicing_variant"] == "splicing_variant"
    assert mapping_dict["splice_site_variant"] == "splice_site_variant"


def test_build_mapping_dict_nan_summary_term(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    value = mapping_dict["protein_altering_variant"]
    # pandas represents a missing "SO summary term" as None (or NaN, depending on
    # version) after groupby().first(); simplify_mutational_consequences renders
    # either as a blank segment.
    assert value is None or value != value


# ---------------------------------------------------------------------------
# get_simplified_consequence
# ---------------------------------------------------------------------------


def test_get_simplified_consequence_returns_first_match(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    assert get_simplified_consequence(
        ["missense_variant", "splice_acceptor_variant"], mapping_dict
    ) == "missense_variant"


def test_get_simplified_consequence_no_match_returns_none(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    assert get_simplified_consequence(["some_unknown_term"], mapping_dict) is None


def test_get_simplified_consequence_strips_whitespace(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    assert get_simplified_consequence([" missense_variant "], mapping_dict) == "missense_variant"


# ---------------------------------------------------------------------------
# simplify_mutational_consequences
# ---------------------------------------------------------------------------


def test_simplify_mutational_consequences_single_set_per_transcript(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    assert simplify_mutational_consequences(
        "missense_variant|missense_variant", mapping_dict
    ) == "missense_variant|missense_variant"


def test_simplify_mutational_consequences_multi_term_sets(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    value = "missense_variant|missense_variant^splice_acceptor_variant"
    assert simplify_mutational_consequences(value, mapping_dict) == "missense_variant|missense_variant"


def test_simplify_mutational_consequences_no_match_yields_blank_segment(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    value = "missense_variant|some_unknown_term"
    assert simplify_mutational_consequences(value, mapping_dict) == "missense_variant|"


def test_simplify_mutational_consequences_nan_summary_term_yields_blank(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    assert simplify_mutational_consequences("protein_altering_variant", mapping_dict) == ""


def test_simplify_mutational_consequences_blank_input(tmp_path):
    mapping_dict = build_mapping_dict(str(_write_map_file(tmp_path)))
    assert simplify_mutational_consequences("", mapping_dict) == ""


# ---------------------------------------------------------------------------
# compute_condensed_consequence
# ---------------------------------------------------------------------------


def test_compute_condensed_consequence_blank_input():
    assert compute_condensed_consequence("") == ""


def test_compute_condensed_consequence_single_value():
    assert compute_condensed_consequence("missense_variant") == "missense_variant"


def test_compute_condensed_consequence_all_identical():
    assert compute_condensed_consequence("missense_variant|missense_variant|missense_variant") == "missense_variant"


def test_compute_condensed_consequence_conflicting_values():
    assert compute_condensed_consequence("missense_variant|synonymous_variant") == CONFLICT_CONSEQUENCE


def test_compute_condensed_consequence_blank_segment_conflicts_with_value():
    assert compute_condensed_consequence("missense_variant|") == CONFLICT_CONSEQUENCE


# ---------------------------------------------------------------------------
# compute_splice_variant
# ---------------------------------------------------------------------------


def test_compute_splice_variant_score_above_threshold():
    assert compute_splice_variant("0.5", "0.0", "0.0", "0.0", "missense_variant") == "Yes"


def test_compute_splice_variant_score_below_threshold_and_not_splice_consequence():
    assert compute_splice_variant("0.1", "0.1", "0.1", "0.1", "missense_variant") == "No"


def test_compute_splice_variant_splice_site_consequence_below_threshold():
    assert compute_splice_variant("0.0", "0.0", "0.0", "0.0", "splice_site_variant") == "Yes"


def test_compute_splice_variant_threshold_is_inclusive():
    assert compute_splice_variant("0.2", "", "", "", "missense_variant") == "Yes"


def test_compute_splice_variant_custom_threshold():
    assert compute_splice_variant("0.3", "", "", "", "missense_variant", threshold=0.5) == "No"


def test_compute_splice_variant_pipe_aligned_positions():
    ds_ag = "0.5|0.0"
    ds_al = "0.0|0.0"
    ds_dg = "0.0|0.0"
    ds_dl = "0.0|0.0"
    cons = "missense_variant|splice_site_variant"
    assert compute_splice_variant(ds_ag, ds_al, ds_dg, ds_dl, cons) == "Yes|Yes"


def test_compute_splice_variant_no_data_returns_blank():
    assert compute_splice_variant("", "", "", "", "") == ""


def test_compute_splice_variant_no_evidence_position_is_blank():
    assert compute_splice_variant("||", "||", "||", "||", "||") == "||"
    assert compute_splice_variant("0.0||", "0.0||", "0.0||", "0.0||", "|missense_variant|") == "No|No|"


def test_compute_splice_variant_scores_without_consequence_are_evidence():
    assert compute_splice_variant("0.0", "0.0", "0.0", "0.0", "") == "No"
    assert compute_splice_variant("0.5", "0.0", "0.0", "0.0", "") == "Yes"


def test_compute_splice_variant_missing_scores_falls_back_to_consequence():
    # spliceai columns entirely blank (not computed), consequence still present.
    assert compute_splice_variant("", "", "", "", "splice_site_variant") == "Yes"
    assert compute_splice_variant("", "", "", "", "missense_variant") == "No"


# ---------------------------------------------------------------------------
# compute_splice_var_amino
# ---------------------------------------------------------------------------


def test_compute_splice_var_amino_blank_input():
    assert compute_splice_var_amino("") == ""


def test_compute_splice_var_amino_all_blank_segments():
    assert compute_splice_var_amino("|") == ""


def test_compute_splice_var_amino_single_unique_value():
    assert compute_splice_var_amino("No|No|No") == "No"
    assert compute_splice_var_amino("Yes") == "Yes"


def test_compute_splice_var_amino_multiple_unique_values():
    assert compute_splice_var_amino("Yes|No") == "Yes"


def test_compute_splice_var_amino_ignores_blank_segments():
    assert compute_splice_var_amino("No||No") == "No"


# ---------------------------------------------------------------------------
# main (end-to-end)
# ---------------------------------------------------------------------------


def test_main_adds_simplified_consequence_column(tmp_path):
    map_path = _write_map_file(tmp_path)
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "vep.mutational_consequences"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({
            "variant_urn": "v1",
            "vep.mutational_consequences": "missense_variant|missense_variant^splice_acceptor_variant",
        })
        writer.writerow({"variant_urn": "v2", "vep.mutational_consequences": ""})
        writer.writerow({
            "variant_urn": "v3",
            "vep.mutational_consequences": "missense_variant|synonymous_variant",
        })

    mod.main([str(in_path), str(out_path), "--consequence-map-file", str(map_path)])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert rows[0]["variant_urn"] == "v1"
    assert rows[0]["simplified_consequence"] == "missense_variant|missense_variant"
    assert rows[0]["condensed_consequence"] == "missense_variant"

    assert rows[1]["variant_urn"] == "v2"
    assert rows[1]["simplified_consequence"] == ""
    assert rows[1]["condensed_consequence"] == ""

    assert rows[2]["variant_urn"] == "v3"
    assert rows[2]["simplified_consequence"] == "missense_variant|synonymous_variant"
    assert rows[2]["condensed_consequence"] == CONFLICT_CONSEQUENCE

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    assert header.index("condensed_consequence") == header.index("simplified_consequence") + 1


def test_main_adds_splice_variant_and_splice_var_amino_columns(tmp_path):
    map_path = _write_map_file(tmp_path)
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    fieldnames = [
        "variant_urn",
        "vep.mutational_consequences",
        "spliceai.ds_ag",
        "spliceai.ds_al",
        "spliceai.ds_dg",
        "spliceai.ds_dl",
    ]
    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        # Two transcripts: one high-scoring splice donor hit, one plain missense.
        writer.writerow({
            "variant_urn": "v1",
            "vep.mutational_consequences": "missense_variant^splice_acceptor_variant|missense_variant",
            "spliceai.ds_ag": "0.5|0.0",
            "spliceai.ds_al": "0.0|0.0",
            "spliceai.ds_dg": "0.0|0.0",
            "spliceai.ds_dl": "0.0|0.0",
        })
        # SpliceAI not computed at all; falls back to simplified_consequence alone.
        writer.writerow({
            "variant_urn": "v2",
            "vep.mutational_consequences": "missense_variant",
            "spliceai.ds_ag": "",
            "spliceai.ds_al": "",
            "spliceai.ds_dg": "",
            "spliceai.ds_dl": "",
        })

    mod.main([str(in_path), str(out_path), "--consequence-map-file", str(map_path)])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert rows[0]["simplified_consequence"] == "missense_variant|missense_variant"
    assert rows[0]["condensed_consequence"] == "missense_variant"
    assert rows[0]["splice_variant"] == "Yes|No"
    assert rows[0]["splice_var_amino"] == "Yes"

    assert rows[1]["simplified_consequence"] == "missense_variant"
    assert rows[1]["condensed_consequence"] == "missense_variant"
    assert rows[1]["splice_variant"] == "No"
    assert rows[1]["splice_var_amino"] == "No"
