"""
Unit tests for flatten_dna_variants module.

Tests the core functionality of expanding pipe-delimited DNA variants
to one row per DNA candidate.
"""

import io
import tempfile
from pathlib import Path
from typing import List

import pandas as pd
import pytest

from src.flatten_dna_variants import (
    flatten_dna_variants,
    get_dna_variant_columns,
    has_dna_variants,
)


class TestGetDnaVariantColumns:
    """Tests for auto-detecting DNA variant columns."""

    def test_detects_hardcoded_columns_when_present(self):
        """Should detect hard-coded column names."""
        df = pd.DataFrame({
            "mapped_hgvs_g": ["g.1A>T"],
            "mapped_hgvs_c": ["c.1A>T"],
            "dna_clingen_allele_id": ["CA123"],
            "gene_symbol": ["TP53"],
        })
        cols = get_dna_variant_columns(df)
        assert set(cols) == {"mapped_hgvs_g", "mapped_hgvs_c", "dna_clingen_allele_id"}

    def test_detects_spliceai_columns(self):
        """Should detect columns with spliceai. prefix."""
        df = pd.DataFrame({
            "mapped_hgvs_g": ["g.1A>T"],
            "spliceai.ds_ag": ["0.5"],
            "spliceai.max_delta_score": ["0.5"],
            "other_col": ["value"],
        })
        cols = get_dna_variant_columns(df)
        assert "mapped_hgvs_g" in cols
        assert "spliceai.ds_ag" in cols
        assert "spliceai.max_delta_score" in cols
        assert "other_col" not in cols

    def test_detects_clinvar_columns(self):
        """Should detect columns with clinvar. prefix."""
        df = pd.DataFrame({
            "mapped_hgvs_g": ["g.1A>T"],
            "clinvar.202601.clinical_significance": ["Pathogenic"],
            "other_col": ["value"],
        })
        cols = get_dna_variant_columns(df)
        assert "mapped_hgvs_g" in cols
        assert "clinvar.202601.clinical_significance" in cols
        assert "other_col" not in cols

    def test_detects_gnomad_columns(self):
        """Should detect columns with gnomad. prefix."""
        df = pd.DataFrame({
            "mapped_hgvs_g": ["g.1A>T"],
            "gnomad.v4.1.allele_frequency": ["0.01"],
            "other_col": ["value"],
        })
        cols = get_dna_variant_columns(df)
        assert "mapped_hgvs_g" in cols
        assert "gnomad.v4.1.allele_frequency" in cols
        assert "other_col" not in cols

    def test_detects_parsed_position_allele_columns(self):
        """Should detect columns for parsed positions and alleles."""
        df = pd.DataFrame({
            "mapped_hgvs_g": ["g.1A>T"],
            "mapped_hgvs_g_start": ["100"],
            "mapped_hgvs_g_stop": ["100"],
            "mapped_hgvs_g_ref": ["A"],
            "mapped_hgvs_g_alt": ["T"],
            "mapped_hgvs_c_start": ["10"],
            "mapped_hgvs_c_stop": ["10"],
            "mapped_hgvs_c_ref": ["A"],
            "mapped_hgvs_c_alt": ["T"],
            "other_col": ["value"],
        })
        cols = get_dna_variant_columns(df)
        assert "mapped_hgvs_g" in cols
        assert "mapped_hgvs_g_start" in cols
        assert "mapped_hgvs_g_stop" in cols
        assert "mapped_hgvs_g_ref" in cols
        assert "mapped_hgvs_g_alt" in cols
        assert "mapped_hgvs_c_start" in cols
        assert "other_col" not in cols

    def test_detects_intronic_region_columns(self):
        """Should detect intronic region flag columns."""
        df = pd.DataFrame({
            "mapped_hgvs_g": ["g.1A>T"],
            "touches_intronic_region": ["True"],
            "spans_intron": ["False"],
            "other_col": ["value"],
        })
        cols = get_dna_variant_columns(df)
        assert "touches_intronic_region" in cols
        assert "spans_intron" in cols
        assert "other_col" not in cols

    def test_returns_empty_list_for_no_dna_columns(self):
        """Should return empty list if no DNA columns present."""
        df = pd.DataFrame({
            "gene_symbol": ["TP53"],
            "other_col": ["value"],
        })
        cols = get_dna_variant_columns(df)
        assert cols == []


class TestHasDnaVariants:
    """Tests for checking if a row has DNA variants."""

    def test_returns_true_when_any_dna_column_has_value(self):
        """Should return True if any DNA variant column has a value."""
        row = pd.Series({
            "mapped_hgvs_g": "g.1A>T",
            "mapped_hgvs_c": "",
            "gene_symbol": "TP53",
        })
        dna_cols = ["mapped_hgvs_g", "mapped_hgvs_c"]
        assert has_dna_variants(row, dna_cols) is True

    def test_returns_false_when_all_dna_columns_empty(self):
        """Should return False if all DNA variant columns are empty."""
        row = pd.Series({
            "mapped_hgvs_g": "",
            "mapped_hgvs_c": "",
            "gene_symbol": "TP53",
        })
        dna_cols = ["mapped_hgvs_g", "mapped_hgvs_c"]
        assert has_dna_variants(row, dna_cols) is False

    def test_ignores_whitespace_in_columns(self):
        """Should treat whitespace-only values as empty."""
        row = pd.Series({
            "mapped_hgvs_g": "   ",
            "mapped_hgvs_c": "\t",
            "gene_symbol": "TP53",
        })
        dna_cols = ["mapped_hgvs_g", "mapped_hgvs_c"]
        assert has_dna_variants(row, dna_cols) is False


class TestFlattenDnaVariants:
    """Tests for the main flatten_dna_variants function."""

    def test_expands_single_candidate_row_unchanged(self):
        """Should return single row when there's only one candidate."""
        input_df = pd.DataFrame({
            "gene_symbol": ["TP53"],
            "mapped_hgvs_g": ["g.1A>T"],
            "mapped_hgvs_c": ["c.1A>T"],
            "dna_clingen_allele_id": ["CA123"],
        })

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_in:
            input_df.to_csv(f_in, sep="\t", index=False)
            input_path = Path(f_in.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_out:
            output_path = Path(f_out.name)

        try:
            flatten_dna_variants(input_path, output_path)
            output_df = pd.read_csv(output_path, sep="\t", dtype=str)

            assert len(output_df) == 1
            assert output_df.iloc[0]["mapped_hgvs_g"] == "g.1A>T"
            assert output_df.iloc[0]["gene_symbol"] == "TP53"
        finally:
            input_path.unlink()
            output_path.unlink()

    def test_expands_multi_candidate_row(self):
        """Should expand pipe-delimited candidates into separate rows."""
        input_df = pd.DataFrame({
            "gene_symbol": ["TP53"],
            "mapped_hgvs_g": ["g.1A>T|g.2A>T|g.3A>T"],
            "mapped_hgvs_c": ["c.1A>T|c.2A>T|c.3A>T"],
            "dna_clingen_allele_id": ["CA123|CA124|CA125"],
        })

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_in:
            input_df.to_csv(f_in, sep="\t", index=False)
            input_path = Path(f_in.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_out:
            output_path = Path(f_out.name)

        try:
            flatten_dna_variants(input_path, output_path)
            output_df = pd.read_csv(output_path, sep="\t", dtype=str)

            assert len(output_df) == 3
            assert output_df.iloc[0]["mapped_hgvs_g"] == "g.1A>T"
            assert output_df.iloc[1]["mapped_hgvs_g"] == "g.2A>T"
            assert output_df.iloc[2]["mapped_hgvs_g"] == "g.3A>T"
            # Non-list columns should be repeated
            assert all(output_df["gene_symbol"] == "TP53")
        finally:
            input_path.unlink()
            output_path.unlink()

    def test_preserves_empty_slots_in_misaligned_candidates(self):
        """Should preserve empty slots when candidates have different lengths."""
        input_df = pd.DataFrame({
            "gene_symbol": ["TP53"],
            "mapped_hgvs_g": ["g.1A>T|g.2A>T"],
            "mapped_hgvs_c": ["c.1A>T||c.3A>T"],  # 3 candidates, 2nd is empty
            "dna_clingen_allele_id": ["CA123||CA125"],
        })

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_in:
            input_df.to_csv(f_in, sep="\t", index=False)
            input_path = Path(f_in.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_out:
            output_path = Path(f_out.name)

        try:
            flatten_dna_variants(input_path, output_path)
            output_df = pd.read_csv(output_path, sep="\t", dtype=str)

            # Should expand to 3 rows (max candidates)
            assert len(output_df) == 3
            # Row 0: g.1A>T, c.1A>T, CA123
            assert output_df.iloc[0]["mapped_hgvs_g"] == "g.1A>T"
            assert output_df.iloc[0]["mapped_hgvs_c"] == "c.1A>T"
            assert output_df.iloc[0]["dna_clingen_allele_id"] == "CA123"
            # Row 1: g.2A>T, (empty), (empty)
            assert output_df.iloc[1]["mapped_hgvs_g"] == "g.2A>T"
            assert pd.isna(output_df.iloc[1]["mapped_hgvs_c"]) or output_df.iloc[1]["mapped_hgvs_c"] == ""
            assert pd.isna(output_df.iloc[1]["dna_clingen_allele_id"]) or output_df.iloc[1]["dna_clingen_allele_id"] == ""
            # Row 2: (beyond g, so empty), c.3A>T, CA125
            assert pd.isna(output_df.iloc[2]["mapped_hgvs_g"]) or output_df.iloc[2]["mapped_hgvs_g"] == ""
            assert output_df.iloc[2]["mapped_hgvs_c"] == "c.3A>T"
            assert output_df.iloc[2]["dna_clingen_allele_id"] == "CA125"
        finally:
            input_path.unlink()
            output_path.unlink()

    def test_drops_protein_only_rows_without_dna_variants(self):
        """Should drop rows with no DNA variant information."""
        input_df = pd.DataFrame({
            "gene_symbol": ["TP53", "BRCA1"],
            "mapped_hgvs_g": ["g.1A>T", ""],  # BRCA1 has no DNA variants
            "mapped_hgvs_c": ["c.1A>T", ""],
            "mapped_hgvs_p": ["p.Arg175His", "p.Gln18*"],
        })

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_in:
            input_df.to_csv(f_in, sep="\t", index=False)
            input_path = Path(f_in.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_out:
            output_path = Path(f_out.name)

        try:
            flatten_dna_variants(input_path, output_path)
            output_df = pd.read_csv(output_path, sep="\t", dtype=str)

            # Should only have TP53 row (BRCA1 dropped)
            assert len(output_df) == 1
            assert output_df.iloc[0]["gene_symbol"] == "TP53"
        finally:
            input_path.unlink()
            output_path.unlink()

    def test_handles_annotation_columns(self):
        """Should expand annotation columns like spliceai.* and clinvar.*."""
        input_df = pd.DataFrame({
            "gene_symbol": ["TP53"],
            "mapped_hgvs_g": ["g.1A>T|g.2A>T"],
            "spliceai.ds_ag": ["0.5|0.7"],
            "spliceai.max_delta_score": ["0.5|0.8"],
            "clinvar.202601.clinical_significance": ["Pathogenic|Benign"],
        })

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_in:
            input_df.to_csv(f_in, sep="\t", index=False)
            input_path = Path(f_in.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_out:
            output_path = Path(f_out.name)

        try:
            flatten_dna_variants(input_path, output_path)
            output_df = pd.read_csv(output_path, sep="\t", dtype=str)

            assert len(output_df) == 2
            # Row 0
            assert output_df.iloc[0]["mapped_hgvs_g"] == "g.1A>T"
            assert output_df.iloc[0]["spliceai.ds_ag"] == "0.5"
            assert output_df.iloc[0]["spliceai.max_delta_score"] == "0.5"
            assert output_df.iloc[0]["clinvar.202601.clinical_significance"] == "Pathogenic"
            # Row 1
            assert output_df.iloc[1]["mapped_hgvs_g"] == "g.2A>T"
            assert output_df.iloc[1]["spliceai.ds_ag"] == "0.7"
            assert output_df.iloc[1]["spliceai.max_delta_score"] == "0.8"
            assert output_df.iloc[1]["clinvar.202601.clinical_significance"] == "Benign"
        finally:
            input_path.unlink()
            output_path.unlink()

    def test_raises_on_nonexistent_input_file(self):
        """Should raise FileNotFoundError if input file doesn't exist."""
        input_path = Path("/nonexistent/file.tsv")
        output_path = Path("/tmp/output.tsv")

        with pytest.raises(FileNotFoundError):
            flatten_dna_variants(input_path, output_path)

    def test_raises_on_empty_input_file(self):
        """Should raise ValueError for empty input file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_in:
            f_in.write("")  # Empty file
            input_path = Path(f_in.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_out:
            output_path = Path(f_out.name)

        try:
            with pytest.raises(ValueError, match="empty"):
                flatten_dna_variants(input_path, output_path)
        finally:
            input_path.unlink()
            if output_path.exists():
                output_path.unlink()

    def test_raises_when_all_rows_have_no_dna_variants(self):
        """Should raise ValueError when all rows are protein-only."""
        input_df = pd.DataFrame({
            "gene_symbol": ["TP53", "BRCA1"],
            "mapped_hgvs_g": ["", ""],
            "mapped_hgvs_c": ["", ""],
            "mapped_hgvs_p": ["p.Arg175His", "p.Gln18*"],
        })

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_in:
            input_df.to_csv(f_in, sep="\t", index=False)
            input_path = Path(f_in.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_out:
            output_path = Path(f_out.name)

        try:
            with pytest.raises(ValueError, match="No rows with DNA variants"):
                flatten_dna_variants(input_path, output_path)
        finally:
            input_path.unlink()
            if output_path.exists():
                output_path.unlink()

    def test_handles_explicit_column_specification(self):
        """Should only expand explicitly specified columns."""
        input_df = pd.DataFrame({
            "gene_symbol": ["TP53"],
            "mapped_hgvs_g": ["g.1A>T|g.2A>T"],
            "mapped_hgvs_c": ["c.1A>T|c.2A>T"],
            "spliceai.ds_ag": ["0.5|0.7"],  # Should NOT expand (not in DNA cols)
        })

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_in:
            input_df.to_csv(f_in, sep="\t", index=False)
            input_path = Path(f_in.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f_out:
            output_path = Path(f_out.name)

        try:
            # Only expand mapped_hgvs_g and mapped_hgvs_c
            flatten_dna_variants(
                input_path,
                output_path,
                dna_variant_columns=["mapped_hgvs_g", "mapped_hgvs_c"],
            )
            output_df = pd.read_csv(output_path, sep="\t", dtype=str)

            assert len(output_df) == 2
            # spliceai.ds_ag should NOT be expanded (still pipe-delimited)
            assert output_df.iloc[0]["spliceai.ds_ag"] == "0.5|0.7"
        finally:
            input_path.unlink()
            output_path.unlink()
