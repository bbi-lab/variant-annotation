#!/usr/bin/env python3
"""
Flatten pipe-delimited DNA variants to one row per DNA candidate.

This script takes annotated variant data with pipe-delimited DNA candidates
and expands it so each DNA candidate gets its own row. Protein-only variants
without DNA reverse translations are dropped.

Columns that can contain multiple pipe-delimited values (DNA variants,
ClinGen IDs, annotation scores) are split across rows. Non-list columns
are repeated as-is for each DNA candidate.

Example:
    Input row with 2 DNA candidates:
        mapped_hgvs_g    mapped_hgvs_c    spliceai.ds_ag
        c.1A>T|c.2A>T    g.100A>T|g.101A>T    0.5|0.6

    Output (2 rows):
        mapped_hgvs_g    mapped_hgvs_c    spliceai.ds_ag
        c.1A>T           g.100A>T         0.5
        c.2A>T           g.101A>T         0.6
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, List, Set

import pandas as pd


def get_dna_variant_columns(df: pd.DataFrame) -> List[str]:
    """
    Auto-detect columns that contain DNA variant information.

    Includes:
    - Hard-coded DNA variant columns (mapped_hgvs_g, mapped_hgvs_c, dna_clingen_allele_id)
    - Parsed position/allele columns (mapped_hgvs_*_{start,stop,ref,alt})
    - Intronic region flags (touches_intronic_region, spans_intron)
    - Any columns with known annotation prefixes (spliceai., clinvar., gnomad.)

    Args:
        df: Input DataFrame.

    Returns:
        List of column names that may contain pipe-delimited DNA variants.
    """
    # Hard-coded columns
    dna_cols = [
        "mapped_hgvs_g",
        "mapped_hgvs_c",
        "dna_clingen_allele_id",
        "touches_intronic_region",
        "spans_intron",
    ]

    # Annotation prefixes that get pipe-delimited for multi-candidate rows
    annotation_prefixes = ("spliceai.", "clinvar.", "gnomad.")

    # Add any columns with annotation prefixes
    for col in df.columns:
        if col.startswith(annotation_prefixes):
            dna_cols.append(col)
    
    # Add parsed position/allele columns (mapped_hgvs_*_start, etc.)
    for col in df.columns:
        if col.startswith("mapped_hgvs_") and col.endswith(("_start", "_stop", "_ref", "_alt")):
            dna_cols.append(col)

    # Return only columns that actually exist in the DataFrame
    return [c for c in dna_cols if c in df.columns]


def has_dna_variants(row: pd.Series, dna_cols: List[str]) -> bool:
    """
    Check if a row has any DNA variant information.

    Returns True if at least one DNA variant column has a non-empty value.
    This filters out protein-only variants without reverse translations.

    Args:
        row: A DataFrame row (Series).
        dna_cols: List of DNA variant column names.

    Returns:
        True if row has DNA variants, False otherwise.
    """
    return any(str(row[col]).strip() for col in dna_cols)


def flatten_dna_variants(
    input_file: Path,
    output_file: Path,
    dna_variant_columns: Optional[List[str]] = None,
) -> None:
    """
    Expand pipe-delimited DNA variants to one row per DNA candidate.

    Reads input TSV, identifies rows with DNA variants (drops protein-only rows),
    splits pipe-delimited columns, and writes one row per DNA candidate.

    Args:
        input_file: Input TSV file path.
        output_file: Output TSV file path.
        dna_variant_columns: Explicit list of columns to expand. If None, auto-detects.

    Raises:
        FileNotFoundError: If input file does not exist.
        ValueError: If input file is empty or has no data rows.
    """
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Read input TSV
    try:
        df = pd.read_csv(input_file, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        raise ValueError("Input file is empty or has no columns")

    if df.empty:
        raise ValueError("Input file is empty or has no data rows")

    # Get columns to expand
    if dna_variant_columns is None:
        dna_variant_columns = get_dna_variant_columns(df)
    else:
        # Filter to only columns that exist
        dna_variant_columns = [c for c in dna_variant_columns if c in df.columns]

    if not dna_variant_columns:
        raise ValueError("No DNA variant columns found in input file")

    # Filter to rows with DNA variants (drop protein-only without reverse translation)
    df = df[df.apply(lambda row: has_dna_variants(row, dna_variant_columns), axis=1)]

    if df.empty:
        raise ValueError(
            "No rows with DNA variants found. "
            "All rows appear to be protein-only without reverse translations."
        )

    # Expand rows
    expanded_rows = []
    for _, row in df.iterrows():
        # Find max number of candidates across all DNA variant columns
        max_candidates = 1
        for col in dna_variant_columns:
            cell_value = str(row[col]).strip()
            if cell_value:
                candidates = [c.strip() for c in cell_value.split("|")]
                max_candidates = max(max_candidates, len(candidates))

        # Create one expanded row per candidate
        for i in range(max_candidates):
            new_row = row.copy()

            # Split DNA variant columns
            for col in dna_variant_columns:
                cell_value = str(row[col]).strip()
                if cell_value:
                    candidates = [c.strip() for c in cell_value.split("|")]
                    # Use candidate i if it exists, otherwise empty string
                    new_row[col] = candidates[i] if i < len(candidates) else ""
                else:
                    new_row[col] = ""

            expanded_rows.append(new_row)

    expanded_df = pd.DataFrame(expanded_rows)
    expanded_df.to_csv(output_file, sep="\t", index=False)


def main() -> int:
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        type=Path,
        help="Input TSV file (from annotation pipeline)",
    )

    parser.add_argument(
        "output_file",
        type=Path,
        help="Output TSV file (one row per DNA variant)",
    )

    parser.add_argument(
        "--dna-variant-columns",
        type=str,
        default=None,
        help=(
            "Comma-separated list of columns to expand on pipes. "
            "If not provided, auto-detects based on column names."
        ),
    )

    args = parser.parse_args()

    dna_cols = None
    if args.dna_variant_columns:
        dna_cols = [c.strip() for c in args.dna_variant_columns.split(",")]

    try:
        flatten_dna_variants(args.input_file, args.output_file, dna_cols)
        print(f"✓ Flattened variants written to {args.output_file}")
        return 0
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
