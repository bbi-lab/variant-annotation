#!/usr/bin/env bash
# Usage: src/scripts/generate_mavedb_and_igvf_files.sh [data-dir]
#
# Regenerates both the IGVF TSV and the MaveDB variant/scores/counts/functional-
# classes CSVs for every Excel/CSV/TSV source file in <data-dir> (default: data/rich).
#
# Shared steps (both outputs derive from the same corrected coding-variant data):
#
#   0. fix_missing_ref_allele.py    - reconstruct the occasional manual-entry typo where a
#                                      substitution is missing its reference letter (e.g.
#                                      "c.865-4>G" instead of "c.865-4C>G"), read directly
#                                      from the row's own seq/reference columns.
#   1. expand_unchanged_hgvs_c.py   - expand bare "c.=" (unchanged/WT control) cells
#                                      into "c.<start>_<end>=" by inferring the window
#                                      position from sibling rows' resolved substitutions.
#   2. translate_hgvs_accessions.py - normalize the coding-variant accession's version
#                                      separator to a dot (e.g. NM_032415_7 -> NM_032415.7),
#                                      keeping the RefSeq accession itself unchanged.
#   3. correct_reference_alleles.py - correct substitutions (coding, UTR, and intronic-offset
#                                      positions, including each component of a multi-variant
#                                      haplotype "c.[...;...]") whose reported reference allele
#                                      doesn't match the true RefSeq transcript/genomic sequence
#                                      (fixed edits baked into the assay target), e.g.
#                                      "c.462A>C" -> "c.462=" when the true reference at c.462
#                                      is C, not A. Runs via Docker (needs UTA). This pipeline
#                                      tabulates corrections by category (coding/intronic/UTR)
#                                      per dataset and writes three reports:
#                                        <data-dir>/correction_report_variants.tsv    variants corrected, by category
#                                        <data-dir>/correction_report_positions.tsv   distinct positions corrected, by category
#                                        <data-dir>/correction_details.tsv            dataset, original variant, corrected variant (one row per correction)
#                                      translate_hgvs_accessions.py also normalizes malformed
#                                      unbracketed multi-variant HGVS ("c.2523A>C;2533A>T") to
#                                      the correct bracketed form ("c.[2523A>C;2533A>T]").
#
# IGVF branch (from the shared, corrected data):
#
#   4. translate_hgvs_accessions.py - translate coding variant to Ensembl (MANE), keep
#                                      genomic variant as RefSeq (normalizing its version
#                                      separator too), and prepend an Ensembl protein
#                                      accession to the prot column.
#   5. filter_columns.py            - drop the seq and reference columns.
#      -> <stem>.igvf.tsv                  same columns as the source, coding/prot
#                                           translated to Ensembl
#
# MaveDB branch (from the shared, corrected data):
#
#   4. filter_columns.py            - drop the prot and genomic variant columns.
#   5. inversions_to_delins.py      - MAVE-HGVS doesn't support "inv" notation; rewrite any
#                                      inversion as an equivalent delins (e.g. "c.2687_2688inv"
#                                      -> "c.2687_2688delinsAA") by fetching the true reference
#                                      from UTA and reverse-complementing it. Runs via Docker.
#   6. dedupe_mavedb_variants.py    - MaveDB requires one row per variant; correction (or a
#                                      pre-existing source issue, e.g. overlapping tiled
#                                      windows) can produce duplicate coding-variant values.
#                                      For each duplicate group, ranks rows (not-haplotype-
#                                      and-not-corrected > not-haplotype > any) and keeps the
#                                      best, ties broken by file order; the rest go to
#                                      <stem>.mavedb.duplicates_removed.tsv.
#      -> <stem>.mavedb.csv
#   7. rename-columns (x3, via src.utilities) - split into scores/counts/functional_classes,
#      each renaming "coding variant" -> hgvs_nt and dropping seq/reference/variant_type/region:
#        <stem>.mavedb.scores.csv             hgvs_nt, Rep1_FC, Rep2_FC, Rep3_FC, score
#        <stem>.mavedb.counts.csv             hgvs_nt, every column containing "counts"
#        <stem>.mavedb.functional_classes.csv hgvs_nt, class_name

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

DATA_DIR="${1:-data/rich}"
MANE_FILE="data/MANE.GRCh38.v1.5.summary.txt"

shopt -s nullglob
source_files=("$DATA_DIR"/*.xlsx "$DATA_DIR"/*.csv "$DATA_DIR"/*.tsv)
shopt -u nullglob

# Exclude Excel lock files (~$foo.xlsx) and this pipeline's own generated output/
# report/intermediate files, so a rerun over a directory that already has output
# in it doesn't try to reprocess that output as a new source file.
filtered=()
for f in "${source_files[@]}"; do
  base="$(basename "$f")"
  case "$base" in
    '~$'*) continue ;;
    *.igvf.tsv|*.igvf.raw.tsv) continue ;;
    *.mavedb.csv|*.mavedb.*.csv|*.mavedb.*.tsv) continue ;;
    *.step0.csv|*.step1.csv|*.step2.csv|*.step3.csv) continue ;;
    *.correction_summary.txt|*.correction_details.tsv) continue ;;
    correction_report_variants.tsv|correction_report_positions.tsv|correction_details.tsv) continue ;;
  esac
  filtered+=("$f")
done
source_files=("${filtered[@]}")

if [[ ${#source_files[@]} -eq 0 ]]; then
  echo "error: no .xlsx/.csv/.tsv source files found in $DATA_DIR" >&2
  exit 1
fi

total_checked=0
total_corrected_rows=0
total_corrected_positions=0
total_corrected_rows_coding=0
total_corrected_positions_coding=0
total_corrected_rows_intronic=0
total_corrected_positions_intronic=0
total_corrected_rows_utr=0
total_corrected_positions_utr=0
summary_rows=()
variant_report_rows=()
position_report_rows=()

correction_details_report="$DATA_DIR/correction_details.tsv"
printf 'dataset\toriginal\tcorrected\n' > "$correction_details_report"

for src in "${source_files[@]}"; do
  stem="${src%.*}"
  step0="${stem}.step0.csv"
  step1="${stem}.step1.csv"
  step2="${stem}.step2.csv"
  step3="${stem}.step3.csv"
  correction_summary="${stem}.correction_summary.txt"
  correction_details="${stem}.correction_details.tsv"

  echo "=== $src ===" >&2

  # --- shared steps ---

  # Fix rare manual-entry typos where the reference letter is missing from a
  # substitution (e.g. "c.865-4>G" instead of "c.865-4C>G"), reconstructed from
  # the row's own seq/reference columns, before anything else touches the cell.
  python3 src/fix_missing_ref_allele.py "$src" "$step0"

  python3 src/expand_unchanged_hgvs_c.py "$step0" "$step1"
  rm -f "$step0"

  python3 src/translate_hgvs_accessions.py "$step1" "$step2" \
    --hgvs-c-col "coding variant" --hgvs-c-mode keep --normalize-version-separator

  src/scripts/run_correct_reference_alleles.sh "$step2" "$step3" \
    --coding-col "coding variant" \
    --summary-file "$correction_summary" --corrections-file "$correction_details"

  checked=0
  corrected_rows=0
  corrected_positions=0
  checked_coding=0; corrected_rows_coding=0; corrected_positions_coding=0
  checked_intronic=0; corrected_rows_intronic=0; corrected_positions_intronic=0
  checked_utr=0; corrected_rows_utr=0; corrected_positions_utr=0
  if [[ -f "$correction_summary" ]]; then
    # shellcheck disable=SC1090
    source "$correction_summary"
    rm -f "$correction_summary"
  fi
  total_checked=$((total_checked + checked))
  total_corrected_rows=$((total_corrected_rows + corrected_rows))
  total_corrected_positions=$((total_corrected_positions + corrected_positions))
  total_corrected_rows_coding=$((total_corrected_rows_coding + corrected_rows_coding))
  total_corrected_positions_coding=$((total_corrected_positions_coding + corrected_positions_coding))
  total_corrected_rows_intronic=$((total_corrected_rows_intronic + corrected_rows_intronic))
  total_corrected_positions_intronic=$((total_corrected_positions_intronic + corrected_positions_intronic))
  total_corrected_rows_utr=$((total_corrected_rows_utr + corrected_rows_utr))
  total_corrected_positions_utr=$((total_corrected_positions_utr + corrected_positions_utr))

  dataset_name="$(basename "$stem")"
  summary_rows+=("$dataset_name|$checked|$corrected_rows|$corrected_positions")
  variant_report_rows+=("$dataset_name|$corrected_rows_coding|$corrected_rows_intronic|$corrected_rows_utr|$corrected_rows")
  position_report_rows+=("$dataset_name|$corrected_positions_coding|$corrected_positions_intronic|$corrected_positions_utr|$corrected_positions")

  if [[ -f "$correction_details" ]]; then
    tail -n +2 "$correction_details" | while IFS=$'\t' read -r original corrected; do
      printf '%s\t%s\t%s\n' "$dataset_name" "$original" "$corrected"
    done >> "$correction_details_report"
    rm -f "$correction_details"
  fi

  rm -f "$step1"

  # --- IGVF branch ---

  igvf_raw="${stem}.igvf.raw.tsv"

  python3 src/translate_hgvs_accessions.py "$step3" "$igvf_raw" \
    --mane-file "$MANE_FILE" \
    --direction to-ensembl \
    --hgvs-c-col "coding variant" --hgvs-c-mode translate \
    --hgvs-g-col "genomic variant" --hgvs-g-mode keep \
    --hgvs-p-col prot --hgvs-p-mode add --hgvs-p-accession ensembl \
    --normalize-version-separator --drop-blank-columns

  python3 src/filter_columns.py "$igvf_raw" "${stem}.igvf.tsv" --omit-col seq,reference
  rm -f "$igvf_raw"

  # --- MaveDB branch ---

  mavedb_stage1="${stem}.mavedb.stage1.csv"
  mavedb_stage2="${stem}.mavedb.stage2.csv"
  mavedb_raw="${stem}.mavedb.raw.csv"
  mavedb_final="${stem}.mavedb.csv"
  mavedb_duplicates="${stem}.mavedb.duplicates_removed.tsv"

  python3 src/filter_columns.py "$step3" "$mavedb_stage1" --omit-col prot,"genomic variant"
  rm -f "$step3"

  # MAVE-HGVS doesn't support inversion ("inv") notation; rewrite as an
  # equivalent delins (needs UTA, so runs via Docker like the correction step).
  src/scripts/run_inversions_to_delins.sh "$mavedb_stage1" "$mavedb_stage2" \
    --hgvs-col "coding variant"
  rm -f "$mavedb_stage1"

  python3 src/filter_columns.py "$mavedb_stage2" "$mavedb_raw" --omit-col "orig_coding variant"
  rm -f "$mavedb_stage2"

  python3 src/dedupe_mavedb_variants.py "$mavedb_raw" "$step2" "$mavedb_final" "$mavedb_duplicates" \
    --coding-col "coding variant"
  rm -f "$mavedb_raw" "$step2"

  counts_cols=$(python3 -c "
import csv, sys
with open(sys.argv[1], newline='') as fh:
    header = next(csv.reader(fh, delimiter=','))
print(','.join(c for c in header if 'counts' in c))
" "$mavedb_final")

  python3 -m src.utilities rename-columns "$mavedb_final" "${stem}.mavedb.scores.csv" \
    --keep-col "coding variant:hgvs_nt" \
    --keep-col Rep1_FC --keep-col Rep2_FC --keep-col Rep3_FC \
    --keep-col "Functional Score:score" --reorder

  python3 -m src.utilities rename-columns "$mavedb_final" "${stem}.mavedb.counts.csv" \
    --keep-col "coding variant:hgvs_nt" --keep-col "$counts_cols" --reorder

  python3 -m src.utilities rename-columns "$mavedb_final" "${stem}.mavedb.functional_classes.csv" \
    --keep-col "coding variant:hgvs_nt" --keep-col "Classification:class_name" --reorder
done

echo >&2
echo "=== Reference-allele correction summary ===" >&2
printf '%-30s %10s %12s %12s\n' "Dataset" "Checked" "Corrected" "Positions" >&2
for row in "${summary_rows[@]}"; do
  IFS='|' read -r name row_checked row_corrected row_positions <<< "$row"
  printf '%-30s %10s %12s %12s\n' "$name" "$row_checked" "$row_corrected" "$row_positions" >&2
done
printf '%-30s %10s %12s %12s\n' "TOTAL" "$total_checked" "$total_corrected_rows" "$total_corrected_positions" >&2

# Two reports, by category (coding/intronic/utr): number of variants corrected,
# and the number of distinct positions that had a fixed-edit reference mismatch.
variant_report="$DATA_DIR/correction_report_variants.tsv"
position_report="$DATA_DIR/correction_report_positions.tsv"

{
  printf 'dataset\tcoding\tintronic\tutr\ttotal\n'
  for row in "${variant_report_rows[@]}"; do
    IFS='|' read -r name coding intronic utr total <<< "$row"
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$coding" "$intronic" "$utr" "$total"
  done
  printf 'TOTAL\t%s\t%s\t%s\t%s\n' "$total_corrected_rows_coding" "$total_corrected_rows_intronic" "$total_corrected_rows_utr" "$total_corrected_rows"
} > "$variant_report"

{
  printf 'dataset\tcoding\tintronic\tutr\ttotal\n'
  for row in "${position_report_rows[@]}"; do
    IFS='|' read -r name coding intronic utr total <<< "$row"
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$coding" "$intronic" "$utr" "$total"
  done
  printf 'TOTAL\t%s\t%s\t%s\t%s\n' "$total_corrected_positions_coding" "$total_corrected_positions_intronic" "$total_corrected_positions_utr" "$total_corrected_positions"
} > "$position_report"

echo >&2
echo "=== Corrected variants by category (also written to $variant_report) ===" >&2
printf '%-30s %10s %10s %10s %10s\n' "Dataset" "Coding" "Intronic" "UTR" "Total" >&2
for row in "${variant_report_rows[@]}"; do
  IFS='|' read -r name coding intronic utr total <<< "$row"
  printf '%-30s %10s %10s %10s %10s\n' "$name" "$coding" "$intronic" "$utr" "$total" >&2
done
printf '%-30s %10s %10s %10s %10s\n' "TOTAL" "$total_corrected_rows_coding" "$total_corrected_rows_intronic" "$total_corrected_rows_utr" "$total_corrected_rows" >&2

echo >&2
echo "=== Corrected positions by category (also written to $position_report) ===" >&2
printf '%-30s %10s %10s %10s %10s\n' "Dataset" "Coding" "Intronic" "UTR" "Total" >&2
for row in "${position_report_rows[@]}"; do
  IFS='|' read -r name coding intronic utr total <<< "$row"
  printf '%-30s %10s %10s %10s %10s\n' "$name" "$coding" "$intronic" "$utr" "$total" >&2
done
printf '%-30s %10s %10s %10s %10s\n' "TOTAL" "$total_corrected_positions_coding" "$total_corrected_positions_intronic" "$total_corrected_positions_utr" "$total_corrected_positions" >&2

echo >&2
echo "=== Per-variant correction details written to $correction_details_report ===" >&2
