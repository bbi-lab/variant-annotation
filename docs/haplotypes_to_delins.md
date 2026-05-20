# haplotypes_to_delins

`src/normalize_intra_codon_haplotypes.py` — rewrite intra-codon c.-haplotype HGVS strings to delins notation.

Some variant datasets represent compound substitutions within a single codon as a haplotype expression: `c.[123A>G;125T>C]`. HGVS parsers and downstream tools generally expect the canonical delins form instead: `c.123_125delinsGCC`. This script rewrites the `raw_hgvs_nt` column in-place for rows that match this pattern, preserving the original value in a new `orig_raw_hgvs_nt` column.

Only intra-codon haplotypes (all substitutions within the same codon) are rewritten. Haplotypes that span codon boundaries, contain non-substitution operations (insertions, deletions), or use UTR/intronic coordinates are left unchanged.

Runs in Docker.

---

## What changes

For rewritten rows:

- `raw_hgvs_nt` — replaced with `c.<start>_<end>delins<ALT>` (or accession-prefixed if the original had one, e.g. `NM_000001.1:c.123_125delinsGCC`)
- `orig_raw_hgvs_nt` — new column inserted immediately before `raw_hgvs_nt`; contains the original haplotype expression

All other rows and columns are written unchanged.

---

## Usage

```bash
src/scripts/run_haplotypes_to_delins.sh input.tsv output.tsv
src/scripts/run_haplotypes_to_delins.sh input.tsv output.tsv --raw-hgvs-nt-col raw_hgvs_nt
```

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `--raw-hgvs-nt-col COL` | `raw_hgvs_nt` | Input column containing the haplotype HGVS strings to rewrite |
| `--target-sequence-col COL` | `target_sequence` | Column containing the coding reference sequence (used to fill in non-mutated positions within the codon) |
| `--orig-raw-hgvs-nt-col COL` | `orig_raw_hgvs_nt` | Name of the new column that will store the original value |
| `--skip N` | `0` | Skip first N data rows |
| `--limit N` | no limit | Process at most N rows |
| `--log-level` | `INFO` | Logging verbosity |
| `--csv-field-size-limit BYTES` | system default | Increase for large fields |
