# diff_vcf_coords

`src/diff_vcf_coords.py` — compare the VCF-style coordinate columns between two versions of a variants TSV file.

Joins old and new files on `variant_urn` and reports every coordinate component that changed across three groups: genomic (GRCh38), transcript (coding), and protein. Each changed index position within a pipe-delimited multi-candidate row produces one output row.

Runs locally using the project's Python virtual environment (not Docker).

---

## Output

One TSV row per changed group/index position:

| Column | Description |
|---|---|
| `variant_urn` | Variant identifier |
| `group` | `genome`, `transcript`, or `protein` |
| `index` | 0-based position within the pipe-delimited value |
| `old_hgvs` / `new_hgvs` | HGVS string for the group at this position |
| `old_chromosome` / `new_chromosome` | Chromosome or transcript accession |
| `old_start` / `new_start` | Start coordinate |
| `old_stop` / `new_stop` | Stop coordinate |
| `old_ref` / `new_ref` | Reference allele |
| `old_alt` / `new_alt` | Alternate allele |

Variants present in only one file appear with the missing side's values left blank.

Output defaults to `<new_tsv_stem>.vcf_diffs.tsv` in the current directory.

---

## Columns compared

| Group | Columns |
|---|---|
| `genome` | `mapped_hgvs_g`, `mapped_hgvs_g_chromosome`, `mapped_hgvs_g_start`, `mapped_hgvs_g_stop`, `mapped_hgvs_g_ref`, `mapped_hgvs_g_alt` |
| `transcript` | `mapped_hgvs_c`, `mapped_hgvs_c_transcript`, `mapped_hgvs_c_start`, `mapped_hgvs_c_stop`, `mapped_hgvs_c_ref`, `mapped_hgvs_c_alt` |
| `protein` | `mapped_hgvs_p`, `mapped_hgvs_p_chromosome`, `mapped_hgvs_p_start`, `mapped_hgvs_p_stop`, `mapped_hgvs_p_ref`, `mapped_hgvs_p_alt` |

Only columns present in both files are compared. Missing groups are silently skipped.

---

## Usage

```bash
src/scripts/run_diff_vcf_coords.sh old.tsv new.tsv
src/scripts/run_diff_vcf_coords.sh old.tsv new.tsv output.tsv
```

| Argument | Description |
|---|---|
| `old_tsv` | Baseline TSV file |
| `new_tsv` | Updated TSV file |
| `output.tsv` | Optional output path (third positional argument) |
