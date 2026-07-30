# compare_cvfg_datasets

`src/compare_cvfg_datasets.py` — a development utility for cross-validating a CVFG pipeline output file against a manually-curated reference dataset.

This tool is not part of the annotation pipeline proper. It is used during pipeline development to detect variants present in one dataset but missing from the other, and to inspect where classifications or annotations agree or diverge.

---

## Inputs

| Positional | Description | Default |
|---|---|---|
| `file_a` | Integrated variant effect dataset (condensed TSV) | `data/cvfg/integrated_variant_effect_dataset.tsv` |
| `file_b` | CVFG variants flat TSV (pipeline output, after `flatten_dna_variants`) | `data/cvfg/v5/cvfg_variants.3.flat.tsv` |
| `file_c` | Dataset-to-score-set URN mapping TSV (maps `score_set_urn` → `dataset_name`) | `data/cvfg/dataset_to_score_set.tsv` |

All positional arguments are optional; Python defaults are used when omitted. Paths are passed through the run script which maps them to the container's `/work` directory.

---

## Join logic

Rows in A and B are matched using a compound genomic key:

| A column | B column |
|---|---|
| `Dataset` | `dataset_name` (added from C via `score_set_urn`) |
| `hg38_start` | `mapped_hgvs_g_start` |
| `hg38_end` | `mapped_hgvs_g_stop` |
| `ref_allele` | `mapped_hgvs_g_ref` |
| `alt_allele` | `mapped_hgvs_g_alt` |

Position columns are coerced to nullable integers before joining (A stores them as floats).

### Pre-filter: VCF/HGVS-c length mismatch

Before the join, A rows where the allele lengths implied by the VCF `ref_allele`/`alt_allele` fields are inconsistent with the lengths in the `hgvs_c` substitution notation are written to `compare_errors_in_a.tsv` and excluded from the join entirely.

### Relaxed fallback joins

Rows unmatched after the strict join are re-attempted with three successive relaxed strategies, applied to the residual unmatched set at each stage:

1. **Transcript-alt discrepancy** — For codon changes that span an intron, A stores only the substituted codon bases while B stores the full genomic alt. Rows are matched by position + ref only, then the match is confirmed by comparing A's `alt_allele` to B's `mapped_hgvs_c_alt`. On minus-strand genes (`strand == -1`) A's genomic alt is the reverse complement of the transcript alt.

2. **VCF-anchored deletion discrepancy** — A stores deletions with an extra anchor nucleotide included in both `ref_allele` and `alt_allele` (VCF convention), while B stores only the deleted bases with a blank alt. Three position conventions and both left- and right-anchor orientations are tried.

3. **Off-by-one start coordinate** — Some A rows have `len(ref_allele) != hg38_end - hg38_start + 1`. Trying `hg38_start + 1` and matching against B's transcript alt (with strand-aware reverse complement) resolves these.

### Post-filter: coordinate and ref-allele errors

Among the unmatched A rows remaining after all relaxed joins, rows with any of the following are written to `compare_only_in_a_coord_errors.tsv` separately from the clean unmatched set:
- `len(ref_allele) != hg38_end - hg38_start + 1`
- `transcript_alt` disagrees with what `hgvs_c` implies
- `ref_allele` does not match the GRCh38 reference sequence at the given coordinates (requires UTA; see `--check-ref-alleles`)

---

## Output files

All files are written to `--output-dir` (default: `data/cvfg`).

| File | Content |
|---|---|
| `compare_errors_in_a.tsv` | A rows excluded before the join due to VCF/HGVS-c length mismatch |
| `compare_only_in_a.tsv` | Clean A rows with no match in B (all A columns) |
| `compare_only_in_a_coord_errors.tsv` | Unmatched A rows with coordinate or ref-allele errors |
| `compare_only_in_b.tsv` | B rows with no match in A (all B columns) |
| `compare_joined.tsv` | Matched rows with selected columns from both A and B |

### Joined output columns

The joined file draws from both sources:

**From A:** `ID`, `clinvar_sig_2018`, `clinvar_sig_2025`, `gnomad_MAF`, `has_spliceai`, `Uuid_ClinGen_repo`, `auth_reported_func_class`, `consequence`, `AM_score`, `MutPred2`, `REVEL`

**From B:** `dataset_name`, `gene_symbol`, `variant_urn`, `raw_hgvs_nt`, `raw_hgvs_pro`, `mapped_hgvs_g`, `mapped_hgvs_c`, `mapped_hgvs_p`, `mapping_error`, `reverse_translation_error`, `clingen_allele_id`, `dna_clingen_allele_id`, `score`, ClinVar columns (2018/2025/2026), `gnomad.v4_1.minor_allele_frequency`, `spliceai.max_delta_score`, `clingen_evidence_repository.Uuid`, `mavedb.primary_calibration.functional_class_label`, `vep.mutational_consequence`, `alphamissense.pathogenicity`, `mutpred2.score`, `revel.score`

`has_spliceai` is a derived boolean added to A: `True` if any of the four SpliceAI DS columns (`spliceAI_DS_AG`, `spliceAI_DS_AL`, `spliceAI_DS_DG`, `spliceAI_DS_DL`) is non-empty.

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `--output-dir DIR` | `data/cvfg` | Directory for all output TSVs |
| `--check-ref-alleles` / `--no-check-ref-alleles` | enabled | Verify ref alleles against GRCh38 via UTA (requires `UTA_DB_URL`) |

---

## Usage

```bash
src/scripts/run_compare_cvfg_datasets.sh \
  data/cvfg/integrated_variant_effect_dataset.tsv \
  data/cvfg/v7/cvfg_variants.13.flat.tsv \
  data/cvfg/dataset_to_score_set.tsv

# Skip the ref-allele check (faster; no UTA required)
src/scripts/run_compare_cvfg_datasets.sh \
  data/cvfg/integrated_variant_effect_dataset.tsv \
  data/cvfg/v7/cvfg_variants.13.flat.tsv \
  data/cvfg/dataset_to_score_set.tsv \
  --no-check-ref-alleles
```

---

## Dependencies

- **`pandas`** — included in project dependencies.
- **`hgvs`** + `UTA_DB_URL` — optional; required only for the genomic ref-allele check (`--check-ref-alleles`). If UTA is unavailable the check is skipped with a warning.
