# annotate_erepo

`src/annotate_erepo.py` — Step 8 of the variant-annotation pipeline (optional).

Annotates each variant row with expert-panel variant classifications from the [ClinGen Evidence Repository](https://erepo.clinicalgenome.org/) (erepo). The full classification table is downloaded once from the erepo API, cached locally, and then joined against input rows using up to three keys.

---

## Output columns

All output columns are prefixed `clingen_evidence_repository.` and are pipe-delimited, aligned to the `mapped_hgvs_c` candidate positions.

| Column | Description |
|---|---|
| `clingen_evidence_repository.ClinVar Variation Id` | ClinVar Variation ID from the erepo record |
| `clingen_evidence_repository.Allele Registry Id` | ClinGen Allele Registry ID (CAID) from the erepo record |
| `clingen_evidence_repository.Disease Mondo Id` | Mondo disease ontology ID |
| `clingen_evidence_repository.Mode of Inheritance` | Inheritance pattern (e.g. `Autosomal dominant`) |
| `clingen_evidence_repository.Assertion` | ACMG/AMP classification (e.g. `Pathogenic`, `Likely Benign`) |
| `clingen_evidence_repository.Applied Evidence Codes (Met)` | ACMG evidence codes that were met |
| `clingen_evidence_repository.Applied Evidence Codes (Not Met)` | ACMG evidence codes that were not met |
| `clingen_evidence_repository.Summary of interpretation` | Free-text interpretation summary |
| `clingen_evidence_repository.PubMed Articles` | Supporting PubMed IDs |
| `clingen_evidence_repository.Expert Panel` | Name of the ClinGen Expert Panel that issued the classification |
| `clingen_evidence_repository.Guideline` | Guideline version applied |
| `clingen_evidence_repository.Approval Date` | Date the classification was approved |
| `clingen_evidence_repository.Published Date` | Date the classification was published |
| `clingen_evidence_repository.Retracted` | Whether the classification has been retracted |
| `clingen_evidence_repository.Evidence Repo Link` | URL to the erepo classification record |
| `clingen_evidence_repository.Uuid` | Unique classification identifier in erepo |
| `clingen_evidence_repository.warnings` | Cross-key discrepancy notes (see [Join key diagnostics](#join-key-diagnostics)) |

When a single variant candidate matches multiple erepo records (e.g. the same variant classified by two different expert panels), the values from all matching records are merged within each field, separated by ` | `.

---

## Join keys

The script supports three independent join keys, all used by default:

| Key | Source column | Erepo column | Notes |
|---|---|---|---|
| `hgvs` | `mapped_hgvs_c` (one pipe-delimited candidate at a time) | `HGVS Expressions` (comma-separated) | Most reliable key; erepo HGVS expressions are pre-normalised across multiple representations |
| `clinvar` | `--clinvar-variation-id-col` (default: `clinvar.202601.variation_id`) | `ClinVar Variation Id` | Requires a ClinVar annotation step before erepo |
| `caid` | `--caid-col` (default: `dna_clingen_allele_id`) | `Allele Registry Id` | Requires step 3 (add_dna_clingen_allele_ids) before erepo |

All three keys are tried for each candidate. Results are deduplicated by erepo UUID; the same record is never emitted twice for the same candidate even if it is found via multiple keys.

Use `--join-keys hgvs,caid` (or any subset) to restrict which keys are used.

### HGVS normalisation

The `HGVS Expressions` column in the erepo TSV does not include gene-symbol annotations. Input `mapped_hgvs_c` values that carry a parenthesised gene symbol (e.g. `NM_004333.5(BRAF):c.740T>C`) are stripped to `NM_004333.5:c.740T>C` before lookup.

### Join key diagnostics

When a candidate is found via one key but not another (and the missing key's column is populated), a diagnostic note is written to `clingen_evidence_repository.warnings`. When two keys return non-overlapping records, both record sets are merged and the discrepancy is noted. These warnings are informational; they do not suppress annotation.

---

## Cache

The full erepo classification TSV is downloaded once from:

```
https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download
```

and cached in `--cache-dir` (default: `$EREPO_CACHE_DIR` or `/tmp/erepo_cache`) as `clingen_erepo_classifications.tsv`. All three in-memory indexes (by HGVS, by ClinVar Variation ID, by CAID) are built from this file at startup.

Use `--refresh-cache` to re-download the file on the next run even if a cached copy exists. Because the erepo is updated periodically, refresh the cache when running against new data or after significant erepo releases.

---

## Usage

```bash
src/scripts/run_annotate_erepo.sh input.tsv output.tsv
```

With explicit join-key selection (HGVS and CAID only):
```bash
src/scripts/run_annotate_erepo.sh input.tsv output.tsv \
  --join-keys hgvs,caid
```

Refresh the cached erepo TSV:
```bash
src/scripts/run_annotate_erepo.sh input.tsv output.tsv \
  --refresh-cache
```

---

## CLI options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--cache-dir DIR` | `EREPO_CACHE_DIR` | `/tmp/erepo_cache` | Directory for the cached erepo TSV |
| `--refresh-cache` | — | off | Re-download the erepo TSV even if cached |
| `--join-keys KEY,...` | — | `hgvs,clinvar,caid` | Comma-separated join keys (allowed: `hgvs`, `clinvar`, `caid`) |
| `--mapped-hgvs-c-col COL` | — | `mapped_hgvs_c` | Input column with pipe-delimited transcript HGVS strings |
| `--clinvar-variation-id-col COL` | — | `clinvar.202601.variation_id` | Input column with ClinVar Variation IDs |
| `--caid-col COL` | — | `dna_clingen_allele_id` | Input column with ClinGen Allele Registry IDs |
| `--skip N` | — | `0` | Skip first N data rows |
| `--limit N` | — | no limit | Stop after N rows |
| `--log-level` | — | `INFO` | Logging verbosity |
| `--csv-field-size-limit BYTES` | — | system default | Increase for large HGVS fields |
| `--delimiter` | — | `\t` | Input/output delimiter |

---

## Dependencies

- **`requests`** — used to download the erepo TSV (included in project dependencies).
- **Network access** — required on the first run (or with `--refresh-cache`) to reach `erepo.clinicalgenome.org`. Subsequent runs use the cached file.

---

## Troubleshooting

**All output columns are empty despite populated `mapped_hgvs_c`**

- Check that the HGVS strings in the input match the pre-normalised form used by erepo. Gene-symbol annotations are stripped automatically (`NM_004333.5(BRAF):c.740T>C` → `NM_004333.5:c.740T>C`), but other formatting differences (e.g. spacing, version mismatches) are not corrected.
- The erepo only covers variants that have received a formal expert-panel classification. Many variants will legitimately have no match.

**`clinvar` join key yields no results**

- Verify that `--clinvar-variation-id-col` points to a populated column. The default column (`clinvar.202601.variation_id`) requires step 5 (annotate_clinvar) to have been run. Use `--clinvar-variation-id-col clinvar.202501.variation_id` if an earlier ClinVar version was used.

**Warnings column shows "found via 'hgvs' but not 'caid'"**

- This is expected when a variant has no ClinGen Allele Registry ID (e.g. a novel variant, or one not yet submitted). The HGVS-based match is still used. The warning is informational only.

**Download fails or times out**

- The erepo API has a 120-second request timeout. On slow connections or when the erepo service is under load, the download may fail. Retry or use `--refresh-cache` on the next run.
- The TSV download URL does not require authentication.
