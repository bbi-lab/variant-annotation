# filter_dna_variants

`src/filter_dna_variants.py` — optional filtering step for DNA-candidate lists.

This script keeps the number of rows unchanged, but removes DNA candidates that are not in the retained classes. It is intended for post-processing rows that already contain pipe-delimited DNA candidates from reverse translation and downstream annotation.

By default the script keeps:
- SNVs
- wild-type / unchanged candidates

Rows that lose all DNA candidates are still written. If such a row also has no `mapped_hgvs_p`, the script logs a warning count at the end.

---

## Behaviour

### Candidate selection

The filter uses `mapped_hgvs_c` when that column contains any non-empty candidate. If `mapped_hgvs_c` is absent or completely blank for the row, the script falls back to `mapped_hgvs_g`.

Each candidate is classified independently, and candidates not in the retained classes are removed from every aligned list field at the same position.

### Columns that are compacted

The script uses the same DNA-list column set as `flatten_dna_variants.py`:

- `mapped_hgvs_g`
- `mapped_hgvs_c`
- `reverse_translation_warnings`
- `mapped_hgvs_g_chromosome`
- `mapped_hgvs_g_start`
- `mapped_hgvs_g_stop`
- `mapped_hgvs_g_ref`
- `mapped_hgvs_g_alt`
- `mapped_hgvs_c_transcript`
- `mapped_hgvs_c_start`
- `mapped_hgvs_c_stop`
- `mapped_hgvs_c_ref`
- `mapped_hgvs_c_alt`
- `touches_intronic_region`
- `spans_intron`
- `dna_clingen_allele_id`
- Annotation columns beginning with `alphamissense.`, `clingen_evidence_repository.`, `clinvar.`, `gnomad.`, `mutpred2.`, `spliceai.`, `revel.`, or `vep.`

Non-list columns are copied through unchanged.

### Retained classes

Supported classes are:
- `snv`
- `wt`
- `del`
- `ins`
- `dup`
- `inv`
- `delins`
- `complex`
- `other`
- `all`

`all` disables filtering and keeps every candidate.

---

## Usage

```bash
# Default: keep SNVs and wild-type / unchanged candidates
src/scripts/run_filter_dna_variants.sh annotated.tsv filtered.tsv

# Keep only SNVs
src/scripts/run_filter_dna_variants.sh annotated.tsv filtered.tsv --keep-class snv

# Keep SNVs and WT explicitly
src/scripts/run_filter_dna_variants.sh annotated.tsv filtered.tsv \
  --keep-class snv --keep-class wt
```

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `--keep-class CLASS` | `snv,wt` | Candidate class to retain. May be repeated and/or comma-separated. |
| `--delimiter` | `\t` | Input/output delimiter |
| `--skip N` | `0` | Skip first N data rows |
| `--limit N` | no limit | Process at most N rows (after skip) |
| `--csv-field-size-limit BYTES` | system default | Increase for large HGVS fields |
| `--log-level` | `INFO` | Logging verbosity |

---

## Troubleshooting

**Rows end up with no DNA candidates**

This usually means none of the candidates matched the retained classes. If the row also has no `mapped_hgvs_p`, the script reports a warning count at the end.

**Too many or too few candidates are retained**

The filter class is inferred from the HGVS string itself, prioritizing `mapped_hgvs_c` and then falling back to `mapped_hgvs_g`. Use `--keep-class all` to disable filtering if you want to inspect the candidate lists first.
