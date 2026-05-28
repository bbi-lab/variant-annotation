# remap_transcript_ids

`src/remap_transcript_ids.py` — optional utility step.

Replaces Ensembl transcript or protein accessions with their RefSeq equivalents (or vice versa) in the HGVS columns and derived accession columns produced by the pipeline.  The mapping is restricted to identifiers that are guaranteed to be sequence-identical across namespaces: by default this means MANE Select (and optionally MANE Plus Clinical) transcripts, but a fully custom accession mapping file may also be supplied.

---

## What it does

For each row, the following columns are updated when their accession appears in the mapping:

| Column | How it is updated |
|---|---|
| `mapped_hgvs_g` | Accession prefix before `:` replaced |
| `mapped_hgvs_c` | Accession prefix before `:` replaced |
| `mapped_hgvs_p` | Accession prefix before `:` replaced |
| `mapped_hgvs_c_transcript` | Bare accession value replaced |
| `mapped_hgvs_p_protein` | Bare accession value replaced |

All other columns are written unchanged.  Pipe-delimited cells (multiple candidates per cell, as produced by `reverse_translate_protein_variants`) are handled element-by-element, so candidate cardinality is preserved.

Accession matching uses the full versioned string first (e.g. `NM_007294.4`).  If no exact match is found the version suffix is stripped and the base identifier is tried; the target version from the mapping table is then used in the output (e.g. `ENST00000357654.9`).

---

## Mapping sources

### MANE Select file (recommended)

Download from NCBI:

```
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.*.summary.txt.gz
```

The file contains both `MANE Select` and `MANE Plus Clinical` transcripts, all of which are sequence-identical between RefSeq and Ensembl.  Pass it with `--mane-file` and specify `--direction`.

Expected column headers (tab-delimited):

| Column | Description |
|---|---|
| `RefSeq_nuc` | RefSeq transcript accession (e.g. `NM_007294.4`) |
| `RefSeq_prot` | RefSeq protein accession (e.g. `NP_009225.1`) |
| `Ensembl_nuc` | Ensembl transcript accession (e.g. `ENST00000357654.9`) |
| `Ensembl_prot` | Ensembl protein accession (e.g. `ENSP00000350283.3`) |
| `MANE_status` | `MANE Select` or `MANE Plus Clinical` |

### Custom mapping file

A two-column TSV/CSV with headers `source_id` and `target_id`.  Direction is implied by which namespace appears in the `source_id` column.

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `input` | — | Input TSV/CSV file |
| `output` | — | Output TSV/CSV file |
| `--mane-file FILE` | — | NCBI MANE summary file (mutually exclusive with `--mapping-file`) |
| `--mapping-file FILE` | — | Custom two-column mapping file (mutually exclusive with `--mane-file`) |
| `--direction {to-refseq,to-ensembl}` | `to-refseq` | Swap direction (only used with `--mane-file`) |
| `--hgvs-col COL` | *(see above)* | Additional or override HGVS column; may be repeated |
| `--accession-col COL` | *(see above)* | Additional or override plain-accession column; may be repeated |
| `--csv-field-size-limit N` | `131072` | Maximum CSV field size |
| `--log-level LEVEL` | `WARNING` | Logging verbosity |

---

## Usage

```bash
# Remap to RefSeq (default direction)
src/scripts/run_remap_transcript_ids.sh input.tsv output.tsv \
  --mane-file data/MANE.GRCh38.v1.3.summary.txt.gz \
  --direction to-refseq

# Remap to Ensembl
src/scripts/run_remap_transcript_ids.sh input.tsv output.tsv \
  --mane-file data/MANE.GRCh38.v1.3.summary.txt.gz \
  --direction to-ensembl

# Custom mapping file
src/scripts/run_remap_transcript_ids.sh input.tsv output.tsv \
  --mapping-file data/custom_accession_map.tsv
```

---

## Notes

- The genomic column (`mapped_hgvs_g`) uses `NC_` RefSeq accessions that have no Ensembl equivalent; in practice those entries will never be remapped by a MANE file.  The column is processed anyway in case a custom mapping is supplied.
- If a column listed in `--hgvs-col` or `--accession-col` is absent from the input file it is silently skipped.
- The script is idempotent: running it twice in the same direction produces the same output.
