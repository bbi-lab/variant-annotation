# annotate_spliceai

`src/annotate_spliceai.py` — Step 7 of the variant-annotation pipeline (optional).

Annotates each variant row with SpliceAI splice-impact delta scores. Supports two execution modes: `precomputed` (default) — fast lookup against pre-scored tabix-indexed VCFs — and `compute` — local on-the-fly scoring using the SpliceAI deep-learning model.

---

## Output columns

| Column | Description |
|---|---|
| `spliceai.ds_ag` | Delta score — acceptor gain |
| `spliceai.ds_al` | Delta score — acceptor loss |
| `spliceai.ds_dg` | Delta score — donor gain |
| `spliceai.ds_dl` | Delta score — donor loss |
| `spliceai.dp_ag` | Delta position — acceptor gain |
| `spliceai.dp_al` | Delta position — acceptor loss |
| `spliceai.dp_dg` | Delta position — donor gain |
| `spliceai.dp_dl` | Delta position — donor loss |
| `spliceai.max_delta_score` | `max(DS_AG, DS_AL, DS_DG, DS_DL)` |

All scores are formatted to four decimal places. All nine columns are pipe-delimited and candidate-aligned when `mapped_hgvs_g` contains multiple pipe-separated candidates (from step 2 reverse translation).

When a gnomAD record has multiple SpliceAI `INFO` entries (e.g. multiple ALT alleles), the maximum value is taken per score field, filtered to the matching ALT allele when possible.

---

## Execution modes

### Precomputed mode (default, recommended)

Looks up each variant in user-provided bgzipped, tabix-indexed SpliceAI score VCFs (the Illumina SpliceAI pre-scored data release). Tabix queries are performed per variant using thread-local `pysam.TabixFile` handles when pysam is available, falling back to subprocess `tabix` otherwise.

**Two-pass strategy:** Pass 1 scans the input file to collect all unique HGVS strings. Pass 2 performs tabix lookups (optionally in parallel via `--max-workers`), then streams the annotated output.

```bash
src/scripts/run_annotate_spliceai.sh input.tsv output.tsv \
  --mode precomputed \
  --precomputed-snv-vcf spliceai_scores.masked.snv.hg38.vcf.gz \
  --precomputed-indel-vcf spliceai_scores.masked.indel.hg38.vcf.gz
```

#### Cache preparation (optional one-time step)

The source VCFs are copied into the cache directory and indexed with `tabix -p vcf` on first use. Run `--prepare-cache-only` before annotation to do this separately:

```bash
src/scripts/run_annotate_spliceai.sh /dev/null /dev/null \
  --mode precomputed \
  --precomputed-snv-vcf spliceai_scores.masked.snv.hg38.vcf.gz \
  --precomputed-indel-vcf spliceai_scores.masked.indel.hg38.vcf.gz \
  --prepare-cache-only
```

To force re-copy and re-index of already-cached files, add `--refresh-cache`.

#### Multiple VCF files

The SNV and indel VCFs are searched in order; the first file that returns a match wins. Additional VCFs can be added with `--precomputed-vcf PATH` (can be specified multiple times). Files are searched in the order: SNV VCF → indel VCF → any `--precomputed-vcf` entries.

#### HGVS parsing and indel normalization

Variants are looked up by tabix region query using coordinates parsed from `NC_XXXXXX.XX:g.<expr>` HGVS notation. Supported variant classes:

| HGVS expression | Example | FASTA required? |
|---|---|---|
| SNV (`pos REF>ALT`) | `g.12345A>T` | No |
| Single-base deletion (`pos del`) | `g.12345del` | Yes — for anchor base |
| Range deletion (`start_end del`) | `g.12345_12347del` | Yes — for anchor base |
| Insertion (`pos_pos ins SEQ`) | `g.12345_12346insACG` | Yes — for ref anchor |
| Deletion-insertion (`pos delins SEQ`) | `g.12345delinsAC` | Yes — for ref allele |
| Duplication (`pos dup`, `start_end dup`) | `g.12345dup` | Yes — for anchor base |

For SNVs and tabix-only queries (no `--genome`), indel lookups fall back to positional matching by indel length/type — ref/alt bases are verified where possible but the anchor base is not checked. Supplying `--genome` enables exact VCF-normalised coordinate computation for all variant classes.

Both GRCh37 and GRCh38 NC accession mappings are built in. Set `--annotation grch37` to activate the GRCh37 mapping (also affects the `-A` argument passed to the SpliceAI CLI in compute mode).

### Compute mode

Runs the `spliceai` CLI subprocess to score variants on the fly using the deep-learning model. Much slower than precomputed mode and requires a local SpliceAI and TensorFlow installation. Suitable for variants not covered by the Illumina pre-scored files (e.g. uncommon or novel indels).

```bash
src/scripts/run_annotate_spliceai.sh input.tsv output.tsv \
  --mode compute \
  --genome /path/to/hg38.fa.gz \
  --annotation grch38
```

The script builds a temporary VCF from the unique variants in the input, invokes `spliceai -I ... -O ... -R ... -A ...`, parses the output VCF, and maps scores back to input rows. Temporary files are cleaned up after each run.

`--prepare-cache-only` is not supported in compute mode. `--genome` is required.

---

## CLI options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--mode {precomputed,compute}` | — | `precomputed` | Execution mode |
| `--annotation {grch38,grch37}` | — | `grch38` | Genome build for NC accession mapping and SpliceAI `-A` flag |
| `--hgvs-g-col COL` | — | `mapped_hgvs_g` | Input column with genomic HGVS values |
| `--cache-dir DIR` | `SPLICEAI_CACHE_DIR` | `/tmp/spliceai_cache` | Cache directory for VCFs and tabix indices |
| `--precomputed-snv-vcf PATH` | `SPLICEAI_PRECOMPUTED_SNV_VCF` | — | Precomputed SpliceAI SNV VCF (`.vcf.gz`) |
| `--precomputed-indel-vcf PATH` | `SPLICEAI_PRECOMPUTED_INDEL_VCF` | — | Precomputed SpliceAI indel VCF (`.vcf.gz`) |
| `--precomputed-vcf PATH` | — | — | Additional precomputed VCFs; can be specified multiple times |
| `--genome PATH` | — | — | Reference FASTA; required for compute mode; enables exact indel normalization in precomputed mode |
| `--prepare-cache-only` | — | off | Copy and index precomputed VCFs; do not annotate (precomputed mode only) |
| `--refresh-cache` | — | off | Re-copy source VCFs and rebuild tabix indices (precomputed mode only) |
| `--max-workers N` | — | `1` | Parallel tabix lookup threads (precomputed mode only) |
| `--skip N` | — | `0` | Skip first N data rows |
| `--limit N` | — | no limit | Stop after N rows |
| `--progress-every N` | — | `1000` | Log progress every N unique variants (0 disables) |
| `--log-level` | — | `INFO` | Logging verbosity |
| `--csv-field-size-limit BYTES` | — | system default | Increase for large HGVS fields |
| `--delimiter` | — | `\t` | Input/output delimiter |

---

## Dependencies

### Precomputed mode

- **`tabix`** — must be on `PATH` (from htslib). Used for indexing and as a fallback lookup mechanism.
- **`pysam`** (optional but recommended) — enables faster thread-local tabix handles, avoiding subprocess overhead per lookup. If pysam is not installed, every tabix query spawns a subprocess.
- **`pyfaidx`** (optional) — required for exact indel VCF normalisation (`--genome`). Without it, indel lookups use positional matching only.

### Compute mode

- **`spliceai`** Python package and **TensorFlow** — install separately; not included in the default project dependencies.
- **Reference FASTA** — must be a bgzipped, faidx-indexed `.fa.gz`.

---

## Troubleshooting

**Many empty `spliceai.*` columns for indels**

The Illumina pre-scored SNV and indel VCFs have coverage limitations. Some indel classes (e.g. insertions above a certain length, duplications) may not appear in the pre-scored files at all. Use `--mode compute` or inspect the Illumina SpliceAI release notes for coverage details.

**`tabix failed` error on first run**

The source VCF was not bgzipped, or `tabix` is not on `PATH`. Confirm the file ends in `.vcf.gz` (block-gzip compressed) and that `tabix --version` works. Use `bgzip` to compress plain VCF files before use.

**Slow precomputed lookups**

Increase `--max-workers` to parallelize tabix queries. Ensure `pysam` is installed to use thread-local handles instead of per-query subprocess invocations. For very large inputs, confirm the VCF files are in the Docker volume (not the overlay filesystem) to avoid I/O bottlenecks.

**`pyfaidx not installed` warning with indels**

Install `pyfaidx` and pass `--genome` to enable exact VCF-normalised indel coordinate computation. Without it, indel tabix queries use positional matching (chromosome + position) and verify length/type but not the exact ref/alt bases from the FASTA anchor.

**Compute mode: `spliceai failed`**

Ensure SpliceAI and TensorFlow are installed in the Python environment, and that the reference FASTA is bgzipped and faidx-indexed. Check that `--annotation` matches the genome build of your FASTA (`grch38` for hg38, `grch37` for hg19).
