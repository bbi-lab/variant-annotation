# liftover

`src/liftover.py` — lift genomic positions between genome assemblies (e.g. hg19 -> hg38) using a UCSC chain file.

Not part of the annotation pipeline; this is a standalone utility for preparing or cross-checking input data that lives on a different assembly than the pipeline expects.

Runs locally using the project's Python virtual environment (not Docker). Requires the `pyliftover` package, an optional extra not installed by default:

```bash
pip install variant-annotation[liftover]
# or: pip install pyliftover
```

Chain files are downloaded automatically from UCSC on first use and cached in `~/.pyliftover/`; pass `--chain-file` to supply a pre-downloaded file instead.

---

## Modes

Two modes are available and may be combined on the same input file.

### HGVS mode (`--input-col`)

Lifts genomic HGVS `g.` strings such as `NC_000001.10:g.925952C>T`. Supported variant types: substitutions, deletions, insertions, delins, duplications, and inversions.

Input strings may use RefSeq accession identifiers (`NC_000001.10`) or plain chromosome names (`chr1`); the output format mirrors the input. Pipe-delimited multi-candidate values are lifted independently and results are pipe-aligned to the input.

A position that maps to the opposite strand, or a range whose start and end map to different chromosomes, is reported as an error rather than silently lifted — these require manual review.

### Coordinate mode (`--chrom-col`)

Lifts raw chromosome + position columns such as those produced by a BED or region file. Positions are 1-based inclusive. `--end-col` is optional; omit it for single-position rows.

---

## Output

For each mode, new columns are appended to the input (originals are preserved):

| Mode | New columns |
|---|---|
| HGVS | `<input-col>_<to-assembly>` (or `--output-col`), `liftover_error` (or `--error-col`) |
| Coordinate | `<chrom-col>_<to-assembly>`, `<start-col>_<to-assembly>`, `<end-col>_<to-assembly>` (if `--end-col` given), `coord_liftover_error` (or `--coord-error-col`) |

A row with a non-blank error column failed to lift; the corresponding lifted column(s) are left blank for that row.

---

## Usage

```bash
# HGVS mode
src/scripts/run_liftover.sh input.tsv output.tsv \
  --from-assembly hg19 --to-assembly hg38 \
  --input-col hgvs_g

# Coordinate mode
src/scripts/run_liftover.sh regions.tsv lifted.tsv \
  --from-assembly hg19 --to-assembly hg38 \
  --chrom-col chrom --start-col start --end-col end

# Both modes on the same file
src/scripts/run_liftover.sh variants.tsv lifted.tsv \
  --from-assembly hg19 --to-assembly hg38 \
  --input-col hgvs_g --output-col hgvs_g_hg38 \
  --chrom-col chrom --start-col start --end-col end

# CSV, custom assembly names, pre-downloaded chain file
src/scripts/run_liftover.sh variants.csv lifted.csv \
  --delimiter , \
  --from-assembly GRCh37 --to-assembly GRCh38 \
  --chain-file /data/hg19ToHg38.over.chain.gz \
  --input-col mapped_hgvs_g --output-col hgvs_g_hg38
```

| Option | Description |
|---|---|
| `--from-assembly ASM` | Source assembly (`hg19`/`GRCh37` or `hg38`/`GRCh38`). Default: `hg19` |
| `--to-assembly ASM` | Target assembly. Default: `hg38` |
| `--chain-file PATH` | UCSC `.over.chain`/`.over.chain.gz` file; auto-downloaded if omitted |
| `--input-col COL` | HGVS mode: column containing HGVS g. strings to lift |
| `--output-col COL` | HGVS mode: output column name (default: `<input-col>_<to-assembly>`) |
| `--error-col COL` | HGVS mode: error column name (default: `liftover_error`) |
| `--chrom-col COL` | Coordinate mode: column containing chromosome names (`chr1` or `1`) |
| `--start-col COL` | Coordinate mode: column containing 1-based start positions (required with `--chrom-col`) |
| `--end-col COL` | Coordinate mode: column containing 1-based end positions (optional) |
| `--lifted-chrom-col COL` | Coordinate mode: output chromosome column (default: `<chrom-col>_<to-assembly>`) |
| `--lifted-start-col COL` | Coordinate mode: output start column (default: `<start-col>_<to-assembly>`) |
| `--lifted-end-col COL` | Coordinate mode: output end column (default: `<end-col>_<to-assembly>`) |
| `--coord-error-col COL` | Coordinate mode: error column name (default: `coord_liftover_error`) |
| `--skip N` | Skip first N data rows |
| `--limit N` | Process at most N data rows |
| `--delimiter DELIM` | Field delimiter (default: TAB; use `,` for CSV) |
| `--csv-field-size-limit BYTES` | Increase for large fields |
| `--log-level LEVEL` | `DEBUG`, `INFO` (default), `WARNING`, or `ERROR` |

At least one of `--input-col` or `--chrom-col` must be given.
