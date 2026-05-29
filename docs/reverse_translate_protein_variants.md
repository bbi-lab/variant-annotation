# reverse_translate_protein_variants

`src/reverse_translate_protein_variants.py` — Step 2 of the variant-annotation pipeline.

Converts protein-only mapped variants into DNA (c./g.) HGVS candidates by exhaustively enumerating every codon substitution that produces the observed amino acid change. This is necessary because a single protein change can arise from two or three synonymous DNA codons, and downstream annotation steps (ClinVar, gnomAD, VEP) operate on DNA variants.

---

## What it does

For each input row where `mapped_hgvs_p` is non-empty and both `mapped_hgvs_c` and `mapped_hgvs_g` are blank, the script:

1. Resolves a transcript accession to use for backtranslation (see [Transcript resolution](#transcript-resolution)).
2. Calls the external `reverse-translate-variants` CLI (from the `variant-translation` package) in batches, passing `(transcript, hgvs_p)` pairs.
3. Writes the returned pipe-delimited `hgvs_c` and `hgvs_g` candidates back to the same `mapped_hgvs_c` / `mapped_hgvs_g` columns.
4. Populates derived position/allele columns (`_start`, `_stop`, `_ref`, `_alt`, `touches_intronic_region`, `spans_intron`) for every candidate.
5. Sets `assayed_variant_level` to `"protein"` for these rows.

When `--wt-codon-mode` is enabled, rows with synonymous protein HGVS changes (for example `p.Met1Met` or `p.=`) can also receive an additional WT codon candidate in `mapped_hgvs_c` and `mapped_hgvs_g`:

- `unambiguous`: only for amino acids with exactly one codon (Met = `ATG`, Trp = `TGG`)
- `all`: for all synonymous rows, using UTA to look up the transcript's actual codon for multi-codon amino acids
- `none` (default): no extra WT codon candidate is appended

Any appended WT candidate is deduplicated against candidates already returned by the reverse-translation CLI.

Rows that already have `mapped_hgvs_c` or `mapped_hgvs_g` (DNA variants from step 1) pass through unchanged; their derived columns are also populated, and `assayed_variant_level` is set to `"dna"`.

---

## Transcript resolution

For each protein-only row the script determines the transcript accession in this order:

1. **Accession prefix in `mapped_hgvs_p`** — if the protein HGVS string already has a transcript-like prefix (e.g. `NM_007294.3:p.Glu23fs`), that accession is used directly.
2. **RefSeq protein accession in `mapped_hgvs_p`** — if the prefix is an NP_/XP_/YP_/WP_ protein accession, the script queries UTA (`uta_20241220.associated_accessions`) to find the corresponding transcript. When multiple transcripts are found, preference is given to NM_ > XM_ > ENST, then highest version number.
3. **Transcript fallback columns** — any columns listed in `--transcript-fallback-column` are scanned for transcript-prefixed HGVS values (e.g. `raw_hgvs_nt`), and the accession from the first match is used.

If no transcript accession can be resolved, the row is skipped with `reverse_translation_error = "Unable to resolve transcript accession for protein reverse translation"`.

---

## Input and output columns

### Required input columns

| Column | Description |
|---|---|
| `mapped_hgvs_p` | Protein-level HGVS string from `map_variants` (e.g. `NP_000268.1:p.Ala406Thr`) |
| `mapped_hgvs_c` | Must be blank for a row to be reverse-translated |
| `mapped_hgvs_g` | Must be blank for a row to be reverse-translated |

### Output columns (added or updated)

| Column | Description |
|---|---|
| `mapped_hgvs_c` | Pipe-delimited transcript (c.) HGVS candidates, e.g. `NM_000277.3:c.1216G>A\|NM_000277.3:c.1217C>A` |
| `mapped_hgvs_g` | Pipe-delimited genomic (g.) HGVS candidates aligned to the same positions |
| `assayed_variant_level` | `"protein"` for protein-origin rows; `"dna"` for DNA-origin rows |
| `reverse_translation_error` | Error message if reverse translation failed; blank on success |
| `reverse_translation_warnings` | Any warnings from HGVS parsing (blank on success) |
| `mapped_hgvs_c_transcript` | Pipe-delimited transcript accessions extracted from `mapped_hgvs_c` candidates |
| `mapped_hgvs_c_start` | Pipe-delimited c. start positions |
| `mapped_hgvs_c_stop` | Pipe-delimited c. stop positions |
| `mapped_hgvs_c_ref` | Pipe-delimited c. reference alleles |
| `mapped_hgvs_c_alt` | Pipe-delimited c. alternate alleles |
| `mapped_hgvs_g_chromosome` | Pipe-delimited chromosome identifiers extracted from g. candidates |
| `mapped_hgvs_g_start` | Pipe-delimited g. start positions |
| `mapped_hgvs_g_stop` | Pipe-delimited g. stop positions |
| `mapped_hgvs_g_ref` | Pipe-delimited g. reference alleles |
| `mapped_hgvs_g_alt` | Pipe-delimited g. alternate alleles |
| `touches_intronic_region` | Pipe-delimited boolean flags (`"true"`/`"false"`) aligned to `mapped_hgvs_c` candidates |
| `spans_intron` | Pipe-delimited boolean flags (`"true"`/`"false"`) aligned to `mapped_hgvs_c` candidates |

All pipe-delimited output columns are position-aligned: index 0 of `mapped_hgvs_c`, `mapped_hgvs_g`, `mapped_hgvs_c_start`, … all refer to the same candidate. An empty slot (e.g. `||` in a 3-candidate row) means that candidate could not be resolved.

---

## Usage

```bash
src/scripts/run_reverse_translate_protein_variants.sh input.tsv output.tsv [options]
```

Or directly:

```bash
docker compose run --rm reverse-translate-protein-variants \
    /work/input.tsv /work/output.tsv [options]
```

---

## CLI options

### Candidate generation

| Option | Default | Description |
|---|---|---|
| `--include-indels` | off | Include small insertion/deletion candidates in addition to substitutions and delins |
| `--max-indel-size N` | `3` | Maximum indel size (nt) when `--include-indels` is set |
| `--wt-codon-mode {none,unambiguous,all}` | `none` | Append WT codon `c.`/`g.` candidates for synonymous protein variants (`ref AA == alt AA`). `unambiguous` adds only Met/Trp; `all` adds all and queries UTA for multi-codon residues. Requires `--include-indels`. |
| `--no-strict-ref-aa` | off | Disable reference amino-acid validation against the resolved transcript sequence |
| `--use-inv-notation` | off | Emit inversions using `inv` HGVS notation instead of `delins` |
| `--substitutions-and-delins-only` | off | For premature-stop changes, suppress length-changing insertion/deletion candidates |
| `--assembly ASSEMBLY` | `GRCh38` | Reference genome assembly |

### Column names

All default column names match the defaults from `map_variants`. Override only if your input uses different names.

| Option | Default |
|---|---|
| `--mapped-hgvs-g` | `mapped_hgvs_g` |
| `--mapped-hgvs-c` | `mapped_hgvs_c` |
| `--mapped-hgvs-p` | `mapped_hgvs_p` |
| `--reverse-translation-error` | `reverse_translation_error` |
| `--reverse-translation-warnings` | `reverse_translation_warnings` |
| `--assayed-variant-level` | `assayed_variant_level` |
| `--mapped-hgvs-g-chromosome` | `mapped_hgvs_g_chromosome` |
| `--mapped-hgvs-c-transcript` | `mapped_hgvs_c_transcript` |
| `--mapped-hgvs-g-start` | `mapped_hgvs_g_start` |
| `--mapped-hgvs-g-stop` | `mapped_hgvs_g_stop` |
| `--mapped-hgvs-g-ref` | `mapped_hgvs_g_ref` |
| `--mapped-hgvs-g-alt` | `mapped_hgvs_g_alt` |
| `--mapped-hgvs-c-start` | `mapped_hgvs_c_start` |
| `--mapped-hgvs-c-stop` | `mapped_hgvs_c_stop` |
| `--mapped-hgvs-c-ref` | `mapped_hgvs_c_ref` |
| `--mapped-hgvs-c-alt` | `mapped_hgvs_c_alt` |
| `--touches-intronic-region` | `touches_intronic_region` |
| `--spans-intron` | `spans_intron` |

### Transcript lookup

| Option | Description |
|---|---|
| `--transcript-fallback-column COL` | Additional column to scan for transcript accessions when none can be found in `mapped_hgvs_p`. Repeatable. Useful when `raw_hgvs_nt` carries the original transcript (e.g. `NM_000277.3:c.1218G>A`) for rows that were later reclassified as protein-only. |
| `--no-resolve-missing-ref-alleles` | Skip UTA ref-allele resolution for del/delins-like HGVS strings in derived columns. |

### Run control

| Option | Default | Description |
|---|---|---|
| `--skip N` | `0` | Skip the first N input rows |
| `--limit N` | `0` (no limit) | Stop after processing N rows |
| `--log-level` | `INFO` | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`) |
| `--csv-field-size-limit BYTES` | system default | Increase if large HGVS fields cause CSV parse errors |

---

## Dependencies

### `reverse-translate-variants` CLI

This script delegates all backtranslation to the `reverse-translate-variants` command-line tool from the `variant-translation` package. The tool must be on `PATH`; the Docker image installs it automatically. The CLI is invoked with `--one-row-per-input` and `--join-delimiter |` so that every input row produces exactly one output row with pipe-delimited candidates.

### UTA database

Transcript accession lookup (NP_ → NM_) queries the `uta_20241220.associated_accessions` table via PostgreSQL. Set `UTA_DB_URL` to your UTA connection string, e.g.:

```
postgresql://uta_admin:uta@uta:5432/uta
```

The Docker Compose service starts UTA automatically. If UTA is unreachable the script raises `RuntimeError` immediately.

When `--wt-codon-mode all` is used, UTA is also queried for codon sequence lookup at the relevant CDS position for synonymous protein variants.

---

## Batching strategy

Protein-only rows that appear consecutively in the input are accumulated into a contiguous block and flushed to the `reverse-translate-variants` CLI as a single batch call. When a non-protein row interrupts the block, the accumulated block is flushed first, then the non-protein row is written immediately. This minimises subprocess overhead while preserving strict output row order.

---

## Troubleshooting

**No candidates returned for a row**

Check that:
- `mapped_hgvs_p` uses three-letter amino acid codes (e.g. `p.Ala406Thr`, not `p.A406T`). `map_variants` converts one-letter codes automatically, but externally supplied HGVS may not be normalised.
- A valid transcript accession is resolvable (check logs for "Unable to resolve transcript accession").
- The amino acid change is valid for the transcript (`--no-strict-ref-aa` relaxes this check).
- For synonymous (silent) variants (`p.=` or `p.Ala406=`), default behavior (`--wt-codon-mode none`) may return no candidate because the protein is unchanged. Enable `--wt-codon-mode unambiguous` or `--wt-codon-mode all` (with `--include-indels`) to emit WT codon delins candidates.

**Many rows have `reverse_translation_error`**

- Confirm UTA is running: `docker compose logs uta`
- If the error is from `reverse-translate-variants` itself (exit code non-zero), the whole block of protein rows in that batch will fail together. Re-run with `--skip N` to resume after the failing block.

**Indel candidates are missing**

By default, only substitution and delins candidates are generated. Pass `--include-indels --max-indel-size 3` to include small insertions and deletions.
