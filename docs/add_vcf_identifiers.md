# add_vcf_identifiers

`src/add_vcf_identifiers.py` — Step 4 of the variant-annotation pipeline (optional).

Parses three mapped HGVS columns (`mapped_hgvs_g`, `mapped_hgvs_c`, `mapped_hgvs_p`) into VCF-style component fields: chromosome/transcript, start position, stop position, reference allele, and alternate allele. Deletions and insertions are written using the standard VCF left-anchor convention. Inversions are converted to the reverse-complement sequence. Boolean intronic-region flags are derived from the transcript HGVS coordinates.

When HGVS columns are pipe-delimited (as produced by `reverse_translate_protein_variants` for protein-origin rows), every candidate is parsed independently and all output columns are similarly pipe-delimited, preserving alignment with `dna_clingen_allele_id`.

---

## What it does

For each input row the script:

1. Splits each of the three HGVS input columns on `|` to get a list of candidates (one element for DNA rows; potentially many for reverse-translated protein rows).
2. Parses each candidate HGVS string into `(start, stop, ref, alt)` using HGVS edit notation rules.
3. Applies the VCF left-anchor convention: deletions shift the `start` position one base to the left and prepend the anchor base; insertions prepend the anchor base to both `ref` and `alt`.
4. Converts inversions to the reverse-complement sequence.
5. Writes pipe-delimited strings of the extracted values as new output columns.
6. Resolves missing `ref` alleles for `del`/`delins`/`dup`/`inv` variants where no sequence was supplied in the HGVS string, by querying UTA via the `hgvs` Python library.
7. Sets per-row boolean flags `touches_intronic_region` and `spans_intron` — both are `"true"` if **any** transcript candidate is intronic or spans an intron boundary.

---

## Relationship to reverse_translate_protein_variants

`reverse_translate_protein_variants` (step 2) imports `_parse_hgvs` and `_apply_vcf_anchor` from this module and writes the same derived columns for every row as part of its processing. If you run step 2, the derived position columns will already exist in the output.

Step 4 re-computes these columns. Two differences to be aware of:

- **`touches_intronic_region` / `spans_intron`:** step 2 writes these as *pipe-delimited* per-candidate flags. Step 4 collapses them to a single `"true"` / `"false"` per row (True if any candidate fires). Running step 4 after step 2 overwrites the per-candidate flags with the row-level flags.
- **`--resolve-missing-ref-alleles`:** step 4 always enables this from the CLI (default True). Step 2 also enables it by default. Results should be identical but the code paths are separate.

If you have no protein-only rows, step 2 was skipped, and step 4 is the only place where derived position columns are added.

---

## Input and output columns

### Required input columns

| Column | Description |
|---|---|
| `mapped_hgvs_g` | Genomic HGVS, e.g. `NC_000012.12:g.102917130G>A`; may be pipe-delimited |
| `mapped_hgvs_c` | Transcript HGVS, e.g. `NM_000277.3:c.1218G>A`; may be pipe-delimited |
| `mapped_hgvs_p` | Protein HGVS, e.g. `NP_000268.1:p.Ala406Thr` |

Blank or absent columns are silently skipped.

### Output columns (added or updated)

| Column | Genomic | Transcript | Protein | Notes |
|---|---|---|---|---|
| `mapped_hgvs_g_chromosome` | ✓ | — | — | Chromosome number/letter (1–22, X, Y, M) extracted from NC_ accession |
| `mapped_hgvs_c_transcript` | — | ✓ | — | Transcript accession extracted from the HGVS prefix |
| `mapped_hgvs_p_chromosome` | — | — | ✓ | Protein accession extracted from the HGVS prefix |
| `mapped_hgvs_*_start` | ✓ | ✓ | ✓ | Start position (1-based); left-shifted by 1 for deletions |
| `mapped_hgvs_*_stop` | ✓ | ✓ | ✓ | Stop position (1-based; equals start for SNVs) |
| `mapped_hgvs_*_ref` | ✓ | ✓ | ✓ | Reference allele; one-letter amino acid codes for protein |
| `mapped_hgvs_*_alt` | ✓ | ✓ | ✓ | Alternate allele; one-letter amino acid codes for protein |
| `touches_intronic_region` | — | ✓ | — | `"true"` if any c. candidate coordinate includes an intronic offset (`+`/`-`) |
| `spans_intron` | — | ✓ | — | `"true"` if any c. candidate spans both sides of an intron boundary |

All g. and c. output columns are pipe-delimited when the input is pipe-delimited; p. columns are never pipe-delimited (at most one protein HGVS per row).

---

## HGVS parsing rules

### Nucleotide edits (g., c., n.)

| HGVS edit | `ref` | `alt` |
|---|---|---|
| `A>T` (substitution) | `A` | `T` |
| `delAT` (deletion with sequence) | `AT` (pre-anchor) | `` (empty) |
| `del` (deletion, no sequence) | resolved from UTA | `` (empty) |
| `insGC` (insertion) | `` (empty, pre-anchor) | `GC` (pre-anchor) |
| `delinsGC` (delins) | deleted sequence | `GC` |
| `dupT` (duplication with sequence) | `T` | `TT` |
| `dup` (duplication, no sequence) | resolved from UTA | `{ref}{ref}` |
| `invATG` (inversion) | `ATG` | `CAT` (reverse complement) |
| `=` (synonymous) | `` | `` |

### VCF left-anchor convention

For deletions and insertions the VCF format requires the anchor base (the last unchanged base to the left). This script:

- **Deletions:** decrements `start` by 1 and prepends the anchor base (fetched from UTA) to both `ref` and `alt`.
- **Insertions:** prepends the anchor base (the base already at `start`) to both `ref` and `alt`.

When UTA is unavailable or the anchor cannot be determined, the values are written without anchoring.

### Genomic haplotypes

`g.[comp1;comp2;…]` expressions (multi-component haplotypes) are collapsed into a single allele covering the full span. Bases between component variants are filled from the reference sequence via UTA. If UTA is unavailable and there are gaps between components, the haplotype cannot be resolved and all fields are left blank.

### Protein edits (p.)

Protein HGVS is parsed to 1-letter amino acid codes:

| HGVS edit | `ref` | `alt` |
|---|---|---|
| `p.Ala406Thr` (missense) | `A` | `T` |
| `p.Ala406=` (synonymous) | `A` | `A` |
| `p.Ala406del` (deletion) | `A` | `` (empty) |
| `p.Ala406Ter` or `p.Ala406*` (nonsense) | `A` | `*` |
| `p.Ala406fs` (frameshift) | `A` | `fs` |
| `p.Ala406_Arg410del` (range deletion) | `A_R` | `` (empty) |
| `p.Ala406_Arg410delinsGly` | `A_R` | `G` |
| `p.Ala405_Ala406insThr` (insertion) | `A` (anchor) | `AT` (anchor + inserted) |

---

## Missing ref allele resolution

For `del`, `dup`, `delins`, and `inv` edits where the deleted/duplicated sequence is not encoded in the HGVS string (e.g. `NM_000277.3:c.1218del` rather than `c.1218delG`), the script uses the `hgvs` Python library to normalise the variant against UTA and extract the reference sequence. The resolved sequence is then used as `ref`.

This requires UTA. When UTA is unavailable the ref field is left blank for unresolved variants.

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `--mapped-hgvs-g COL` | `mapped_hgvs_g` | Genomic HGVS input column |
| `--mapped-hgvs-c COL` | `mapped_hgvs_c` | Transcript HGVS input column |
| `--mapped-hgvs-p COL` | `mapped_hgvs_p` | Protein HGVS input column |
| `--touches-intronic-region-column COL` | `touches_intronic_region` | Output boolean column name |
| `--spans-intron-column COL` | `spans_intron` | Output boolean column name |
| `--resolve-missing-ref-alleles` | on (CLI default) | Resolve missing ref alleles via UTA/HGVS normalization |
| `--max-workers N` | `8` | Number of concurrent worker threads |
| `--skip N` | `0` | Skip the first N data rows |
| `--limit N` | no limit | Stop after N rows |
| `--log-level` | `INFO` | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`) |
| `--csv-field-size-limit BYTES` | system default | Increase for large HGVS fields |

---

## Dependencies

### UTA database

Used for:
- Resolving missing ref alleles (`hgvs` library connects to UTA)
- Fetching anchor bases for indel VCF anchoring
- Translating `c.` positions to transcript sequence positions for anchor lookup

Set `UTA_DB_URL` to your UTA connection string, e.g. `postgresql://uta_admin:uta@uta:5432/uta`. When UTA is unavailable, the script continues — fields that require UTA are left blank rather than failing.

---

## Concurrency

Rows are processed concurrently using a `ThreadPoolExecutor` (default 8 workers). Results are written in input order. The in-flight cap is `max_workers × 4` to bound memory usage. HGVS and UTA lookups use `lru_cache` to avoid redundant queries across workers.

---

## Troubleshooting

**Many blank `_ref` / `_alt` columns for deletions**

UTA is likely unavailable. Confirm it is running (`docker compose logs uta`) and that `UTA_DB_URL` is set correctly.

**`_start` / `_stop` / `_ref` / `_alt` all blank for a row**

The HGVS string could not be parsed. Common causes:
- Non-standard edit notation (e.g. `?` positions, predicted protein changes in parentheses like `p.(Ala406Thr)` — these are parsed with the parentheses stripped)
- Accession-only strings with no posedit (e.g. `NM_000277.3:`)
- Malformed pipe-delimited candidates from an earlier step

**`touches_intronic_region` is always `"false"`**

This column reflects the transcript (c.) HGVS column. If the input only has genomic (g.) HGVS strings and no c. strings, the flag will always be false.
