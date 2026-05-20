# annotate_vep

`src/annotate_vep.py` — Step 9 of the variant-annotation pipeline (optional).

Annotates each variant row with a mutational consequence term from the Ensembl VEP REST API. The core design goal is to produce the most transcript-specific consequence possible: for transcript HGVS inputs (`NM_`, `NR_`, `ENST`) the script selects the consequence for the specific referenced transcript rather than the global worst across all overlapping transcripts.

---

## Output columns

| Column | Description |
|---|---|
| `vep.mutational_consequence` | Sequence Ontology consequence term (e.g. `missense_variant`, `synonymous_variant`). Empty when candidates disagree or VEP returned nothing. |
| `vep.consequence_source` | `transcript` — consequence taken from the matched `transcript_consequences` entry; `most_severe` — global fallback |
| `vep.access_date` | ISO date the Ensembl API was queried (or the date from a cache hit) |
| `vep.error` | Non-empty when resolved candidates disagree; contains a pipe-delimited consequence list aligned to input candidates |

The namespace prefix defaults to `vep` and can be changed with `--vep-namespace`.

---

## Transcript selection

Ensembl VEP annotates every transcript that overlaps the genomic position, not only the one referenced in the input HGVS. The top-level `most_severe_consequence` field is therefore the worst outcome across the entire locus.

This script selects the consequence for the specific input transcript when the input HGVS makes that possible:

| Input HGVS type | API flag sent | How the consequence is selected | `vep.consequence_source` |
|---|---|---|---|
| RefSeq transcript (`NM_`, `NR_`) | `refseq=1` | Matched by versionless `transcript_id` in `transcript_consequences` | `transcript` |
| Ensembl transcript (`ENST`) | _(default)_ | Matched by versionless `transcript_id` in `transcript_consequences` | `transcript` |
| Genomic (`NC_`, `chr:g.`) | _(default)_ | `most_severe_consequence` used — no specific transcript can be identified | `most_severe` |
| Protein (`NP_`) / Recoder fallback | _(default)_ | `most_severe_consequence` used | `most_severe` |

When no matching transcript entry is found (e.g. the transcript version is absent from Ensembl's RefSeq set), the script falls back to `most_severe_consequence` and records `most_severe` in `vep.consequence_source`.

**Column priority:** `mapped_hgvs_c` should appear before `mapped_hgvs_g` in `--hgvs-cols` (the default) because transcript HGVS strings yield transcript-specific consequences while genomic HGVS strings always yield `most_severe`.

---

## API call strategy

For each batch of input rows:

1. **VEP POST** (`/vep/human/hgvs`) — RefSeq HGVS strings (`NM_`/`NR_`) are sent with `refseq=1` in separate batches from Ensembl/genomic HGVS strings. Batches are dispatched concurrently via a `ThreadPoolExecutor`.
2. **Variant Recoder POST** (`/variant_recoder/human`) — For HGVS strings that returned no result from the first VEP call, the Recoder API resolves them to genomic HGVS equivalents (`NC_...`).
3. **Second VEP POST** — The recoded genomic HGVS strings are queried in a second concurrent VEP pass. When multiple recoded equivalents exist, the most severe consequence across them is chosen.
4. **API errors** — When an entire HTTP request fails (timeout, non-200), the affected HGVS strings are marked `source=api_error` and are **not cached**, so they will be retried on the next run.

### Consequence severity ordering

When multiple consequence terms are present (e.g. from `transcript_consequences` entries), the most severe is selected using the canonical VEP severity ordering: `transcript_ablation` > `splice_acceptor_variant` > `splice_donor_variant` > `stop_gained` > `frameshift_variant` > … > `intergenic_variant`. This list mirrors MaveDB worker logic.

---

## Multi-candidate rows

For rows with pipe-delimited HGVS candidates (from step 2 reverse translation):

- Each candidate is resolved independently.
- If all resolved candidates agree (or only some candidates failed to resolve), `vep.mutational_consequence` contains the single agreed consequence.
- If resolved candidates return **different** consequences, `vep.mutational_consequence` is left empty and `vep.error` records the discrepancy as a pipe-delimited list aligned to the candidate positions.

---

## Redis caching

VEP API responses are cached in Redis as `(consequence, source)` pairs per HGVS string. Misses (VEP returned nothing) are stored under a sentinel so repeated no-hit queries don't re-query the API. API errors are not cached and will be retried.

| Variable | Default | Description |
|---|---|---|
| `VEP_CACHE_ENABLED` | `true` | Set to `0` / `false` to disable caching entirely |
| `VEP_CACHE_REDIS_URL` | `redis://redis:6379/0` | Redis connection URL (also falls back to `REDIS_URL`) |
| `VEP_CACHE_PREFIX` | `vep:v1` | Key namespace; bump the version suffix to invalidate after a significant Ensembl release |
| `VEP_CACHE_TTL_SECONDS` | `86400` | TTL for hits (seconds) |
| `VEP_CACHE_MISS_TTL_SECONDS` | `86400` | TTL for misses (seconds) |

---

## CLI options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--vep-namespace NS` | — | `vep` | Output column name prefix |
| `--hgvs-cols COLS` | — | `mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p` | Comma-separated columns tried in priority order; first non-blank value used per row |
| `--vep-api-url URL` | `ENSEMBL_API_URL` | `https://rest.ensembl.org` | Ensembl REST API base URL |
| `--vep-batch-size N` | `VEP_BATCH_SIZE` | `500` | HGVS values per VEP POST request |
| `--vep-workers N` | `VEP_WORKERS` | `8` | Concurrent VEP/Recoder batch requests |
| `--row-batch-size N` | `VEP_ROW_BATCH_SIZE` | `1000` | Input rows per lookup/write batch |
| `--vep-timeout-seconds N` | `VEP_TIMEOUT_SECONDS` | `60` | HTTP timeout per API request |
| `--keep-existing` | — | off | Skip rows already annotated (non-empty `vep.mutational_consequence`); only annotate blank rows |
| `--vep-cache-file FILE` | — | — | Path to a pre-computed TSV fallback cache (columns: `hgvs`, `vep.mutational_consequence`, `vep.access_date`, `vep.error`); rows with non-empty `vep.error` are ignored |
| `--skip N` | — | `0` | Skip first N data rows |
| `--limit N` | — | no limit | Stop after N rows |
| `--log-level` | — | `INFO` | Logging verbosity |
| `--csv-field-size-limit BYTES` | — | system default | Increase for large HGVS fields |
| `--delimiter` | — | `\t` | Input/output delimiter |

---

## Usage

```bash
src/scripts/run_annotate_vep.sh input.tsv output.tsv
```

Resume an interrupted run (re-use existing annotations, only query blank rows):
```bash
src/scripts/run_annotate_vep.sh input.tsv output.tsv --keep-existing
```

Use a local Ensembl mirror:
```bash
src/scripts/run_annotate_vep.sh input.tsv output.tsv \
  --vep-api-url http://grch38.ensembl.org
```

---

## Dependencies

- **`requests`** — included in project dependencies.
- **Ensembl REST API** — public endpoint `https://rest.ensembl.org`; no authentication required. Rate limits apply; the default `--vep-workers 8` is designed to stay within the public API's burst tolerance.
- **Redis** (optional) — greatly reduces repeat API calls across runs. Fails gracefully when unavailable.

---

## Troubleshooting

**Many rows with empty `vep.mutational_consequence` and non-empty `vep.error`**

Pipe-delimited candidates disagree on consequence. This is expected for protein-derived rows where one synonymous codon causes a `synonymous_variant` and another a `splice_region_variant`. Check the discrepant consequences in `vep.error` and decide whether to post-filter.

**`api_error` entries on every run**

- Ensembl REST API may be under maintenance or rate-limiting. Check `https://rest.ensembl.org` directly.
- Reduce `--vep-workers` or increase `--vep-timeout-seconds`.
- `api_error` entries are never cached, so they are always retried. Use `--keep-existing` to skip already-annotated rows and only retry the blanks.

**Consequence is `most_severe` for all rows despite transcript HGVS input**

- Confirm `mapped_hgvs_c` is listed first in `--hgvs-cols`.
- The transcript version in the input HGVS may not be present in Ensembl's current RefSeq set, triggering the `most_severe` fallback. This is common for older transcript versions. The `vep.consequence_source` column records when this happens.

**`VEP_CACHE_PREFIX` after an Ensembl release**

Bump the prefix suffix (e.g. `vep:v1` → `vep:v2`) to force all cached entries to be re-queried against the new Ensembl release. Old keys expire naturally via TTL.
