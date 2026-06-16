# annotate_vep

`src/annotate_vep.py` — Step 9 of the variant-annotation pipeline (optional).

Annotates each variant row with a mutational consequence term from the Ensembl VEP REST API. The core design goal is to produce the most transcript-specific consequence possible: for transcript HGVS inputs (`NM_`, `NR_`, `ENST`) the script selects the consequence for the specific referenced transcript rather than the global worst across all overlapping transcripts.

---

## Output columns

| Column | Description |
|---|---|
| `vep.mutational_consequences` | `^`-delimited list of all consequence terms when the result came from a matched transcript entry; single most-severe term otherwise. Pipe-delimited across candidates. Empty string for a candidate with an API error. |
| `vep.most_severe_mutational_consequence` | Single most-severe consequence term per candidate. Pipe-delimited across candidates. Empty string for a candidate with an API error. |
| `vep.consequence_source` | `transcript`, `most_severe`, or `no_change` per candidate; empty string when the candidate had an API error and no fallback applied. Pipe-delimited across candidates. |
| `vep.access_date` | ISO access date per candidate. Pipe-delimited across candidates, aligned to the input candidate positions (empty slot for empty candidates). |
| `vep.error` | Per-candidate API error message (`api_error` or `api_error:<sanitized detail>`) when the candidate's request failed; empty string otherwise. Pipe-delimited across candidates. When `no_change` fallback is confirmed, this field is cleared. |

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

### Special `no_change` handling for unchanged transcript `c.delins`

Some transcript HGVS strings represent no actual sequence change (for example, replacing a codon with the same codon) but can still trigger a VEP parse error. To avoid surfacing these as API failures, `annotate_vep` applies a special fallback:

- Pattern: transcript-level `c.start_stopdelinsALT` (for example `NM_000133.4:c.1_3delinsATG`).
- Verification: query UTA for the transcript CDS reference sequence at `start..stop` and compare to `ALT`.
- If identical (`ref == alt`):
  - `vep.mutational_consequences = no_change`
  - `vep.most_severe_mutational_consequence = no_change`
  - `vep.consequence_source = no_change`
  - `vep.error = ""`

This fallback is used when VEP returns an API error, returns no consequence, or returns no cache entry for that candidate.

Prerequisite: `UTA_DB_URL` must be set and reachable from the runtime container/environment.

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

For rows with pipe-delimited HGVS candidates (from step 2 reverse translation), each candidate is resolved independently. All output columns are pipe-delimited with one value per candidate position in the same order as the input.

- `vep.mutational_consequences`: `^`-delimited consequence terms for transcript HGVS candidates (`source == "transcript"`), or the single most-severe term for genomic/protein inputs. Empty string for a candidate with an API error.
- `vep.most_severe_mutational_consequence`: single most-severe term per candidate; empty string on API error.
- `vep.access_date`: ISO access date per candidate; pipe-delimited and position-aligned to candidates.
- `vep.consequence_source`: `transcript`, `most_severe`, or `no_change` per candidate; empty string on API error when no fallback applies.
- `vep.error`: `api_error` or `api_error:<sanitized message>` for candidates whose API request failed; empty string otherwise. Cleared when `no_change` fallback applies.

---

## File-based consequence cache

A pre-built TSV of VEP consequences can be loaded at startup with `--vep-file-cache FILE`. Entries from the file populate the consequence cache before Redis is consulted, so file-cache hits bypass both Redis and the VEP API entirely. This is useful for offline use, reproducibility, or pre-seeding a fresh run from a prior annotated output.

### File format

Tab-separated, with a header row. Column names use the same namespace prefix as the current run (`--vep-namespace`, default `vep`):

| Column | Description |
|---|---|
| `hgvs` | Input HGVS string (cache key) |
| `{prefix}.most_severe_mutational_consequence` | Most-severe consequence term; blank for a confirmed miss |
| `{prefix}.mutational_consequences` | `^`-delimited list of all consequence terms; may be blank |
| `{prefix}.consequence_source` | `transcript`, `most_severe`, `no_change`, or an `api_error:...` value |
| `{prefix}.access_date` | ISO date string from the original annotation run; propagated to output |
| `{prefix}.error` | Original error message; re-encoded into `consequence_source` on load |

A row with a blank `most_severe_mutational_consequence` is treated as a confirmed miss and will not be re-queried. A non-blank `{prefix}.error` value is used to reconstruct the `api_error:...` source token so `annotate_row` emits it correctly.

The file format matches the output produced by a prior `annotate_vep` run, so an existing annotated TSV can be cut, de-duplicated by `hgvs`, and fed directly back in.

### Access dates from the file cache

When a candidate is resolved from the file cache, the `{prefix}.access_date` value from the file is written to the output rather than today's date. This preserves provenance when replaying annotations from an older run.

---

## Redis caching

VEP API responses are cached in Redis as `(most_severe, all_consequences, source)` triples per HGVS string. Misses (VEP returned nothing) are stored under a sentinel so repeated no-hit queries don't re-query the API. API errors are not cached and will be retried. Cache entries from prior versions that lack the `all_consequences` field are silently discarded on read and re-queried.

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
| `--vep-file-cache FILE` | — | — | TSV of pre-computed VEP consequences; entries pre-populate the consequence cache and take precedence over Redis and the VEP API |
| `--keep-existing` | — | off | Skip rows already annotated (non-empty `vep.most_severe_mutational_consequence`); only annotate blank rows |
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

**`api_error` entries on every run**

- The `vep.error` column now includes the HTTP status and response excerpt (e.g. `api_error:VEP HTTP 503: ...`) to identify the specific failure.
- Ensembl REST API may be under maintenance or rate-limiting. Check `https://rest.ensembl.org` directly.
- Reduce `--vep-workers` or increase `--vep-timeout-seconds`.
- `api_error` entries are never cached, so they are always retried. Use `--keep-existing` to skip already-annotated rows and only retry the blanks.

**Consequence is `most_severe` for all rows despite transcript HGVS input**

- Confirm `mapped_hgvs_c` is listed first in `--hgvs-cols`.
- The transcript version in the input HGVS may not be present in Ensembl's current RefSeq set, triggering the `most_severe` fallback. This is common for older transcript versions. The `vep.consequence_source` column records when this happens.

**`VEP_CACHE_PREFIX` after an Ensembl release**

Bump the prefix suffix (e.g. `vep:v1` → `vep:v2`) to force all cached entries to be re-queried against the new Ensembl release. Old keys expire naturally via TTL.
