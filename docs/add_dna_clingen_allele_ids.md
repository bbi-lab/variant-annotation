# add_dna_clingen_allele_ids

`src/add_dna_clingen_allele_ids.py` — Step 3 of the variant-annotation pipeline.

Resolves a DNA-level ClinGen Allele Registry (CAR) identifier for each variant and writes it to a new `dna_clingen_allele_id` column. For protein-origin rows that were reverse-translated in step 2, each pipe-delimited DNA candidate gets its own ID, producing a pipe-delimited output aligned to the candidates in `mapped_hgvs_c` / `mapped_hgvs_g`. For DNA-origin rows the existing `clingen_allele_id` from step 1 is reused when it is already a CA-prefixed identifier.

`dna_clingen_allele_id` is the primary key used by all downstream annotation steps (ClinVar, gnomAD, VEP).

---

## What it does

For each row:

1. Determines whether the row is **DNA-origin** (`raw_hgvs_nt` is non-blank) or **protein-origin** (only `raw_hgvs_pro` present, no `raw_hgvs_nt`).
2. Builds a list of `(hgvs_c, hgvs_g)` candidate pairs from the pipe-delimited `mapped_hgvs_c` and `mapped_hgvs_g` columns.
3. For each candidate pair, resolves a ClinGen allele ID using c. first, falling back to g. if c. is absent or returns no result.
4. **DNA-origin rows with a single non-pipe-delimited candidate:** if `clingen_allele_id` already contains a CA-prefixed ID (from step 1), it is reused directly without a network request.
5. Joins the per-candidate IDs with `|` and writes the result to `dna_clingen_allele_id`. Empty slots are preserved: `CA101||CA103` means the middle candidate had no ClinGen match.

Rows are processed concurrently and written in input order.

---

## Row classification

| Condition | Row type | `dna_clingen_allele_id` behaviour |
|---|---|---|
| `raw_hgvs_nt` non-blank, `clingen_allele_id` starts with `CA`, single candidate | DNA (reuse) | Existing `clingen_allele_id` value copied as-is |
| `raw_hgvs_nt` non-blank, any other case | DNA (query) | ClinGen queried for each c./g. candidate |
| `raw_hgvs_nt` blank | Protein-origin | ClinGen queried for each reverse-translated candidate |

ClinGen allele IDs beginning with `_:` (blank-node placeholders that ClinGen emits for unregistered alleles) are treated as absent — the slot is left empty.

IDs beginning with `PA` (protein alleles) on a DNA-origin row, or `CA` on a protein-origin row, raise a `ValueError` (misconfigured columns indicate a data pipeline issue).

---

## Lookup strategy per candidate

For each `(hgvs_c, hgvs_g)` pair:

1. Try `hgvs_c`. If it returns a CA-prefixed allele ID, use it.
2. If `hgvs_c` is blank or returns no ID, try `hgvs_g`.
3. If both fail, leave the slot empty.

ClinGen HTTP responses are cached in Redis (see [Caching](#caching)). An in-process `dict` further deduplicates identical HGVS strings within the same run without a second Redis round-trip.

---

## Caching

All ClinGen lookups go through `src/lib/clingen.py`, which uses Redis:

- **Hits** (ClinGen returns a valid allele ID): stored under `{prefix}:{hgvs}` → allele ID, and under `{prefix}:{allele_id}` → full JSON response. TTL: 1 day (configurable via `CLINGEN_CACHE_TTL_SECONDS`).
- **Misses** (404 or 200 with no allele ID): stored under the HGVS map key as `__MISS__` sentinel. TTL: 1 day (configurable via `CLINGEN_CACHE_MISS_TTL_SECONDS`). Prevents redundant network requests on re-runs.
- **Network errors** (timeouts, 5xx): not cached; will be retried on the next run.

Cache prefix defaults to `clingen:v1`. Override with `CLINGEN_CACHE_PREFIX`.

### `--known-misses-file`

An optional TSV file listing HGVS strings that are already known not to have a ClinGen record (one string per row under a `hgvs` header column). Matching strings are skipped without querying Redis or the ClinGen API. This is useful when re-running on a dataset where a prior run already established that certain HGVS strings yield no result — load the miss list from the prior run's output rather than waiting for Redis to warm up.

The file can be generated from a prior pipeline output with a short script (see [Generating a known-misses file](#generating-a-known-misses-file)).

---

## Input and output columns

### Required input columns

| Column | Description |
|---|---|
| `mapped_hgvs_c` | Pipe-delimited transcript HGVS candidates (from steps 1–2) |
| `mapped_hgvs_g` | Pipe-delimited genomic HGVS candidates (aligned to `mapped_hgvs_c`) |
| `clingen_allele_id` | ClinGen allele ID from step 1 (used for DNA-row reuse; may be blank) |
| `raw_hgvs_nt` | Original DNA HGVS input (used to classify the row type; may be blank) |
| `raw_hgvs_pro` | Original protein HGVS input (used to classify the row type; may be blank) |

### Output columns

| Column | Description |
|---|---|
| `dna_clingen_allele_id` | Pipe-delimited CA-prefixed ClinGen allele IDs, aligned to `mapped_hgvs_c`/`_g` candidates. Empty slots indicate no match for that candidate. |

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `--output-col COL` | `dna_clingen_allele_id` | Name of the output column |
| `--clingen-allele-id-col COL` | `clingen_allele_id` | Column with existing step-1 ClinGen ID (for DNA-row reuse) |
| `--hgvs-c-col COL` | `mapped_hgvs_c` | Transcript HGVS candidates column |
| `--hgvs-g-col COL` | `mapped_hgvs_g` | Genomic HGVS candidates column |
| `--raw-hgvs-nt-col COL` | `raw_hgvs_nt` | Column used to detect DNA-origin rows |
| `--raw-hgvs-pro-col COL` | `raw_hgvs_pro` | Column used to detect protein-origin rows |
| `--max-retries N` | `3` | Retries per ClinGen HTTP request |
| `--max-workers N` | `8` | Concurrent worker threads |
| `--known-misses-file FILE` | none | TSV of HGVS strings to skip without querying ClinGen |
| `--skip N` | `0` | Skip the first N data rows |
| `--limit N` | no limit | Stop after processing N rows |
| `--log-level` | `INFO` | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`) |
| `--csv-field-size-limit BYTES` | system default | Increase for large HGVS fields |

### Environment variables (from `src/lib/clingen.py`)

| Variable | Default | Description |
|---|---|---|
| `CLINGEN_CACHE_ENABLED` | `true` | Set to `false` to disable Redis caching entirely |
| `CLINGEN_CACHE_REDIS_URL` | `redis://redis:6379/0` | Redis connection URL (also falls back to `REDIS_URL`) |
| `CLINGEN_CACHE_PREFIX` | `clingen:v1` | Redis key prefix |
| `CLINGEN_CACHE_TTL_SECONDS` | `86400` | TTL for cache hits (seconds) |
| `CLINGEN_CACHE_MISS_TTL_SECONDS` | `86400` | TTL for cache misses (seconds) |
| `CLINGEN_API_URL` | `https://reg.genome.network/allele` | ClinGen Allele Registry endpoint |

---

## Concurrency

Rows are submitted to a `ThreadPoolExecutor` (default 8 workers) and written to the output file in input order. An in-process `dict` guarded by a `threading.Lock` deduplicates concurrent requests for the same HGVS string: the first thread to claim a key performs the lookup and stores the result; subsequent threads wait for the result rather than issuing duplicate requests.

---

## Generating a known-misses file

After running `add_dna_clingen_allele_ids`, you can extract HGVS strings that still have no allele ID for use as `--known-misses-file` on a subsequent run:

```python
import csv

csv.field_size_limit(10_000_000)

in_path  = "output_clingen.tsv"
out_path = "known_misses.tsv"

no_id: list[str] = []
seen: set[str] = set()

def add_miss(hgvs: str) -> None:
    h = hgvs.strip()
    if h and h not in seen:
        seen.add(h)
        no_id.append(h)

with open(in_path, newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        dna_ca_parts = (row.get("dna_clingen_allele_id") or "").split("|")
        for col in ("mapped_hgvs_c", "mapped_hgvs_g"):
            for hgvs, ca in zip((row.get(col) or "").split("|"), dna_ca_parts):
                if hgvs.strip() and not ca.strip():
                    add_miss(hgvs)

with open(out_path, "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["hgvs"])
    for h in no_id:
        writer.writerow([h])
```

Pass the result as `--known-misses-file known_misses.tsv` on the next run.

---

## Troubleshooting

**Many empty `dna_clingen_allele_id` slots**

- Some variants genuinely have no ClinGen record (especially rare or novel variants). This is expected.
- If the rate is unexpectedly high, check that `mapped_hgvs_c` / `mapped_hgvs_g` are well-formed HGVS strings from prior steps.
- Check Redis connectivity: `docker compose logs redis`.
- If the prior run hit many network errors, those were not cached. Re-run; Redis will now cache the misses from the current run.

**ClinGen rate-limit errors (HTTP 429)**

The script retries with exponential backoff automatically. For large datasets reduce `--max-workers` to lower concurrency.

**`ValueError: Row N: Invalid ClinGen allele ID prefix`**

A `clingen_allele_id` from step 1 starts with `PA` (a protein allele) for a row that `add_dna_clingen_allele_ids` classifies as a DNA-origin row, or vice versa. Check that `raw_hgvs_nt` and `raw_hgvs_pro` columns are correct for the flagged row.
