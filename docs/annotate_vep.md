# annotate_vep

`src/annotate_vep.py` — Step 9 of the variant-annotation pipeline (optional).

Annotates each variant row with a molecular consequence term from the Ensembl VEP REST API. The core design goal is to produce the most transcript-specific consequence possible: for transcript HGVS inputs (`NM_`, `NR_`, `XM_`, `XR_`, `ENST`) the resolution selects the consequence for the specific referenced transcript rather than the global worst across all overlapping transcripts.

The script is a **composition root**: it owns CSV streaming, the Redis cache, and CLI wiring. Every decision about *which* consequence applies to a variant lives in `variant_annotation/lib/vep/`, which mavedb-api imports directly — so this pipeline and the API cannot drift apart. See [architecture.md](architecture.md) for the layering and `variant_annotation/lib/sequence_ontology.py` for the term vocabulary.

---

## Output columns

| Column | Description |
|---|---|
| `vep.mutational_consequences` | `^`-delimited list of all consequence terms when the result came from a matched transcript entry; single term otherwise. Pipe-delimited across candidates. Empty for a candidate with no consequence or a failed request. |
| `vep.most_severe_mutational_consequence` | Single most-severe consequence term per candidate. Pipe-delimited across candidates. Empty for a candidate with no consequence or a failed request. |
| `vep.consequence_source` | `transcript`, `most_severe`, or `reference_identical` per candidate; empty when no consequence was determined. Pipe-delimited across candidates. |
| `vep.access_date` | ISO access date per candidate, set whenever a lookup was attempted. Pipe-delimited and position-aligned to the input candidates (empty slot for empty candidates). |
| `vep.error` | Per-candidate failure description when the candidate's request failed; empty otherwise, **including** when VEP answered with no consequence. Pipe-delimited across candidates. |

The namespace prefix defaults to `vep` and can be changed with `--vep-namespace`.

**A blank consequence with a blank error means VEP was asked and had no answer** — a settled negative. A blank consequence with a populated error means the request failed and the answer is unknown; rerun those.

---

## Transcript selection

Ensembl VEP annotates every transcript that overlaps the genomic position, not only the one referenced in the input HGVS. The top-level `most_severe_consequence` field is therefore the worst outcome across the entire locus — a single BRCA1 coding variant returns 368 `transcript_consequences` entries, and the headline term routinely describes a transcript with a different reading frame than the one that was assayed.

Resolution selects the consequence for the input's own transcript whenever one can be identified:

| Input | API flag sent | How the consequence is selected | `vep.consequence_source` |
|---|---|---|---|
| RefSeq transcript (`NM_`, `NR_`, `XM_`, `XR_`) | `refseq=1` | Matched by versionless `transcript_id` in `transcript_consequences` | `transcript` |
| Ensembl transcript (`ENST`, `LRG_`) | _(default)_ | Matched by versionless `transcript_id` in `transcript_consequences` | `transcript` |
| Genomic (`NC_`) **with** a transcript supplied | `refseq=1` when the transcript is RefSeq | Matched against the supplied transcript | `transcript` |
| Genomic (`NC_`) with no transcript | _(default)_ | `most_severe_consequence` — no transcript can be identified | `most_severe` |
| Protein (`NP_`) → Recoder fallback | _(default)_ | `most_severe_consequence` across the recoded genomic forms | `most_severe` |
| Reference-identical (`c.123=`, `p.Met4=`, unchanged `delins`) | _(none)_ | No VEP call; labelled `no_change` | `reference_identical` |

When no matching transcript entry is found (e.g. the transcript is absent from Ensembl's RefSeq set), resolution falls back to `most_severe_consequence` and records `most_severe`.

> The refseq flag is keyed on the **effective transcript**, not on the HGVS prefix. A genomic input carrying a RefSeq transcript hint is therefore sent with `refseq=1` so the hint can actually match. This CLI has no transcript-hint column — it always derives the transcript from the HGVS string — but the library supports one, and mavedb-api uses it to resolve genomic-level alleles against their paired coding transcript.

### Consequence severity ordering

Severity is Ensembl's published ranking of its 41 calculated variant consequences, reproduced in `variant_annotation/lib/sequence_ontology.py` with SO accessions so ranking is deterministic and offline. A term Ensembl adds after that list was refreshed ranks last but is **never discarded** — an unrecognised term is still a real answer, so it is reported rather than turned into a spurious miss.

Refresh the list from `https://rest.ensembl.org/info/variation/consequence_types?rank=1` when a release adds or reorders terms, and bump `RESOLVER_VERSION` in `variant_annotation/lib/vep/__init__.py` so cached answers computed under the old ranking are invalidated.

### Reference-identical (`no_change`) handling

Some HGVS expressions describe no sequence change. VEP rejects them, reasonably — they are not sequence variants — but they are ordinary in functional-assay data (a wild-type control), so they get an explicit label instead of a null that cannot be told apart from "not yet annotated".

Two forms are recognised:

- **Explicit HGVS `=` notation** (`NM_000049.4:c.123=`, `NC_000003.12:g.10141916C=`, `NP_000537.3:p.Met4=`). The notation is itself the claim, so no reference lookup is needed. This covers residues with no codon redundancy — `Met` and `Trp` — which cannot be recoded to an alternate synonymous codon at all.
- **A coding `delins` that restates the reference** (`NM_000133.4:c.1_3delinsATG` where the transcript already reads `ATG` there). Verified by querying UTA for the transcript CDS reference sequence over `start..stop` and comparing to the replacement bases.

Both yield `no_change` in the consequence columns with `vep.consequence_source = reference_identical` and an empty `vep.error`. `no_change` is deliberately **not** an SO term: `synonymous_variant` means "a change that leaves the amino acid intact", which is a different claim from "no change was described".

Prerequisite for the `delins` form: `UTA_DB_URL` must be set and reachable. Without it the explicit `=` form still works, and unchanged `delins` rows are reported as having no consequence.

**Column priority:** `mapped_hgvs_c` should appear before `mapped_hgvs_g` in `--hgvs-cols` (the default) because transcript HGVS strings yield transcript-specific consequences while bare genomic HGVS strings can only yield `most_severe`.

---

## API call strategy

For each batch of input rows:

1. **VEP POST** (`/vep/human/hgvs`) — inputs are partitioned by which transcript set they need; RefSeq inputs are sent with `refseq=1` in separate batches from Ensembl/genomic ones, because `refseq` is a per-request flag. Batches dispatch concurrently via a `ThreadPoolExecutor`.
2. **Variant Recoder POST** (`/variant_recoder/human`) — inputs VEP neither answered nor errored on are recoded to genomic HGVS equivalents (`NC_…`). Inputs whose *request failed* are **not** recoded: the answer is unknown, so a fallback would be answering a question that was never asked.
3. **Second VEP POST** — the recoded genomic HGVS are queried in a second concurrent pass. When several recoded equivalents exist, the most severe consequence across them is chosen and the source is necessarily `most_severe`.
4. **Reference-identical check** — inputs still unresolved are tested for no-change form (see above).
5. **Failures** — when a request fails after retries, only that batch's inputs are affected; every other batch's results stand. Failed inputs are reported with a populated `vep.error` and are **not cached**, so they retry on the next run.

Retries cover timeouts, connection resets, 5xx, and 429 (honouring `Retry-After`). A 4xx other than 429 is **not** retried — VEP returning 400 for a protein HGVS it cannot parse is a settled answer about the input, and retrying burns quota for the same rejection.

---

## Multi-candidate rows

For rows with pipe-delimited HGVS candidates (from step 2 reverse translation), each candidate is resolved independently. All output columns are pipe-delimited with one value per candidate position, in input order. A blank candidate yields a blank in every column rather than shifting the alignment.

---

## Redis caching

Resolved and confirmed-absent answers are cached per `(resolver version, Ensembl release, transcript, HGVS)`. Absent answers are stored under a sentinel so repeated no-hit queries don't re-query the API. **Failed requests are never cached** — caching an outage would turn a transient failure into a persistent wrong result.

The cache key carries both versioning axes, so a stored answer is reused only when it was computed under the same rules *and* the same upstream data: the library's `RESOLVER_VERSION` (a change to the resolution rule or severity ranking) and the Ensembl release (`/info/software`, fetched once per run). Either bump invalidates the affected answers automatically, rather than relying on an operator to remember to bump the prefix. The transcript is in the key because the same HGVS resolved against two transcripts is two different questions. If the release cannot be determined, the key uses `eunknown` and those entries simply never collide with release-keyed ones.

| Variable | Default | Description |
|---|---|---|
| `VEP_CACHE_ENABLED` | `true` | Set to `0` / `false` to disable caching entirely |
| `VEP_CACHE_REDIS_URL` | `redis://redis:6379/0` | Redis connection URL (also falls back to `REDIS_URL`) |
| `VEP_CACHE_PREFIX` | `vep:v3` | Key namespace; the version suffix tracks the cached *value shape* |
| `VEP_CACHE_TTL_SECONDS` | `86400` | TTL for hits (seconds) |
| `VEP_CACHE_MISS_TTL_SECONDS` | `86400` | TTL for confirmed-absent answers (seconds) |

Entries written before `vep:v3` used a different value shape and are simply never consulted under the current prefix.

---

## CLI options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--vep-namespace NS` | — | `vep` | Output column name prefix |
| `--hgvs-cols COLS` | — | `mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p` | Comma-separated columns tried in priority order; first non-blank value used per row |
| `--vep-api-url URL` | `ENSEMBL_API_URL` | `https://rest.ensembl.org` | Ensembl REST API base URL |
| `--vep-batch-size N` | `VEP_BATCH_SIZE` | `200` | HGVS values per VEP POST request; **200 is Ensembl's documented maximum** and a larger value is rejected |
| `--vep-workers N` | `VEP_WORKERS` | `8` | Concurrent VEP/Recoder batch requests |
| `--row-batch-size N` | `VEP_ROW_BATCH_SIZE` | `1000` | Input rows per lookup/write batch |
| `--vep-timeout-seconds N` | `VEP_TIMEOUT_SECONDS` | `60` | HTTP timeout per VEP request (Recoder gets its own longer budget) |
| `--no-recoder` | — | off | Skip the Variant Recoder fallback; report unresolvable inputs as having no consequence rather than accepting a cross-transcript answer |
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
- **UTA** (optional, via `UTA_DB_URL`) — only for reference-identical `delins` detection. One connection is opened for the whole run.
- **Redis** (optional) — greatly reduces repeat API calls across runs. Fails gracefully when unavailable.

---

## Troubleshooting

**`vep.error` populated on every run**

- The message includes the HTTP status and a response excerpt to identify the specific failure.
- Ensembl REST API may be under maintenance or rate-limiting. Check `https://rest.ensembl.org` directly.
- Reduce `--vep-workers` or increase `--vep-timeout-seconds`.
- Failed requests are never cached, so they always retry. Use `--keep-existing` to skip already-annotated rows and only retry the blanks.

**Consequence is `most_severe` for all rows despite transcript HGVS input**

- Confirm `mapped_hgvs_c` is listed first in `--hgvs-cols`.
- The transcript may be absent from Ensembl's current RefSeq set, triggering the fallback. This is common for older transcript versions. `vep.consequence_source` records when it happens.

**A consequence looks wrong for the gene that was assayed**

Check `vep.consequence_source`. A `most_severe` value is the worst call across every overlapping transcript and may not describe the assayed transcript at all. A `transcript` value is specific to the input's own transcript and is the trustworthy case.

**After an Ensembl release**

The Ensembl release is part of the cache key, so a release bump invalidates the affected answers automatically — a run against the new release re-resolves rather than serving a prior-release answer, with no prefix or TTL action needed. If the release *adds or reorders* consequence terms, also refresh `ENSEMBL_CONSEQUENCE_RANKING` and bump `RESOLVER_VERSION` so answers computed under the old ranking (including any still cached under the prior release key) are invalidated too.
