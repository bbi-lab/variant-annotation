# Differences from MaveDB 2026.1.2

This document describes how each script in this pipeline diverges from the equivalent logic in the MaveDB production codebase (tag 2026.1.2). The intent is to make it easy to evaluate whether a change in MaveDB should be ported here, or vice versa.

---

## map_variants

**MaveDB equivalents:** `map_variants_for_score_set` (`worker/jobs/variant_processing/mapping.py`) and `populate_hgvs_for_score_set` (`worker/jobs/external_services/hgvs.py`).

### Architecture

| | MaveDB | This pipeline |
|---|---|---|
| Execution model | ARQ background job, triggered asynchronously by the API when a score set is published | Standalone CLI / Docker container, invoked directly on a TSV file |
| Storage | Results stored in PostgreSQL (`MappedVariant`, `TargetGeneMapping` tables) | Results written as columns in the output TSV |
| dcd_mapping access | HTTP API service (`vrs_mapper().map_score_set(urn)`) | Python library imported directly |
| Progress tracking | `job_manager.update_progress()` writes to the job queue | Log lines every 1 000 rows |

### Transcript accession selection

This is the most significant semantic difference.

**MaveDB:** One transcript accession is determined per *score set*, not per row. `get_target_coding_info()` checks the score set's `target_accession` (for accession-based targets) or the `cdna.sequence_accessions` field in `post_mapped_metadata` (for sequence-based targets). All variants in the score set use that single accession when extracting `hgvs_c` from ClinGen. If `post_mapped_metadata` contains more than one cDNA accession, the job raises `ValueError`. Multi-target score sets are explicitly not supported (`NotImplementedError`) — `populate_hgvs_for_score_set` skips them entirely.

**This pipeline:** The transcript accession is extracted per row from `raw_hgvs_nt` (e.g. `NM_000277.3` from `NM_000277.3:c.1218G>A`). Different rows may use different transcripts. Multi-target inputs (multiple gene groups with different sequences/accessions) are naturally supported.

### VRS metadata persisted

**MaveDB:** Stores full `pre_mapped` and `post_mapped` VRS JSON objects in `MappedVariant`. Also stores per-target alignment QC statistics (`TargetGeneMapping`: alignment score, percent identity, mismatch count, gap count, etc.) and the `mapping_api_version` (dcd_mapping tool version).

**This pipeline:** Stores only the VRS digest string derived from `post_mapped.digest` (or `post_mapped.id`). No alignment QC statistics are retained. The columns `dna_vrs_digest` and `protein_vrs_digest` exist in the output but carry less information than MaveDB's full VRS payloads.

### ClinGen HGVS population (two-step vs. one-step)

**MaveDB:** The mapping and HGVS population steps are separate background jobs. `map_variants_for_score_set` runs VRS mapping and stores `hgvs_assay_level` from the post-mapped VRS expression. `populate_hgvs_for_score_set` then runs separately, querying ClinGen *by stored ClinGen allele ID* to populate `hgvs_g`, `hgvs_c`, `hgvs_p`. Variants without a `clingen_allele_id` are skipped. Variants with comma-separated (multi-variant) allele IDs are also skipped.

**This pipeline:** VRS mapping and ClinGen HGVS population happen in a single pass. ClinGen is queried by the assay-level HGVS string (not by allele ID) immediately after VRS mapping. There is no concept of deferring the ClinGen lookup to a separate job. Multi-variant comma-separated allele IDs are not a concern because we query by HGVS string, not by stored allele ID.

### Concurrency model

**MaveDB:** `populate_hgvs_for_score_set` queries ClinGen sequentially per variant (one `await get_clingen_allele_data(clingen_id)` per loop iteration). `map_variants_for_score_set` delegates to the dcd-mapping API service which handles its own parallelism.

**This pipeline:** Case-1 rows are batched and queried concurrently via `asyncio` + `asyncio.Semaphore` up to `--max-clingen-concurrency` (default 5). Sequence-based groups are processed one at a time (BLAT is single-threaded per run).

### Multi-variant haplotype normalization

**MaveDB:** No special handling in the mapping or HGVS population jobs. Intra-codon c.-haplotypes (e.g. `c.[1A>G;3G>T]`) would be passed as-is to dcd_mapping.

**This pipeline:** Intra-codon c.-haplotypes are detected and normalized to a `delins` expression (e.g. `c.1_3delinsGA`) before VRS mapping. The normalization is attempted only when all component substitutions fall within the same codon.

### Protein HGVS normalization

**MaveDB:** Variants arrive at the mapping job already validated and normalized by the `create_variants_for_score_set` job upstream. One-letter protein HGVS codes are not expected at this stage.

**This pipeline:** Protein HGVS strings in one-letter amino acid code (e.g. `p.A406T`) are automatically converted to three-letter code (`p.Ala406Thr`) before processing, to match the format expected by dcd_mapping and ClinGen.

### Ensembl ENST → RefSeq resolution

**MaveDB:** `get_target_coding_info` accepts ENST accessions directly and passes them to ClinGen HGVS extraction. If no MANE transcript is found matching the ENST accession, the MANE entry from ClinGen is used as a fallback.

**This pipeline:** If `raw_hgvs_nt` contains an ENST-prefixed accession, `_normalize_transcript_accession()` first calls the Ensembl REST API (`/xrefs/id/{enst_base}?external_db=RefSeq_mRNA`) to convert the ENST to a RefSeq NM_ accession. ClinGen is then queried with the NM_ accession. This extra lookup is performed because ClinGen transcript allele matching is more reliable against NM_ accessions than ENST IDs.

### Failure modes and output

**MaveDB:** Per-job failure sets `score_set.mapping_state` to `failed` or `incomplete`. Variant-level failures are recorded in `AnnotationStatus` rows (with `AnnotationFailureCategory`). Failed variants still produce `MappedVariant` rows (with `error_message` and `pre_mapped`/`post_mapped` as null).

**This pipeline:** Per-row failures are written to the `mapping_error` column in the output TSV. Rows with errors still appear in the output with blank mapped HGVS columns. No job-level state machine.

### Resume / partial-run support

**MaveDB:** Jobs are retried automatically by ARQ up to `job.max_retries` times on failure.

**This pipeline:** Supports `--skip N` to resume from a specific row, and `--merge-existing prior.tsv` to reuse results from a partial prior run without reprocessing matched rows.

---

## reverse_translate_protein_variants

**MaveDB equivalent:** None — MaveDB does not perform reverse translation.

MaveDB relates DNA and protein variants through the ClinGen Allele Registry: `populate_hgvs_for_score_set` retrieves `hgvs_c` and `hgvs_p` for a variant by querying ClinGen with the stored `clingen_allele_id`. This means protein-level HGVS (`hgvs_p`) is derived from whatever ClinGen returns for a DNA allele ID, not by enumerating synonymous codons.

**Scope limitation of the MaveDB approach:** ClinGen's coverage is concentrated on SNVs and alleles that have been explicitly submitted by submitters or imported from curated databases (ClinVar, dbSNP, LOVD, etc.). Variants that have never been submitted — especially rare or novel missense changes — often have no ClinGen record at all. For these variants `populate_hgvs_for_score_set` simply skips the row (no `hgvs_p` is populated). There is no mechanism to recover DNA candidates for a protein variant that lacks a pre-existing ClinGen allele.

**This pipeline:** `reverse_translate_protein_variants` performs an exhaustive reverse translation using the `reverse-translate-variants` CLI. For a given `(transcript, p.AminoChange)` pair it enumerates every synonymous codon in the genetic code and emits a pipe-delimited list of `hgvs_c` / `hgvs_g` candidates. This works regardless of whether the variant has ever been submitted to ClinGen. The resulting candidates are then queried against ClinGen in step 3 (`add_dna_clingen_allele_ids`) — but even candidates that have no ClinGen record still appear in the output with an empty allele-ID slot, preserving them for downstream annotation steps that don't require a ClinGen ID (e.g. VEP, gnomAD lookups by position).

| | MaveDB | This pipeline |
|---|---|---|
| How protein ↔ DNA link is made | Via stored ClinGen allele ID | Exhaustive codon-level reverse translation |
| Coverage | Only variants with a ClinGen record | All theoretically possible DNA changes for the amino acid substitution |
| Indels | Not generated; depends on what was submitted to ClinGen | Optional via `--include-indels --max-indel-size N` |
| Novel / unsubmitted variants | Skipped | Handled (candidates produced, ClinGen lookup attempted separately) |
| Output cardinality | One `hgvs_c` / `hgvs_p` per variant | One or more pipe-delimited candidates per protein row |

---

## `annotate_gnomad` (step 6)

**MaveDB** (`mavedb-api/src/mavedb/lib/gnomad.py`): queries gnomAD using AWS Athena exclusively. The query selects only `caid`, `joint.freq.all.ac`, `joint.freq.all.an`, `joint.fafmax.faf95_max_gen_anc`, and `joint.fafmax.faf95_max`. No QC filter columns are fetched. Lookups are CAID-only. No local cache is used; every run queries Athena.

**This pipeline:** `annotate_gnomad` offers two execution backends (Hail and Athena) and two lookup strategies (coordinate-based or CAID-based), and the Athena backend's query structure deliberately matches MaveDB's for compatibility.

Key differences:

| | MaveDB | This pipeline |
|---|---|---|
| Execution backend | Athena only | Hail (default) or Athena |
| Local cache | None | Hail mode builds an indexed local cache (one-time write, reused on subsequent runs) |
| Lookup strategy | CAID only | `coordinates` (default) or `caid`; coordinate lookup covers variants that have no CAID in gnomAD |
| QC filtering | Not available | `--require-pass` (combined filters) and `--callset-pass-filter any/all` (per-callset); Hail mode only |
| QC filter columns in output | Not produced | `filters`, `exome_filters`, `genome_filters` populated in Hail mode; always empty in Athena mode |
| Additional output fields | Not produced | `minor_allele_frequency`, `gene_symbols` |
| Athena query fields | `caid`, AC, AN, faf95_max, faf95_max_gen_anc | Identical when using `--execution-mode athena` |

### Why filtering matters

MaveDB returns raw frequency data regardless of gnomAD QC status. A variant that failed gnomAD quality control (e.g. `AC0`, `AS_VQSR`) will be annotated with potentially unreliable frequency values in MaveDB. This pipeline allows callers to exclude such variants from annotation via `--require-pass` and `--callset-pass-filter`, which leave the frequency columns empty for QC-failed variants rather than propagating potentially misleading values.

Filtering is **only available in Hail mode**. The Athena query (like MaveDB's) does not fetch filter columns, so `--require-pass` and `--callset-pass-filter` are silently no-ops when `--execution-mode athena` is used.

### Gene-level cache filtering

Hail mode supports `--genes BRCA1,BRCA2` to restrict the local cache to specific gene symbols (from `vep.worst_csq_by_gene_canonical`). This has no MaveDB equivalent.

---

## `annotate_spliceai` (step 7)

**MaveDB:** SpliceAI annotation has not been implemented in MaveDB.

**This pipeline:** `annotate_spliceai` adds nine SpliceAI delta-score columns (`spliceai.ds_ag`, `spliceai.ds_al`, `spliceai.ds_dg`, `spliceai.ds_dl`, `spliceai.dp_ag`, `spliceai.dp_al`, `spliceai.dp_dg`, `spliceai.dp_dl`, `spliceai.max_delta_score`) to each row. Two execution modes are available: `precomputed` (default) — tabix lookup against the Illumina pre-scored SpliceAI VCFs — and `compute` — local on-the-fly scoring using the SpliceAI deep-learning model. For rows with multiple pipe-delimited DNA candidates (from step 2), all output columns are pipe-delimited and candidate-aligned.

---

## `annotate_erepo` (step 8)

**MaveDB:** ClinGen Evidence Repository annotation has not been implemented in MaveDB.

**This pipeline:** `annotate_erepo` downloads the full ClinGen erepo expert-panel classification TSV and joins it against each variant candidate using up to three keys: HGVS expression, ClinVar Variation ID, and CAID. Sixteen classification columns are added per candidate (prefixed `clingen_evidence_repository.`), including `Assertion`, `Expert Panel`, `Disease Mondo Id`, `Mode of Inheritance`, applied ACMG evidence codes, and supporting metadata. A `warnings` column records cross-key discrepancies.

---

## `annotate_vep` (step 9)

**MaveDB** (`mavedb-api/src/mavedb/lib/vep.py`): queries Ensembl VEP via `/vep/human/hgvs` and falls back to Variant Recoder → second VEP pass for unresolved inputs. Always uses the top-level `most_severe_consequence` field, which reflects the worst outcome across **all** transcripts overlapping the variant's genomic position. A single genomic HGVS or `hgvs_nt` value is queried per variant.

**This pipeline:** `annotate_vep` follows the same two-step API strategy (VEP → Recoder → VEP) but adds transcript-specific consequence selection and handles pipe-delimited multi-candidate rows.

Key differences:

| | MaveDB | This pipeline |
|---|---|---|
| Consequence selection | Always `most_severe_consequence` (global worst across all transcripts at the locus) | Transcript-specific when input is `NM_`/`NR_`/`ENST` HGVS; falls back to `most_severe_consequence` otherwise |
| Output consequence columns | Single consequence term (stored in the database, no named TSV column) | `vep.most_severe_mutational_consequence` (single most-severe term per candidate) + `vep.mutational_consequences` (`^`-delimited full list of terms from the matched transcript entry; single most-severe term otherwise) |
| RefSeq transcript handling | Not distinguished from genomic HGVS | Sends `refseq=1` flag for `NM_`/`NR_` inputs so `transcript_consequences` uses RefSeq IDs |
| `vep.consequence_source` output column | Not produced | `transcript` when a matched transcript entry was used; `most_severe` for global fallback |
| HGVS input priority | Single HGVS string per variant | Configurable column list (`--hgvs-cols`); default `mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p` — transcript HGVS tried first |
| Multi-candidate rows | Not applicable | Each pipe-delimited candidate resolved independently; all output columns pipe-delimited per candidate; `vep.error` contains per-candidate API error messages (e.g. `api_error:VEP HTTP 503: …`) |
| Redis cache structure | Consequence string only | `(most_severe, all_consequences, source)` triple; misses stored as explicit sentinel; API errors skipped (never cached); entries from the prior `(consequence, source)` format (missing `all_consequences`) are silently discarded on read and re-queried |
| Resume support | Not built-in | `--keep-existing` skips rows with an existing non-empty annotation |

### Why transcript specificity matters

When a gene has many overlapping transcripts at a locus, `most_severe_consequence` may reflect a consequence on a minor isoform rather than the canonical transcript of interest. For example, a variant in a coding exon of the canonical transcript may be annotated as `missense_variant` on that transcript but `intron_variant` on a longer overlapping non-coding transcript — and vice versa. Using `most_severe_consequence` would pick whichever is most severe globally, which may not be the biologically relevant consequence for the experiment's transcript. This pipeline preferentially selects the consequence on the specific transcript named in the input HGVS (`mapped_hgvs_c`).

---

## `annotate_mavedb` (step 10)

**MaveDB:** stores all calibrations for a score set and exposes all of them (primary, investigator-provided, research-use-only, alternatives) through the API and UI. The MaveDB worker applies calibrations at ingest time using whichever calibrations are associated with a score set.

**This pipeline:** `annotate_mavedb` selects specifically the **primary** calibration and the **investigator-provided** calibration, annotating each variant with the functional classification label under those two calibrations only.

Key differences:

| | MaveDB | This pipeline |
|---|---|---|
| Calibrations applied | All associated calibrations | Primary + investigator-provided only |
| Additional calibrations (alternative thresholds, RUO) | Exposed in the UI and stored | Fetched but not applied |
| When primary = investigator-provided | One result | Both column groups populated with the same values |
| Calibration type support | Range-based and class-based | Range-based (local score comparison) and class-based (API variant-to-class lookup) |
| Caching | Persisted in the MaveDB database | In-memory per run only; no on-disk cache |

### Calibration selection rationale

The primary calibration is the one the MaveDB curators have designated as the recommended interpretation threshold. The investigator-provided calibration reflects the thresholds determined by the original study authors. Exposing both independently allows downstream users to compare the curator-recommended and author-recommended classifications, which sometimes differ. Research-use-only and alternative calibrations are omitted because they are not intended for clinical or variant-interpretation use.

---

## `annotate_predictors` (step 11)

**MaveDB:** does not incorporate in-silico predictor scores (REVEL, AlphaMissense, MutPred2, etc.). MaveDB focuses on experimentally-derived functional scores from multiplexed assays (MAVE data), not computational predictors.

**This pipeline:** `annotate_predictors` adds pre-computed computational pathogenicity predictions from REVEL, AlphaMissense, and/or MutPred2 as supplementary columns. These are entirely separate from MaveDB scores and calibrations.
