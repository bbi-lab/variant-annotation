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
