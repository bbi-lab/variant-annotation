# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [1.2.0-beta.1] - 2026-10-15

### Changed

- Every coding↔genomic projection goes through `variant_annotation.lib.translation.coding_to_genomic` / `genomic_to_coding`, so reverse translation's candidates, its WT-codon candidate and consumers' own projections follow one set of rules. The `reverse-translate-variants` CLI's genomic column is no longer read.
- `CoordinateTranslator` gains `g_to_c_literal` and `c_to_n`; `HgvsMapper` implements both.
- A genomic projection that hgvs cannot make for a reason other than the variant (an unreachable data provider) is raised from `construct_equivalent_variants` instead of being reported as a candidate with no genomic form.

### Fixed

- A coding change across an exon junction (a codon split by an intron, an insertion between two exons' bases, or the WT codon at a split codon) no longer projects to a genomic change that deletes the intron. It has no genomic projection.

## [1.1.0] - 2026-10-12

### Changed

- Bumped the `variant-translation` dependency to 1.2.0, picking up improved reverse-translation behavior in `reverse_translate_protein_variants`: non-adjacent delins candidates can now be expressed as genomic multivariant alleles (`g.[...]`), and protein multivariant (`p.[...]`) input is recognized.

## [1.0.0] - 2026-10-12

First stable release. No functional changes since v1.0.0-beta.2.

## [1.0.0-beta.2] - 2026-10-12

### Added

- Author ORCID iD in `CITATION.cff`.

## [1.0.0-beta.1] - 2026-10-12

First beta release of the variant annotation pipeline.

### Added

- Core annotation pipeline:
  - `map_variants` — map input variants to HGVS g./c./p. representations.
  - `reverse_translate_protein_variants` — generate DNA candidates for protein-level variants.
  - `add_dna_clingen_allele_ids` — resolve DNA-level ClinGen Allele Registry IDs.
  - `add_vcf_identifiers` (optional) — parse VCF-style chrom/pos/ref/alt fields.
  - `annotate_clinvar` (optional) — ClinVar clinical significance and star ratings.
  - `annotate_gnomad` (optional) — gnomAD allele frequencies.
  - `annotate_spliceai` (optional) — SpliceAI delta scores.
  - `annotate_erepo` (optional) — ClinGen Expert Panel classifications.
  - `annotate_vep` (optional) — Ensembl VEP mutational consequence.
  - `annotate_mavedb` (optional) — MaveDB functional classifications and calibrations.
  - `annotate_predictors` (optional) — REVEL, AlphaMissense, and MutPred2 scores.
  - `flatten_dna_variants` (optional) — collapse to one row per DNA variant.
- Bundled utilities for TSV/CSV manipulation, dataset comparison, cache maintenance, and genome assembly liftover (`utilities`, `diff_hgvs`, `diff_tsv`, `diff_vcf_coords`, `compare_cvfg_datasets`, `filter_dna_variants`, `haplotypes_to_delins`, `backfill_clingen_allele_ids`, `clear_clingen_cache`, `liftover`, `gnomad_hail_shell`).
- Redis-backed response caching across API-calling pipeline stages, with graceful degradation when Redis is unreachable.
- Docker Compose-based development and execution environment, including a dedicated Hail-enabled service for gnomAD lookups.

[Unreleased]: https://github.com/bbi-lab/variant-annotation/compare/v1.2.0-beta.1...HEAD
[1.2.0-beta.1]: https://github.com/bbi-lab/variant-annotation/compare/v1.1.0...v1.2.0-beta.1
[1.1.0]: https://github.com/bbi-lab/variant-annotation/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/bbi-lab/variant-annotation/compare/v1.0.0-beta.2...v1.0.0
[1.0.0-beta.2]: https://github.com/bbi-lab/variant-annotation/compare/v1.0.0-beta.1...v1.0.0-beta.2
[1.0.0-beta.1]: https://github.com/bbi-lab/variant-annotation/releases/tag/v1.0.0-beta.1
