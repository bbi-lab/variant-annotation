# mavedb_and_igvf_data_prep

This pipeline prepares MAVE (Multiplexed Assay of Variant Effect — a family of assays, including saturation genome editing (SGE), that measure the functional effect of many variants at once) data from IGVF-affiliated labs for submission to **[MaveDB](https://mavedb.org)** and the **[IGVF DataPortal](https://data.igvf.org)**. The real lab datasets it runs against are not included in this repository; what's here is a small, runnable **example** that exercises the same code end to end, so the pipeline can be understood and tried out.

`generate_mavedb_and_igvf_files.sh` in this directory turns a raw MAVE/SGE dataset (Excel, CSV, or TSV — one row per variant, with local `seq`/`reference` windows) into two kinds of output: a MaveDB submission bundle (using RefSeq transcript accession numbers; split into scores/counts/functional-classification CSVs) and an IGVF TSV (Ensembl-translated, same shape as the source).

This directory contains:

- `generate_mavedb_and_igvf_files.sh` — a generic copy of [`src/scripts/generate_mavedb_and_igvf_files.sh`](../../src/scripts/generate_mavedb_and_igvf_files.sh). It calls the exact same underlying tools in `src/`; the only difference is its default data directory points here instead of the lab data directory.
- `data/CARD11_example.tsv` — a small, hand-picked input dataset (17 rows) that exercises every special case the pipeline handles (see below).
- `data/*` — the generated output, already run once so you can inspect it without needing Docker/UTA running.

## Running it

From the repository root, with the Docker services up (`docker compose up -d uta seqrepo`, or just let `docker compose run` start them on demand):

```bash
examples/mavedb_and_igvf_data_prep/generate_mavedb_and_igvf_files.sh
```

Or point it at any other directory of MAVE/SGE source files:

```bash
examples/mavedb_and_igvf_data_prep/generate_mavedb_and_igvf_files.sh path/to/your/data
```

---

## Source data shape

Each input file is one MAVE/SGE experiment: one row per variant (or per wild-type control read), with columns including `seq` (observed local read), `reference` (local target window), `coding variant` (RefSeq `c.` HGVS, using a RefSeq transcript that is **not necessarily the true reference** — see [Fixed-edit reference correction](#fixed-edit-reference-correction) below), `prot`, `genomic variant`, `variant_type`, `region` (a tiling-window ID), replicate counts, fold-change scores, and a functional classification.

The data this example is modeled on comes from **SGE (saturation genome editing)**, a MAVE assay that edits variants directly into their endogenous genomic locus, rather than into a synthetic construct. That has two consequences for the shape of the data this pipeline handles:

- **Variants are specified at nucleotide resolution.** SGE variants are DNA edits, not amino-acid substitutions, so there's a full nucleotide-level HGVS description available for every variant. This pipeline's MaveDB `hgvs_nt` column carries a transcript-qualified HGVS c. string, and its IGVF output carries HGVS c. (transcript), g. (genomic), and p. (protein) strings side by side. Other MAVE assay types operate at amino-acid resolution instead (e.g. a deep mutational scan built from a synthesized variant-protein library) — for those, both outputs might carry only HGVS p. strings, since there's no nucleotide-level ground truth to derive c./g. from.
- **Because the edits land in the endogenous genome, MaveDB submissions use fully-qualified HGVS.** MAVE data from libraries built by editing an endogenous locus is uploaded to MaveDB as HGVS c. or g. strings that carry an accession (`NM_032415.7:c.462A>C`, not bare `c.462A>C`) — the accession itself supplies the reference, so no extra sequence needs to accompany the upload. This is why every coding-variant value this pipeline produces keeps its accession prefix throughout (RefSeq for MaveDB, Ensembl for IGVF). By contrast, a MAVE conducted against a synthetic target sequence — not tied to a specific genomic locus — is uploaded to MaveDB as *unqualified* MAVE-HGVS strings (bare `c.462A>C`, no accession) plus the target sequence itself; MaveDB then maps that sequence onto the genome to establish coordinates.

Several idiosyncrasies of this source format, all first noticed in the real lab data, drive most of the pipeline's design — and all of them are reproduced in `data/CARD11_example.tsv` (see [About the example data](#about-the-example-data)):

- Accession versions use `_` instead of `.` (`NM_032415_7`, not `NM_032415.7`).
- Wild-type/"unchanged" control reads are recorded as a bare `ACCESSION:c.=`, with no position — but different control rows cover different windows.
- The assay's local `reference` window is sometimes a deliberately recoded (fixed-edit) sequence, not the true transcript sequence, at specific positions.
- A few rows encode two co-occurring variants as `ACCESSION:c.2523A>C;2533A>T` — malformed HGVS (missing the bracket syntax for a haplotype).
- Some data contains inversions (`c.2687_2688inv`), which MAVE-HGVS does not support.
- Overlapping tiled sequencing windows can report the same coding position twice, with different local reference sequences.
- Manual data entry occasionally drops a reference letter from a substitution call (`c.865-4>G` instead of `c.865-4C>G`).

**`seq`/`reference` aren't guaranteed.** Steps 0 and 1 (`fix_missing_ref_allele.py`, `expand_unchanged_hgvs_c.py`) both work by comparing a row's `seq` against its `reference` — that's how a missing reference letter or an unpositioned `c.=` control gets reconstructed. Not every lab's export includes those two columns; some report only the HGVS calls and counts. Against a dataset without them, steps 0 and 1 don't apply — trying to run them would fail with a missing-column error, so they'd need to be left out of the invocation (or their column names remapped via `--seq-col`/`--reference-col`, if the data has equivalent columns under different names). This has no effect on steps 2 onward: `correct_reference_alleles.py` reads the true reference from UTA regardless of what's in the source file, so fixed-edit correction, haplotype handling, `inv`→`delins`, and deduplication all work the same either way — only the *typo-repair* and *WT-control-expansion* steps depend on `seq`/`reference` being present.

**`region` is never actually read.** SGE datasets tile a target into overlapping windows and label each row with a `region` ID for that window — useful for a human skimming the file, and used in prose throughout this document to talk about "which window" a row belongs to. But no script in this pipeline parses its value: step 1 groups rows by matching `reference` *text* (identical, or a substring of a longer window), not by `region` equality, and nothing else touches it at all. It's carried through as an ordinary passthrough column — present in `<stem>.igvf.tsv`, present in the intermediate `<stem>.mavedb.csv`, and then dropped only because the final `rename-columns` split never lists it among the columns to keep. A dataset without a `region` column (or an equivalent) would run through this pipeline exactly the same way.

## Pipeline overview

Both outputs are derived from the **same corrected data** — steps run once per dataset, then the pipeline branches:

```
source (.xlsx/.csv/.tsv)
  │
  ▼
0. fix_missing_ref_allele.py      reconstruct a missing reference letter from seq/reference
  ▼
1. expand_unchanged_hgvs_c.py     expand bare "c.=" controls into position ranges
  ▼
2. translate_hgvs_accessions.py   normalize accession version separator (_ → .)
  ▼
3. correct_reference_alleles.py   fix reference-allele mismatches from fixed edits (Docker/UTA)
  │
  ├─────────────────────────────┐
  ▼ IGVF branch                 ▼ MaveDB branch
4a. translate_hgvs_accessions.py 4b. filter_columns.py         drop prot, genomic variant
    translate to Ensembl (MANE)  ▼
    keep genomic as RefSeq       5b. inversions_to_delins.py   inv → delins (Docker/UTA)
    add Ensembl protein ID       ▼
  ▼                              6b. dedupe_mavedb_variants.py deduplicate (see below)
5a. filter_columns.py            ▼
    drop seq, reference          7b. rename-columns × 3        split scores/counts/functional_classes
  ▼                              ▼
<stem>.igvf.tsv                  <stem>.mavedb.csv
                                  <stem>.mavedb.{scores,counts,functional_classes}.csv
                                  <stem>.mavedb.duplicates_removed.tsv
```

Two steps (3 and 5b) need [UTA](https://github.com/biocommons/uta) (the true transcript/genomic reference sequence) and run via Docker, reusing existing services in `compose.yaml`: step 3 uses the `correct-reference-alleles` service, and step 5b's `inversions_to_delins.py` piggybacks on the `map-variants` service. Everything else runs directly with `python3`.

## Step-by-step

### 0. `fix_missing_ref_allele.py` — reconstruct a missing reference letter

Occasionally a manual-entry typo drops the reference letter from a substitution call entirely, e.g. `c.865-4>G` instead of `c.865-4C>G`. This is reconstructed directly from the row's own `seq`/`reference` columns: the single position where they differ gives the missing reference letter (and is cross-checked against the malformed cell's own alternate letter before accepting it), before anything else in the pipeline touches the cell.

### 1. `expand_unchanged_hgvs_c.py` — expand wild-type controls

A wild-type control row is reported as `ACCESSION:c.=` with no position. Different control rows cover different tiled windows, distinguishable from each other only by their `reference` text (not by `region` — see the note above), so a bare `c.=` alone doesn't say *which* window is unchanged.

For each `c.=` row, the script finds another row elsewhere in the file whose `reference` text is either identical, or a substring of a longer resolved window, **and** which has a resolvable single-nucleotide substitution (e.g. `c.208G>A`). That substitution's known position, minus its offset within the shared window, pins down the window's absolute start; the row's own reference length gives the end.

```
NM_003921.5:c.=  →  NM_003921.5:c.208_243=
```

Rows that can't be anchored raise an error by default (`--on-unresolved leave` to pass them through unchanged with a warning instead).

### 2. `translate_hgvs_accessions.py` (normalize) — fix the version separator

```bash
translate_hgvs_accessions.py IN OUT --hgvs-c-col "coding variant" --hgvs-c-mode keep --normalize-version-separator
```

`--hgvs-c-mode keep` means: don't translate the accession's namespace, just normalize `NM_032415_7` → `NM_032415.7`. This same call also normalizes any malformed multi-variant HGVS it encounters — see [Malformed haplotype syntax](#malformed-haplotype-syntax).

### 3. `correct_reference_alleles.py` — fix fixed-edit reference mismatches

See [Fixed-edit reference correction](#fixed-edit-reference-correction) below — this is the most involved step. Runs via `src/scripts/run_correct_reference_alleles.sh` (Docker, needs UTA + SeqRepo).

### 4a/5a. IGVF branch

`translate_hgvs_accessions.py` again, now `--hgvs-c-mode translate --direction to-ensembl` (MANE-based transcript translation; errors if the accession isn't a MANE transcript), `--hgvs-g-mode keep` (genomic variant stays RefSeq, just version-normalized), `--hgvs-p-mode add --hgvs-p-accession ensembl` (prepends an Ensembl protein accession to the bare `prot` cell, inferred from the transcript in the same row), and `--drop-blank-columns` (drops stray empty columns some source spreadsheets have). Then `filter_columns.py` drops `seq`/`reference`. Result: same columns as the source, coding/protein translated to Ensembl.

Correcting reference alleles *before* translating to Ensembl (rather than after) is deliberate: MANE Select transcripts are sequence-identical between RefSeq and Ensembl by definition, so correcting against the RefSeq reference and translating the accession afterward is valid, and avoids needing to know whether UTA resolves Ensembl accessions directly.

### 4b–7b. MaveDB branch

1. `filter_columns.py --omit-col prot,"genomic variant"` — MaveDB doesn't want these.
2. `inversions_to_delins.py` — MAVE-HGVS incompatibility fix (below).
3. `dedupe_mavedb_variants.py` — MaveDB requires one row per variant (below).
4. Three `rename-columns` calls (via `src.utilities`) split the result:

   | Output | Columns |
   |---|---|
   | `<stem>.mavedb.scores.csv` | `hgvs_nt`, `Rep1_FC`, `Rep2_FC`, `Rep3_FC`, `score` (renamed from `Functional Score`) |
   | `<stem>.mavedb.counts.csv` | `hgvs_nt`, every column whose name contains `counts` |
   | `<stem>.mavedb.functional_classes.csv` | `hgvs_nt`, `class_name` (renamed from `Classification`) |

   `coding variant` is renamed to `hgvs_nt` in all three. `seq`, `reference`, `variant_type`, `region` are dropped implicitly (keep-mode never lists them).

---

## Fixed-edit reference correction

### The problem

MAVE/SGE target constructs sometimes contain deliberately engineered edits — e.g. silent recoding to remove a restriction site — baked into the assay's local `reference` sequence. These fixed edits are often excluded from analysis entirely, with any variant that falls on a fixed-edit site simply omitted. In this dataset, though, mutations were deliberately introduced at those sites too, and the investigators chose to keep them, describing them as variants relative to the MANE Select transcript rather than relative to the actual, edited background genome. So when the saturation-mutagenesis scan also mutates one of those already-edited positions, the variant call inherits the *engineered* reference base instead of the transcript's true one. Confirmed against the real `NM_032415.7` (CARD11) sequence, and reproduced in `data/CARD11_example.tsv`:

- True reference at c.462 is **C**. The engineered target has **A** there instead.
- The scan reports `c.462A>C` (which actually just restores the true base — not a real variant) and `c.462A>T` (a real variant, but mislabeled relative to the wrong reference).
- Correct descriptions: `c.462=` (HGVS "predicted no change") and `c.462C>T`.

### How it's detected and fixed

`correct_reference_alleles.py` looks up the transcript's *true* reference base directly from UTA (not from the assay's `reference` column) for every single-nucleotide-substitution position, and rewrites:

| True reference vs. reported | Result |
|---|---|
| Matches reported reference | No fixed edit here — unchanged |
| Matches reported alternate | `ACCESSION:c.<pos>=` |
| Matches neither | `ACCESSION:c.<pos><true_ref>><alt>` |

Three position kinds are each read from the appropriate source:

- **Coding** (`c.462`) and **UTR** (`c.-45`, `c.*12`) — both exonic, read directly out of the transcript's spliced cDNA sequence (fetched once per transcript, cached, then indexed — UTR positions are offset from `cds_start_i`/`cds_end_i`).
- **Intronic** (`c.359-4`, `c.-45+3`, `c.*12-5`) — introns aren't in the spliced transcript sequence, so the position is mapped to a genomic coordinate via `hgvs`'s `VariantMapper.c_to_g` (which — usefully — doesn't validate the reference for intronic positions, unlike exonic ones, so it tolerates a mismatched input ref), then read from the genomic reference sequence and complemented back to the transcript strand if needed. The genomic sequence is fetched once per transcript (the full exon-spanning window, not per position) and cached, to avoid a network round trip per row.

Only single-nucleotide substitutions are checked; ranges, `c.=` cells, and indels are left alone.

The UTR path has no coverage in the real lab data on hand — none of it has a UTR-region variant — so `data/CARD11_example.tsv` includes a fabricated (but UTA-verified) UTR trio specifically to exercise that code path for real; see [About the example data](#about-the-example-data).

### Multi-variant haplotypes

Correctly bracketed multi-variant HGVS (`ACCESSION:c.[2523A>C;2533A>T]`) is unpacked and each component checked independently:

- A component that turns out to be a true no-change is dropped from the variant list (it's not really a co-occurring variant).
- If real variants remain, they're re-bracketed (or left bracket-free if only one remains).
- If a component collapses but the position matters (e.g. avoiding a formatting collision, see below), the position-qualified `<pos>=` form is kept rather than an ambiguous bare `c.=`.

Real example from the lab data (also in `data/CARD11_example.tsv`):

```
NM_032415.7:c.[2523A>C;2533A>T]  →  NM_032415.7:c.2533A>T
```

(`2523A>C` restores the true reference and is dropped; `2533A>T` is real and the brackets come off since only one component remains.)

### Malformed haplotype syntax

Some source rows encode a haplotype without the bracket syntax — `ACCESSION:c.2523A>C;2533A>T` instead of `ACCESSION:c.[2523A>C;2533A>T]`. `translate_hgvs_accessions.py`'s HGVS-splitting logic normalizes this unconditionally (there's only one correct interpretation) for coding, genomic, and protein columns alike, before any other processing happens.

### Why `c.<pos>=`, not bare `c.=`

An earlier version of this script collapsed a true no-change substitution to a bare `ACCESSION:c.=`, discarding the position. That's wrong for two reasons: it loses information, and — more concretely — it causes **collisions**: any two different positions that both happen to be fixed-edit artifacts collapse to the exact same string and look like duplicate rows. Keeping the position (`c.462=`) makes every no-change result unique to its position — which is also how `dedupe_mavedb_variants.py` (below) tells a real collision from a formatting artifact.

### Reporting

The pipeline prints and writes (to `<data-dir>/`):

| File | Contents |
|---|---|
| `correction_report_variants.tsv` | Corrected-variant counts per dataset, by category (coding/intronic/UTR) |
| `correction_report_positions.tsv` | Distinct corrected-position counts per dataset, by category |
| `correction_details.tsv` | `dataset`, `original`, `corrected` — one row per corrected cell (a haplotype's bracketed group counts as one row, since the whole cell is what changed) |

---

## MAVE-HGVS `inv` incompatibility

MAVE-HGVS doesn't support HGVS inversion notation (`c.2687_2688inv`). `inversions_to_delins.py` rewrites it to the equivalent `delins`, by fetching the true reference from UTA and reverse-complementing it:

```
NM_032415.7:c.2687_2688inv  →  NM_032415.7:c.2687_2688delinsAA
```

(True reference at c.2687–2688 is `TT`; reverse complement is `AA`.) This runs only in the MaveDB branch — IGVF keeps standard `inv` notation.

---

## MaveDB deduplication

MaveDB requires exactly one row per variant. Three distinct things can produce duplicate `coding variant` values by the time the MaveDB branch reaches this point:

1. **Haplotype collapse** — a haplotype's real component(s) get dropped during correction and it reduces to a plain variant that already exists as its own row (the `c.2533A>T` example above).
2. **`c.=` collisions** — mitigated by the position-qualified `c.<pos>=` format, but not eliminated: if a *haplotype* fully collapses to no real components, or two genuinely identical positions collide, this can still happen.
3. **Pre-existing source duplicates**, unrelated to correction — e.g. one real lab dataset has two overlapping tiled sequencing windows (different `region` IDs) that both cover the same coding positions, so the same nominal variant is reported twice in the *raw* source file, each with a different local `reference`. (Extracting those two windows' reference sequences at the shared positions showed them to genuinely differ. Cross-checking both against the true transcript sequence showed one window matched and the other didn't.)

`dedupe_mavedb_variants.py` resolves each duplicate group by ranking rows into tiers and keeping the best (ties broken by first occurrence in file order):

1. Not haplotype-derived **and** did not need reference correction (i.e. its reported reference already matched the truth).
2. Not haplotype-derived (regardless of correction status).
3. Any row (all candidates were haplotype-derived).

Both "was a haplotype" and "needed correction" are determined by comparing each row against a row-order-aligned pre-correction snapshot of the same file (the pipeline passes the step-2 intermediate through for this purpose — it isn't deleted until after dedup runs).

Applied to the overlapping-window case described above: at the positions where one window needed correction and the other didn't, the uncorrected window's row wins (tier 1 beats tier 2) — matching the biologically correct choice. At positions where *neither* window needed correction (both tier 1), the tie falls through to first-occurrence-in-file-order — not tied to which window it came from.

Dropped rows are written to `<stem>.mavedb.duplicates_removed.tsv` (all original columns preserved) rather than silently discarded.

---

## About the example data

`CARD11_example.tsv` is **not synthetic biology** — 14 of its 17 rows are copied verbatim (sequences, counts, scores, everything) from real CARD11 MAVE data, specifically chosen because each one hits a different special case described above. The remaining 3 rows (the 3′ UTR trio) are fabricated, since no real dataset on hand has a UTR-region variant — but they're not arbitrary: the reference window was built from the *true* `NM_032415.7` sequence (fetched from UTA) with one deliberately wrong base, the same shape as the real fixed-edit rows, so the correction math runs for real rather than being hand-waved.

All coding-variant accessions use the underscore-versioned form (`NM_032415_7`, not `NM_032415.7`) — the same as the real source spreadsheets — so version normalization is exercised on every single row, not called out separately below.

### Special characteristics of this example

This file is built to teach the pipeline, not to look like a typical input, and a few things about it don't generalize:

- **It has `seq`/`reference` columns.** As noted in [Source data shape](#source-data-shape), not every real dataset does. This example needed them, because two of the special cases it demonstrates (steps 0 and 1) are reconstructed *from* those columns.
- **It's unrealistically dense with edge cases.** Every one of the 17 rows hits a distinct special case (or is a deliberate baseline for contrast). In real lab data these are rare: across five real CARD11/IRF4/BCL10 datasets this pipeline has processed, well under 10% of rows needed any reference correction at all, inversions appeared exactly once per affected file, and malformed haplotypes only in two of the five files (six rows each) — not the one-in-every-few-rows density this example uses for teaching purposes.
- **It's one gene, one strand.** Everything here is `NM_032415.7` (CARD11), which is on the minus strand — the transcript-to-genomic complementing this exercises always runs in that direction. A plus-strand gene isn't represented in this file (though the pipeline's own strand handling is generic, not CARD11-specific).
- **The 3′ UTR trio (`region` 901) is fabricated, not observed.** No real dataset on hand has a UTR variant, so these three rows — sequence, counts, and scores alike — were constructed to exercise that code path rather than copied from an experiment. `region` 901 itself is an arbitrary placeholder, not a real tiling-window ID; it was picked simply to avoid colliding with CARD11's real region numbers.
- **It's small.** 17 rows, versus thousands in a real MAVE/SGE dataset — chosen so the whole file, and this table, stay readable.

| Rows | Region | Special case | What happens |
|---|---|---|---|
| `c.=` | 259 | **Unchanged/WT control expansion** | No position on its own; `expand_unchanged_hgvs_c.py` infers it from the sibling `c.460`/`c.462` rows sharing the same `reference` window → `c.460_495=` |
| `c.460G>A` | 259 | *(baseline, unaffected)* | A normal missense call — included for contrast; nothing about it needs fixing |
| `c.462A>C`, `c.462A>G`, `c.462A>T` | 259 | **Fixed-edit reference correction (coding)** | The assay's local reference reports `A` at c.462; the true `NM_032415.7` reference is `C`. `A>C` is really no change (`c.462=`); `A>G`/`A>T` get their reference letter corrected (`c.462C>G`, `c.462C>T`) |
| `c.1017+4A>C`, `c.1017+4A>G`, `c.1017+4A>T` | 273 | **Fixed-edit reference correction (intronic)** | Same idea, but the position is intronic — the true reference (`G`) is read via UTA's genomic mapping, not the transcript's own spliced sequence. `A>G` collapses to `c.1017+4=`; the others become `c.1017+4G>C` / `c.1017+4G>T` |
| `c.*8A>C`, `c.*8A>G`, `c.*8A>T` | 901 (synthetic) | **Fixed-edit reference correction (3′ UTR)** | Same pattern again, this time in the 3′ UTR (true reference `G`, reported `A`). `A>G` → `c.*8=`; the others → `c.*8G>C` / `c.*8G>T`. This is the only place this code path has run against real pipeline output so far — no lab dataset on hand has a UTR variant |
| `c.2533A>T` | 453 | *(baseline, unaffected — also a dedup target)* | A normal nonsense call, kept as-is. It's also the row that a later collision lands on (see below) |
| `c.2523A>C;2533A>T` | 453 | **Malformed multi-variant syntax + haplotype correction + dedup** | Missing the HGVS bracket syntax for a two-variant haplotype; `translate_hgvs_accessions.py` normalizes it to `c.[2523A>C;2533A>T]`. `correct_reference_alleles.py` then finds `2523A>C` is itself a fixed-edit artifact (true ref `C`) and drops it, leaving just `c.2533A>T` — which collides with the row above. `dedupe_mavedb_variants.py` keeps the original plain row and drops this one (it "was a haplotype"), writing it to `CARD11_example.mavedb.duplicates_removed.tsv` |
| `c.2533A>T;2535A>G` | 453 | **Malformed multi-variant syntax, no correction needed** | Also missing brackets, normalized to `c.[2533A>T;2535A>G]` — but neither component needs correction, so it passes through as a well-formed two-variant haplotype, untouched by dedup (its value is unique) |
| `c.2687_2688inv` | 458 | **`inv` → `delins` (MaveDB only)** | MAVE-HGVS doesn't support inversion notation. `inversions_to_delins.py` fetches the true reference (`TT`) and reverse-complements it: `c.2687_2688delinsAA`. Only in the MaveDB output — the IGVF output keeps standard `inv` notation |
| `c.865-4C>A` | 300 | *(baseline, unaffected)* | A normal, well-formed intronic call — included for contrast with the row below |
| `c.865-4>G` | 300 | **Missing reference-letter typo** | A manual-entry typo — missing the reference letter entirely. `fix_missing_ref_allele.py` reconstructs it from this row's own `seq`/`reference` columns (which differ at exactly one position) as `c.865-4C>G`, *before* anything else in the pipeline touches it |

## What comes out

Running the pipeline against `CARD11_example.tsv` produces:

- `CARD11_example.igvf.tsv` — 17 rows (all of them; IGVF doesn't require uniqueness), Ensembl-translated (`ENST00000396946.9`, `ENSP00000380150.4`)
- `CARD11_example.mavedb.csv` — 16 rows (one dropped by dedup), RefSeq (`NM_032415.7`)
- `CARD11_example.mavedb.scores.csv` / `.counts.csv` / `.functional_classes.csv` — the submission-ready split, 16 rows each
- `CARD11_example.mavedb.duplicates_removed.tsv` — the 1 row dropped by dedup, for audit
- `correction_report_variants.tsv` / `correction_report_positions.tsv` — for this dataset: **10 variants corrected at 4 positions** (4/2 coding, 3/1 intronic, 3/1 UTR)
- `correction_details.tsv` — all 10 corrections, original → corrected, e.g.:

  ```
  NM_032415.7:c.462A>C            -> NM_032415.7:c.462=
  NM_032415.7:c.1017+4A>G         -> NM_032415.7:c.1017+4=
  NM_032415.7:c.*8A>G             -> NM_032415.7:c.*8=
  NM_032415.7:c.[2523A>C;2533A>T] -> NM_032415.7:c.2533A>T
  ```

Every intermediate file the pipeline creates along the way (`*.step*.csv`, `*.mavedb.stage*.csv`, `*.mavedb.raw.csv`, `*.correction_summary.txt`) is cleaned up automatically; only the files listed above should remain after a run.

## Known limitations / possible follow-ups

- **Only single-nucleotide substitutions are checked** for reference correction and appear in the correction reports. Indel components — including those inside a haplotype (`c.2559_2560insT`, `c.2562del`) — are passed through untouched, even if their reference context is also affected by a fixed edit.
- **`variant_type`/`Classification` labels aren't reconciled** after a substitution collapses to a true no-change. A row originally labeled `synonymous` that becomes `c.462=` still carries that label, even though it's now more like the explicit wild-type control rows than a variant.
- **The "first occurrence in file order" fallback is arbitrary** wherever it applies (tier 2/3 duplicates, or tier-1 ties). It's deterministic and auditable (via `duplicates_removed.tsv`), but not chosen on any biological basis beyond what's specified above.
- **Pre-existing source duplicates are only caught incidentally.** `dedupe_mavedb_variants.py` doesn't know anything about tiled windows or `region` — it only ever compares the final `hgvs_nt` text. Two overlapping windows describing the *same* position get caught because they converge on identical text after correction, as in the case above. But if two windows describing the same position, for some reason, didn't converge on identical post-correction text, they'd silently coexist as two separate, uncaught rows — there's nothing in this pipeline that would notice or flag them for review.
