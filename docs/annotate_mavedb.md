# annotate_mavedb

`src/annotate_mavedb.py` — Step 10 of the variant-annotation pipeline (optional).

Annotates each variant row with functional classification labels from MaveDB calibrations. For each variant, the script fetches the calibrations associated with its score set from the MaveDB API, then classifies the variant under two selected calibrations — the **primary** calibration and the **investigator-provided** calibration. An optional third calibration — the **requested** calibration — can be added via `--requested-calibrations-file`.

---

## Output columns

| Column | Description |
|---|---|
| `mavedb.primary_calibration.urn` | URN of the primary calibration for the score set |
| `mavedb.primary_calibration.name` | Title of the primary calibration |
| `mavedb.primary_calibration.url` | URL to the score set page on MaveDB |
| `mavedb.primary_calibration.functional_class_label` | Functional classification label for this variant under the primary calibration |
| `mavedb.primary_calibration.functional_classification` | The `functionalClassification` value (`normal`, `abnormal`, or `not_specified`) reported by the MaveDB API for the matching classification under the primary calibration |
| `mavedb.investigator_provided_calibration.urn` | URN of the investigator-provided calibration |
| `mavedb.investigator_provided_calibration.name` | Title of the investigator-provided calibration |
| `mavedb.investigator_provided_calibration.url` | URL to the score set page on MaveDB |
| `mavedb.investigator_provided_calibration.functional_class_label` | Functional classification label under the investigator-provided calibration |
| `mavedb.investigator_provided_calibration.functional_classification` | The `functionalClassification` value reported by the MaveDB API for the matching classification under the investigator-provided calibration |

All ten columns are always written, defaulting to empty string when no applicable calibration exists or the variant is unclassified.

When the primary and investigator-provided calibrations are the same object (which is common when an expert panel both provides and designates a single calibration), both column groups will contain identical values.

### Requested calibration columns (optional)

When `--requested-calibrations-file` is given, five more columns are written:

| Column | Description |
|---|---|
| `mavedb.requested_calibration.urn` | URN of the requested calibration for this score set, as given in the requested-calibrations file |
| `mavedb.requested_calibration.name` | Title of the requested calibration |
| `mavedb.requested_calibration.url` | URL to the score set page on MaveDB |
| `mavedb.requested_calibration.functional_class_label` | Functional classification label for this variant under the requested calibration |
| `mavedb.requested_calibration.functional_classification` | The `functionalClassification` value reported by the MaveDB API for the matching classification under the requested calibration |

These columns are left blank for score sets that are absent from the requested-calibrations file, that have a blank `requested_calibration_urn`, or whose `requested_calibration_urn` is not among the calibrations returned for that score set (e.g. a typo, or a URN that belongs to a different score set). **If `--requested-calibrations-file` is not given, none of these five columns appear in the output at all.**

---

## Calibration selection

A MaveDB score set may have multiple calibrations. The script selects at most two automatically, plus an optional third supplied explicitly:

### Primary calibration

The calibration with `primary == true`. If no calibration is explicitly marked primary, the script falls back in order to:
1. The investigator-provided calibration (if one exists)
2. The first non-research-use-only calibration in the list

### Investigator-provided calibration

The calibration with `investigatorProvided == true`. Populated independently of the primary calibration; may be the same object.

### Requested calibration

The requested calibration is selected the same way as the primary and investigator-provided calibrations — by looking it up (by URN) within the same calibration list already fetched for the score set (`GET /api/v1/score-calibrations/score-set/{score_set_urn}`), using the `requested_calibration_urn` given for that score set in the requested-calibrations file. It therefore need not be the primary or investigator-provided calibration, and can be any calibration on the score set (including a research-use-only one) — but it must actually belong to that score set.

Because the lookup is restricted to the score set's own calibration list, a `requested_calibration_urn` that doesn't match any calibration for that score set (a typo, a deleted calibration, or a URN that belongs to a *different* score set) simply isn't found: no extra API call is made to fetch it, a warning is logged, and the `mavedb.requested_calibration.*` columns are left blank for that row.

Other calibrations (e.g. alternative thresholds, research-use-only calibrations) are fetched but not applied unless requested.

---

## Requested-calibrations file

Passed via `--requested-calibrations-file PATH`. A CSV/TSV file with two required columns:

| Column | Description |
|---|---|
| `score_set_urn` | MaveDB score set URN |
| `requested_calibration_urn` | URN of the specific calibration to fetch and classify for that score set |

A convenient starting point is a score-set list like `data/cvfg/score_sets.tsv` with a `requested_calibration_urn` column added. Other columns (e.g. `dataset_name`) are ignored. Rows with a blank `score_set_urn` or `requested_calibration_urn` are treated as "no requested calibration" for that score set — they do not need to be removed from the file.

Looking up the requested calibration makes no additional API calls: it reuses the calibration list already fetched (and cached) for the primary/investigator-provided calibrations.

---

## Classification strategies

Each calibration uses one of two classification strategies, detected from its structure:

### Range-based

The calibration's `functionalClassifications` entries each have a `range` field. The variant's numeric score is evaluated against each range in order; the first matching range's `label` and `functionalClassification` are returned.

Bounds are evaluated per-entry using the `inclusiveLowerBound` and `inclusiveUpperBound` flags (defaults: inclusive lower, exclusive upper). Unbounded ends are represented as `null` (treated as ±∞).

The numeric score is taken from `--score-col`; a missing or non-numeric score results in empty `functional_class_label` and `functional_classification` values.

### Class-based

The calibration's `functionalClassifications` entries have no `range` field. The API endpoint `GET /api/v1/score-calibrations/{urn}/variants` is called to retrieve the full variant-to-class mapping for that calibration. The variant is looked up by its URN; if not present, `functional_class_label` and `functional_classification` are left empty.

Class-based variant mappings are cached in memory per calibration URN for the duration of the run.

---

## Caching

All data is fetched live from the MaveDB API. In-memory caches are maintained per run:

- **Calibration list** — keyed by score set URN. Fetched once per unique score set regardless of how many rows reference it. The requested calibration (when `--requested-calibrations-file` is given) is looked up within this same cached list, so it adds no extra API calls.
- **Class-based variant assignments** — keyed by calibration URN. Fetched once per class-based calibration encountered.

There is no persistent on-disk cache. Large runs that span many score sets will make one API call per unique score set at startup.

---

## Usage

```bash
src/scripts/run_annotate_mavedb.sh input.tsv output.tsv
```

With custom column names:
```bash
src/scripts/run_annotate_mavedb.sh input.tsv output.tsv \
  --variant-urn-col variant_urn \
  --score-col score
```

With a requested calibration per score set:
```bash
src/scripts/run_annotate_mavedb.sh input.tsv output.tsv \
  --requested-calibrations-file data/cvfg/score_sets_with_requested_calibrations.tsv
```

---

## CLI options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--mavedb-api-url URL` | `MAVEDB_API_URL` | `https://api.mavedb.org` | MaveDB REST API base URL |
| `--variant-urn-col COL` | — | `variant_urn` | Input column containing MaveDB variant URNs |
| `--score-col COL` | — | `score` | Input column containing the numeric variant score |
| `--requested-calibrations-file PATH` | — | none | CSV/TSV file with `score_set_urn` and `requested_calibration_urn` columns; adds the `mavedb.requested_calibration.*` output columns, looked up by URN within each score set's own calibration list |
| `--skip N` | — | `0` | Skip first N data rows |
| `--limit N` | — | no limit | Stop after N rows |
| `--log-level` | — | `INFO` | Logging verbosity |
| `--csv-field-size-limit BYTES` | — | system default | Increase for large fields |

The `url` output columns use `MAVEDB_FRONTEND_URL` (default: `https://mavedb.org`) as the base URL for score-set page links.

The input delimiter is auto-detected from the file extension (`.tsv`/`.txt` → tab; otherwise comma). Output delimiter matches.

---

## Dependencies

- **`requests`** — included in project dependencies. Used for all MaveDB API calls (no authentication required for public score sets).

---

## Troubleshooting

**All `functional_class_label` / `functional_classification` columns are empty**

- Confirm `--variant-urn-col` points to a column containing valid MaveDB variant URNs (`urn:mavedb:...<score-set-urn>#<id>` format). The script extracts the score set URN as the portion before the final `#`.
- The score set may have no calibrations yet, or only research-use-only calibrations. Check the score set page on `https://mavedb.org`.

**Range-based classification returns empty for a scored variant**

- The variant's score may fall outside all defined ranges, or between ranges with a gap. This is a genuine "unclassified" result per the calibration's thresholds.
- Confirm `--score-col` points to the correct numeric score column.

**Warnings about failed calibration fetch**

- Transient network errors are logged as warnings; the row is written with empty annotation columns. Retry the run using `--skip` to resume from the failed position, or re-run the full file (the calibration cache resets each run).

**`mavedb.requested_calibration.*` columns are empty for a score set that has a requested calibration**

- Check the logs for a "not found among calibrations for score set ..." warning. This means `requested_calibration_urn` doesn't match any calibration returned for that score set — either it's mistyped/deleted, or it actually belongs to a *different* score set. Fix the URN in the requested-calibrations file.
- Confirm `score_set_urn` in the requested-calibrations file exactly matches the score set URN embedded in `variant_urn` (the portion before the final `#`).
- The score set may have no calibrations yet, or only research-use-only calibrations that were nonetheless excluded — check the score set page on `https://mavedb.org` to confirm the calibration URN is listed there.
