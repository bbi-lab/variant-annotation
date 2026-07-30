# annotate_mavedb

`src/annotate_mavedb.py` — Step 10 of the variant-annotation pipeline (optional).

Annotates each variant row with functional classification labels from MaveDB calibrations. For each variant, the script fetches the calibrations associated with its score set from the MaveDB API, then classifies the variant under two selected calibrations — the **primary** calibration and the **investigator-provided** calibration.

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

---

## Calibration selection

A MaveDB score set may have multiple calibrations. The script selects at most two:

### Primary calibration

The calibration with `primary == true`. If no calibration is explicitly marked primary, the script falls back in order to:
1. The investigator-provided calibration (if one exists)
2. The first non-research-use-only calibration in the list

### Investigator-provided calibration

The calibration with `investigatorProvided == true`. Populated independently of the primary calibration; may be the same object.

Other calibrations (e.g. alternative thresholds, research-use-only calibrations) are fetched but not applied.

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

All data is fetched live from the MaveDB API. Two in-memory caches are maintained per run:

- **Calibration list** — keyed by score set URN. Fetched once per unique score set regardless of how many rows reference it.
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

---

## CLI options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--mavedb-api-url URL` | `MAVEDB_API_URL` | `https://api.mavedb.org` | MaveDB REST API base URL |
| `--variant-urn-col COL` | — | `variant_urn` | Input column containing MaveDB variant URNs |
| `--score-col COL` | — | `score` | Input column containing the numeric variant score |
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
