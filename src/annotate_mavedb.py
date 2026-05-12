"""Annotate variant rows with MaveDB functional classification labels.

For each variant row, the script:

1. Extracts the score set URN from the variant URN (the portion before ``#``).
2. Fetches all readable calibrations for that score set from the MaveDB API:
   ``GET /api/v1/score-calibrations/score-set/{score_set_urn}``
3. Identifies two calibrations of interest (each row may have at most one of
   each type):

   - **Primary calibration** – the calibration whose ``primary`` flag is
     ``True``.
   - **Investigator-provided calibration** – the calibration whose
     ``investigatorProvided`` flag is ``True``.

4. Classifies the variant against each applicable calibration:

   - *Range-based*: finds the functional classification whose numeric range
     contains the variant's score (bounds are evaluated using the per-range
     ``inclusiveLowerBound`` / ``inclusiveUpperBound`` flags).
   - *Class-based*: fetches all variant-to-class assignments for that
     calibration (``GET /api/v1/score-calibrations/{urn}/variants``) and
     looks up the current variant by URN.

5. Writes six annotation columns per row.

Output columns
--------------

  ``mavedb.primary_calibration.urn``
      URN of the primary calibration for the score set (empty if none).

  ``mavedb.primary_calibration.name``
      Title of the primary calibration (empty if none).

  ``mavedb.primary_calibration.functional_class``
      Label of the functional classification that matches the variant's score
      under the primary calibration (empty when unclassified or no calibration).

  ``mavedb.investigator_provided_calibration.urn``
      URN of the investigator-provided calibration (empty if none).

  ``mavedb.investigator_provided_calibration.name``
      Title of the investigator-provided calibration (empty if none).

  ``mavedb.investigator_provided_calibration.functional_class``
      Label of the matching classification under the investigator-provided
      calibration (empty when unclassified or no calibration).

When no calibration is explicitly marked ``primary``, the script falls back
first to the investigator-provided calibration, then to the first
non-research-use-only calibration available for the score set.

The ``*.url`` columns point to the score set page on the MaveDB website
(``https://mavedb.org/score-sets/{score_set_urn}``), where calibrations are
displayed.  Set the ``MAVEDB_FRONTEND_URL`` environment variable to override
the base URL (default: ``https://mavedb.org``).

Usage::

    python -m src.annotate_mavedb input.tsv output.tsv [OPTIONS]

Example::

    python -m src.annotate_mavedb variants.tsv annotated.tsv \\
        --mavedb-api-url https://api.mavedb.org \\
        --variant-urn-col variant_urn \\
        --score-col score
"""

from __future__ import annotations

import argparse
import csv
from itertools import islice
import logging
import os
from pathlib import Path
from typing import Any, Optional

import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_MAVEDB_API_URL = "https://api.mavedb.org"
MAVEDB_FRONTEND_URL = os.environ.get("MAVEDB_FRONTEND_URL", "https://mavedb.org")

OUTPUT_COLS = [
    "mavedb.primary_calibration.urn",
    "mavedb.primary_calibration.name",
    "mavedb.primary_calibration.url",
    "mavedb.primary_calibration.functional_class",
    "mavedb.investigator_provided_calibration.urn",
    "mavedb.investigator_provided_calibration.name",
    "mavedb.investigator_provided_calibration.url",
    "mavedb.investigator_provided_calibration.functional_class",
]


# ---------------------------------------------------------------------------
# URN helpers
# ---------------------------------------------------------------------------


def score_set_urn_from_variant_urn(variant_urn: str) -> Optional[str]:
    """Return the score set URN embedded in *variant_urn* (the part before ``#``).

    MaveDB variant URNs have the form ``urn:mavedb:<accession>#<id>``.
    Returns ``None`` if the URN contains no ``#`` separator.
    """
    idx = variant_urn.rfind("#")
    if idx < 0:
        return None
    return variant_urn[:idx]


# ---------------------------------------------------------------------------
# MaveDB API fetchers
# ---------------------------------------------------------------------------


def fetch_calibrations(
    api_url: str,
    score_set_urn: str,
    session: requests.Session,
) -> list[dict[str, Any]]:
    """Return all readable calibrations for *score_set_urn* from the MaveDB API.

    Returns an empty list when the score set is not found (HTTP 404) or has no
    calibrations.  Any other HTTP error is re-raised.
    """
    url = f"{api_url}/api/v1/score-calibrations/score-set/{score_set_urn}"
    resp = session.get(url, timeout=30)
    if resp.status_code == 404:
        return []
    resp.raise_for_status()
    return resp.json()


def fetch_calibration_variant_class_ids(
    api_url: str,
    calibration_urn: str,
    session: requests.Session,
) -> dict[str, int]:
    """Return a mapping of ``{variant_urn: functional_classification_id}`` for a
    class-based calibration.

    Calls ``GET /api/v1/score-calibrations/{urn}/variants``, which returns a
    list of ``{functionalClassificationId: int, variants: [...]}`` objects.
    """
    url = f"{api_url}/api/v1/score-calibrations/{calibration_urn}/variants"
    resp = session.get(url, timeout=60)
    resp.raise_for_status()
    data: list[dict[str, Any]] = resp.json()
    mapping: dict[str, int] = {}
    for fc_entry in data:
        fc_id: Optional[int] = fc_entry.get("functionalClassificationId")
        if fc_id is None:
            continue
        for variant in fc_entry.get("variants", []):
            v_urn: Optional[str] = variant.get("urn")
            if v_urn:
                mapping[v_urn] = fc_id
    return mapping


# ---------------------------------------------------------------------------
# Classification logic
# ---------------------------------------------------------------------------


def classify_score_range(
    score: float,
    functional_classifications: list[dict[str, Any]],
) -> str:
    """Return the label of the first range-based functional classification that
    contains *score*, or an empty string if no range matches.

    Bounds are evaluated using the per-entry ``inclusiveLowerBound`` and
    ``inclusiveUpperBound`` flags (defaulting to ``True`` and ``False``
    respectively, matching the MaveDB server-side defaults).
    """
    for fc in functional_classifications:
        rng = fc.get("range")
        if rng is None:
            continue  # class-based entry – skip

        lower_raw, upper_raw = rng[0], rng[1]
        lo: float = float("-inf") if lower_raw is None else float(lower_raw)
        hi: float = float("inf") if upper_raw is None else float(upper_raw)
        inc_lower: bool = fc.get("inclusiveLowerBound", True)
        inc_upper: bool = fc.get("inclusiveUpperBound", False)

        lower_ok = (score >= lo) if inc_lower else (score > lo)
        upper_ok = (score <= hi) if inc_upper else (score < hi)

        if lower_ok and upper_ok:
            return fc.get("label", "")
    return ""


def classify_variant(
    variant_urn: str,
    score_str: str,
    calibration: dict[str, Any],
    api_url: str,
    session: requests.Session,
    class_id_cache: dict[str, dict[str, int]],
) -> tuple[str, str, str]:
    """Return ``(calibration_urn, calibration_name, functional_class_label)``
    for *variant_urn* under *calibration*.

    For range-based calibrations the classification is computed locally using
    *score_str*.  For class-based calibrations the variant-to-class mapping is
    fetched from the API (cached in *class_id_cache* keyed by calibration URN).
    """
    cal_urn: str = calibration.get("urn", "")
    cal_name: str = calibration.get("title", "")
    fcs: list[dict[str, Any]] = calibration.get("functionalClassifications") or []

    if not fcs:
        return cal_urn, cal_name, ""

    first_fc = fcs[0]
    if first_fc.get("range") is not None:
        # Range-based calibration – classify locally.
        if not score_str:
            return cal_urn, cal_name, ""
        try:
            score = float(score_str)
        except ValueError:
            logger.warning(
                "Non-numeric score %r for variant %s; cannot classify.", score_str, variant_urn
            )
            return cal_urn, cal_name, ""
        label = classify_score_range(score, fcs)
        return cal_urn, cal_name, label
    else:
        # Class-based calibration – look up variant assignment from API.
        if cal_urn not in class_id_cache:
            try:
                class_id_cache[cal_urn] = fetch_calibration_variant_class_ids(
                    api_url, cal_urn, session
                )
            except requests.RequestException as exc:
                logger.warning(
                    "Failed to fetch variant classes for calibration %s: %s", cal_urn, exc
                )
                class_id_cache[cal_urn] = {}

        fc_id = class_id_cache[cal_urn].get(variant_urn)
        if fc_id is None:
            return cal_urn, cal_name, ""
        for fc in fcs:
            if fc.get("id") == fc_id:
                return cal_urn, cal_name, fc.get("label", "")
        return cal_urn, cal_name, ""


# ---------------------------------------------------------------------------
# Per-row annotation
# ---------------------------------------------------------------------------


def annotate_row(
    row: dict[str, str],
    *,
    api_url: str,
    variant_urn_col: str,
    score_col: str,
    session: requests.Session,
    calibration_cache: dict[str, list[dict[str, Any]]],
    class_id_cache: dict[str, dict[str, int]],
) -> dict[str, str]:
    """Return the six MaveDB annotation columns for one input *row*.

    All output values default to empty strings.  Caches calibration lists in
    *calibration_cache* (keyed by score set URN) and class-based variant
    assignments in *class_id_cache* (keyed by calibration URN) so each
    distinct score set or calibration is fetched at most once.
    """
    out: dict[str, str] = {col: "" for col in OUTPUT_COLS}

    variant_urn = row.get(variant_urn_col, "").strip()
    score_str = row.get(score_col, "").strip()

    if not variant_urn:
        return out

    ss_urn = score_set_urn_from_variant_urn(variant_urn)
    if ss_urn is None:
        logger.debug("No score set URN found in variant_urn=%r; skipping row.", variant_urn)
        return out

    if ss_urn not in calibration_cache:
        try:
            calibration_cache[ss_urn] = fetch_calibrations(api_url, ss_urn, session)
        except requests.RequestException as exc:
            logger.warning("Failed to fetch calibrations for %s: %s", ss_urn, exc)
            calibration_cache[ss_urn] = []

    calibrations = calibration_cache[ss_urn]

    primary_cal = next((c for c in calibrations if c.get("primary")), None)
    inv_cal = next((c for c in calibrations if c.get("investigatorProvided")), None)
    if primary_cal is None and inv_cal is not None:
        primary_cal = inv_cal
    if primary_cal is None:
        primary_cal = next(
            (c for c in calibrations if not c.get("researchUseOnly", False)), None
        )

    cal_url = f"{MAVEDB_FRONTEND_URL}/score-sets/{ss_urn}"

    if primary_cal is not None:
        urn, name, fc_label = classify_variant(
            variant_urn, score_str, primary_cal, api_url, session, class_id_cache
        )
        out["mavedb.primary_calibration.urn"] = urn
        out["mavedb.primary_calibration.name"] = name
        out["mavedb.primary_calibration.url"] = cal_url
        out["mavedb.primary_calibration.functional_class"] = fc_label

    if inv_cal is not None:
        urn, name, fc_label = classify_variant(
            variant_urn, score_str, inv_cal, api_url, session, class_id_cache
        )
        out["mavedb.investigator_provided_calibration.urn"] = urn
        out["mavedb.investigator_provided_calibration.name"] = name
        out["mavedb.investigator_provided_calibration.url"] = cal_url
        out["mavedb.investigator_provided_calibration.functional_class"] = fc_label

    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Annotate variant rows with MaveDB functional classification labels."
    )
    p.add_argument("input_file", help="Input CSV/TSV file")
    p.add_argument("output_file", help="Output CSV/TSV file")
    p.add_argument(
        "--mavedb-api-url",
        default=os.environ.get("MAVEDB_API_URL", DEFAULT_MAVEDB_API_URL),
        metavar="URL",
        help=(
            "Base URL for the MaveDB REST API. "
            "Defaults to the MAVEDB_API_URL env var, or https://api.mavedb.org."
        ),
    )
    p.add_argument(
        "--variant-urn-col",
        default="variant_urn",
        metavar="COLUMN",
        help="Input column containing MaveDB variant URNs (default: variant_urn)",
    )
    p.add_argument(
        "--score-col",
        default="score",
        metavar="COLUMN",
        help="Input column containing the numeric variant score (default: score)",
    )
    p.add_argument(
        "--skip",
        type=int,
        default=0,
        help="Number of data rows to skip before annotation (default: 0)",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of data rows to annotate",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO)",
    )
    p.add_argument(
        "--csv-field-size-limit",
        type=int,
        default=csv.field_size_limit(),
        metavar="BYTES",
        help="Maximum per-field character length for CSV/TSV parsing (default: %(default)s)",
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    csv.field_size_limit(args.csv_field_size_limit)

    input_path = Path(args.input_file)
    output_path = Path(args.output_file)
    delim = "\t" if input_path.suffix.lower() in (".tsv", ".txt") else ","
    out_delim = "\t" if output_path.suffix.lower() in (".tsv", ".txt") else ","

    api_url = args.mavedb_api_url.rstrip("/")
    session = requests.Session()
    calibration_cache: dict[str, list[dict[str, Any]]] = {}
    class_id_cache: dict[str, dict[str, int]] = {}

    rows_written = 0
    rows_skipped = 0

    with open(input_path, newline="", encoding="utf-8") as inf:
        reader = csv.DictReader(inf, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears to be empty.")
            raise SystemExit(1)

        output_fieldnames = list(reader.fieldnames) + OUTPUT_COLS
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", newline="", encoding="utf-8") as outf:
            writer = csv.DictWriter(outf, fieldnames=output_fieldnames, delimiter=out_delim)
            writer.writeheader()

            for row in reader:
                if rows_skipped < args.skip:
                    rows_skipped += 1
                    continue
                if args.limit is not None and rows_written >= args.limit:
                    break

                annotations = annotate_row(
                    row,
                    api_url=api_url,
                    variant_urn_col=args.variant_urn_col,
                    score_col=args.score_col,
                    session=session,
                    calibration_cache=calibration_cache,
                    class_id_cache=class_id_cache,
                )
                row.update(annotations)
                writer.writerow(row)
                rows_written += 1

                if rows_written % 500 == 0:
                    logger.info(
                        "Annotated %d rows; %d unique score sets fetched.",
                        rows_written,
                        len(calibration_cache),
                    )

    logger.info(
        "Done. %d rows written; %d unique score sets fetched.",
        rows_written,
        len(calibration_cache),
    )


if __name__ == "__main__":
    main()
