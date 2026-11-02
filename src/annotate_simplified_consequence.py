"""Add simplified_consequence, condensed_consequence, splice_variant, and splice_var_amino columns.

``vep.mutational_consequences`` (see ``src/annotate_vep.py``) is a
``|``-delimited list of per-transcript consequence sets, each itself a
``^``-delimited list of raw VEP terms, e.g.::

    missense_variant|missense_variant^splice_region_variant

For each ``|``-delimited set, its ``^``-delimited terms are looked up (in
order) in a VEP-output-term -> SO-summary-term table
(--consequence-map-file, e.g. data/cvfg/extended_ensembl_consequence.csv.gz)
and the summary term for the first matching VEP term is kept — mirroring the
priority-ordered ``get_simplified_consequence`` mapping used elsewhere in the
CVFG pipeline, except the inner split is on ``^`` (per-transcript VEP terms)
rather than ``,``. The per-set results are re-joined with ``|`` to produce
``simplified_consequence``.

``condensed_consequence`` collapses ``simplified_consequence`` to a single
value: that value if all of its pipe-delimited entries are identical, else
"conflict_consequence".

``splice_variant`` is derived, position-by-position, from the pipe-aligned
SpliceAI delta-score columns (``spliceai.ds_ag``/``ds_al``/``ds_dg``/``ds_dl``,
see ``src/annotate_spliceai.py``) and ``simplified_consequence``: a position is
"Yes" if any of its four delta scores is >= 0.2 (--splice-score-threshold) or
its ``simplified_consequence`` entry is "splice_site_variant", blank if it has
no delta scores and no ``simplified_consequence`` entry, else "No".

``splice_var_amino`` collapses ``splice_variant`` to a single value: blank if
no non-blank entries, that value if all non-blank entries agree, else "Yes".

Usage::

    python -m src.annotate_simplified_consequence input.tsv output.tsv \\
        --consequence-map-file data/cvfg/extended_ensembl_consequence.csv.gz
"""

from __future__ import annotations

import argparse
import csv
import logging
import time
from itertools import islice
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# Terms not present in the VEP output term table but already at the desired
# summary-term granularity, so they map to themselves.
EXTRA_TERMS = ["UTR_variant", "splicing_variant", "splice_site_variant", "no_change"]

# condensed_consequence value used when simplified_consequence's pipe-delimited
# entries disagree.
CONFLICT_CONSEQUENCE = "conflict_consequence"


def build_mapping_dict(consequence_map_file: str) -> dict[str, object]:
    """Build a VEP-output-term -> SO-summary-term lookup, prioritized by Importance."""
    consequence_df = pd.read_csv(consequence_map_file)
    mapping_dict = (
        consequence_df.sort_values(by="Importance", ascending=False)
        .groupby("VEP output term")[["SO summary term", "Importance"]]
        .first()
        .to_dict()["SO summary term"]
    )
    mapping_dict.update({term: term for term in EXTRA_TERMS})
    return mapping_dict


def get_simplified_consequence(terms: list[str], mapping_dict: dict[str, object]) -> Optional[str]:
    """Return the mapped summary term for the first *terms* entry found in mapping_dict, else None."""
    matched_terms = [
        mapping_dict[t] for t in (raw.strip() for raw in terms) if t in mapping_dict
    ]
    return matched_terms[0] if matched_terms else None


def _term_to_str(value: object) -> str:
    """Render a mapping_dict value (str, None, or NaN float) as a TSV-safe string."""
    if value is None:
        return ""
    if isinstance(value, float) and value != value:  # NaN
        return ""
    return str(value)


def simplify_mutational_consequences(value: str, mapping_dict: dict[str, object]) -> str:
    """Map a full ``|``/``^``-delimited vep.mutational_consequences value to simplified_consequence."""
    if not value:
        return ""
    simplified = [
        _term_to_str(get_simplified_consequence(part.split("^"), mapping_dict))
        for part in value.split("|")
    ]
    return "|".join(simplified)


def compute_condensed_consequence(simplified_consequence: str) -> str:
    """Collapse a pipe-delimited simplified_consequence value to a single value.

    "" if blank; that value if all pipe-delimited entries are identical; else
    CONFLICT_CONSEQUENCE.
    """
    if not simplified_consequence:
        return ""
    values = simplified_consequence.split("|")
    unique = set(values)
    if len(unique) == 1:
        return next(iter(unique))
    return CONFLICT_CONSEQUENCE


def _split_pipe(value: str) -> list[str]:
    """Split a pipe-delimited value, treating a blank value as zero entries."""
    return value.split("|") if value else []


def _meets_splice_threshold(value: str, threshold: float) -> bool:
    """Return True if *value* parses as a float >= threshold."""
    if not value:
        return False
    try:
        return float(value) >= threshold
    except ValueError:
        return False


def compute_splice_variant(
    ds_ag: str,
    ds_al: str,
    ds_dg: str,
    ds_dl: str,
    simplified_consequence: str,
    threshold: float = 0.2,
) -> str:
    """Return a pipe-delimited Yes/No column aligned to the SpliceAI/consequence positions.

    A position is "Yes" if any of its four delta scores is >= threshold, or its
    simplified_consequence entry is "splice_site_variant". It is blank if it has
    no delta scores and no simplified_consequence entry, since there is no
    evidence either way; otherwise it is "No".
    Returns "" if none of the inputs have any positions.
    """
    ag_list = _split_pipe(ds_ag)
    al_list = _split_pipe(ds_al)
    dg_list = _split_pipe(ds_dg)
    dl_list = _split_pipe(ds_dl)
    cons_list = _split_pipe(simplified_consequence)

    n = max(len(ag_list), len(al_list), len(dg_list), len(dl_list), len(cons_list))
    if n == 0:
        return ""

    results = []
    for i in range(n):
        scores = (
            ag_list[i] if i < len(ag_list) else "",
            al_list[i] if i < len(al_list) else "",
            dg_list[i] if i < len(dg_list) else "",
            dl_list[i] if i < len(dl_list) else "",
        )
        consequence = cons_list[i] if i < len(cons_list) else ""
        if not consequence and not any(scores):
            results.append("")
            continue
        is_splice = (
            any(_meets_splice_threshold(score, threshold) for score in scores)
            or consequence == "splice_site_variant"
        )
        results.append("Yes" if is_splice else "No")
    return "|".join(results)


def compute_splice_var_amino(splice_variant: str) -> str:
    """Collapse a pipe-delimited splice_variant value to a single value.

    Blank entries are dropped; "" if none remain, that value if all remaining
    entries agree, else "Yes".
    """
    values = {v for v in splice_variant.split("|") if v} if splice_variant else set()
    if not values:
        return ""
    if len(values) == 1:
        return next(iter(values))
    return "Yes"


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Add a simplified_consequence column derived from a pipe/caret-delimited "
            "VEP mutational-consequences column."
        )
    )
    p.add_argument("input_file", help="Input CSV/TSV file")
    p.add_argument("output_file", help="Output CSV/TSV file")
    p.add_argument(
        "--consequence-map-file",
        default="data/cvfg/extended_ensembl_consequence.csv.gz",
        metavar="PATH",
        help=(
            "CSV(.gz) with 'VEP output term', 'SO summary term', and 'Importance' "
            "columns (default: %(default)s)"
        ),
    )
    p.add_argument(
        "--consequences-col",
        default="vep.mutational_consequences",
        help="Input column with pipe/caret-delimited VEP terms (default: %(default)s)",
    )
    p.add_argument(
        "--output-col",
        default="simplified_consequence",
        help="Output column to write (default: %(default)s)",
    )
    p.add_argument(
        "--condensed-consequence-col",
        default="condensed_consequence",
        help="Output column for the collapsed single-value consequence (default: %(default)s)",
    )
    p.add_argument(
        "--spliceai-ds-ag-col",
        default="spliceai.ds_ag",
        help="Input column with pipe-delimited SpliceAI DS_AG scores (default: %(default)s)",
    )
    p.add_argument(
        "--spliceai-ds-al-col",
        default="spliceai.ds_al",
        help="Input column with pipe-delimited SpliceAI DS_AL scores (default: %(default)s)",
    )
    p.add_argument(
        "--spliceai-ds-dg-col",
        default="spliceai.ds_dg",
        help="Input column with pipe-delimited SpliceAI DS_DG scores (default: %(default)s)",
    )
    p.add_argument(
        "--spliceai-ds-dl-col",
        default="spliceai.ds_dl",
        help="Input column with pipe-delimited SpliceAI DS_DL scores (default: %(default)s)",
    )
    p.add_argument(
        "--splice-score-threshold",
        type=float,
        default=0.2,
        help="Minimum SpliceAI delta score to count as a splice variant (default: %(default)s)",
    )
    p.add_argument(
        "--splice-variant-col",
        default="splice_variant",
        help="Output column for the pipe-delimited splice-variant Yes/No flags (default: %(default)s)",
    )
    p.add_argument(
        "--splice-var-amino-col",
        default="splice_var_amino",
        help="Output column for the collapsed single-value splice-variant flag (default: %(default)s)",
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

    map_path = Path(args.consequence_map_file)
    if not map_path.exists():
        logger.error("Consequence map file not found: %s", map_path)
        raise SystemExit(1)
    mapping_dict = build_mapping_dict(str(map_path))
    logger.info("Loaded %d VEP output term mappings from %s", len(mapping_dict), map_path)

    input_path = Path(args.input_file)
    output_path = Path(args.output_file)
    delim = "\t" if input_path.suffix.lower() in (".tsv", ".txt") else ","

    processed = 0
    simplified_count = 0
    condensed_count = 0
    splice_variant_count = 0
    splice_var_amino_count = 0
    started = time.monotonic()

    with input_path.open("r", encoding="utf-8", newline="") as in_fh, \
         output_path.open("w", encoding="utf-8", newline="") as out_fh:

        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears empty: %s", input_path)
            raise SystemExit(1)

        fieldnames = list(reader.fieldnames)
        new_cols = [
            args.output_col,
            args.condensed_consequence_col,
            args.splice_variant_col,
            args.splice_var_amino_col,
        ]
        out_fieldnames = fieldnames + [c for c in new_cols if c not in fieldnames]
        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=delim,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        selected = islice(
            reader,
            args.skip,
            None if args.limit is None else args.skip + args.limit,
        )

        for row in selected:
            simplified = simplify_mutational_consequences(row.get(args.consequences_col, ""), mapping_dict)
            row[args.output_col] = simplified

            condensed = compute_condensed_consequence(simplified)
            row[args.condensed_consequence_col] = condensed

            splice_variant = compute_splice_variant(
                row.get(args.spliceai_ds_ag_col, ""),
                row.get(args.spliceai_ds_al_col, ""),
                row.get(args.spliceai_ds_dg_col, ""),
                row.get(args.spliceai_ds_dl_col, ""),
                simplified,
                threshold=args.splice_score_threshold,
            )
            row[args.splice_variant_col] = splice_variant

            splice_var_amino = compute_splice_var_amino(splice_variant)
            row[args.splice_var_amino_col] = splice_var_amino

            writer.writerow(row)

            processed += 1
            if simplified:
                simplified_count += 1
            if condensed:
                condensed_count += 1
            if splice_variant:
                splice_variant_count += 1
            if splice_var_amino:
                splice_var_amino_count += 1

            if processed % 10000 == 0:
                elapsed = max(time.monotonic() - started, 1e-9)
                logger.info(
                    "Progress: %d rows processed (%.1f rows/s)",
                    processed,
                    processed / elapsed,
                )

    elapsed = max(time.monotonic() - started, 1e-9)
    logger.info(
        "Done. %d/%d rows given a non-blank %s, %d a non-blank %s, %d a non-blank %s, "
        "%d a non-blank %s. %d rows written to %s (%.1f rows/s)",
        simplified_count, processed, args.output_col,
        condensed_count, args.condensed_consequence_col,
        splice_variant_count, args.splice_variant_col,
        splice_var_amino_count, args.splice_var_amino_col,
        processed, output_path, processed / elapsed,
    )


if __name__ == "__main__":
    main()
