"""Annotate variants with VEP mutational consequences.

This script annotates each input row with a single mutational consequence derived
from one or more pipe-delimited genomic HGVS candidates (default column:
``mapped_hgvs_g``).

For rows with multiple DNA candidates, each candidate is checked independently.
If all resolved candidate consequences agree (treating unresolved candidates as
equal to the shared resolved consequence), that value is emitted in
``vep.mutational_consequence``.

If resolved candidates disagree, ``vep.mutational_consequence`` is left blank and
``vep.error`` records the discrepancy including a pipe-delimited consequence list
aligned to the input candidate positions.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from itertools import islice
import logging
import os
from pathlib import Path
import sys
from typing import Optional, TextIO

import requests  # type: ignore[import-untyped]

logger = logging.getLogger(__name__)

ENSEMBL_API_URL_DEFAULT = os.environ.get("ENSEMBL_API_URL", "https://rest.ensembl.org")

# Ordered from most to least severe (mirrors MaveDB worker logic).
VEP_CONSEQUENCES = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
]


def _split_pipe_preserve_positions(value: str) -> list[str]:
    raw = value or ""
    if "|" not in raw:
        return [raw.strip()]
    return [part.strip() for part in raw.split("|")]


def run_variant_recoder(
    hgvs_strings: list[str],
    *,
    api_url: str,
    timeout_seconds: int,
) -> dict[str, list[str]]:
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    out: dict[str, list[str]] = {}
    try:
        response = requests.post(
            f"{api_url}/variant_recoder/human",
            headers=headers,
            json={"ids": hgvs_strings},
            timeout=timeout_seconds,
        )
    except requests.RequestException as exc:
        logger.error("Variant Recoder request failed: %s", exc)
        return out

    if response.status_code != 200:
        logger.error("Variant Recoder request failed: %s %s", response.status_code, response.text)
        return out

    data = response.json()
    for entry in data:
        for _, variant_data in entry.items():
            # The API returns either a single result dict or a list of result dicts.
            records = variant_data if isinstance(variant_data, list) else [variant_data]
            for record in records:
                input_hgvs = record.get("input")
                if not input_hgvs:
                    continue
                genomic_hgvs_list: list[str] = []
                for genomic_hgvs in record.get("hgvsg") or []:
                    if genomic_hgvs.startswith("NC_"):
                        genomic_hgvs_list.append(genomic_hgvs)
                if genomic_hgvs_list:
                    out.setdefault(input_hgvs, []).extend(genomic_hgvs_list)
    return out


def _vep_lookup_batch(
    hgvs_strings: list[str],
    *,
    api_url: str,
    timeout_seconds: int,
) -> dict[str, Optional[str]]:
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    out: dict[str, Optional[str]] = {}
    try:
        response = requests.post(
            f"{api_url}/vep/human/hgvs",
            headers=headers,
            json={"hgvs_notations": hgvs_strings},
            timeout=timeout_seconds,
        )
    except requests.RequestException as exc:
        logger.error("VEP request failed: %s", exc)
        return out

    if response.status_code != 200:
        logger.error("VEP request failed: %s %s", response.status_code, response.text)
        return out

    for entry in response.json():
        hgvs = entry.get("input")
        if not hgvs:
            continue
        out[hgvs] = entry.get("most_severe_consequence")
    return out


def _run_batches_concurrent(
    hgvs_list: list[str],
    *,
    api_url: str,
    timeout_seconds: int,
    batch_size: int,
    max_workers: int,
) -> dict[str, Optional[str]]:
    """Submit all VEP batches concurrently and merge results."""
    batches = [hgvs_list[i : i + batch_size] for i in range(0, len(hgvs_list), batch_size)]
    result: dict[str, Optional[str]] = {}
    with ThreadPoolExecutor(max_workers=min(max_workers, len(batches))) as pool:
        futures = {
            pool.submit(_vep_lookup_batch, batch, api_url=api_url, timeout_seconds=timeout_seconds): batch
            for batch in batches
        }
        for fut in as_completed(futures):
            result.update(fut.result())
    return result


def get_functional_consequence(
    hgvs_strings: list[str],
    *,
    api_url: str,
    timeout_seconds: int,
    batch_size: int,
    max_workers: int = 1,
) -> dict[str, Optional[str]]:
    """Return HGVS -> most severe VEP consequence with recoder fallback.

    Mirrors MaveDB worker behavior:
    1) VEP lookup for input HGVS (batches run concurrently when max_workers > 1)
    2) For unresolved entries, Variant Recoder to genomic HGVS
    3) VEP lookup for recoded genomic HGVS (also concurrent)
    4) Choose the most severe consequence by fixed priority order
    """
    result: dict[str, Optional[str]] = {}
    if not hgvs_strings:
        return result

    # Preserve caller order while de-duplicating.
    unique_hgvs = list(dict.fromkeys(hgvs_strings))

    result.update(
        _run_batches_concurrent(
            unique_hgvs,
            api_url=api_url,
            timeout_seconds=timeout_seconds,
            batch_size=batch_size,
            max_workers=max_workers,
        )
    )

    missing_hgvs = [h for h in unique_hgvs if h not in result]
    if not missing_hgvs:
        return result

    recoded = run_variant_recoder(missing_hgvs, api_url=api_url, timeout_seconds=timeout_seconds)
    for missing in missing_hgvs:
        if missing not in recoded:
            result[missing] = None

    all_recoded_hgvs: list[str] = []
    for input_hgvs in missing_hgvs:
        all_recoded_hgvs.extend(recoded.get(input_hgvs, []))

    recoded_results: dict[str, Optional[str]] = {}
    if all_recoded_hgvs:
        recoded_results.update(
            _run_batches_concurrent(
                all_recoded_hgvs,
                api_url=api_url,
                timeout_seconds=timeout_seconds,
                batch_size=batch_size,
                max_workers=max_workers,
            )
        )

    for input_hgvs, recoded_hgvs_list in recoded.items():
        consequences = [recoded_results.get(recoded_hgvs) for recoded_hgvs in recoded_hgvs_list]
        chosen: Optional[str] = None
        for consequence in VEP_CONSEQUENCES:
            if consequence in consequences:
                chosen = consequence
                break
        result[input_hgvs] = chosen

    return result


def _get_hgvs_for_row(row: dict[str, str], hgvs_cols: list[str]) -> str:
    """Return the first non-blank HGVS value from the priority-ordered column list."""
    for col in hgvs_cols:
        val = (row.get(col) or "").strip()
        if val:
            return val
    return ""


def annotate_row(
    row: dict[str, str],
    consequence_cache: dict[str, Optional[str]],
    *,
    col_prefix: str,
    hgvs_cols: list[str],
    access_date: str,
) -> dict[str, str]:
    consequence_col = f"{col_prefix}.mutational_consequence"
    access_col = f"{col_prefix}.access_date"
    error_col = f"{col_prefix}.error"

    out = {
        consequence_col: "",
        access_col: access_date,
        error_col: "",
    }

    candidates = _split_pipe_preserve_positions(_get_hgvs_for_row(row, hgvs_cols))
    if not candidates or all(c == "" for c in candidates):
        return out

    candidate_consequences: list[str] = []
    known_non_empty: list[str] = []
    for hgvs in candidates:
        if not hgvs:
            candidate_consequences.append("")
            continue
        consequence = consequence_cache.get(hgvs)
        value = "" if not consequence else str(consequence)
        candidate_consequences.append(value)
        if value:
            known_non_empty.append(value)

    unique_known = set(known_non_empty)

    if len(unique_known) > 1:
        out[consequence_col] = ""
        out[error_col] = (
            "Discrepant VEP mutational consequences across DNA candidates: "
            + "|".join(candidate_consequences)
        )
        return out

    if len(unique_known) == 1:
        out[consequence_col] = next(iter(unique_known))
        return out

    # No resolved consequence for any candidate.
    return out


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Annotate rows with VEP mutational consequence using genomic HGVS candidates"
        )
    )
    p.add_argument("input_file", help="Input CSV/TSV file")
    p.add_argument("output_file", help="Output CSV/TSV file")
    p.add_argument(
        "--vep-namespace",
        default="vep",
        help="Namespace for output columns (default: vep)",
    )
    p.add_argument(
        "--hgvs-cols",
        default="mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p",
        metavar="COLS",
        help=(
            "Comma-separated list of columns to try for HGVS input, in priority order. "
            "The first non-blank value found is used for each row. "
            "(default: mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p)"
        ),
    )
    p.add_argument(
        "--vep-api-url",
        default=ENSEMBL_API_URL_DEFAULT,
        help="Base URL for Ensembl REST API (default from ENSEMBL_API_URL env var)",
    )
    p.add_argument(
        "--vep-batch-size",
        type=int,
        default=int(os.environ.get("VEP_BATCH_SIZE", "500")),
        help="Number of HGVS values per VEP/recoder request batch (default: 500)",
    )
    p.add_argument(
        "--vep-workers",
        type=int,
        default=int(os.environ.get("VEP_WORKERS", "8")),
        help="Number of concurrent VEP batch requests (default: 8)",
    )
    p.add_argument(
        "--row-batch-size",
        type=int,
        default=int(os.environ.get("VEP_ROW_BATCH_SIZE", "1000")),
        help="Number of input rows per lookup/write batch (default: 1000)",
    )
    p.add_argument(
        "--vep-timeout-seconds",
        type=int,
        default=int(os.environ.get("VEP_TIMEOUT_SECONDS", "60")),
        help="HTTP timeout in seconds for VEP API calls (default: 60)",
    )
    p.add_argument("--delimiter", default="\t", help="Input/output delimiter (default TAB)")
    p.add_argument("--skip", type=int, default=0, help="Number of data rows to skip before annotation")
    p.add_argument("--limit", type=int, default=None, help="Maximum number of data rows to annotate")
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    p.add_argument(
        "--csv-field-size-limit",
        type=int,
        default=csv.field_size_limit(),
        metavar="BYTES",
        help="Maximum per-field character length for CSV/TSV parsing (default: %(default)s).",
    )
    p.add_argument(
        "--keep-existing",
        action="store_true",
        default=False,
        help=(
            "Preserve rows that already have a non-empty annotation value and only "
            "annotate rows with a blank consequence column. Reports how many rows "
            "were newly annotated."
        ),
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    csv.field_size_limit(args.csv_field_size_limit)

    if args.skip < 0:
        logger.error("--skip must be >= 0, got: %d", args.skip)
        sys.exit(1)
    if args.limit is not None and args.limit < 1:
        logger.error("--limit must be >= 1 when provided, got: %d", args.limit)
        sys.exit(1)
    if args.vep_batch_size < 1:
        logger.error("--vep-batch-size must be >= 1, got: %d", args.vep_batch_size)
        sys.exit(1)
    if args.vep_workers < 1:
        logger.error("--vep-workers must be >= 1, got: %d", args.vep_workers)
        sys.exit(1)
    if args.row_batch_size < 1:
        logger.error("--row-batch-size must be >= 1, got: %d", args.row_batch_size)
        sys.exit(1)
    if args.vep_timeout_seconds < 1:
        logger.error("--vep-timeout-seconds must be >= 1, got: %d", args.vep_timeout_seconds)
        sys.exit(1)

    delim = "\t" if args.delimiter == "\\t" else args.delimiter
    input_path = Path(args.input_file)
    output_path = Path(args.output_file)

    hgvs_cols = [c.strip() for c in args.hgvs_cols.split(",") if c.strip()]
    if not hgvs_cols:
        logger.error("--hgvs-cols must contain at least one column name")
        sys.exit(1)

    prefix = args.vep_namespace
    access_date = date.today().isoformat()
    ann_cols = [
        f"{prefix}.mutational_consequence",
        f"{prefix}.access_date",
        f"{prefix}.error",
    ]

    consequence_cache: dict[str, Optional[str]] = {}
    total_rows = 0
    kept_rows = 0
    newly_resolved_rows = 0
    discrepancy_rows = 0

    with input_path.open("r", encoding="utf-8", newline="") as in_fh, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears empty: %s", input_path)
            sys.exit(1)
        fieldnames = list(reader.fieldnames)
        out_fieldnames = fieldnames + [c for c in ann_cols if c not in fieldnames]

        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=delim,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        selected_reader = islice(
            reader,
            args.skip,
            None if args.limit is None else args.skip + args.limit,
        )

        batch_rows: list[dict[str, str]] = []
        for row in selected_reader:
            batch_rows.append(row)
            if len(batch_rows) < args.row_batch_size:
                continue

            batch_total, batch_kept, batch_resolved, batch_discrepancies = _process_batch(
                batch_rows,
                writer,
                out_fh,
                consequence_cache,
                hgvs_cols=hgvs_cols,
                col_prefix=prefix,
                access_date=access_date,
                api_url=args.vep_api_url,
                timeout_seconds=args.vep_timeout_seconds,
                vep_batch_size=args.vep_batch_size,
                vep_workers=args.vep_workers,
                keep_existing=args.keep_existing,
            )
            total_rows += batch_total
            kept_rows += batch_kept
            newly_resolved_rows += batch_resolved
            discrepancy_rows += batch_discrepancies
            batch_rows = []

        if batch_rows:
            batch_total, batch_kept, batch_resolved, batch_discrepancies = _process_batch(
                batch_rows,
                writer,
                out_fh,
                consequence_cache,
                hgvs_cols=hgvs_cols,
                col_prefix=prefix,
                access_date=access_date,
                api_url=args.vep_api_url,
                timeout_seconds=args.vep_timeout_seconds,
                vep_batch_size=args.vep_batch_size,
                vep_workers=args.vep_workers,
                keep_existing=args.keep_existing,
            )
            total_rows += batch_total
            kept_rows += batch_kept
            newly_resolved_rows += batch_resolved
            discrepancy_rows += batch_discrepancies

    logger.info(
        "Done. %d rows processed; %d kept (existing); %d rows newly annotated; %d discrepancy rows -> %s",
        total_rows,
        kept_rows,
        newly_resolved_rows,
        discrepancy_rows,
        output_path,
    )


def _process_batch(
    rows: list[dict[str, str]],
    writer: csv.DictWriter,
    out_fh: TextIO,
    consequence_cache: dict[str, Optional[str]],
    *,
    hgvs_cols: list[str],
    col_prefix: str,
    access_date: str,
    api_url: str,
    timeout_seconds: int,
    vep_batch_size: int,
    vep_workers: int = 1,
    keep_existing: bool = False,
) -> tuple[int, int, int, int]:
    """Process a batch of rows.

    Returns ``(total, kept, newly_resolved, discrepancies)``.
    *kept* is the number of rows whose existing annotation was preserved
    (only non-zero when *keep_existing* is True).
    """
    consequence_col = f"{col_prefix}.mutational_consequence"
    error_col = f"{col_prefix}.error"

    batch_hgvs: list[str] = []
    for row in rows:
        if keep_existing and (row.get(consequence_col) or "").strip():
            continue
        for hgvs in _split_pipe_preserve_positions(_get_hgvs_for_row(row, hgvs_cols)):
            if hgvs and hgvs not in consequence_cache:
                batch_hgvs.append(hgvs)

    if batch_hgvs:
        looked_up = get_functional_consequence(
            batch_hgvs,
            api_url=api_url,
            timeout_seconds=timeout_seconds,
            batch_size=vep_batch_size,
            max_workers=vep_workers,
        )
        consequence_cache.update(looked_up)

    kept = 0
    newly_resolved = 0
    discrepancies = 0
    for row in rows:
        if keep_existing and (row.get(consequence_col) or "").strip():
            kept += 1
            writer.writerow(row)
            continue
        ann = annotate_row(
            row,
            consequence_cache,
            col_prefix=col_prefix,
            hgvs_cols=hgvs_cols,
            access_date=access_date,
        )
        row.update(ann)
        writer.writerow(row)
        if ann[consequence_col]:
            newly_resolved += 1
        elif ann[error_col]:
            discrepancies += 1

    out_fh.flush()
    return len(rows), kept, newly_resolved, discrepancies


if __name__ == "__main__":
    main()
