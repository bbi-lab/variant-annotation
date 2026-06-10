"""Construct the set of DNA variants protein-equivalent to each input variant.

CLI entry point. Handles CSV streaming; all translation and column logic lives
in src.lib.translation and src.lib.pipeline.reverse_translate_step.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import time
from pathlib import Path
from typing import Any, Optional

from dotenv import load_dotenv

from variant_annotation.lib.clients.coordinates import HgvsMapper
from variant_annotation.lib.clients.uta import UtaClient, connect_uta
from variant_annotation.lib.pipeline.reverse_translate_step import ColumnConfig, process_rows
from variant_annotation.lib.translation.types import TranslationConfig, WtCodonMode

load_dotenv()

logger = logging.getLogger(__name__)

PROGRESS_EVERY_ROWS = 1000


def _detect_separator(file_path: str) -> str:
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def reverse_translate_protein_variants(
    input_file: str,
    output_file: str,
    config: TranslationConfig = TranslationConfig(),
    columns: ColumnConfig = ColumnConfig(),
    *,
    resolve_missing_ref_alleles: bool = True,
    skip: int = 0,
    limit: int = 0,
) -> None:
    in_sep = _detect_separator(input_file)
    out_sep = _detect_separator(output_file)
    started = time.monotonic()

    uta_db_url = (os.environ.get("UTA_DB_URL") or "").strip()
    if not uta_db_url:
        raise RuntimeError("UTA_DB_URL must be set for protein reverse translation.")

    with connect_uta(uta_db_url) as uta_conn:
        transcripts = UtaClient(uta_conn)
        coordinates = HgvsMapper.from_url(uta_db_url, assembly=config.assembly)

        def flush_block(
            writer: csv.DictWriter,
            out_handle: Any,
            block_rows: list[dict[str, str]],
            block_start_idx: int,
        ) -> None:
            if not block_rows:
                return
            logger.info(
                "Translating protein block: rows %d-%d (%d row(s)).",
                block_start_idx,
                block_start_idx + len(block_rows) - 1,
                len(block_rows),
            )
            process_rows(
                block_rows,
                transcripts,
                coordinates,
                config,
                columns,
                derive_coordinate_fields=True,
                resolve_missing_ref_alleles=resolve_missing_ref_alleles,
            )
            for row in block_rows:
                writer.writerow(row)
            out_handle.flush()

        with (
            open(input_file, newline="", encoding="utf-8") as in_fh,
            open(output_file, "w", newline="", encoding="utf-8") as out_fh,
        ):
            reader = csv.DictReader(in_fh, delimiter=in_sep)
            in_fieldnames = list(reader.fieldnames or [])
            if not in_fieldnames:
                raise ValueError(f"Input file {input_file!r} is empty or missing a header row.")

            out_fieldnames = list(in_fieldnames)
            all_output_cols = (
                columns.mapped_hgvs_g,
                columns.mapped_hgvs_c,
                columns.mapped_hgvs_p,
                columns.assayed_variant_level,
                columns.mapped_hgvs_g_chromosome,
                columns.mapped_hgvs_c_transcript,
                columns.mapped_hgvs_g_start,
                columns.mapped_hgvs_g_stop,
                columns.mapped_hgvs_g_ref,
                columns.mapped_hgvs_g_alt,
                columns.mapped_hgvs_c_start,
                columns.mapped_hgvs_c_stop,
                columns.mapped_hgvs_c_ref,
                columns.mapped_hgvs_c_alt,
                columns.touches_intronic_region,
                columns.spans_intron,
                columns.reverse_translation_error,
                columns.reverse_translation_warnings,
            )
            for col in all_output_cols:
                if col not in out_fieldnames:
                    out_fieldnames.append(col)

            writer = csv.DictWriter(out_fh, fieldnames=out_fieldnames, delimiter=out_sep, extrasaction="ignore")
            writer.writeheader()
            out_fh.flush()

            pending_block: list[dict[str, str]] = []
            pending_block_start_idx = 0
            n_rows = n_protein_origin = n_translated = 0

            for idx, row in enumerate(reader):
                if idx < skip:
                    continue
                if limit > 0 and n_rows >= limit:
                    break

                protein_origin, target_for_translation = columns.assign_level(row)
                if protein_origin:
                    n_protein_origin += 1

                if target_for_translation:
                    if not pending_block:
                        pending_block_start_idx = idx
                    pending_block.append(row)
                    n_translated += 1
                else:
                    if pending_block:
                        flush_block(writer, out_fh, pending_block, pending_block_start_idx)
                        pending_block.clear()
                    # Pass-through row: derive coordinate fields from existing hgvs values.
                    columns.derive_passthrough_fields(row, resolve_missing_ref_alleles=resolve_missing_ref_alleles)
                    writer.writerow(row)
                    out_fh.flush()

                n_rows += 1
                if n_rows % PROGRESS_EVERY_ROWS == 0:
                    elapsed = max(time.monotonic() - started, 1e-9)
                    logger.info(
                        "Progress: %d rows processed (%d protein-origin, %d pending, %.1f rows/s).",
                        n_rows,
                        n_protein_origin,
                        len(pending_block),
                        n_rows / elapsed,
                    )

            if pending_block:
                flush_block(writer, out_fh, pending_block, pending_block_start_idx)

    elapsed = max(time.monotonic() - started, 1e-9)
    logger.info(
        "Done. %d rows processed, %d protein-origin, %d translated (%.1f rows/s).",
        n_rows,
        n_protein_origin,
        n_translated,
        n_rows / elapsed,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Construct the set of DNA variants protein-equivalent to each input variant.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_file")
    parser.add_argument("output_file")

    col = parser.add_argument_group("column names")
    col.add_argument("--mapped-hgvs-g", dest="mapped_hgvs_g", default="mapped_hgvs_g")
    col.add_argument("--mapped-hgvs-c", dest="mapped_hgvs_c", default="mapped_hgvs_c")
    col.add_argument("--mapped-hgvs-p", dest="mapped_hgvs_p", default="mapped_hgvs_p")
    col.add_argument(
        "--reverse-translation-error", dest="reverse_translation_error", default="reverse_translation_error"
    )
    col.add_argument(
        "--reverse-translation-warnings", dest="reverse_translation_warnings", default="reverse_translation_warnings"
    )
    col.add_argument("--assayed-variant-level", dest="assayed_variant_level", default="assayed_variant_level")
    col.add_argument("--mapped-hgvs-g-chromosome", dest="mapped_hgvs_g_chromosome", default="mapped_hgvs_g_chromosome")
    col.add_argument("--mapped-hgvs-c-transcript", dest="mapped_hgvs_c_transcript", default="mapped_hgvs_c_transcript")
    col.add_argument("--mapped-hgvs-g-start", dest="mapped_hgvs_g_start", default="mapped_hgvs_g_start")
    col.add_argument("--mapped-hgvs-g-stop", dest="mapped_hgvs_g_stop", default="mapped_hgvs_g_stop")
    col.add_argument("--mapped-hgvs-g-ref", dest="mapped_hgvs_g_ref", default="mapped_hgvs_g_ref")
    col.add_argument("--mapped-hgvs-g-alt", dest="mapped_hgvs_g_alt", default="mapped_hgvs_g_alt")
    col.add_argument("--mapped-hgvs-c-start", dest="mapped_hgvs_c_start", default="mapped_hgvs_c_start")
    col.add_argument("--mapped-hgvs-c-stop", dest="mapped_hgvs_c_stop", default="mapped_hgvs_c_stop")
    col.add_argument("--mapped-hgvs-c-ref", dest="mapped_hgvs_c_ref", default="mapped_hgvs_c_ref")
    col.add_argument("--mapped-hgvs-c-alt", dest="mapped_hgvs_c_alt", default="mapped_hgvs_c_alt")
    col.add_argument("--touches-intronic-region", dest="touches_intronic_region", default="touches_intronic_region")
    col.add_argument("--spans-intron", dest="spans_intron", default="spans_intron")
    col.add_argument(
        "--transcript-fallback-column",
        dest="transcript_fallback_columns",
        action="append",
        default=[],
        metavar="COLUMN",
    )

    tx = parser.add_argument_group("translation")
    tx.add_argument("--assembly", default="GRCh38")
    tx.add_argument("--include-indels", action="store_true", default=False)
    tx.add_argument("--max-indel-size", type=int, default=3)
    tx.add_argument("--no-strict-ref-aa", dest="strict_ref_aa", action="store_false", default=True)
    tx.add_argument("--use-inv-notation", action="store_true", default=False)
    tx.add_argument(
        "--substitutions-and-delins-only",
        dest="allow_length_changing_stop_candidates",
        action="store_false",
        default=True,
    )
    tx.add_argument("--wt-codon-mode", choices=["none", "unambiguous", "all"], default="none")

    parser.add_argument(
        "--resolve-missing-ref-alleles",
        dest="resolve_missing_ref_alleles",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--skip", type=int, default=0, metavar="N")
    parser.add_argument("--limit", type=int, default=0, metavar="N")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("--csv-field-size-limit", type=int, default=csv.field_size_limit(), metavar="BYTES")
    return parser


def main(argv: Optional[list[str]] = None) -> None:
    args = _build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    csv.field_size_limit(args.csv_field_size_limit)

    reverse_translate_protein_variants(
        args.input_file,
        args.output_file,
        TranslationConfig(
            wt_codon_mode=WtCodonMode(args.wt_codon_mode),
            assembly=args.assembly,
            include_indels=args.include_indels,
            max_indel_size=args.max_indel_size,
            strict_ref_aa=args.strict_ref_aa,
            use_inv_notation=args.use_inv_notation,
            allow_length_changing_stop_candidates=args.allow_length_changing_stop_candidates,
        ),
        ColumnConfig(
            mapped_hgvs_p=args.mapped_hgvs_p,
            mapped_hgvs_c=args.mapped_hgvs_c,
            mapped_hgvs_g=args.mapped_hgvs_g,
            transcript_fallback_columns=tuple(args.transcript_fallback_columns),
            assayed_variant_level=args.assayed_variant_level,
            reverse_translation_error=args.reverse_translation_error,
            reverse_translation_warnings=args.reverse_translation_warnings,
            mapped_hgvs_g_chromosome=args.mapped_hgvs_g_chromosome,
            mapped_hgvs_g_start=args.mapped_hgvs_g_start,
            mapped_hgvs_g_stop=args.mapped_hgvs_g_stop,
            mapped_hgvs_g_ref=args.mapped_hgvs_g_ref,
            mapped_hgvs_g_alt=args.mapped_hgvs_g_alt,
            mapped_hgvs_c_transcript=args.mapped_hgvs_c_transcript,
            mapped_hgvs_c_start=args.mapped_hgvs_c_start,
            mapped_hgvs_c_stop=args.mapped_hgvs_c_stop,
            mapped_hgvs_c_ref=args.mapped_hgvs_c_ref,
            mapped_hgvs_c_alt=args.mapped_hgvs_c_alt,
            touches_intronic_region=args.touches_intronic_region,
            spans_intron=args.spans_intron,
        ),
        resolve_missing_ref_alleles=args.resolve_missing_ref_alleles,
        skip=args.skip,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
