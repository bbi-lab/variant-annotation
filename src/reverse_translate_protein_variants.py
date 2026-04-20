"""Reverse-translate mapped protein-only variants into HGVS c./g. candidates.

This is the second step of the variant-annotation pipeline. It only processes
rows that represent assay-level protein variants, defined as rows with a
non-blank protein HGVS value and blank c./g. HGVS values. Reverse-translated
candidate HGVS c./g. strings are joined with ``|`` and written back to the
existing mapped columns while preserving input row order.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import subprocess
import tempfile
import time
from urllib.parse import urlsplit
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional

from dotenv import load_dotenv
import psycopg2  # type: ignore[import-untyped]

load_dotenv()

logger = logging.getLogger(__name__)

PROGRESS_EVERY_ROWS = 1000
TRANSCRIPT_ACCESSION_PREFIXES = ("NM_", "NR_", "XM_", "XR_", "ENST", "LRG_")
REFSEQ_PROTEIN_ACCESSION_PREFIXES = ("NP_", "XP_", "YP_", "WP_")


def _detect_separator(file_path: str) -> str:
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _extract_accession(hgvs_string: str) -> str:
    token = (hgvs_string or "").strip()
    if ":" not in token:
        return ""
    return token.split(":", 1)[0].strip()


def _looks_like_transcript_accession(accession: str) -> bool:
    return accession.startswith(TRANSCRIPT_ACCESSION_PREFIXES)


def _looks_like_refseq_protein_accession(accession: str) -> bool:
    return accession.startswith(REFSEQ_PROTEIN_ACCESSION_PREFIXES)


def _transcript_sort_key(accession: str) -> tuple[int, int, str]:
    if accession.startswith(("NM_", "NR_")):
        prefix_rank = 0
    elif accession.startswith(("XM_", "XR_")):
        prefix_rank = 1
    elif accession.startswith(("ENST", "LRG_")):
        prefix_rank = 2
    else:
        prefix_rank = 3

    _, _, version = accession.partition(".")
    version_number = int(version) if version.isdigit() else -1
    return prefix_rank, -version_number, accession


def _append_error(existing_error: str, new_error: str) -> str:
    existing = existing_error.strip()
    if not existing:
        return new_error
    if new_error in existing:
        return existing
    return f"{existing}; {new_error}"


def _find_reverse_translate_cli() -> str:
    cli_path = shutil.which("reverse-translate-variants")
    if cli_path is None:
        raise RuntimeError(
            "reverse-translate-variants was not found on PATH. Ensure variant-translation is installed in the image."
        )
    return cli_path


def _build_pg_connection_kwargs(database_url: str) -> dict[str, Any]:
    parsed = urlsplit(database_url)
    path_parts = [part for part in parsed.path.split("/") if part]
    if not path_parts:
        raise RuntimeError(f"UTA_DB_URL is missing a database name: {database_url}")

    kwargs: dict[str, Any] = {"dbname": path_parts[0]}
    if parsed.hostname:
        kwargs["host"] = parsed.hostname
    if parsed.port:
        kwargs["port"] = parsed.port
    if parsed.username:
        kwargs["user"] = parsed.username
    if parsed.password:
        kwargs["password"] = parsed.password
    return kwargs


def _resolve_transcript_from_refseq_protein_id(connection: Any, protein_accession: str) -> str:
    with connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT tx_ac
            FROM uta_20241220.associated_accessions
            WHERE pro_ac = %s
            """,
            (protein_accession,),
        )
        transcript_accessions = sorted(
            {
                row[0]
                for row in cursor.fetchall()
                if row and row[0] and _looks_like_transcript_accession(row[0])
            },
            key=_transcript_sort_key,
        )
    return transcript_accessions[0] if transcript_accessions else ""


def _run_reverse_translate_batch(
    cli_path: str,
    rows: list[dict[str, str]],
    *,
    assembly: str,
    include_indels: bool,
    max_indel_size: int,
    strict_ref_aa: bool,
    use_inv_notation: bool,
    allow_length_changing_stop_candidates: bool,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    with tempfile.TemporaryDirectory(prefix="reverse_translate_block_") as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        input_path = temp_dir / "input.tsv"
        output_path = temp_dir / "output.tsv"
        errors_path = temp_dir / "errors.tsv"

        with open(input_path, "w", newline="", encoding="utf-8") as input_handle:
            writer = csv.DictWriter(input_handle, fieldnames=["transcript", "hgvs_p"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        command = [
            cli_path,
            "--input",
            str(input_path),
            "--input-format",
            "tsv",
            "--transcript-column",
            "transcript",
            "--hgvs-p-column",
            "hgvs_p",
            "--assembly",
            assembly,
            "--max-indel-size",
            str(max_indel_size),
            "--one-row-per-input",
            "--join-delimiter",
            "|",
            "--output",
            str(output_path),
            "--errors",
            str(errors_path),
        ]
        if include_indels:
            command.append("--include-indels")
        if not strict_ref_aa:
            command.append("--no-strict-ref-aa")
        if use_inv_notation:
            command.append("--use-inv-notation")
        if not allow_length_changing_stop_candidates:
            command.append("--substitutions-and-delins-only")

        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            stderr = (completed.stderr or completed.stdout or "").strip()
            raise RuntimeError(
                "reverse-translate-variants failed"
                + (f": {stderr}" if stderr else ".")
            )

        translated_rows: list[dict[str, str]] = []
        if output_path.is_file():
            with open(output_path, newline="", encoding="utf-8") as output_handle:
                translated_rows = list(csv.DictReader(output_handle, delimiter="\t"))

        error_rows: list[dict[str, str]] = []
        if errors_path.is_file() and errors_path.stat().st_size > 0:
            with open(errors_path, newline="", encoding="utf-8") as error_handle:
                error_rows = list(csv.DictReader(error_handle, delimiter="\t"))

        return translated_rows, error_rows


def _classify_row(
    row: dict[str, str],
    *,
    mapped_hgvs_g_col: str,
    mapped_hgvs_c_col: str,
    mapped_hgvs_p_col: str,
    assayed_variant_level_col: str,
) -> tuple[bool, bool]:
    existing_level = (row.get(assayed_variant_level_col) or "").strip().lower()
    has_g = bool((row.get(mapped_hgvs_g_col) or "").strip())
    has_c = bool((row.get(mapped_hgvs_c_col) or "").strip())
    has_p = bool((row.get(mapped_hgvs_p_col) or "").strip())

    protein_origin = existing_level == "protein" or (existing_level == "" and has_p and not has_c and not has_g)
    target_for_reverse_translation = protein_origin and not has_c and not has_g and has_p
    return protein_origin, target_for_reverse_translation


def _resolve_transcript_accession(
    row: dict[str, str],
    *,
    mapped_hgvs_p_col: str,
    transcript_fallback_columns: tuple[str, ...],
    uta_connection: Any,
) -> str:
    protein_hgvs = (row.get(mapped_hgvs_p_col) or "").strip()
    protein_accession = _extract_accession(protein_hgvs)

    if protein_accession:
        if _looks_like_transcript_accession(protein_accession):
            return protein_accession

        if _looks_like_refseq_protein_accession(protein_accession):
            transcript_accession = _resolve_transcript_from_refseq_protein_id(uta_connection, protein_accession)
            if transcript_accession:
                return transcript_accession

    for col in transcript_fallback_columns:
        value = (row.get(col) or "").strip()
        accession = _extract_accession(value)
        if accession and _looks_like_transcript_accession(accession):
            return accession

    return ""


def reverse_translate_protein_variants(
    input_file: str,
    output_file: str,
    *,
    mapped_hgvs_g_col: str = "mapped_hgvs_g",
    mapped_hgvs_c_col: str = "mapped_hgvs_c",
    mapped_hgvs_p_col: str = "mapped_hgvs_p",
    mapping_error_col: str = "mapping_error",
    assayed_variant_level_col: str = "assayed_variant_level",
    transcript_fallback_columns: tuple[str, ...] = (),
    assembly: str = "GRCh38",
    include_indels: bool = False,
    max_indel_size: int = 3,
    strict_ref_aa: bool = True,
    use_inv_notation: bool = False,
    allow_length_changing_stop_candidates: bool = True,
) -> None:
    in_sep = _detect_separator(input_file)
    out_sep = _detect_separator(output_file)
    started = time.monotonic()
    cli_path = _find_reverse_translate_cli()
    uta_db_url = (os.environ.get("UTA_DB_URL") or "").strip()
    if not uta_db_url:
        raise RuntimeError("UTA_DB_URL must be set for protein reverse translation.")

    def flush_block(
        writer: csv.DictWriter,
        out_handle,
        block_rows: list[dict[str, str]],
        block_start_idx: int,
    ) -> None:
        if not block_rows:
            return

        logger.info(
            "Reverse-translating contiguous protein block: rows %d-%d (%d row(s)).",
            block_start_idx,
            block_start_idx + len(block_rows) - 1,
            len(block_rows),
        )

        input_rows: list[dict[str, str]] = []
        block_errors: list[str] = [""] * len(block_rows)
        translated_target_indexes: list[int] = []

        for offset, row in enumerate(block_rows):
            row[assayed_variant_level_col] = "protein"
            transcript_accession = _resolve_transcript_accession(
                row,
                mapped_hgvs_p_col=mapped_hgvs_p_col,
                transcript_fallback_columns=transcript_fallback_columns,
                uta_connection=uta_connection,
            )
            if not transcript_accession:
                block_errors[offset] = "Unable to resolve transcript accession for protein reverse translation"
                continue

            input_rows.append(
                {
                    "transcript": transcript_accession,
                    "hgvs_p": (row.get(mapped_hgvs_p_col) or "").strip(),
                }
            )
            translated_target_indexes.append(offset)

        translated_rows: list[dict[str, str]] = []
        error_rows: list[dict[str, str]] = []
        if input_rows:
            try:
                translated_rows, error_rows = _run_reverse_translate_batch(
                    cli_path,
                    input_rows,
                    assembly=assembly,
                    include_indels=include_indels,
                    max_indel_size=max_indel_size,
                    strict_ref_aa=strict_ref_aa,
                    use_inv_notation=use_inv_notation,
                    allow_length_changing_stop_candidates=allow_length_changing_stop_candidates,
                )
            except RuntimeError as exc:
                for offset in translated_target_indexes:
                    block_errors[offset] = str(exc)
                translated_rows = []
                error_rows = []

        if input_rows and translated_rows and len(translated_rows) != len(input_rows):
            mismatch_error = (
                "reverse-translate-variants returned an unexpected number of rows "
                f"({len(translated_rows)} for {len(input_rows)} inputs)"
            )
            for offset in translated_target_indexes:
                if not block_errors[offset]:
                    block_errors[offset] = mismatch_error
            translated_rows = []
            error_rows = []

        error_messages_by_key: dict[tuple[str, str], list[str]] = defaultdict(list)
        for error_row in error_rows:
            error_key = (
                (error_row.get("transcript") or "").strip(),
                (error_row.get("hgvs_p") or "").strip(),
            )
            error_message = (error_row.get("error") or "").strip()
            if error_message:
                error_messages_by_key[error_key].append(error_message)

        for input_row, translated_row, offset in zip(input_rows, translated_rows, translated_target_indexes):
            row = block_rows[offset]
            hgvs_c = (translated_row.get("hgvs_c") or "").strip()
            hgvs_g = (translated_row.get("hgvs_g") or "").strip()
            row[mapped_hgvs_c_col] = hgvs_c
            row[mapped_hgvs_g_col] = hgvs_g
            error_key = (input_row["transcript"], input_row["hgvs_p"])
            if error_messages_by_key[error_key]:
                block_errors[offset] = error_messages_by_key[error_key].pop(0)
            if not hgvs_c and not hgvs_g and not block_errors[offset]:
                block_errors[offset] = "Reverse translation returned no candidate DNA variants"

        for offset, row in enumerate(block_rows):
            if block_errors[offset]:
                row[mapping_error_col] = _append_error(
                    row.get(mapping_error_col, ""),
                    block_errors[offset],
                )
            writer.writerow(row)
            out_handle.flush()

    with psycopg2.connect(**_build_pg_connection_kwargs(uta_db_url)) as uta_connection:
        with open(input_file, newline="", encoding="utf-8") as in_fh, open(
            output_file, "w", newline="", encoding="utf-8"
        ) as out_fh:
            reader = csv.DictReader(in_fh, delimiter=in_sep)
            in_fieldnames = list(reader.fieldnames or [])
            if not in_fieldnames:
                raise ValueError(f"Input file {input_file!r} is empty or missing a header row.")

            out_fieldnames = list(in_fieldnames)
            for col in (mapped_hgvs_g_col, mapped_hgvs_c_col, mapped_hgvs_p_col, mapping_error_col, assayed_variant_level_col):
                if col not in out_fieldnames:
                    out_fieldnames.append(col)

            writer = csv.DictWriter(
                out_fh,
                fieldnames=out_fieldnames,
                delimiter=out_sep,
                extrasaction="ignore",
            )
            writer.writeheader()
            out_fh.flush()

            pending_block: list[dict[str, str]] = []
            pending_block_start_idx = 0
            n_rows = 0
            n_protein_origin = 0
            n_reverse_translated = 0

            for idx, row in enumerate(reader):
                protein_origin, target_for_reverse_translation = _classify_row(
                    row,
                    mapped_hgvs_g_col=mapped_hgvs_g_col,
                    mapped_hgvs_c_col=mapped_hgvs_c_col,
                    mapped_hgvs_p_col=mapped_hgvs_p_col,
                    assayed_variant_level_col=assayed_variant_level_col,
                )
                if protein_origin:
                    row[assayed_variant_level_col] = "protein"
                    n_protein_origin += 1
                elif (row.get(mapped_hgvs_c_col) or "").strip() or (row.get(mapped_hgvs_g_col) or "").strip():
                    row[assayed_variant_level_col] = "dna"

                if target_for_reverse_translation:
                    if not pending_block:
                        pending_block_start_idx = idx
                    pending_block.append(row)
                    n_reverse_translated += 1
                else:
                    if pending_block:
                        flush_block(writer, out_fh, pending_block, pending_block_start_idx)
                        pending_block.clear()
                    writer.writerow(row)
                    out_fh.flush()

                n_rows += 1
                if n_rows % PROGRESS_EVERY_ROWS == 0:
                    elapsed = max(time.monotonic() - started, 1e-9)
                    logger.info(
                        "Progress: %d rows processed (%d protein-origin, %d pending reverse translation, %.1f rows/s).",
                        n_rows,
                        n_protein_origin,
                        len(pending_block),
                        n_rows / elapsed,
                    )

            if pending_block:
                flush_block(writer, out_fh, pending_block, pending_block_start_idx)

    elapsed = max(time.monotonic() - started, 1e-9)
    logger.info(
        "Done. %d rows processed, %d protein-origin rows identified, %d protein-only rows reverse-translated (%.1f rows/s).",
        n_rows,
        n_protein_origin,
        n_reverse_translated,
        n_rows / elapsed,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reverse-translate protein-only mapped variants into joined HGVS c./g. candidates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_file", help="Input CSV/TSV file path.")
    parser.add_argument("output_file", help="Output CSV/TSV file path.")
    parser.add_argument("--mapped-hgvs-g", dest="mapped_hgvs_g_col", default="mapped_hgvs_g")
    parser.add_argument("--mapped-hgvs-c", dest="mapped_hgvs_c_col", default="mapped_hgvs_c")
    parser.add_argument("--mapped-hgvs-p", dest="mapped_hgvs_p_col", default="mapped_hgvs_p")
    parser.add_argument("--mapping-error", dest="mapping_error_col", default="mapping_error")
    parser.add_argument(
        "--assayed-variant-level",
        dest="assayed_variant_level_col",
        default="assayed_variant_level",
        help="Output column name distinguishing dna vs protein assay level.",
    )
    parser.add_argument(
        "--transcript-fallback-column",
        dest="transcript_fallback_columns",
        action="append",
        default=[],
        metavar="COLUMN",
        help="Additional column containing transcript-anchored HGVS values to mine for transcript accessions.",
    )
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument("--include-indels", action="store_true", default=False)
    parser.add_argument("--max-indel-size", type=int, default=3)
    parser.add_argument(
        "--no-strict-ref-aa",
        dest="strict_ref_aa",
        action="store_false",
        default=True,
        help="Disable reference amino-acid validation against the resolved transcript.",
    )
    parser.add_argument("--use-inv-notation", action="store_true", default=False)
    parser.add_argument(
        "--substitutions-and-delins-only",
        dest="allow_length_changing_stop_candidates",
        action="store_false",
        default=True,
        help="For premature stops, suppress length-changing insertion/deletion candidates.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return parser


def main(argv: Optional[list[str]] = None) -> None:
    args = _build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    reverse_translate_protein_variants(
        input_file=args.input_file,
        output_file=args.output_file,
        mapped_hgvs_g_col=args.mapped_hgvs_g_col,
        mapped_hgvs_c_col=args.mapped_hgvs_c_col,
        mapped_hgvs_p_col=args.mapped_hgvs_p_col,
        mapping_error_col=args.mapping_error_col,
        assayed_variant_level_col=args.assayed_variant_level_col,
        transcript_fallback_columns=tuple(args.transcript_fallback_columns),
        assembly=args.assembly,
        include_indels=args.include_indels,
        max_indel_size=args.max_indel_size,
        strict_ref_aa=args.strict_ref_aa,
        use_inv_notation=args.use_inv_notation,
        allow_length_changing_stop_candidates=args.allow_length_changing_stop_candidates,
    )


if __name__ == "__main__":
    main()