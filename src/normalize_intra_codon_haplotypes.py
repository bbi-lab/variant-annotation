"""Rewrite intra-codon c.-haplotypes into delins notation.

This script scans a tabular variant file and rewrites entries in ``raw_hgvs_nt``
when they are multivariant/haplotype c.-notation strings that describe multiple
single-nucleotide substitutions within a single codon.

For rewritten rows:
- ``orig_raw_hgvs_nt`` stores the original value
- ``raw_hgvs_nt`` is replaced by ``c.<start>_<end>delins<ALT>`` (or accession-prefixed)

Rows that cannot be handled safely are left unchanged.
"""

from __future__ import annotations

import asyncio
import csv
import logging
import re
from pathlib import Path
from typing import Callable

import click

from src import map_variants as mv


logger = logging.getLogger(__name__)


_HAPLOTYPE_C_RE = re.compile(r"^(?:(?P<accession>[^:]+):)?c\.\[(?P<body>[^\]]+)\]$")
_C_SUB_RE = re.compile(r"^(?P<coord>\d+)(?P<ref>[ACGTN])>(?P<alt>[ACGTN])$")
_MAPPED_C_SUB_RE = re.compile(
    r"^(?:(?P<accession>[^:]+):)?c\.(?P<coord>\d+)(?P<ref>[ACGTN])>(?P<alt>[ACGTN])$"
)


def _detect_separator(file_path: str) -> str:
    """Return a delimiter based on the output suffix."""
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _is_blank(value: str | None) -> bool:
    return not (value or "").strip()


def _parse_haplotype_c(raw_hgvs_nt: str) -> tuple[str | None, list[str]] | None:
    """Parse accession-prefixed or bare c.-haplotype expressions.

    Returns ``(accession, components)`` for strings like:
    - ``NM_000001.1:c.[1A>G;3G>T]``
    - ``c.[1A>G;3G>T]``
    """
    text = raw_hgvs_nt.strip()
    if ";" not in text:
        return None
    m = _HAPLOTYPE_C_RE.match(text)
    if m is None:
        return None

    body = m.group("body")
    components = [part.strip() for part in body.split(";") if part.strip()]
    if len(components) < 2:
        return None

    return m.group("accession"), components


def _parse_component_substitution(component: str) -> tuple[int, str, str] | None:
    """Parse a coding substitution component like ``123A>G``.

    UTR or intronic coordinates (for example ``*12A>G``, ``-3A>G``, ``76+1G>A``)
    are intentionally excluded because they do not match this pattern.
    """
    m = _C_SUB_RE.match(component.strip().upper())
    if m is None:
        return None

    return int(m.group("coord")), m.group("ref"), m.group("alt")


def _parse_mapped_c_substitution(mapped_hgvs_c: str) -> tuple[int, str, str] | None:
    """Parse a mapped transcript substitution HGVS string."""
    m = _MAPPED_C_SUB_RE.match(mapped_hgvs_c.strip())
    if m is None:
        return None
    return int(m.group("coord")), m.group("ref").upper(), m.group("alt").upper()


def _compose_delins(
    parsed_subs: list[tuple[int, str, str]],
    target_sequence: str | None,
    accession: str | None,
) -> str | None:
    """Build a delins HGVS if all substitutions land in one codon."""
    if len(parsed_subs) < 2:
        return None

    codons = {(coord - 1) // 3 for coord, _, _ in parsed_subs}
    if len(codons) != 1:
        return None

    by_pos: dict[int, tuple[str, str]] = {coord: (ref, alt) for coord, ref, alt in parsed_subs}
    start = min(by_pos)
    end = max(by_pos)

    seq = (target_sequence or "").strip().upper()
    ref_parts: list[str] = []
    alt_parts: list[str] = []

    for pos in range(start, end + 1):
        if pos in by_pos:
            ref_base, alt_base = by_pos[pos]
        else:
            if not seq or pos > len(seq):
                return None
            ref_base = seq[pos - 1]
            alt_base = ref_base
            if ref_base not in {"A", "C", "G", "T", "N"}:
                return None

        ref_parts.append(ref_base)
        alt_parts.append(alt_base)

    alt_seq = "".join(alt_parts)
    ref_seq = "".join(ref_parts)
    if ref_seq == alt_seq:
        return None

    if start == end:
        variant_part = f"c.{start}delins{alt_seq}"
    else:
        variant_part = f"c.{start}_{end}delins{alt_seq}"

    if accession:
        return f"{accession}:{variant_part}"
    return variant_part


def _map_components_to_mapped_c(
    accession: str | None,
    components: list[str],
    target_sequence: str | None,
    row_label: str,
    dcd: dict | None,
) -> list[str] | None:
    """Map component substitutions and return mapped c. substitutions.

    Returns None if mapping fails for any component.
    """
    if accession:
        mapped: list[str] = []
        for component in components:
            full_hgvs = f"{accession}:c.{component}"
            hgvs_c, _, _, error, _ = mv._process_case1(full_hgvs, dcd)
            if error or not hgvs_c:
                return None
            mapped.append(hgvs_c)
        return mapped

    if _is_blank(target_sequence):
        return None
    if dcd is None:
        dcd = mv._try_import_dcd_mapping()

    row_entries = [(idx, f"c.{component}", "", 2) for idx, component in enumerate(components)]

    loop = asyncio.new_event_loop()
    try:
        asyncio.set_event_loop(loop)
        mapped_rows, transcript_nm = loop.run_until_complete(
            mv._run_dcd_mapping_pipeline(
                group_name=f"haplotype-row-{row_label}",
                target_sequence=(target_sequence or "").strip(),
                row_entries=row_entries,
                dcd=dcd,
            )
        )
    finally:
        asyncio.set_event_loop(None)
        loop.close()

    mapped_assay_by_idx: dict[int, str] = {}
    for idx, assay_hgvs, error in mapped_rows:
        if error or not assay_hgvs:
            return None
        mapped_assay_by_idx[idx] = assay_hgvs

    if len(mapped_assay_by_idx) != len(components):
        return None

    mapped_c: list[str] = []
    for idx in range(len(components)):
        data = mv._query_clingen_by_hgvs(mapped_assay_by_idx[idx])
        if data is None:
            return None
        _, hgvs_c, _ = mv._extract_hgvs_from_clingen(data, transcript_nm)
        if not hgvs_c:
            return None
        mapped_c.append(hgvs_c)

    return mapped_c


def normalize_haplotype_to_delins(
    raw_hgvs_nt: str,
    target_sequence: str | None,
    row_label: str,
    dcd: dict | None,
    mapper: Callable[[str | None, list[str], str | None, str, dict | None], list[str] | None] | None = None,
) -> str | None:
    """Rewrite one raw_hgvs_nt value if it is an intra-codon c.-haplotype.

    Returns rewritten HGVS on success, else None.
    """
    if mapper is None:
        mapper = _map_components_to_mapped_c

    parsed = _parse_haplotype_c(raw_hgvs_nt)
    if parsed is None:
        return None

    accession, components = parsed

    for component in components:
        if _parse_component_substitution(component) is None:
            return None

    try:
        mapped_c_parts = mapper(accession, components, target_sequence, row_label, dcd)
    except Exception as exc:  # pragma: no cover - defensive logging
        logger.warning("Skipping row %s due to mapping failure: %s", row_label, exc)
        return None

    if not mapped_c_parts or len(mapped_c_parts) != len(components):
        return None

    parsed_mapped: list[tuple[int, str, str]] = []
    for mapped_c in mapped_c_parts:
        parsed_sub = _parse_mapped_c_substitution(mapped_c)
        if parsed_sub is None:
            return None
        parsed_mapped.append(parsed_sub)

    return _compose_delins(parsed_mapped, target_sequence=target_sequence, accession=accession)


def haplotypes_to_delins(
    input_file: str,
    output_file: str,
    raw_hgvs_nt_col: str = "raw_hgvs_nt",
    target_sequence_col: str = "target_sequence",
    orig_raw_hgvs_nt_col: str = "orig_raw_hgvs_nt",
) -> dict[str, int]:
    """Rewrite intra-codon multivariants and preserve originals in a new column."""
    in_sep = _detect_separator(input_file)
    out_sep = _detect_separator(output_file)

    dcd: dict | None = None
    totals = {
        "rows": 0,
        "rewritten": 0,
    }

    with open(input_file, newline="", encoding="utf-8") as in_fh:
        reader = csv.DictReader(in_fh, delimiter=in_sep)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"Input file {input_file!r} is empty or missing a header")
        if raw_hgvs_nt_col not in fieldnames:
            raise ValueError(f"Column {raw_hgvs_nt_col!r} not found in input file")

        if orig_raw_hgvs_nt_col not in fieldnames:
            raw_idx = fieldnames.index(raw_hgvs_nt_col)
            fieldnames.insert(raw_idx, orig_raw_hgvs_nt_col)

        with open(output_file, "w", newline="", encoding="utf-8") as out_fh:
            writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter=out_sep)
            writer.writeheader()

            for row_idx, row in enumerate(reader, start=1):
                totals["rows"] += 1
                raw_hgvs_nt = row.get(raw_hgvs_nt_col, "")
                target_sequence = row.get(target_sequence_col, "")

                row[orig_raw_hgvs_nt_col] = raw_hgvs_nt

                rewritten = normalize_haplotype_to_delins(
                    raw_hgvs_nt=raw_hgvs_nt,
                    target_sequence=target_sequence,
                    row_label=str(row_idx),
                    dcd=dcd,
                )
                if rewritten:
                    row[raw_hgvs_nt_col] = rewritten
                    totals["rewritten"] += 1

                writer.writerow(row)

    return totals


@click.command()
@click.argument("input_file")
@click.argument("output_file")
@click.option("--raw-hgvs-nt-col", default="raw_hgvs_nt", show_default=True)
@click.option("--target-sequence-col", default="target_sequence", show_default=True)
@click.option("--orig-raw-hgvs-nt-col", default="orig_raw_hgvs_nt", show_default=True)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    default="INFO",
    show_default=True,
)
def main(
    input_file: str,
    output_file: str,
    raw_hgvs_nt_col: str,
    target_sequence_col: str,
    orig_raw_hgvs_nt_col: str,
    log_level: str,
) -> None:
    """Rewrite intra-codon c.-haplotypes to delins in INPUT_FILE and write OUTPUT_FILE."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    stats = haplotypes_to_delins(
        input_file=input_file,
        output_file=output_file,
        raw_hgvs_nt_col=raw_hgvs_nt_col,
        target_sequence_col=target_sequence_col,
        orig_raw_hgvs_nt_col=orig_raw_hgvs_nt_col,
    )

    click.echo(
        f"Wrote {stats['rows']} rows to {output_file}; rewritten={stats['rewritten']}"
    )


if __name__ == "__main__":
    main()
