"""CSV composition root for the equivalence-class construction step.

ColumnConfig is the schema adapter — it owns all row read/write operations,
keeping column-name concerns in one place. process_rows wires VariantInputs
from rows and merges TranslationResults back using only ColumnConfig methods.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from ..accessions import extract_accession, looks_like_transcript_accession
from ..hgvs.fields import HgvsFields, derive_fields
from ..translation import (
    CoordinateTranslator,
    TranslationConfig,
    TranslationError,
    TranslationResult,
    TranscriptSource,
    VariantInput,
    construct_equivalent_variants,
)

logger = logging.getLogger(__name__)


@dataclass
class ColumnConfig:
    """CSV column name mapping and row adapter for the construction step.

    Read methods extract typed values from row dicts; write methods merge
    results back. All column-name knowledge lives here.
    """

    # Input
    mapped_hgvs_p: str = "mapped_hgvs_p"
    mapped_hgvs_c: str = "mapped_hgvs_c"
    mapped_hgvs_g: str = "mapped_hgvs_g"
    transcript_fallback_columns: tuple[str, ...] = field(default_factory=tuple)

    # Output
    assayed_variant_level: str = "assayed_variant_level"
    reverse_translation_error: str = "reverse_translation_error"
    reverse_translation_warnings: str = "reverse_translation_warnings"

    # Derived genomic coordinates
    mapped_hgvs_g_chromosome: str = "mapped_hgvs_g_chromosome"
    mapped_hgvs_g_start: str = "mapped_hgvs_g_start"
    mapped_hgvs_g_stop: str = "mapped_hgvs_g_stop"
    mapped_hgvs_g_ref: str = "mapped_hgvs_g_ref"
    mapped_hgvs_g_alt: str = "mapped_hgvs_g_alt"

    # Derived coding coordinates
    mapped_hgvs_c_transcript: str = "mapped_hgvs_c_transcript"
    mapped_hgvs_c_start: str = "mapped_hgvs_c_start"
    mapped_hgvs_c_stop: str = "mapped_hgvs_c_stop"
    mapped_hgvs_c_ref: str = "mapped_hgvs_c_ref"
    mapped_hgvs_c_alt: str = "mapped_hgvs_c_alt"

    # Intronic flags
    touches_intronic_region: str = "touches_intronic_region"
    spans_intron: str = "spans_intron"

    # --- read methods ---

    def read_hgvs(self, row: dict) -> str | None:
        """Return the first non-blank HGVS value from the row (p. > c. > g.)."""
        for col in (self.mapped_hgvs_p, self.mapped_hgvs_c, self.mapped_hgvs_g):
            val = (row.get(col) or "").strip()
            if val:
                return val
        return None

    def read_transcript_hint(self, row: dict) -> str | None:
        """Return a transcript accession from fallback columns, or None."""
        for col in self.transcript_fallback_columns:
            val = (row.get(col) or "").strip()
            acc = extract_accession(val)
            if acc and looks_like_transcript_accession(acc):
                return acc
        return None

    def classify_row(self, row: dict) -> tuple[bool, bool]:
        """Return (protein_origin, target_for_translation) for a CSV row."""
        existing_level = (row.get(self.assayed_variant_level) or "").strip().lower()
        has_g = bool((row.get(self.mapped_hgvs_g) or "").strip())
        has_c = bool((row.get(self.mapped_hgvs_c) or "").strip())
        has_p = bool((row.get(self.mapped_hgvs_p) or "").strip())
        protein_origin = existing_level == "protein" or (existing_level == "" and has_p and not has_c and not has_g)
        target_for_translation = protein_origin and not has_c and not has_g and has_p
        return protein_origin, target_for_translation

    def _read_candidates(self, row: dict, column: str) -> list[str]:
        """Split a pipe-joined candidate column into a list. The CSV-side half of
        the pipe-join boundary that write_result owns on the way out."""
        text = (row.get(column) or "").strip()
        return [p.strip() for p in text.split("|") if p.strip()] if text else []

    # --- write methods ---

    def write_result(self, row: dict, result: TranslationResult) -> None:
        """Merge a successful TranslationResult back into a row dict."""
        row[self.assayed_variant_level] = "protein"
        row[self.mapped_hgvs_c] = "|".join(result.hgvs_c_candidates)
        row[self.mapped_hgvs_g] = "|".join(result.hgvs_g_candidates)
        if result.hgvs_p:
            row[self.mapped_hgvs_p] = result.hgvs_p

    def write_error(self, row: dict, error: TranslationError) -> None:
        """Append a TranslationError message into the error column."""
        row[self.assayed_variant_level] = "protein"
        existing = (row.get(self.reverse_translation_error) or "").strip()
        row[self.reverse_translation_error] = f"{existing}; {error.error}" if existing else error.error

    def write_fields(self, row: dict, fields: HgvsFields) -> None:
        """Write derived coordinate fields into a row dict."""
        row[self.mapped_hgvs_g_chromosome] = fields.g_chromosome
        row[self.mapped_hgvs_g_start] = fields.g_start
        row[self.mapped_hgvs_g_stop] = fields.g_stop
        row[self.mapped_hgvs_g_ref] = fields.g_ref
        row[self.mapped_hgvs_g_alt] = fields.g_alt
        row[self.mapped_hgvs_c_transcript] = fields.c_transcript
        row[self.mapped_hgvs_c_start] = fields.c_start
        row[self.mapped_hgvs_c_stop] = fields.c_stop
        row[self.mapped_hgvs_c_ref] = fields.c_ref
        row[self.mapped_hgvs_c_alt] = fields.c_alt
        row[self.touches_intronic_region] = fields.touches_intronic_region
        row[self.spans_intron] = fields.spans_intron
        row[self.reverse_translation_warnings] = fields.warnings

    def assign_level(self, row: dict) -> tuple[bool, bool]:
        """Classify the row, write its assayed_variant_level, and return
        (protein_origin, target_for_translation)."""
        protein_origin, target_for_translation = self.classify_row(row)
        if protein_origin:
            row[self.assayed_variant_level] = "protein"
        elif (row.get(self.mapped_hgvs_c) or "").strip() or (row.get(self.mapped_hgvs_g) or "").strip():
            row[self.assayed_variant_level] = "dna"
        return protein_origin, target_for_translation

    def derive_passthrough_fields(self, row: dict, *, resolve_missing_ref_alleles: bool = True) -> None:
        """Derive coordinate fields for a pass-through row from its existing
        c./g. candidate columns. No-op when neither column holds candidates."""
        c_candidates = self._read_candidates(row, self.mapped_hgvs_c)
        g_candidates = self._read_candidates(row, self.mapped_hgvs_g)
        if c_candidates or g_candidates:
            self.write_fields(
                row, derive_fields(c_candidates, g_candidates, resolve_missing_ref_alleles=resolve_missing_ref_alleles)
            )


def process_rows(
    rows: list[dict],
    transcripts: TranscriptSource,
    coordinates: CoordinateTranslator,
    config: TranslationConfig = TranslationConfig(),
    columns: ColumnConfig = ColumnConfig(),
    *,
    derive_coordinate_fields: bool = True,
    resolve_missing_ref_alleles: bool = True,
) -> list[dict]:
    """Run equivalence-class construction over a batch of CSV rows.

    Rows without a recognisable HGVS value are returned unchanged.
    """
    indexed_inputs = [
        (i, VariantInput(hgvs=hgvs, transcript=columns.read_transcript_hint(row)))
        for i, row in enumerate(rows)
        if (hgvs := columns.read_hgvs(row))
    ]

    if not indexed_inputs:
        return rows

    indexes, inputs = zip(*indexed_inputs)

    results, errors = construct_equivalent_variants(
        inputs,
        transcripts=transcripts,
        coordinates=coordinates,
        config=config,
    )

    result_map = {id(r.input): r for r in results}
    error_map = {id(e.input): e for e in errors}
    result_by_index = {idx: result_map[id(inp)] for idx, inp in zip(indexes, inputs) if id(inp) in result_map}
    error_by_index = {idx: error_map[id(inp)] for idx, inp in zip(indexes, inputs) if id(inp) in error_map}

    for i, row in enumerate(rows):
        if i in error_by_index:
            columns.write_error(row, error_by_index[i])
        elif i in result_by_index:
            result = result_by_index[i]
            columns.write_result(row, result)
            if derive_coordinate_fields:
                columns.write_fields(
                    row,
                    derive_fields(
                        result.hgvs_c_candidates,
                        result.hgvs_g_candidates,
                        resolve_missing_ref_alleles=resolve_missing_ref_alleles,
                    ),
                )

    return rows
