"""Core equivalence-class construction logic.

Implements the collapse → expand pattern:
  any input kind → ProteinConsequence → full coding/genomic equivalence class.

construct_equivalent_variants batches all inputs into a single subprocess call
(one collapse pass, one batch expand). construct_one is for single-variant
callers that need an immediate result.

This module knows about the ports (TranscriptSource, CoordinateTranslator) and
shared vocabulary. It does not import from clients/ — implementations are
injected by the composition root.
"""

from __future__ import annotations

import csv
import logging
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Sequence

from ..accessions import (
    extract_accession,
    looks_like_refseq_protein_accession,
    looks_like_transcript_accession,
)
from ..consequence import ProteinConsequence
from ._ports import CoordinateTranslator, TranscriptSource
from .types import (
    ProjectionPair,
    TranslationConfig,
    TranslationError,
    TranslationResult,
    VariantInput,
    WtCodonMode,
)

logger = logging.getLogger(__name__)

# Delimiter passed to reverse-translate-variants (--join-delimiter) to pipe-join
# each candidate column, and used to split the columns back apart row-wise.
_JOIN_DELIMITER = "|"

_WT_UNAMBIGUOUS_CODONS: dict[str, str] = {
    "Met": "ATG",
    "Trp": "TGG",
}

# Optional prediction parens: c_to_p renders inferred consequences as p.(Phe2335=), so the WT
# codon path must parse the parenthesized form, not just the bare p.Phe2335=.
_P_HGVS_AA_CHANGE_RE = re.compile(r"(?:^|:)p\.\(?(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|=)\)?$")


class _Kind(Enum):
    PROTEIN = "p"
    CODING = "c"
    GENOMIC = "g"
    UNKNOWN = "unknown"


def _classify_kind(hgvs: str) -> _Kind:
    if ":" not in hgvs:
        return _Kind.UNKNOWN

    body = hgvs.split(":", 1)[1].strip()
    if len(body) < 2 or body[1] != ".":
        return _Kind.UNKNOWN

    prefix = body[0].lower()
    if prefix == "p":
        return _Kind.PROTEIN
    if prefix in ("c", "n"):
        return _Kind.CODING
    if prefix in ("g", "m"):
        return _Kind.GENOMIC

    return _Kind.UNKNOWN


def _find_reverse_translate_cli() -> str:
    cli_path = shutil.which("reverse-translate-variants")
    if cli_path is None:
        raise RuntimeError("reverse-translate-variants was not found on PATH. Ensure variant-translation is installed.")

    return cli_path


@dataclass
class _BatchOutputRow:
    """Parsed output for one consequence from the reverse-translate-variants subprocess.

    projection_pairs holds the consequence's fanned-out equivalence class as
    ProjectionPair objects — one per coding candidate, carrying its genomic
    projection and variant type. The pairing is parsed straight from the CLI's
    position-aligned columns, never re-derived by index.
    """

    projection_pairs: list[ProjectionPair]


@dataclass
class _BatchErrorRow:
    """One per-input error from the reverse-translate-variants subprocess."""

    transcript: str
    hgvs_p: str
    error: str


def _split_aligned(value: str) -> list[str]:
    """Split a join-delimited candidate column into its per-candidate cells.

    Unlike a naive filter-empties split, this preserves empty cells: an empty
    genomic cell is a meaningful position holder (a candidate whose c→g
    projection failed), so the columns must stay index-aligned. An entirely
    blank column means zero candidates and returns an empty list.
    """
    return value.split(_JOIN_DELIMITER) if value.strip() else []


def _parse_projection_pairs(row: dict[str, str]) -> list[ProjectionPair]:
    """Parse one CLI output row into position-aligned ProjectionPair objects.

    The reverse-translate-variants CLI (run with --one-row-per-input) emits the
    equivalence class as parallel pipe-joined columns — variant_type[i],
    hgvs_c[i] and hgvs_g[i] all describe candidate i. We zip them column-wise so
    each candidate's fields land in a single pair. A pair is keyed by its coding
    expression; an empty genomic cell becomes hgvs_g=None (projection failed)
    rather than shifting subsequent candidates out of alignment.
    """
    coding = _split_aligned(row.get("hgvs_c") or "")
    genomic = _split_aligned(row.get("hgvs_g") or "")
    variant_types = _split_aligned(row.get("variant_type") or "")

    pairs: list[ProjectionPair] = []
    for i, hgvs_c in enumerate(coding):
        hgvs_c = hgvs_c.strip()
        # A candidate with no coding expression is not a pair; skip it
        # without disturbing the alignment of the remaining candidates.
        if not hgvs_c:
            continue

        hgvs_g = genomic[i].strip() if i < len(genomic) else ""
        variant_type = variant_types[i].strip() if i < len(variant_types) else ""
        pairs.append(
            ProjectionPair(
                hgvs_c=hgvs_c,
                hgvs_g=hgvs_g or None,
                variant_type=variant_type or None,
            )
        )

    return pairs


def _run_reverse_translate_batch(
    cli_path: str,
    consequences: list[ProteinConsequence],
    *,
    config: TranslationConfig,
) -> tuple[list[_BatchOutputRow], list[_BatchErrorRow]]:
    with tempfile.TemporaryDirectory(prefix="equiv_construct_") as tmp:
        tmp_path = Path(tmp)
        input_path = tmp_path / "input.tsv"
        output_path = tmp_path / "output.tsv"
        errors_path = tmp_path / "errors.tsv"

        with open(input_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=["transcript", "hgvs_p"], delimiter="\t")
            writer.writeheader()
            writer.writerows({"transcript": c.transcript, "hgvs_p": c.hgvs_p} for c in consequences)

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
            config.assembly,
            "--max-indel-size",
            str(config.max_indel_size),
            "--one-row-per-input",
            "--join-delimiter",
            "|",
            "--output",
            str(output_path),
            "--errors",
            str(errors_path),
        ]
        if config.include_indels:
            command.append("--include-indels")
        if not config.strict_ref_aa:
            command.append("--no-strict-ref-aa")
        if config.use_inv_notation:
            command.append("--use-inv-notation")
        if not config.allow_length_changing_stop_candidates:
            command.append("--substitutions-and-delins-only")

        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            stderr = (result.stderr or result.stdout or "").strip()
            raise RuntimeError("reverse-translate-variants failed" + (f": {stderr}" if stderr else "."))

        output_rows: list[_BatchOutputRow] = []
        if output_path.is_file():
            with open(output_path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    output_rows.append(_BatchOutputRow(projection_pairs=_parse_projection_pairs(row)))

        error_rows: list[_BatchErrorRow] = []
        if errors_path.is_file() and errors_path.stat().st_size > 0:
            with open(errors_path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    msg = (row.get("error") or "").strip()
                    if msg:
                        error_rows.append(
                            _BatchErrorRow(
                                transcript=(row.get("transcript") or "").strip(),
                                hgvs_p=(row.get("hgvs_p") or "").strip(),
                                error=msg,
                            )
                        )

        return output_rows, error_rows


def _resolve_consequence(
    inp: VariantInput,
    *,
    transcripts: TranscriptSource,
    coordinates: CoordinateTranslator,
) -> ProteinConsequence | str:
    """Collapse any input kind to ProteinConsequence. Returns error string on failure."""
    kind = _classify_kind(inp.hgvs)

    if kind is _Kind.PROTEIN:
        accession = extract_accession(inp.hgvs)
        transcript = inp.transcript
        if not transcript:
            if looks_like_refseq_protein_accession(accession):
                transcript = transcripts.transcript_for_protein(accession)
            elif looks_like_transcript_accession(accession):
                transcript = accession

        if not transcript:
            return f"Unable to resolve transcript for protein input: {inp.hgvs!r}"
        return ProteinConsequence(hgvs_p=inp.hgvs, transcript=transcript)

    if kind is _Kind.CODING:
        transcript = inp.transcript or extract_accession(inp.hgvs)
        if not transcript:
            return f"Unable to determine transcript for coding input: {inp.hgvs!r}"

        try:
            hgvs_p = coordinates.c_to_p(inp.hgvs)
        except Exception as exc:
            return f"Forward translation failed for {inp.hgvs!r}: {exc}"

        return ProteinConsequence(hgvs_p=hgvs_p, transcript=transcript)

    if kind is _Kind.GENOMIC:
        if not inp.transcript:
            return f"Transcript hint required for genomic input: {inp.hgvs!r}"

        try:
            c_hgvs = coordinates.g_to_c(inp.hgvs, inp.transcript)
            hgvs_p = coordinates.c_to_p(c_hgvs)
        except Exception as exc:
            return f"Forward translation failed for {inp.hgvs!r}: {exc}"

        return ProteinConsequence(hgvs_p=hgvs_p, transcript=inp.transcript)

    # Kind is unknown.
    return f"Unrecognised HGVS kind for input: {inp.hgvs!r}"


def _parse_protein_aa_change(hgvs_p: str) -> tuple[str, int, str] | None:
    m = _P_HGVS_AA_CHANGE_RE.search((hgvs_p or "").strip())
    if not m:
        return None

    ref_aa3 = m.group("ref")
    pos = int(m.group("pos"))
    alt_raw = m.group("alt")
    return ref_aa3, pos, ref_aa3 if alt_raw == "=" else alt_raw


def _build_wt_c_hgvs(transcript: str, aa_position: int, codon: str) -> str:
    c_start = (aa_position - 1) * 3 + 1
    return f"{transcript}:c.{c_start}_{aa_position * 3}delins{codon}"


def _get_wt_codon(
    ref_aa3: str,
    aa_position: int,
    transcript: str,
    mode: WtCodonMode,
    transcripts: TranscriptSource,
) -> str | None:
    if mode is WtCodonMode.NONE:
        return None
    if ref_aa3 in _WT_UNAMBIGUOUS_CODONS:
        return _WT_UNAMBIGUOUS_CODONS[ref_aa3]
    if mode is WtCodonMode.ALL:
        return transcripts.codon_at(transcript, aa_position)

    return None


def _apply_wt_codon(
    pairs: list[ProjectionPair],
    consequence: ProteinConsequence,
    *,
    config: TranslationConfig,
    transcripts: TranscriptSource,
    coordinates: CoordinateTranslator,
) -> list[ProjectionPair]:
    """Append the WT-codon delins projection pair when applicable.

    The WT codon is a synonymous intra-codon delins (the reference codon spelled
    out) that the reverse-translate tool does not itself emit. We add it as a
    whole ProjectionPair — coding expression, its genomic projection (None if
    c→g fails), and variant_type "delins" — so the coding/genomic pairing is
    added atomically. Appending one paired object (rather than pushing onto two
    independent lists) makes the coding/genomic desync impossible by
    construction: a failed genomic projection yields hgvs_g=None on the same
    pair, never a shorter genomic list.
    """
    if config.wt_codon_mode is WtCodonMode.NONE:
        return pairs

    aa_change = _parse_protein_aa_change(consequence.hgvs_p)
    if aa_change is None:
        return pairs

    ref_aa3, aa_pos, alt_aa3 = aa_change
    if ref_aa3 != alt_aa3:
        return pairs

    codon = _get_wt_codon(ref_aa3, aa_pos, consequence.transcript, config.wt_codon_mode, transcripts)
    if codon is None:
        return pairs

    wt_c = _build_wt_c_hgvs(consequence.transcript, aa_pos, codon)
    if any(pair.hgvs_c == wt_c for pair in pairs):
        return pairs

    try:
        wt_g = coordinates.c_to_g(wt_c)
    except Exception:
        wt_g = None

    return pairs + [ProjectionPair(hgvs_c=wt_c, hgvs_g=wt_g or None, variant_type="delins")]


def _build_result(
    inp: VariantInput,
    consequence: ProteinConsequence,
    output_row: _BatchOutputRow,
    error_messages: list[str],
    *,
    config: TranslationConfig,
    transcripts: TranscriptSource,
    coordinates: CoordinateTranslator,
) -> TranslationResult | TranslationError:
    pairs = _apply_wt_codon(
        output_row.projection_pairs,
        consequence,
        config=config,
        transcripts=transcripts,
        coordinates=coordinates,
    )
    error = error_messages[0] if error_messages else None

    if not pairs:
        return TranslationError(
            input=inp,
            error=error or "Reverse translation returned no candidate DNA variants",
        )

    hgvs_p_out = consequence.hgvs_p if _classify_kind(inp.hgvs) is not _Kind.PROTEIN else None
    return TranslationResult(
        input=inp,
        projection_pairs=pairs,
        hgvs_p=hgvs_p_out,
    )


def construct_one(
    inp: VariantInput,
    *,
    transcripts: TranscriptSource,
    coordinates: CoordinateTranslator,
    config: TranslationConfig = TranslationConfig(),
) -> TranslationResult | TranslationError:
    """Construct the equivalence class for a single VariantInput (one subprocess call)."""
    results, errors = construct_equivalent_variants(
        [inp], transcripts=transcripts, coordinates=coordinates, config=config
    )
    return results[0] if results else errors[0]


def construct_equivalent_variants(
    inputs: Sequence[VariantInput],
    *,
    transcripts: TranscriptSource,
    coordinates: CoordinateTranslator,
    config: TranslationConfig = TranslationConfig(),
) -> tuple[list[TranslationResult], list[TranslationError]]:
    """Construct the full equivalence class for each input variant.

    All inputs are collapsed to ProteinConsequence first, then a single
    subprocess call expands the entire batch. Every input appears in exactly
    one output list.

    The output set includes the input variant itself when it is a c./g. input
    that encodes its own protein consequence — the library makes no exclusions.
    Callers (e.g. mavedb-api) can deduplicate via vrs_digest before DB writes.
    """
    if not inputs:
        return [], []

    cli_path = _find_reverse_translate_cli()

    # Phase 1: collapse all inputs to ProteinConsequence (or error).
    consequences: list[ProteinConsequence | str] = [
        _resolve_consequence(inp, transcripts=transcripts, coordinates=coordinates) for inp in inputs
    ]

    # Phase 2: batch expand — one subprocess call for all valid consequences.
    batch_consequences: list[ProteinConsequence] = []
    batch_positions: list[int] = []
    early_errors: list[TranslationError] = []

    for i, (inp, consequence) in enumerate(zip(inputs, consequences)):
        if isinstance(consequence, str):
            early_errors.append(TranslationError(input=inp, error=consequence))
        else:
            batch_consequences.append(consequence)
            batch_positions.append(i)

    output_rows: list[_BatchOutputRow] = []
    error_rows: list[_BatchErrorRow] = []
    subprocess_error: str | None = None

    if batch_consequences:
        try:
            output_rows, error_rows = _run_reverse_translate_batch(cli_path, batch_consequences, config=config)
        except RuntimeError as exc:
            subprocess_error = str(exc)

    if subprocess_error:
        return [], early_errors + [TranslationError(input=inputs[i], error=subprocess_error) for i in batch_positions]

    if batch_consequences and len(output_rows) != len(batch_consequences):
        mismatch = f"reverse-translate-variants returned {len(output_rows)} rows for {len(batch_consequences)} inputs"
        return [], early_errors + [TranslationError(input=inputs[i], error=mismatch) for i in batch_positions]

    # Build error lookup keyed by (transcript, hgvs_p).
    error_messages_by_key: dict[tuple[str, str], list[str]] = defaultdict(list)
    for err in error_rows:
        error_messages_by_key[(err.transcript, err.hgvs_p)].append(err.error)

    # Phase 3: build results, applying WT codon mode per consequence.
    results: list[TranslationResult] = []
    late_errors: list[TranslationError] = []

    for consequence, output_row, pos in zip(batch_consequences, output_rows, batch_positions):
        inp = inputs[pos]

        key = (consequence.transcript, consequence.hgvs_p)
        error_msgs = error_messages_by_key.pop(key, [])

        outcome = _build_result(
            inp,
            consequence,
            output_row,
            error_msgs,
            config=config,
            transcripts=transcripts,
            coordinates=coordinates,
        )
        if isinstance(outcome, TranslationResult):
            results.append(outcome)
        else:
            late_errors.append(outcome)

    return results, early_errors + late_errors
