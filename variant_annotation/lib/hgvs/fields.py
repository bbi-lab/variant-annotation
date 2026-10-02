"""HgvsFields — per-candidate coordinate derivation from HGVS strings."""

from __future__ import annotations

import io
import logging
from dataclasses import dataclass

from .parse import apply_vcf_anchor, parse_hgvs, reverse_complement

logger = logging.getLogger(__name__)


def _capture_warnings(hgvs_logger: logging.Logger) -> tuple[logging.Handler, io.StringIO]:
    buf = io.StringIO()
    handler = logging.StreamHandler(buf)
    handler.setLevel(logging.WARNING)
    hgvs_logger.addHandler(handler)
    return handler, buf


@dataclass
class HgvsFields:
    """VCF-style coordinate fields for a set of pipe-joined HGVS candidates.

    Each string field is pipe-joined, one entry per candidate, preserving
    candidate cardinality. Boolean fields are pipe-joined "true"/"false" strings
    to match the existing CSV column convention.
    """

    # Genomic fields (from g. candidates)
    g_chromosome: str = ""
    g_start: str = ""
    g_stop: str = ""
    g_ref: str = ""
    g_alt: str = ""

    # Coding fields (from c. candidates)
    c_transcript: str = ""
    c_start: str = ""
    c_stop: str = ""
    c_ref: str = ""
    c_alt: str = ""

    # Quality
    touches_intronic_region: str = ""
    spans_intron: str = ""
    warnings: str = ""


def derive_fields(
    hgvs_c_candidates: list[str],
    hgvs_g_candidates: list[str],
    *,
    resolve_missing_ref_alleles: bool = True,
) -> HgvsFields:
    """Derive coordinate fields from parallel c. and g. candidate lists.

    hgvs_c_candidates and hgvs_g_candidates are parallel — index i in each
    corresponds to the same candidate. Either list may be shorter (or empty)
    if candidates are missing at that level.
    """
    hgvs_logger = logging.getLogger("hgvs")

    g_chromosomes, g_starts, g_stops, g_refs, g_alts = [], [], [], [], []
    for candidate in hgvs_g_candidates:
        handler, buf = _capture_warnings(hgvs_logger)

        try:
            start, stop, ref, alt, _, _, chromosome = parse_hgvs(
                candidate, resolve_missing_ref_alleles=resolve_missing_ref_alleles
            )
        finally:
            hgvs_logger.removeHandler(handler)

        if alt == "inv" and ref:
            alt = reverse_complement(ref)

        start, stop, ref, alt = apply_vcf_anchor(candidate, start, stop, ref, alt)
        g_chromosomes.append(chromosome or "")
        g_starts.append(start or "")
        g_stops.append(stop or "")
        g_refs.append(ref or "")
        g_alts.append(alt or "")

    c_transcripts, c_starts, c_stops, c_refs, c_alts = [], [], [], [], []
    touches_intronic, spans_introns, all_warnings = [], [], []
    for candidate in hgvs_c_candidates:
        handler, buf = _capture_warnings(hgvs_logger)

        try:
            start, stop, ref, alt, touches, spans, transcript = parse_hgvs(
                candidate, resolve_missing_ref_alleles=resolve_missing_ref_alleles
            )
            warning_text = buf.getvalue().strip()
        finally:
            hgvs_logger.removeHandler(handler)

        if alt == "inv" and ref:
            alt = reverse_complement(ref)

        start, stop, ref, alt = apply_vcf_anchor(candidate, start, stop, ref, alt)
        c_transcripts.append(transcript or "")
        c_starts.append(start or "")
        c_stops.append(stop or "")
        c_refs.append(ref or "")
        c_alts.append(alt or "")
        touches_intronic.append("true" if touches else "false")
        spans_introns.append("true" if spans else "false")
        all_warnings.append(warning_text)

    return HgvsFields(
        g_chromosome="|".join(g_chromosomes),
        g_start="|".join(g_starts),
        g_stop="|".join(g_stops),
        g_ref="|".join(g_refs),
        g_alt="|".join(g_alts),
        c_transcript="|".join(c_transcripts),
        c_start="|".join(c_starts),
        c_stop="|".join(c_stops),
        c_ref="|".join(c_refs),
        c_alt="|".join(c_alts),
        touches_intronic_region="|".join(touches_intronic),
        spans_intron="|".join(spans_introns),
        warnings="|".join(all_warnings),
    )
