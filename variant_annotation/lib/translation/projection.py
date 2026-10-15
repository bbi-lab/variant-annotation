# Copyright (C) 2026 Jeremy Stone
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Projection of a variant between coding and genomic coordinates.

Every coding↔genomic projection the library and its consumers make goes through these two functions, so
the rules for writing a projection live in one place:

* A change counts as one genomic change only when its coding and genomic spans hold the same number of
  bases. hgvs maps the two ends of an interval separately, so a coding change across an exon junction (a
  codon split by an intron, an insertion between two exons' bases) otherwise becomes a genomic change that
  also deletes the intron. No single contiguous genomic change describes it, so its projection is ``None``.
* hgvs normalization writes a delins whose inserted bases are the reverse complement of the reference as
  ``inv``, which VRS cannot represent. It is written as the equivalent delins instead.

A projection hgvs rejects as invalid is ``None``. Any other failure is raised, including missing reference
data, which hgvs cannot tell apart from an unreachable data provider, so a caller can tell a change with no
projection from a projection it could not make.
"""

from __future__ import annotations

import re

from ._ports import CoordinateTranslator

# An accession-qualified c./g./n. position or interval with its edit. Positions keep any intronic offset
# (``123+5``) and UTR prefix (``-12``, ``*7``).
_LOCATED_EDIT_RE = re.compile(
    r"^(?P<prefix>[^:]+:[cgn]\.)(?P<start>[-*]?\d+(?:[+-]\d+)?)(?:_(?P<end>[-*]?\d+(?:[+-]\d+)?))?(?P<edit>\D.*)$"
)
_PLAIN_POSITION_RE = re.compile(r"^\d+$")


def coding_to_genomic(hgvs_c: str, coordinates: CoordinateTranslator) -> str | None:
    """The genomic HGVS of a coding change, or ``None`` when no single genomic change describes it.

    :raises Exception: if the projection fails for a reason other than the variant itself
    """
    try:
        literal = coordinates.c_to_g_literal(hgvs_c)
        if not _spans_match(transcript_hgvs=hgvs_c, genomic_hgvs=literal, coordinates=coordinates):
            return None

        projected = coordinates.c_to_g(hgvs_c)
    except Exception as exc:
        if _is_definitive(exc):
            return None
        raise

    return _inversion_as_delins(projected, literal)


def genomic_to_coding(hgvs_g: str, transcript: str, coordinates: CoordinateTranslator) -> str | None:
    """The coding HGVS of a genomic change on ``transcript``, or ``None`` when no single coding change
    describes it.

    :raises Exception: if the projection fails for a reason other than the variant itself
    """
    try:
        literal = coordinates.g_to_c_literal(hgvs_g, transcript)
        if not _spans_match(transcript_hgvs=literal, genomic_hgvs=hgvs_g, coordinates=coordinates):
            return None

        projected = coordinates.g_to_c(hgvs_g, transcript)
    except Exception as exc:
        if _is_definitive(exc):
            return None
        raise

    return _inversion_as_delins(projected, literal)


def _is_definitive(exc: Exception) -> bool:
    """Whether hgvs rejected the variant, rather than failing to reach its reference data.

    A lost UTA connection is a bare ``HGVSError``, and hgvs's sequence fetchers report any fetch failure,
    an outage included, as ``HGVSDataNotAvailableError``, so neither counts.
    """
    # Imported here so the translation package stays importable without hgvs, as the rest of it is.
    import hgvs.exceptions  # type: ignore[import-untyped]

    return isinstance(
        exc,
        (
            hgvs.exceptions.HGVSInvalidIntervalError,
            hgvs.exceptions.HGVSInvalidVariantError,
            hgvs.exceptions.HGVSNormalizationError,
            hgvs.exceptions.HGVSParseError,
            hgvs.exceptions.HGVSUnsupportedOperationError,
        ),
    )


def _spans_match(*, transcript_hgvs: str, genomic_hgvs: str, coordinates: CoordinateTranslator) -> bool:
    """Whether a coding change and its unnormalized genomic projection cover the same number of bases.

    Compared in transcript (``n.``) coordinates, which count the bases between two exonic positions with
    introns excluded. A position with an intronic offset is anchored on the genome already, so a change
    with one is taken as written.
    """
    genomic = _LOCATED_EDIT_RE.match(genomic_hgvs)
    located = _LOCATED_EDIT_RE.match(transcript_hgvs)
    # A single position cannot straddle an exon junction.
    if genomic is None or located is None or located["end"] is None:
        return True

    transcript = _LOCATED_EDIT_RE.match(coordinates.c_to_n(transcript_hgvs))
    if transcript is None:
        return True

    n_start, n_end = transcript["start"], transcript["end"] or transcript["start"]
    if not (_PLAIN_POSITION_RE.match(n_start) and _PLAIN_POSITION_RE.match(n_end)):
        return True

    g_start, g_end = int(genomic["start"]), int(genomic["end"] or genomic["start"])
    return abs(g_end - g_start) == abs(int(n_end) - int(n_start))


def _inversion_as_delins(projected: str, literal: str) -> str:
    """Rewrite an ``inv`` projection as the delins it stands for.

    The normalized ``inv`` carries the trimmed span; the unnormalized projection of the same change carries
    the inserted bases over a span that may include unchanged flanking bases. Cutting those bases to the
    trimmed span gives the delins normalization would have written. Where that cut is not possible (a UTR or
    intronic position), the unnormalized delins describes the same change and is returned as is.
    """
    inv = _LOCATED_EDIT_RE.match(projected)
    if inv is None or inv["edit"] != "inv":
        return projected

    delins = _LOCATED_EDIT_RE.match(literal)
    if delins is None or delins["prefix"] != inv["prefix"] or not delins["edit"].startswith("delins"):
        return literal

    positions = (inv["start"], inv["end"] or inv["start"], delins["start"], delins["end"] or delins["start"])
    if not all(_PLAIN_POSITION_RE.match(position) for position in positions):
        return literal

    start, end, literal_start, literal_end = (int(position) for position in positions)
    alt = delins["edit"].removeprefix("delins")
    if len(alt) != literal_end - literal_start + 1 or not literal_start <= start <= end <= literal_end:
        return literal

    return f"{inv['prefix']}{start}_{end}delins{alt[start - literal_start : end - literal_start + 1]}"
