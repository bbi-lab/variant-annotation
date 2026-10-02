"""ProteinConsequence — the hub type the collapse/expand operation pivots on."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ProteinConsequence:
    """The protein-level representation of a variant, anchored to a transcript.

    hgvs_p is the protein HGVS string (e.g. "NP_000001.1:p.Ala123Val").
    transcript is the coding transcript needed to drive reverse translation
    (e.g. "NM_000001.1") — it is the transcript from which hgvs_p was derived
    or resolved via UTA.
    """

    hgvs_p: str
    transcript: str
