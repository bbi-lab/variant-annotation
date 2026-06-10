"""Consumer-owned narrow protocols for translation's backing services.

These protocols are defined here (in the consumer) so that clients/ and
mavedb-api can satisfy them structurally without importing from translation/.
The arrows point down: _core imports these; clients/ does not.
"""

from __future__ import annotations

from typing import Protocol


class TranscriptSource(Protocol):
    """UTA-backed transcript facts."""

    def transcript_for_protein(self, protein_accession: str) -> str | None:
        """Return the preferred coding transcript for a RefSeq protein accession.

        Returns None when no transcript can be resolved.
        """
        ...

    def codon_at(self, transcript: str, aa_position: int) -> str | None:
        """Return the 3-letter uppercase codon at aa_position (1-based) in transcript.

        Returns None when the codon cannot be determined. Used only by WtCodonMode.ALL.
        """
        ...


class CoordinateTranslator(Protocol):
    """HGVS-mapper-backed coordinate projection."""

    def c_to_p(self, c_hgvs: str) -> str:
        """Project a coding HGVS string to its protein consequence."""
        ...

    def g_to_c(self, g_hgvs: str, transcript: str) -> str:
        """Project a genomic HGVS string to a coding HGVS string on transcript."""
        ...

    def c_to_g(self, c_hgvs: str) -> str:
        """Project a coding HGVS string to its genomic equivalent."""
        ...
