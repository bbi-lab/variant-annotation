"""Pure accession-classification utilities — no I/O, no external dependencies."""

from __future__ import annotations

TRANSCRIPT_ACCESSION_PREFIXES = ("NM_", "NR_", "XM_", "XR_", "ENST", "LRG_")
REFSEQ_TRANSCRIPT_ACCESSION_PREFIXES = ("NM_", "NR_", "XM_", "XR_")
REFSEQ_PROTEIN_ACCESSION_PREFIXES = ("NP_", "XP_", "YP_", "WP_")


def extract_accession(hgvs_string: str) -> str:
    token = (hgvs_string or "").strip()
    if ":" not in token:
        return ""

    return token.split(":", 1)[0].strip()


def looks_like_transcript_accession(accession: str) -> bool:
    return accession.startswith(TRANSCRIPT_ACCESSION_PREFIXES)


def looks_like_refseq_transcript_accession(accession: str) -> bool:
    """Whether ``accession`` is a RefSeq transcript (``NM_``/``NR_``/``XM_``/``XR_``).

    Ensembl serves RefSeq and Ensembl transcript sets separately, so a caller must know which set an
    accession belongs to before it can ask VEP for the matching one (``refseq=1``).
    """
    return accession.startswith(REFSEQ_TRANSCRIPT_ACCESSION_PREFIXES)


def looks_like_refseq_protein_accession(accession: str) -> bool:
    return accession.startswith(REFSEQ_PROTEIN_ACCESSION_PREFIXES)


def strip_accession_version(accession: str) -> str:
    """Return ``accession`` without its trailing ``.<version>`` suffix.

    ``'NM_000049.4'`` -> ``'NM_000049'``; ``'ENST00000366667.8'`` -> ``'ENST00000366667'``.

    Only a trailing all-digit segment is treated as a version, so an accession with no version, or one
    whose final dotted segment is not numeric, is returned unchanged.
    """
    stem, dot, version = (accession or "").rpartition(".")
    if dot and version.isdigit():
        return stem
    return accession or ""


def transcript_sort_key(accession: str) -> tuple[int, int, str]:
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
