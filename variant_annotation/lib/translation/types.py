"""Public types for the translation API."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class WtCodonMode(str, Enum):
    NONE = "none"
    UNAMBIGUOUS = "unambiguous"
    ALL = "all"


@dataclass
class VariantInput:
    """A single variant to construct equivalents for.

    hgvs may be a p., c., or g. HGVS string. transcript is an optional hint
    that short-circuits UTA lookup — required for g. inputs where the target
    transcript is ambiguous, optional for p./c. inputs.
    """

    hgvs: str
    transcript: str | None = None


@dataclass
class ProjectionPair:
    """One member of a protein consequence's coding/genomic equivalence class:
    a coding allele paired with its genomic projection.

    Reverse translation fans a protein consequence out into N coding candidates,
    each with a genomic projection. That coding↔genomic edge is deterministic —
    given ``(assembly, transcript)``, the genomic form of a coding candidate is
    fixed. (The only ambiguity, protein→codon, sits one level up and is carried
    by co-membership in the equivalence class, not here.) A ProjectionPair holds
    that deterministic pairing as one object: the coding allele, its genomic
    projection, and the variant type the reverse-translate-variants CLI assigned
    it. The CLI emits these three columns position-aligned per candidate (with
    ``--one-row-per-input``); this dataclass is that row-wise unit, so the
    pairing survives instead of being split into two lists that downstream must
    re-zip by index.

    Fields:
        hgvs_c: The coding (c.) HGVS expression. Always present — a pair is
            keyed by its coding candidate.
        hgvs_g: The genomic (g.) HGVS projection, or ``None`` when the c→g
            projection failed (e.g. an intronic candidate the mapper cannot
            resolve). ``None`` is a well-formed one-sided pair, NOT a
            list-desync — this is the case the old two-list shape silently
            dropped, misaligning every candidate after it.
        variant_type: The CLI's classification of the candidate ("snv",
            "insertion", "deletion", "delins"), or ``None`` if unreported.

    This is the consumable unit the mavedb-api reverse-translation worker
    iterates to stamp the projection pairing (``projection_group``) onto its
    derived allele links.
    """

    hgvs_c: str
    hgvs_g: str | None = None
    variant_type: str | None = None


@dataclass
class TranslationResult:
    """Successful equivalence-class construction for one VariantInput.

    projection_pairs holds the coding/genomic equivalence class as
    ProjectionPair objects, preserving the deterministic coding↔genomic pairing
    the reverse-translate-variants CLI emits. hgvs_p is the shared protein apex
    (None for protein-assay inputs, where the protein allele is authoritative
    rather than derived).
    """

    input: VariantInput
    projection_pairs: list[ProjectionPair] = field(default_factory=list)
    hgvs_p: str | None = None


@dataclass
class TranslationError:
    """Failed or skipped equivalence-class construction for one VariantInput."""

    input: VariantInput
    error: str


@dataclass
class TranslationConfig:
    """Translation behaviour knobs — everything except the port implementations.

    Passed as a single object through construct_equivalent_variants and
    process_rows so the same defaults live in one place and function signatures
    stay narrow.
    """

    wt_codon_mode: WtCodonMode = WtCodonMode.NONE
    assembly: str = "GRCh38"
    include_indels: bool = False
    max_indel_size: int = 3
    strict_ref_aa: bool = True
    use_inv_notation: bool = False
    allow_length_changing_stop_candidates: bool = True

    def __post_init__(self) -> None:
        if self.wt_codon_mode is not WtCodonMode.NONE and not self.include_indels:
            raise ValueError(
                "wt_codon_mode requires include_indels=True, "
                "because WT codon delins variants are intra-codon insertions/deletions."
            )
