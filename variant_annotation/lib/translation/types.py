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
class TranslationResult:
    """Successful equivalence-class construction for one VariantInput."""

    input: VariantInput
    hgvs_c_candidates: list[str] = field(default_factory=list)
    hgvs_g_candidates: list[str] = field(default_factory=list)
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
