"""variant_annotation.translation — public API.

Composition roots (pipeline/, mavedb-api) import from here.
"""

from ._core import construct_equivalent_variants, construct_one
from ._ports import CoordinateTranslator, TranscriptSource
from .types import (
    TranslationConfig,
    TranslationError,
    TranslationResult,
    VariantInput,
    WtCodonMode,
)

__all__ = [
    "construct_equivalent_variants",
    "construct_one",
    "CoordinateTranslator",
    "TranscriptSource",
    "TranslationConfig",
    "TranslationError",
    "TranslationResult",
    "VariantInput",
    "WtCodonMode",
]
