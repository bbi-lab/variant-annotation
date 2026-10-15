# Copyright (C) 2026 Jeremy Stone
# SPDX-License-Identifier: AGPL-3.0-or-later

"""variant_annotation.translation — public API.

Composition roots (pipeline/, mavedb-api) import from here.
"""

from ._core import construct_equivalent_variants, construct_one
from ._ports import CoordinateTranslator, TranscriptSource
from .projection import coding_to_genomic, genomic_to_coding
from .types import (
    ProjectionPair,
    TranslationConfig,
    TranslationError,
    TranslationErrorReason,
    TranslationResult,
    VariantInput,
    WtCodonMode,
)

__all__ = [
    "coding_to_genomic",
    "construct_equivalent_variants",
    "construct_one",
    "genomic_to_coding",
    "ProjectionPair",
    "CoordinateTranslator",
    "TranscriptSource",
    "TranslationConfig",
    "TranslationError",
    "TranslationErrorReason",
    "TranslationResult",
    "VariantInput",
    "WtCodonMode",
]
