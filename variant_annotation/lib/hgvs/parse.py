"""HGVS parsing utilities — shared vocabulary, no I/O.

reverse_complement lives here permanently.

parse_hgvs and apply_vcf_anchor are re-exported here from add_vcf_identifiers
as a transitional step; they will be migrated into this module when
add_vcf_identifiers.py is refactored to import downward from variant_annotation.
"""

from __future__ import annotations

# TODO: migrate _parse_hgvs and _apply_vcf_anchor bodies here from
# add_vcf_identifiers.py (along with their dependencies). Once done, update
# add_vcf_identifiers.py to import from variant_annotation.hgvs.parse instead.
from add_vcf_identifiers import (
    _apply_vcf_anchor as apply_vcf_anchor,
    _parse_hgvs as parse_hgvs,
)


def reverse_complement(seq: str) -> str:
    complement_table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement_table)[::-1]


__all__ = ["apply_vcf_anchor", "parse_hgvs", "reverse_complement"]
