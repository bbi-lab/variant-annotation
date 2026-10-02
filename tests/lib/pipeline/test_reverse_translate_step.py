import pytest

from variant_annotation.lib.pipeline.reverse_translate_step import ColumnConfig
from variant_annotation.lib.translation.types import (
    ProjectionPair,
    TranslationResult,
    VariantInput,
)

pytestmark = pytest.mark.unit


def test_column_config_defaults():
    cols = ColumnConfig()
    assert cols.mapped_hgvs_p == "mapped_hgvs_p"
    assert cols.mapped_hgvs_c == "mapped_hgvs_c"
    assert cols.reverse_translation_error == "reverse_translation_error"


def test_write_result_keeps_columns_position_aligned():
    # The coding and genomic pipe-joined columns must stay the same length, with
    # an empty genomic cell where a candidate's projection failed (hgvs_g=None),
    # so downstream index-based pairing lines up candidate-for-candidate.
    cols = ColumnConfig()
    result = TranslationResult(
        input=VariantInput(hgvs="NP_1:p.Ala1Val"),
        projection_pairs=[
            ProjectionPair(hgvs_c="NM_1:c.1A>G", hgvs_g="NC_1:g.1A>G"),
            ProjectionPair(hgvs_c="NM_1:c.1A>C", hgvs_g=None),
            ProjectionPair(hgvs_c="NM_1:c.1A>T", hgvs_g="NC_1:g.1A>T"),
        ],
        hgvs_p="NP_1:p.Ala1Val",
    )
    row: dict = {}
    cols.write_result(row, result)

    coding = row[cols.mapped_hgvs_c].split("|")
    genomic = row[cols.mapped_hgvs_g].split("|")
    assert coding == ["NM_1:c.1A>G", "NM_1:c.1A>C", "NM_1:c.1A>T"]
    assert genomic == ["NC_1:g.1A>G", "", "NC_1:g.1A>T"]  # middle cell empty, alignment preserved
    assert len(coding) == len(genomic)
    assert row[cols.mapped_hgvs_p] == "NP_1:p.Ala1Val"
