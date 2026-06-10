import pytest

from variant_annotation.lib.pipeline.reverse_translate_step import ColumnConfig

pytestmark = pytest.mark.unit


def test_column_config_defaults():
    cols = ColumnConfig()
    assert cols.mapped_hgvs_p == "mapped_hgvs_p"
    assert cols.mapped_hgvs_c == "mapped_hgvs_c"
    assert cols.reverse_translation_error == "reverse_translation_error"
