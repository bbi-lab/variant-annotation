import pytest

from variant_annotation.lib.translation.types import TranslationConfig, WtCodonMode

pytestmark = pytest.mark.unit


def test_translation_config_defaults():
    config = TranslationConfig()
    assert config.wt_codon_mode is WtCodonMode.NONE
    assert config.assembly == "GRCh38"
    assert config.include_indels is False
    assert config.strict_ref_aa is True


def test_translation_config_rejects_wt_codon_mode_without_indels():
    with pytest.raises(ValueError, match="include_indels=True"):
        TranslationConfig(wt_codon_mode=WtCodonMode.UNAMBIGUOUS, include_indels=False)
