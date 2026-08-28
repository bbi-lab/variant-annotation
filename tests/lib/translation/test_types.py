import pytest

from variant_annotation.lib.translation.types import (
    TranslationConfig,
    TranslationError,
    TranslationErrorReason,
    VariantInput,
    WtCodonMode,
)

pytestmark = pytest.mark.unit


def test_translation_error_defaults_to_failed():
    # A bare TranslationError is a genuine failure unless a caller says otherwise, so the
    # existing (pre-typed-reason) construction sites keep meaning "failed".
    err = TranslationError(input=VariantInput(hgvs="NP_1:p.Ala1Val"), error="boom")
    assert err.reason is TranslationErrorReason.FAILED


def test_translation_error_reason_members():
    assert {r.value for r in TranslationErrorReason} == {"not_translatable", "failed"}


def test_translation_config_defaults():
    config = TranslationConfig()
    assert config.wt_codon_mode is WtCodonMode.NONE
    assert config.assembly == "GRCh38"
    assert config.include_indels is False
    assert config.strict_ref_aa is True


def test_translation_config_rejects_wt_codon_mode_without_indels():
    with pytest.raises(ValueError, match="include_indels=True"):
        TranslationConfig(wt_codon_mode=WtCodonMode.UNAMBIGUOUS, include_indels=False)
