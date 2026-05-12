"""Tests for src/annotate_mavedb.py."""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
import requests

from src.annotate_mavedb import (
    OUTPUT_COLS,
    annotate_row,
    classify_score_range,
    classify_variant,
    fetch_calibration_variant_class_ids,
    fetch_calibrations,
    score_set_urn_from_variant_urn,
)


# ---------------------------------------------------------------------------
# score_set_urn_from_variant_urn
# ---------------------------------------------------------------------------


def test_score_set_urn_from_variant_urn_basic():
    assert score_set_urn_from_variant_urn("urn:mavedb:00000657-a-1#1") == "urn:mavedb:00000657-a-1"


def test_score_set_urn_from_variant_urn_large_id():
    assert score_set_urn_from_variant_urn("urn:mavedb:00000001-a-1#99999") == "urn:mavedb:00000001-a-1"


def test_score_set_urn_from_variant_urn_no_hash():
    assert score_set_urn_from_variant_urn("urn:mavedb:00000657-a-1") is None


def test_score_set_urn_from_variant_urn_empty():
    assert score_set_urn_from_variant_urn("") is None


def test_score_set_urn_from_variant_urn_hash_only():
    assert score_set_urn_from_variant_urn("#123") == ""


# ---------------------------------------------------------------------------
# classify_score_range
# ---------------------------------------------------------------------------


def _make_fc(label: str, lower, upper, inc_lower=True, inc_upper=False) -> dict[str, Any]:
    return {
        "id": 1,
        "label": label,
        "range": [lower, upper],
        "inclusiveLowerBound": inc_lower,
        "inclusiveUpperBound": inc_upper,
    }


def test_classify_score_range_matches_first():
    fcs = [
        _make_fc("Abnormal", None, 0.5),
        _make_fc("Functional", 0.5, None),
    ]
    assert classify_score_range(0.3, fcs) == "Abnormal"


def test_classify_score_range_matches_second():
    fcs = [
        _make_fc("Abnormal", None, 0.5),
        _make_fc("Functional", 0.5, None),
    ]
    # Default: lower inclusive, upper exclusive → 0.5 belongs to Functional
    assert classify_score_range(0.5, fcs) == "Functional"


def test_classify_score_range_exclusive_upper():
    fcs = [_make_fc("Tier1", 0.0, 1.0, inc_lower=True, inc_upper=False)]
    assert classify_score_range(1.0, fcs) == ""


def test_classify_score_range_inclusive_upper():
    fcs = [_make_fc("Tier1", 0.0, 1.0, inc_lower=True, inc_upper=True)]
    assert classify_score_range(1.0, fcs) == "Tier1"


def test_classify_score_range_no_match():
    fcs = [_make_fc("Band", 0.5, 0.9)]
    assert classify_score_range(0.95, fcs) == ""


def test_classify_score_range_skips_class_based():
    fcs = [
        {"id": 1, "label": "SomeClass", "range": None, "class": "CategoryA"},
        _make_fc("Functional", 0.8, None),
    ]
    assert classify_score_range(0.9, fcs) == "Functional"


def test_classify_score_range_empty_list():
    assert classify_score_range(0.5, []) == ""


# ---------------------------------------------------------------------------
# fetch_calibrations
# ---------------------------------------------------------------------------


def _mock_session(status_code: int, body: Any) -> MagicMock:
    resp = MagicMock()
    resp.status_code = status_code
    resp.json.return_value = body
    resp.raise_for_status = MagicMock()
    session = MagicMock()
    session.get.return_value = resp
    return session


def test_fetch_calibrations_success():
    payload = [{"urn": "urn:mavedb:cal-1", "title": "Cal", "primary": True}]
    session = _mock_session(200, payload)
    result = fetch_calibrations("https://api.mavedb.org", "urn:mavedb:00000001-a-1", session)
    assert result == payload
    session.get.assert_called_once_with(
        "https://api.mavedb.org/api/v1/score-calibrations/score-set/urn:mavedb:00000001-a-1",
        timeout=30,
    )


def test_fetch_calibrations_404():
    session = _mock_session(404, {"detail": "Not found"})
    result = fetch_calibrations("https://api.mavedb.org", "urn:mavedb:no-such", session)
    assert result == []
    session.get.return_value.raise_for_status.assert_not_called()


def test_fetch_calibrations_error_raises():
    resp = MagicMock()
    resp.status_code = 500
    resp.raise_for_status.side_effect = requests.HTTPError("Server error")
    session = MagicMock()
    session.get.return_value = resp
    with pytest.raises(requests.HTTPError):
        fetch_calibrations("https://api.mavedb.org", "urn:mavedb:bad", session)


# ---------------------------------------------------------------------------
# fetch_calibration_variant_class_ids
# ---------------------------------------------------------------------------


def test_fetch_calibration_variant_class_ids_basic():
    payload = [
        {
            "functionalClassificationId": 10,
            "variants": [
                {"urn": "urn:mavedb:00000001-a-1#1"},
                {"urn": "urn:mavedb:00000001-a-1#2"},
            ],
        },
        {
            "functionalClassificationId": 11,
            "variants": [{"urn": "urn:mavedb:00000001-a-1#3"}],
        },
    ]
    session = _mock_session(200, payload)
    result = fetch_calibration_variant_class_ids(
        "https://api.mavedb.org", "urn:mavedb:cal-1", session
    )
    assert result == {
        "urn:mavedb:00000001-a-1#1": 10,
        "urn:mavedb:00000001-a-1#2": 10,
        "urn:mavedb:00000001-a-1#3": 11,
    }


def test_fetch_calibration_variant_class_ids_empty():
    session = _mock_session(200, [])
    result = fetch_calibration_variant_class_ids("https://api.mavedb.org", "urn:mavedb:cal-x", session)
    assert result == {}


def test_fetch_calibration_variant_class_ids_missing_fc_id():
    payload = [{"variants": [{"urn": "urn:mavedb:00000001-a-1#1"}]}]
    session = _mock_session(200, payload)
    result = fetch_calibration_variant_class_ids("https://api.mavedb.org", "urn:mavedb:cal-x", session)
    assert result == {}


# ---------------------------------------------------------------------------
# classify_variant
# ---------------------------------------------------------------------------


def _make_range_calibration(cal_urn="urn:mavedb:cal-1", title="My Cal") -> dict[str, Any]:
    return {
        "urn": cal_urn,
        "title": title,
        "primary": True,
        "investigatorProvided": False,
        "functionalClassifications": [
            {
                "id": 1,
                "label": "Abnormal",
                "range": [None, 0.5],
                "inclusiveLowerBound": True,
                "inclusiveUpperBound": False,
            },
            {
                "id": 2,
                "label": "Functional",
                "range": [0.5, None],
                "inclusiveLowerBound": True,
                "inclusiveUpperBound": False,
            },
        ],
    }


def _make_class_calibration(cal_urn="urn:mavedb:cal-2", title="Class Cal") -> dict[str, Any]:
    return {
        "urn": cal_urn,
        "title": title,
        "primary": False,
        "investigatorProvided": True,
        "functionalClassifications": [
            {"id": 10, "label": "Pathogenic", "range": None, "class": "P"},
            {"id": 11, "label": "Benign", "range": None, "class": "B"},
        ],
    }


def test_classify_variant_range_based_match():
    cal = _make_range_calibration()
    session = MagicMock()
    cache: dict = {}
    urn, name, label = classify_variant(
        "urn:mavedb:00000001-a-1#5", "0.3", cal, "https://api.mavedb.org", session, cache
    )
    assert urn == "urn:mavedb:cal-1"
    assert name == "My Cal"
    assert label == "Abnormal"
    session.get.assert_not_called()


def test_classify_variant_range_based_no_match():
    cal = _make_range_calibration()
    session = MagicMock()
    cache: dict = {}
    _, _, label = classify_variant(
        "urn:mavedb:00000001-a-1#5", "0.75", cal, "https://api.mavedb.org", session, cache
    )
    assert label == "Functional"


def test_classify_variant_range_based_empty_score():
    cal = _make_range_calibration()
    session = MagicMock()
    cache: dict = {}
    _, _, label = classify_variant(
        "urn:mavedb:00000001-a-1#5", "", cal, "https://api.mavedb.org", session, cache
    )
    assert label == ""
    session.get.assert_not_called()


def test_classify_variant_range_based_non_numeric_score():
    cal = _make_range_calibration()
    session = MagicMock()
    cache: dict = {}
    _, _, label = classify_variant(
        "urn:mavedb:00000001-a-1#5", "NA", cal, "https://api.mavedb.org", session, cache
    )
    assert label == ""


def test_classify_variant_class_based_found(monkeypatch):
    cal = _make_class_calibration()
    session = MagicMock()
    cache: dict = {}

    monkeypatch.setattr(
        "src.annotate_mavedb.fetch_calibration_variant_class_ids",
        lambda api_url, cal_urn, sess: {"urn:mavedb:00000001-a-1#5": 10},
    )
    _, _, label = classify_variant(
        "urn:mavedb:00000001-a-1#5", "", cal, "https://api.mavedb.org", session, cache
    )
    assert label == "Pathogenic"


def test_classify_variant_class_based_not_found(monkeypatch):
    cal = _make_class_calibration()
    session = MagicMock()
    cache: dict = {}

    monkeypatch.setattr(
        "src.annotate_mavedb.fetch_calibration_variant_class_ids",
        lambda api_url, cal_urn, sess: {},
    )
    _, _, label = classify_variant(
        "urn:mavedb:00000001-a-1#99", "", cal, "https://api.mavedb.org", session, cache
    )
    assert label == ""


def test_classify_variant_class_based_uses_cache(monkeypatch):
    cal = _make_class_calibration()
    session = MagicMock()
    # Pre-populate the cache so the fetch function is never called.
    cache = {"urn:mavedb:cal-2": {"urn:mavedb:00000001-a-1#5": 11}}

    fetch_called = []
    monkeypatch.setattr(
        "src.annotate_mavedb.fetch_calibration_variant_class_ids",
        lambda *a, **kw: fetch_called.append(True) or {},
    )
    _, _, label = classify_variant(
        "urn:mavedb:00000001-a-1#5", "", cal, "https://api.mavedb.org", session, cache
    )
    assert label == "Benign"
    assert fetch_called == []  # cache hit → no fetch


def test_classify_variant_no_functional_classifications():
    cal = {"urn": "urn:mavedb:cal-empty", "title": "Empty", "functionalClassifications": []}
    session = MagicMock()
    cache: dict = {}
    urn, name, label = classify_variant(
        "urn:mavedb:00000001-a-1#1", "0.5", cal, "https://api.mavedb.org", session, cache
    )
    assert urn == "urn:mavedb:cal-empty"
    assert name == "Empty"
    assert label == ""


def test_classify_variant_api_error_class_based(monkeypatch):
    cal = _make_class_calibration()
    session = MagicMock()
    cache: dict = {}

    monkeypatch.setattr(
        "src.annotate_mavedb.fetch_calibration_variant_class_ids",
        MagicMock(side_effect=requests.RequestException("timeout")),
    )
    _, _, label = classify_variant(
        "urn:mavedb:00000001-a-1#5", "", cal, "https://api.mavedb.org", session, cache
    )
    assert label == ""
    # Cache populated with empty dict so subsequent calls skip the fetch.
    assert cache["urn:mavedb:cal-2"] == {}


# ---------------------------------------------------------------------------
# annotate_row
# ---------------------------------------------------------------------------


def _make_calibrations(*, primary=True, investigator=True) -> list[dict[str, Any]]:
    cals = []
    if primary:
        cals.append(
            {
                "urn": "urn:mavedb:cal-primary",
                "title": "Primary Cal",
                "primary": True,
                "investigatorProvided": False,
                "functionalClassifications": [
                    {
                        "id": 1,
                        "label": "Abnormal",
                        "range": [None, 0.5],
                        "inclusiveLowerBound": True,
                        "inclusiveUpperBound": False,
                    },
                    {
                        "id": 2,
                        "label": "Functional",
                        "range": [0.5, None],
                        "inclusiveLowerBound": True,
                        "inclusiveUpperBound": False,
                    },
                ],
            }
        )
    if investigator:
        cals.append(
            {
                "urn": "urn:mavedb:cal-inv",
                "title": "Investigator Cal",
                "primary": False,
                "investigatorProvided": True,
                "functionalClassifications": [
                    {
                        "id": 10,
                        "label": "Loss of function",
                        "range": [None, 0.3],
                        "inclusiveLowerBound": True,
                        "inclusiveUpperBound": False,
                    },
                    {
                        "id": 11,
                        "label": "Normal function",
                        "range": [0.3, None],
                        "inclusiveLowerBound": True,
                        "inclusiveUpperBound": False,
                    },
                ],
            }
        )
    return cals


def _annotate_row_with_cache(
    row: dict[str, str],
    calibrations: list[dict[str, Any]],
    *,
    variant_urn_col: str = "variant_urn",
    score_col: str = "score",
    class_id_cache: dict | None = None,
) -> dict[str, str]:
    session = MagicMock()
    ss_urn = score_set_urn_from_variant_urn(row.get(variant_urn_col, ""))
    cal_cache = {ss_urn: calibrations} if ss_urn else {}
    return annotate_row(
        row,
        api_url="https://api.mavedb.org",
        variant_urn_col=variant_urn_col,
        score_col=score_col,
        session=session,
        calibration_cache=cal_cache,
        class_id_cache=class_id_cache or {},
    )


_SCORE_SET_URL = "https://mavedb.org/score-sets/urn:mavedb:00000001-a-1"


def test_annotate_row_both_calibrations():
    row = {"variant_urn": "urn:mavedb:00000001-a-1#5", "score": "0.2"}
    result = _annotate_row_with_cache(row, _make_calibrations())
    assert result["mavedb.primary_calibration.urn"] == "urn:mavedb:cal-primary"
    assert result["mavedb.primary_calibration.name"] == "Primary Cal"
    assert result["mavedb.primary_calibration.url"] == _SCORE_SET_URL
    assert result["mavedb.primary_calibration.functional_class"] == "Abnormal"
    assert result["mavedb.investigator_provided_calibration.urn"] == "urn:mavedb:cal-inv"
    assert result["mavedb.investigator_provided_calibration.name"] == "Investigator Cal"
    assert result["mavedb.investigator_provided_calibration.url"] == _SCORE_SET_URL
    assert result["mavedb.investigator_provided_calibration.functional_class"] == "Loss of function"


def test_annotate_row_high_score():
    row = {"variant_urn": "urn:mavedb:00000001-a-1#5", "score": "0.9"}
    result = _annotate_row_with_cache(row, _make_calibrations())
    assert result["mavedb.primary_calibration.functional_class"] == "Functional"
    assert result["mavedb.investigator_provided_calibration.functional_class"] == "Normal function"


def test_annotate_row_no_primary_calibration():
    # When no explicit primary exists, falls back to the investigator-provided calibration.
    row = {"variant_urn": "urn:mavedb:00000001-a-1#5", "score": "0.2"}
    result = _annotate_row_with_cache(row, _make_calibrations(primary=False))
    assert result["mavedb.primary_calibration.urn"] == "urn:mavedb:cal-inv"
    assert result["mavedb.primary_calibration.name"] == "Investigator Cal"
    assert result["mavedb.primary_calibration.url"] == _SCORE_SET_URL
    assert result["mavedb.primary_calibration.functional_class"] == "Loss of function"
    assert result["mavedb.investigator_provided_calibration.url"] == _SCORE_SET_URL
    assert result["mavedb.investigator_provided_calibration.functional_class"] == "Loss of function"


def test_annotate_row_fallback_to_non_research_use_only():
    # When no explicit primary and no investigator-provided calibration, falls back
    # to the first non-research-use-only calibration.
    non_ruo_cal = {
        "urn": "urn:mavedb:cal-standard",
        "title": "Standard Cal",
        "primary": False,
        "investigatorProvided": False,
        "researchUseOnly": False,
        "functionalClassifications": [
            {
                "id": 1,
                "label": "Abnormal",
                "range": [None, 0.5],
                "inclusiveLowerBound": True,
                "inclusiveUpperBound": False,
            },
        ],
    }
    ruo_cal = {
        "urn": "urn:mavedb:cal-ruo",
        "title": "Research Only Cal",
        "primary": False,
        "investigatorProvided": False,
        "researchUseOnly": True,
        "functionalClassifications": [],
    }
    row = {"variant_urn": "urn:mavedb:00000001-a-1#5", "score": "0.2"}
    result = _annotate_row_with_cache(row, [ruo_cal, non_ruo_cal])
    # ruo_cal is skipped; non_ruo_cal is used as the primary fallback.
    assert result["mavedb.primary_calibration.urn"] == "urn:mavedb:cal-standard"
    assert result["mavedb.primary_calibration.name"] == "Standard Cal"
    assert result["mavedb.primary_calibration.url"] == _SCORE_SET_URL
    assert result["mavedb.primary_calibration.functional_class"] == "Abnormal"
    # No investigator-provided calibration.
    assert result["mavedb.investigator_provided_calibration.urn"] == ""
    assert result["mavedb.investigator_provided_calibration.url"] == ""


def test_annotate_row_no_investigator_calibration():
    row = {"variant_urn": "urn:mavedb:00000001-a-1#5", "score": "0.8"}
    result = _annotate_row_with_cache(row, _make_calibrations(investigator=False))
    assert result["mavedb.primary_calibration.url"] == _SCORE_SET_URL
    assert result["mavedb.primary_calibration.functional_class"] == "Functional"
    assert result["mavedb.investigator_provided_calibration.urn"] == ""
    assert result["mavedb.investigator_provided_calibration.name"] == ""
    assert result["mavedb.investigator_provided_calibration.url"] == ""
    assert result["mavedb.investigator_provided_calibration.functional_class"] == ""


def test_annotate_row_empty_variant_urn():
    row = {"variant_urn": "", "score": "0.5"}
    result = _annotate_row_with_cache(row, _make_calibrations())
    assert all(v == "" for v in result.values())


def test_annotate_row_variant_urn_no_hash():
    row = {"variant_urn": "urn:mavedb:00000001-a-1", "score": "0.5"}
    result = _annotate_row_with_cache(row, _make_calibrations())
    assert all(v == "" for v in result.values())


def test_annotate_row_no_calibrations():
    row = {"variant_urn": "urn:mavedb:00000001-a-1#5", "score": "0.5"}
    result = _annotate_row_with_cache(row, [])
    assert all(v == "" for v in result.values())


def test_annotate_row_empty_score_for_range_cal():
    row = {"variant_urn": "urn:mavedb:00000001-a-1#5", "score": ""}
    result = _annotate_row_with_cache(row, _make_calibrations(investigator=False))
    assert result["mavedb.primary_calibration.urn"] == "urn:mavedb:cal-primary"
    assert result["mavedb.primary_calibration.url"] == _SCORE_SET_URL
    assert result["mavedb.primary_calibration.functional_class"] == ""


def test_annotate_row_calibration_cache_hit():
    """Calibration list is fetched at most once per score set URN."""
    session = MagicMock()
    cal_cache: dict = {}
    class_id_cache: dict = {}

    calibrations = _make_calibrations()
    row = {"variant_urn": "urn:mavedb:00000001-a-1#1", "score": "0.4"}

    # First call: cache miss → would need session.get, but we pre-populate.
    cal_cache["urn:mavedb:00000001-a-1"] = calibrations

    annotate_row(
        row,
        api_url="https://api.mavedb.org",
        variant_urn_col="variant_urn",
        score_col="score",
        session=session,
        calibration_cache=cal_cache,
        class_id_cache=class_id_cache,
    )
    annotate_row(
        row,
        api_url="https://api.mavedb.org",
        variant_urn_col="variant_urn",
        score_col="score",
        session=session,
        calibration_cache=cal_cache,
        class_id_cache=class_id_cache,
    )
    session.get.assert_not_called()


def test_annotate_row_api_error_returns_empty():
    resp = MagicMock()
    resp.status_code = 503
    resp.raise_for_status.side_effect = requests.HTTPError("Service unavailable")
    session = MagicMock()
    session.get.return_value = resp

    row = {"variant_urn": "urn:mavedb:00000001-a-1#1", "score": "0.5"}
    result = annotate_row(
        row,
        api_url="https://api.mavedb.org",
        variant_urn_col="variant_urn",
        score_col="score",
        session=session,
        calibration_cache={},
        class_id_cache={},
    )
    assert all(v == "" for v in result.values())


def test_annotate_row_custom_col_names():
    row = {"my_urn": "urn:mavedb:00000001-a-1#5", "my_score": "0.2"}
    result = _annotate_row_with_cache(
        row,
        _make_calibrations(investigator=False),
        variant_urn_col="my_urn",
        score_col="my_score",
    )
    assert result["mavedb.primary_calibration.functional_class"] == "Abnormal"


# ---------------------------------------------------------------------------
# Integration: main() writes correct TSV
# ---------------------------------------------------------------------------


def test_main_writes_output(tmp_path, monkeypatch):
    calibrations = _make_calibrations(investigator=False)

    input_tsv = tmp_path / "input.tsv"
    input_tsv.write_text(
        "variant_urn\tscore\n"
        "urn:mavedb:00000001-a-1#1\t0.2\n"
        "urn:mavedb:00000001-a-1#2\t0.8\n",
        encoding="utf-8",
    )
    output_tsv = tmp_path / "output.tsv"

    monkeypatch.setattr(
        "src.annotate_mavedb.fetch_calibrations",
        lambda api_url, ss_urn, session: calibrations,
    )

    from src.annotate_mavedb import main

    main([str(input_tsv), str(output_tsv)])

    rows = list(csv.DictReader(output_tsv.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["mavedb.primary_calibration.functional_class"] == "Abnormal"
    assert rows[1]["mavedb.primary_calibration.functional_class"] == "Functional"
    assert rows[0]["mavedb.investigator_provided_calibration.functional_class"] == ""


def test_main_skip_and_limit(tmp_path, monkeypatch):
    calibrations = _make_calibrations(investigator=False)

    input_tsv = tmp_path / "input.tsv"
    lines = ["variant_urn\tscore\n"]
    for i in range(10):
        lines.append(f"urn:mavedb:00000001-a-1#{i}\t0.{i}\n")
    input_tsv.write_text("".join(lines), encoding="utf-8")

    output_tsv = tmp_path / "output.tsv"
    monkeypatch.setattr(
        "src.annotate_mavedb.fetch_calibrations",
        lambda api_url, ss_urn, session: calibrations,
    )

    from src.annotate_mavedb import main

    main([str(input_tsv), str(output_tsv), "--skip", "2", "--limit", "3"])

    rows = list(csv.DictReader(output_tsv.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 3
