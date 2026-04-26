import csv

import pytest

from src.filter_columns import filter_columns


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_keep_mode_writes_only_selected_columns(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(
        input_path,
        rows=[
            {"a": "1", "b": "2", "c": "3"},
            {"a": "4", "b": "5", "c": "6"},
        ],
        fieldnames=["a", "b", "c"],
    )

    n = filter_columns(
        input_file=str(input_path),
        output_file=str(output_path),
        keep_cols=["c", "a"],
    )

    assert n == 2
    out = _read_tsv(output_path)
    assert list(out[0].keys()) == ["a", "c"]  # input order preserved
    assert out[0]["a"] == "1"
    assert out[0]["c"] == "3"


def test_omit_mode_removes_selected_columns(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(
        input_path,
        rows=[{"a": "1", "b": "2", "c": "3"}],
        fieldnames=["a", "b", "c"],
    )

    n = filter_columns(
        input_file=str(input_path),
        output_file=str(output_path),
        omit_cols=["b"],
    )

    assert n == 1
    out = _read_tsv(output_path)
    assert list(out[0].keys()) == ["a", "c"]
    assert out[0] == {"a": "1", "c": "3"}


def test_keep_and_omit_mode_together_is_error(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(input_path, rows=[{"a": "1"}], fieldnames=["a"])

    with pytest.raises(ValueError, match="exactly one mode"):
        filter_columns(
            input_file=str(input_path),
            output_file=str(output_path),
            keep_cols=["a"],
            omit_cols=["a"],
        )


def test_missing_keep_column_is_error(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(input_path, rows=[{"a": "1"}], fieldnames=["a"])

    with pytest.raises(ValueError, match="keep columns not found"):
        filter_columns(
            input_file=str(input_path),
            output_file=str(output_path),
            keep_cols=["missing"],
        )


def test_missing_omit_column_is_error(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(input_path, rows=[{"a": "1"}], fieldnames=["a"])

    with pytest.raises(ValueError, match="omit columns not found"):
        filter_columns(
            input_file=str(input_path),
            output_file=str(output_path),
            omit_cols=["missing"],
        )


def test_omit_all_columns_is_error(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(input_path, rows=[{"a": "1"}], fieldnames=["a"])

    with pytest.raises(ValueError, match="remove all columns"):
        filter_columns(
            input_file=str(input_path),
            output_file=str(output_path),
            omit_cols=["a"],
        )
