import csv

import pytest

from src.filter_rows import filter_rows


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_filter_rows_any_mode(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(
        in_path,
        rows=[
            {"id": "1", "a": "x", "b": ""},
            {"id": "2", "a": "", "b": "y"},
            {"id": "3", "a": "", "b": ""},
        ],
        fieldnames=["id", "a", "b"],
    )

    n = filter_rows(str(in_path), str(out_path), ["a", "b"], match="any")

    assert n == 2
    rows = _read_tsv(out_path)
    assert [r["id"] for r in rows] == ["1", "2"]


def test_filter_rows_all_mode(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(
        in_path,
        rows=[
            {"id": "1", "a": "x", "b": "y"},
            {"id": "2", "a": "x", "b": ""},
            {"id": "3", "a": "", "b": "y"},
        ],
        fieldnames=["id", "a", "b"],
    )

    n = filter_rows(str(in_path), str(out_path), ["a", "b"], match="all")

    assert n == 1
    rows = _read_tsv(out_path)
    assert rows == [{"id": "1", "a": "x", "b": "y"}]


def test_filter_rows_missing_column_raises(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, rows=[{"id": "1", "a": "x"}], fieldnames=["id", "a"])

    with pytest.raises(ValueError, match="value columns not found"):
        filter_rows(str(in_path), str(out_path), ["missing"], match="any")


def test_filter_rows_no_columns_raises(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, rows=[{"id": "1"}], fieldnames=["id"])

    with pytest.raises(ValueError, match="At least one --value-col"):
        filter_rows(str(in_path), str(out_path), [], match="any")
