import csv

import pytest

from src.merge_rows import merge_rows


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_merge_rows_replaces_and_appends(tmp_path):
    f1 = tmp_path / "f1.tsv"
    f2 = tmp_path / "f2.tsv"
    out = tmp_path / "out.tsv"
    fieldnames = ["id", "gene", "value"]

    _write_tsv(
        f1,
        [
            {"id": "1", "gene": "A", "value": "old"},
            {"id": "2", "gene": "B", "value": "keep"},
        ],
        fieldnames,
    )
    _write_tsv(
        f2,
        [
            {"id": "1", "gene": "A", "value": "new"},
            {"id": "3", "gene": "C", "value": "append"},
        ],
        fieldnames,
    )

    n = merge_rows([str(f1), str(f2)], str(out), ["id", "gene"])

    assert n == 3
    rows = _read_tsv(out)
    assert rows == [
        {"id": "1", "gene": "A", "value": "new"},
        {"id": "2", "gene": "B", "value": "keep"},
        {"id": "3", "gene": "C", "value": "append"},
    ]


def test_merge_rows_last_file_wins_across_multiple_files(tmp_path):
    f1 = tmp_path / "f1.tsv"
    f2 = tmp_path / "f2.tsv"
    f3 = tmp_path / "f3.tsv"
    out = tmp_path / "out.tsv"
    fieldnames = ["id", "value"]

    _write_tsv(f1, [{"id": "1", "value": "v1"}], fieldnames)
    _write_tsv(f2, [{"id": "1", "value": "v2"}], fieldnames)
    _write_tsv(f3, [{"id": "1", "value": "v3"}], fieldnames)

    n = merge_rows([str(f1), str(f2), str(f3)], str(out), ["id"])

    assert n == 1
    rows = _read_tsv(out)
    assert rows == [{"id": "1", "value": "v3"}]


def test_merge_rows_missing_key_column_raises(tmp_path):
    f1 = tmp_path / "f1.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(f1, [{"id": "1", "value": "x"}], ["id", "value"])

    with pytest.raises(ValueError, match="missing key columns"):
        merge_rows([str(f1)], str(out), ["id", "gene"])


def test_merge_rows_mismatched_columns_raises(tmp_path):
    f1 = tmp_path / "f1.tsv"
    f2 = tmp_path / "f2.tsv"
    out = tmp_path / "out.tsv"

    _write_tsv(f1, [{"id": "1", "value": "x"}], ["id", "value"])
    _write_tsv(f2, [{"id": "1", "other": "y"}], ["id", "other"])

    with pytest.raises(ValueError, match="columns do not match"):
        merge_rows([str(f1), str(f2)], str(out), ["id"])


def test_merge_rows_duplicate_keys_within_one_file_last_row_wins(tmp_path):
    f1 = tmp_path / "f1.tsv"
    out = tmp_path / "out.tsv"
    fieldnames = ["id", "value"]

    _write_tsv(
        f1,
        [
            {"id": "1", "value": "first"},
            {"id": "1", "value": "second"},
        ],
        fieldnames,
    )

    n = merge_rows([str(f1)], str(out), ["id"])

    assert n == 1
    rows = _read_tsv(out)
    assert rows == [{"id": "1", "value": "second"}]
