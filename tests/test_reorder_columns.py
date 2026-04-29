import csv

import pytest

from src.reorder_columns import reorder_columns


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_reorder_columns_moves_requested_columns_first(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    _write_tsv(
        in_path,
        rows=[{"id": "1", "gene": "A", "value": "x"}],
        fieldnames=["id", "gene", "value"],
    )

    n = reorder_columns(str(in_path), str(out_path), ["value", "id"])

    assert n == 1
    rows = _read_tsv(out_path)
    assert list(rows[0].keys()) == ["value", "id", "gene"]


def test_reorder_columns_missing_column_raises(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    _write_tsv(
        in_path,
        rows=[{"id": "1", "gene": "A"}],
        fieldnames=["id", "gene"],
    )

    with pytest.raises(ValueError, match="missing columns"):
        reorder_columns(str(in_path), str(out_path), ["value"])
