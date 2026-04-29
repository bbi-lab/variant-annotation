import csv

import pytest

from src.merge_columns import merge_columns


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_merge_columns_add_selected_columns(tmp_path):
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.tsv"
    out = tmp_path / "out.tsv"

    _write_tsv(base, [{"id": "1", "gene": "A"}, {"id": "2", "gene": "B"}], ["id", "gene"])
    _write_tsv(extra, [{"id": "1", "score": "1.5"}, {"id": "3", "score": "9.9"}], ["id", "score"])

    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["id"],
        add_columns=["score"],
    )

    assert n == 2
    assert _read_tsv(out) == [
        {"id": "1", "gene": "A", "score": "1.5"},
        {"id": "2", "gene": "B", "score": ""},
    ]


def test_merge_columns_add_all_new_non_key_columns(tmp_path):
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.tsv"
    out = tmp_path / "out.tsv"

    _write_tsv(base, [{"id": "1", "gene": "A"}], ["id", "gene"])
    _write_tsv(extra, [{"id": "1", "score": "1.5", "effect": "high"}], ["id", "score", "effect"])

    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["id"],
        add_columns=[],
        add_all_cols_from_extra=True,
    )

    assert n == 1
    assert _read_tsv(out) == [{"id": "1", "gene": "A", "score": "1.5", "effect": "high"}]


def test_merge_columns_missing_key_raises(tmp_path):
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.tsv"
    out = tmp_path / "out.tsv"

    _write_tsv(base, [{"id": "1", "gene": "A"}], ["id", "gene"])
    _write_tsv(extra, [{"id": "1", "score": "1.5"}], ["id", "score"])

    with pytest.raises(ValueError, match="At least one key column"):
        merge_columns(
            base_file=str(base),
            extra_file=str(extra),
            output_file=str(out),
            key_columns=[],
            add_columns=["score"],
        )
