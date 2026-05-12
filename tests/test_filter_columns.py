import csv

import pytest

from src.filter_columns import filter_columns, rename_columns, _parse_keep_col_args


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


# ---------------------------------------------------------------------------
# _parse_keep_col_args
# ---------------------------------------------------------------------------

def test_parse_keep_col_args_plain():
    result = _parse_keep_col_args(("a", "b"))
    assert result == [("a", "a"), ("b", "b")]


def test_parse_keep_col_args_with_rename():
    result = _parse_keep_col_args(("old:new",))
    assert result == [("old", "new")]


def test_parse_keep_col_args_comma_separated():
    result = _parse_keep_col_args(("a:x,b,c:z",))
    assert result == [("a", "x"), ("b", "b"), ("c", "z")]


def test_parse_keep_col_args_deduplicates_src():
    result = _parse_keep_col_args(("a:x", "a:y"))  # second 'a' is dropped
    assert result == [("a", "x")]


# ---------------------------------------------------------------------------
# rename_columns
# ---------------------------------------------------------------------------

def _write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_rename_columns_keep_no_rename(tmp_path):
    """Keep mode without rename preserves input order for selected columns."""
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, [{"a": "1", "b": "2", "c": "3"}], ["a", "b", "c"])

    n = rename_columns(str(in_path), str(out_path), keep_specs=[("c", "c"), ("a", "a")])

    assert n == 1
    rows = _read_tsv(out_path)
    assert list(rows[0].keys()) == ["a", "c"]  # input order, not keep order
    assert rows[0] == {"a": "1", "c": "3"}


def test_rename_columns_keep_with_rename(tmp_path):
    """SRC:DEST renames the column header in the output."""
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, [{"gene": "BRCA1", "score": "0.9"}], ["gene", "score"])

    n = rename_columns(
        str(in_path), str(out_path),
        keep_specs=[("gene", "gene_symbol"), ("score", "revel_score")],
    )

    assert n == 1
    rows = _read_tsv(out_path)
    assert list(rows[0].keys()) == ["gene_symbol", "revel_score"]
    assert rows[0]["gene_symbol"] == "BRCA1"
    assert rows[0]["revel_score"] == "0.9"


def test_rename_columns_keep_reorder(tmp_path):
    """reorder=True outputs columns in keep_specs order, not input order."""
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, [{"a": "1", "b": "2", "c": "3"}], ["a", "b", "c"])

    rename_columns(
        str(in_path), str(out_path),
        keep_specs=[("c", "col_c"), ("a", "col_a")],
        reorder=True,
    )

    rows = _read_tsv(out_path)
    assert list(rows[0].keys()) == ["col_c", "col_a"]
    assert rows[0]["col_c"] == "3"
    assert rows[0]["col_a"] == "1"


def test_rename_columns_omit_mode(tmp_path):
    """Omit mode drops the listed columns, names unchanged."""
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, [{"a": "1", "b": "2", "c": "3"}], ["a", "b", "c"])

    n = rename_columns(str(in_path), str(out_path), omit_cols=["b"])

    assert n == 1
    rows = _read_tsv(out_path)
    assert list(rows[0].keys()) == ["a", "c"]
    assert rows[0] == {"a": "1", "c": "3"}


def test_rename_columns_both_modes_raises(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, [{"a": "1"}], ["a"])

    with pytest.raises(ValueError, match="exactly one mode"):
        rename_columns(
            str(in_path), str(out_path),
            keep_specs=[("a", "a")],
            omit_cols=["a"],
        )


def test_rename_columns_missing_keep_raises(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, [{"a": "1"}], ["a"])

    with pytest.raises(ValueError, match="keep columns not found"):
        rename_columns(str(in_path), str(out_path), keep_specs=[("missing", "missing")])

