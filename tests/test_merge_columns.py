import csv

import pytest

from src.merge_columns import merge_columns, _parse_add_col_args, _parse_key_col_args


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_xlsx(path, rows, fieldnames, sheet_name="Sheet1"):
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet_name
    ws.append(fieldnames)
    for row in rows:
        ws.append([row.get(f, "") for f in fieldnames])
    wb.save(path)


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


# ---------------------------------------------------------------------------
# Column renaming
# ---------------------------------------------------------------------------

def test_parse_add_col_args_no_rename():
    cols, renames = _parse_add_col_args(("score", "effect"))
    assert cols == ["score", "effect"]
    assert renames == {}


def test_parse_add_col_args_with_rename():
    cols, renames = _parse_add_col_args(("Detects Splicing Variants?:splice_measure",))
    assert cols == ["Detects Splicing Variants?"]
    assert renames == {"Detects Splicing Variants?": "splice_measure"}


def test_parse_add_col_args_mixed():
    cols, renames = _parse_add_col_args(("Src Col:dest_col", "plain_col"))
    assert cols == ["Src Col", "plain_col"]
    assert renames == {"Src Col": "dest_col"}


def test_parse_add_col_args_comma_separated():
    cols, renames = _parse_add_col_args(("A:a_out,B,C:c_out",))
    assert cols == ["A", "B", "C"]
    assert renames == {"A": "a_out", "C": "c_out"}


def test_merge_columns_rename_column(tmp_path):
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.tsv"
    out = tmp_path / "out.tsv"

    _write_tsv(base, [{"id": "1", "gene": "A"}, {"id": "2", "gene": "B"}], ["id", "gene"])
    _write_tsv(
        extra,
        [{"id": "1", "Detects Splicing Variants?": "yes"}, {"id": "2", "Detects Splicing Variants?": "no"}],
        ["id", "Detects Splicing Variants?"],
    )

    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["id"],
        add_columns=["Detects Splicing Variants?"],
        column_renames={"Detects Splicing Variants?": "splice_measure"},
    )

    assert n == 2
    rows = _read_tsv(out)
    assert list(rows[0].keys()) == ["id", "gene", "splice_measure"]
    assert rows == [
        {"id": "1", "gene": "A", "splice_measure": "yes"},
        {"id": "2", "gene": "B", "splice_measure": "no"},
    ]


# ---------------------------------------------------------------------------
# Excel extra file
# ---------------------------------------------------------------------------

def test_merge_columns_excel_extra_default_sheet(tmp_path):
    pytest.importorskip("openpyxl")
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.xlsx"
    out = tmp_path / "out.tsv"

    _write_tsv(base, [{"dataset": "DS1", "score": "1.0"}, {"dataset": "DS2", "score": "2.0"}], ["dataset", "score"])
    _write_xlsx(
        extra,
        [{"dataset": "DS1", "Ensembl Transcript ID": "ENST001"}, {"dataset": "DS2", "Ensembl Transcript ID": "ENST002"}],
        ["dataset", "Ensembl Transcript ID"],
    )

    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["dataset"],
        add_columns=["Ensembl Transcript ID"],
    )

    assert n == 2
    rows = _read_tsv(out)
    assert rows == [
        {"dataset": "DS1", "score": "1.0", "Ensembl Transcript ID": "ENST001"},
        {"dataset": "DS2", "score": "2.0", "Ensembl Transcript ID": "ENST002"},
    ]


def test_merge_columns_excel_extra_named_sheet(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.xlsx"
    out = tmp_path / "out.tsv"

    _write_tsv(base, [{"id": "1"}], ["id"])

    # Workbook with two sheets; data is on the second one named "Curation"
    wb = openpyxl.Workbook()
    ws_other = wb.active
    ws_other.title = "Summary"
    ws_other.append(["id", "wrong_col"])
    ws_other.append(["1", "should_not_appear"])
    ws_curation = wb.create_sheet("Curation")
    ws_curation.append(["id", "RefSeq Transcript ID"])
    ws_curation.append(["1", "NM_007294.3"])
    wb.save(extra)

    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["id"],
        add_columns=["RefSeq Transcript ID"],
        extra_worksheet="Curation",
    )

    assert n == 1
    rows = _read_tsv(out)
    assert rows == [{"id": "1", "RefSeq Transcript ID": "NM_007294.3"}]


def test_merge_columns_excel_extra_index_sheet(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.xlsx"
    out = tmp_path / "out.tsv"

    _write_tsv(base, [{"id": "1"}], ["id"])

    wb = openpyxl.Workbook()
    ws1 = wb.active
    ws1.title = "Summary"
    ws1.append(["id", "wrong"])
    ws1.append(["1", "x"])
    ws2 = wb.create_sheet("Curation")
    ws2.append(["id", "value"])
    ws2.append(["1", "v1"])
    wb.save(extra)

    # Select second sheet by 1-based index
    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["id"],
        add_columns=["value"],
        extra_worksheet="2",
    )

    assert n == 1
    assert _read_tsv(out) == [{"id": "1", "value": "v1"}]


def test_merge_columns_excel_rename_and_worksheet(tmp_path):
    """End-to-end: Excel extra file, named worksheet, rename, and missing-key rows."""
    pytest.importorskip("openpyxl")
    base = tmp_path / "base.tsv"
    extra = tmp_path / "meta.xlsx"
    out = tmp_path / "out.tsv"

    _write_tsv(
        base,
        [{"dataset": "DS1", "score": "0.9"}, {"dataset": "DS_UNKNOWN", "score": "0.1"}],
        ["dataset", "score"],
    )
    _write_xlsx(
        extra,
        [{"dataset": "DS1", "Detects Splicing Variants?": "yes", "Ensembl Transcript ID": "ENST001"}],
        ["dataset", "Detects Splicing Variants?", "Ensembl Transcript ID"],
        sheet_name="Curation",
    )

    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["dataset"],
        add_columns=["Detects Splicing Variants?", "Ensembl Transcript ID"],
        column_renames={"Detects Splicing Variants?": "splice_measure"},
        extra_worksheet="Curation",
    )

    assert n == 2
    rows = _read_tsv(out)
    assert rows[0] == {"dataset": "DS1", "score": "0.9", "splice_measure": "yes", "Ensembl Transcript ID": "ENST001"}
    # Row with no match in extra gets empty strings
    assert rows[1] == {"dataset": "DS_UNKNOWN", "score": "0.1", "splice_measure": "", "Ensembl Transcript ID": ""}


# ---------------------------------------------------------------------------
# Per-side key columns
# ---------------------------------------------------------------------------

def test_parse_key_col_args_same_name():
    base_cols, extra_cols = _parse_key_col_args(("id",))
    assert base_cols == ["id"]
    assert extra_cols == ["id"]


def test_parse_key_col_args_different_name():
    base_cols, extra_cols = _parse_key_col_args(("dataset_name:Dataset Name",))
    assert base_cols == ["dataset_name"]
    assert extra_cols == ["Dataset Name"]


def test_parse_key_col_args_mixed_and_comma():
    base_cols, extra_cols = _parse_key_col_args(("a:A", "b,c:C"))
    assert base_cols == ["a", "b", "c"]
    assert extra_cols == ["A", "b", "C"]


def test_merge_columns_different_key_names(tmp_path):
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.tsv"
    out = tmp_path / "out.tsv"

    _write_tsv(
        base,
        [{"dataset_name": "DS1", "score": "1.0"}, {"dataset_name": "DS2", "score": "2.0"}],
        ["dataset_name", "score"],
    )
    _write_tsv(
        extra,
        [{"Dataset Name": "DS1", "transcript": "ENST001"}, {"Dataset Name": "DS2", "transcript": "ENST002"}],
        ["Dataset Name", "transcript"],
    )

    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["dataset_name"],
        extra_key_columns=["Dataset Name"],
        add_columns=["transcript"],
    )

    assert n == 2
    assert _read_tsv(out) == [
        {"dataset_name": "DS1", "score": "1.0", "transcript": "ENST001"},
        {"dataset_name": "DS2", "score": "2.0", "transcript": "ENST002"},
    ]


def test_merge_columns_different_key_no_match(tmp_path):
    """Rows in base with no match in extra get empty strings."""
    base = tmp_path / "base.tsv"
    extra = tmp_path / "extra.tsv"
    out = tmp_path / "out.tsv"

    _write_tsv(base, [{"ds": "X"}, {"ds": "Y"}], ["ds"])
    _write_tsv(extra, [{"name": "X", "val": "v1"}], ["name", "val"])

    n = merge_columns(
        base_file=str(base),
        extra_file=str(extra),
        output_file=str(out),
        key_columns=["ds"],
        extra_key_columns=["name"],
        add_columns=["val"],
    )

    assert n == 2
    assert _read_tsv(out) == [
        {"ds": "X", "val": "v1"},
        {"ds": "Y", "val": ""},
    ]
