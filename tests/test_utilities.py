import csv

from click.testing import CliRunner

from src.utilities import main


def _write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_utilities_filter_columns_command(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(
        in_path,
        rows=[{"a": "1", "b": "2", "c": "3"}],
        fieldnames=["a", "b", "c"],
    )

    runner = CliRunner()
    result = runner.invoke(
        main,
        ["filter-columns", str(in_path), str(out_path), "--keep-col", "a,c"],
    )

    assert result.exit_code == 0, result.output
    rows = _read_tsv(out_path)
    assert rows == [{"a": "1", "c": "3"}]


def test_utilities_merge_rows_command(tmp_path):
    f1 = tmp_path / "f1.tsv"
    f2 = tmp_path / "f2.tsv"
    out_path = tmp_path / "out.tsv"
    fieldnames = ["id", "value"]
    _write_tsv(f1, [{"id": "1", "value": "old"}], fieldnames)
    _write_tsv(f2, [{"id": "1", "value": "new"}, {"id": "2", "value": "x"}], fieldnames)

    runner = CliRunner()
    result = runner.invoke(
        main,
        ["merge-rows", str(out_path), str(f1), str(f2), "--key-col", "id"],
    )

    assert result.exit_code == 0, result.output
    rows = _read_tsv(out_path)
    assert rows == [{"id": "1", "value": "new"}, {"id": "2", "value": "x"}]


def test_utilities_filter_rows_command_any(tmp_path):
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

    runner = CliRunner()
    result = runner.invoke(
        main,
        ["filter-rows", str(in_path), str(out_path), "--value-col", "a,b", "--match", "any"],
    )

    assert result.exit_code == 0, result.output
    rows = _read_tsv(out_path)
    assert [r["id"] for r in rows] == ["1", "2"]


def test_utilities_filter_rows_command_all(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(
        in_path,
        rows=[
            {"id": "1", "a": "x", "b": "y"},
            {"id": "2", "a": "x", "b": ""},
        ],
        fieldnames=["id", "a", "b"],
    )

    runner = CliRunner()
    result = runner.invoke(
        main,
        ["filter-rows", str(in_path), str(out_path), "--value-col", "a,b", "--match", "all"],
    )

    assert result.exit_code == 0, result.output
    rows = _read_tsv(out_path)
    assert rows == [{"id": "1", "a": "x", "b": "y"}]
