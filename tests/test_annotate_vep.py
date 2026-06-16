"""Unit tests for src/annotate_vep.py."""

from __future__ import annotations

import csv
from datetime import date
from unittest import mock

import src.annotate_vep as mod

# Cache tuple shape: (most_severe, all_consequences, source, access_date)

def _hit(most_severe, source="most_severe", all_cons=None, cached_date=""):
    return (most_severe, all_cons, source, cached_date)


def _miss():
    return (None, None, "most_severe", "")


def _api_err(msg="network error"):
    return (None, None, f"api_error:{msg}", "")


def _vep_err(msg="Start (28695710) must be less than or equal to end+1 (28695243)"):
    return (None, None, f"vep_error:{msg}", "")


# ---------------------------------------------------------------------------
# annotate_row — column selection
# ---------------------------------------------------------------------------


def test_annotate_row_uses_first_non_blank_hgvs_col():
    """mapped_hgvs_c takes priority over mapped_hgvs_g when non-blank."""
    row = {
        "mapped_hgvs_c": "NM_000000.1:c.1A>T",
        "mapped_hgvs_g": "NC_000001.11:g.1A>T",
    }
    cache = {
        "NM_000000.1:c.1A>T": _hit("synonymous_variant"),
        "NC_000001.11:g.1A>T": _hit("missense_variant"),
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep",
        hgvs_cols=["mapped_hgvs_c", "mapped_hgvs_g", "mapped_hgvs_p"],
        access_date="2026-04-30",
    )
    assert out["vep.most_severe_mutational_consequence"] == "synonymous_variant"


def test_annotate_row_falls_back_when_first_col_blank():
    """Falls back to mapped_hgvs_g when mapped_hgvs_c is blank."""
    row = {
        "mapped_hgvs_c": "",
        "mapped_hgvs_g": "NC_000001.11:g.1A>T",
    }
    cache = {"NC_000001.11:g.1A>T": _hit("missense_variant")}
    out = mod.annotate_row(
        row, cache, col_prefix="vep",
        hgvs_cols=["mapped_hgvs_c", "mapped_hgvs_g", "mapped_hgvs_p"],
        access_date="2026-04-30",
    )
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant"


# ---------------------------------------------------------------------------
# annotate_row — multiple candidates (pipe-delimited)
# ---------------------------------------------------------------------------


def test_annotate_row_emits_pipe_aligned_values_for_matching_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    cache = {
        "NC_000001.11:g.1A>T": _hit("missense_variant"),
        "NC_000001.11:g.2C>G": _hit("missense_variant"),
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_g"], access_date="2026-04-30",
    )
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant|missense_variant"
    assert out["vep.access_date"] == "2026-04-30|2026-04-30"
    assert out["vep.error"] == "|"


def test_annotate_row_vep_miss_yields_blank_consequence_for_that_position():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    cache = {
        "NC_000001.11:g.1A>T": _hit("synonymous_variant"),
        "NC_000001.11:g.2C>G": _miss(),
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_g"], access_date="2026-04-30",
    )
    assert out["vep.most_severe_mutational_consequence"] == "synonymous_variant|"
    assert out["vep.error"] == "|"


def test_annotate_row_pipe_aligned_multiple_different_consequences():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G|NC_000001.11:g.3G>A"}
    cache = {
        "NC_000001.11:g.1A>T": _hit("missense_variant"),
        "NC_000001.11:g.2C>G": _hit("synonymous_variant"),
        "NC_000001.11:g.3G>A": _miss(),
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_g"], access_date="2026-04-30",
    )
    assert out["vep.most_severe_mutational_consequence"] == "missense_variant|synonymous_variant|"
    assert out["vep.error"] == "||"


# ---------------------------------------------------------------------------
# annotate_row — VEP error handling (part a: retain error)
# ---------------------------------------------------------------------------


def test_annotate_row_vep_error_preserved_in_error_column():
    """A VEP-level error for a single candidate is recorded in the error column."""
    row = {"mapped_hgvs_c": "NM_007194.4:c.1259_1260insAAG"}
    cache = {"NM_007194.4:c.1259_1260insAAG": _vep_err()}
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_c"], access_date="2026-06-15",
    )
    assert out["vep.most_severe_mutational_consequence"] == ""
    assert out["vep.error"].startswith("vep_error:")
    assert out["vep.consequence_source"] == ""


def test_annotate_row_api_error_preserved_in_error_column():
    """An api_error for a single candidate is recorded in the error column."""
    row = {"mapped_hgvs_c": "NM_000001.1:c.1A>T"}
    cache = {"NM_000001.1:c.1A>T": _api_err()}
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_c"], access_date="2026-06-15",
    )
    assert out["vep.most_severe_mutational_consequence"] == ""
    assert out["vep.error"].startswith("api_error:")
    assert out["vep.consequence_source"] == ""


# ---------------------------------------------------------------------------
# annotate_row — VEP error fill-in from siblings (part b)
# ---------------------------------------------------------------------------


def test_annotate_row_vep_error_filled_from_valid_sibling():
    """Errored candidate inherits consequences from the valid siblings."""
    row = {"mapped_hgvs_c": "NM_007194.4:c.1259_1260insAAG|NM_007194.4:c.1259_1260insAAT|NM_007194.4:c.1259_1260insAAC"}
    cache = {
        "NM_007194.4:c.1259_1260insAAG": _vep_err(),
        "NM_007194.4:c.1259_1260insAAT": _hit("inframe_insertion", source="transcript", all_cons=["inframe_insertion"], cached_date="2026-06-10"),
        "NM_007194.4:c.1259_1260insAAC": _hit("inframe_insertion", source="transcript", all_cons=["inframe_insertion"], cached_date="2026-06-10"),
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_c"], access_date="2026-06-15",
    )
    ms = out["vep.most_severe_mutational_consequence"].split("|")
    err = out["vep.error"].split("|")
    src = out["vep.consequence_source"].split("|")

    # Errored position (0) gets filled
    assert ms[0] == "inframe_insertion"
    assert err[0].startswith("vep_error:")   # error preserved
    assert src[0] == "transcript"             # source inherited from sibling

    # Valid siblings unchanged
    assert ms[1] == "inframe_insertion"
    assert err[1] == ""
    assert ms[2] == "inframe_insertion"
    assert err[2] == ""


def test_annotate_row_error_fill_uses_most_frequent_tuple():
    """Fill picks the most common (consequences, most_severe) tuple."""
    row = {"mapped_hgvs_g": "A|B|C|D"}
    cache = {
        "A": _vep_err(),
        "B": _hit("synonymous_variant"),
        "C": _hit("inframe_insertion"),
        "D": _hit("inframe_insertion"),   # inframe_insertion appears twice → wins
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_g"], access_date="2026-06-15",
    )
    ms = out["vep.most_severe_mutational_consequence"].split("|")
    err = out["vep.error"].split("|")

    assert ms[0] == "inframe_insertion"   # most frequent wins
    assert err[0].startswith("vep_error:")


def test_annotate_row_error_fill_tie_breaks_deterministically():
    """On a count tie, the fill is deterministic (alphabetical on the tuple)."""
    row = {"mapped_hgvs_g": "A|B|C"}
    cache = {
        "A": _vep_err(),
        "B": _hit("synonymous_variant"),
        "C": _hit("missense_variant"),   # tie; "missense_variant" < "synonymous_variant"
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_g"], access_date="2026-06-15",
    )
    ms = out["vep.most_severe_mutational_consequence"].split("|")
    assert ms[0] == "missense_variant"    # alphabetically first of the tied tuples


def test_annotate_row_vep_error_no_fill_when_all_errored():
    """When all candidates have errors, fill-in is not possible."""
    row = {"mapped_hgvs_c": "NM_000001.1:c.1_2insA|NM_000001.1:c.3_4insG"}
    cache = {
        "NM_000001.1:c.1_2insA": _vep_err(),
        "NM_000001.1:c.3_4insG": _vep_err(),
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_c"], access_date="2026-06-15",
    )
    ms_parts = out["vep.most_severe_mutational_consequence"].split("|")
    err_parts = out["vep.error"].split("|")
    assert all(p == "" for p in ms_parts)
    assert all(p.startswith("vep_error:") for p in err_parts)


# ---------------------------------------------------------------------------
# _vep_lookup_batch — VEP error detection in API response
# ---------------------------------------------------------------------------


def test_vep_lookup_batch_detects_vep_error_in_response():
    """_vep_lookup_batch stores vep_error:... when VEP returns an error entry."""
    fake_response_data = [
        {
            "input": "NM_007194.4:c.1259_1260insAAG",
            "error": "Start (28695710) must be less than or equal to end+1 (28695243)",
        }
    ]
    mock_response = mock.Mock()
    mock_response.status_code = 200
    mock_response.json.return_value = fake_response_data

    with mock.patch("src.annotate_vep.requests.post", return_value=mock_response):
        results, error = mod._vep_lookup_batch(
            ["NM_007194.4:c.1259_1260insAAG"],
            api_url="https://rest.ensembl.org",
            timeout_seconds=10,
            refseq=True,
        )

    assert error is None
    h = "NM_007194.4:c.1259_1260insAAG"
    assert h in results
    most_severe, all_cons, source, _ = results[h]
    assert most_severe is None
    assert all_cons is None
    assert source.startswith("vep_error:")
    assert "28695710" in source


def test_vep_lookup_batch_normal_entry_not_treated_as_error():
    """A normal response entry without an 'error' key is not treated as vep_error."""
    fake_response_data = [
        {
            "input": "NM_000001.1:c.1A>T",
            "most_severe_consequence": "missense_variant",
        }
    ]
    mock_response = mock.Mock()
    mock_response.status_code = 200
    mock_response.json.return_value = fake_response_data

    with mock.patch("src.annotate_vep.requests.post", return_value=mock_response):
        results, error = mod._vep_lookup_batch(
            ["NM_000001.1:c.1A>T"],
            api_url="https://rest.ensembl.org",
            timeout_seconds=10,
        )

    assert error is None
    assert results["NM_000001.1:c.1A>T"][0] == "missense_variant"
    assert not results["NM_000001.1:c.1A>T"][2].startswith("vep_error:")


# ---------------------------------------------------------------------------
# get_functional_consequence — recoder fallback
# ---------------------------------------------------------------------------


def test_get_functional_consequence_uses_recoder_fallback(monkeypatch):
    def fake_run_batches(hgvs_list, *, api_url, timeout_seconds, batch_size, max_workers):
        results = {}
        for h in hgvs_list:
            if h == "NC_000001.11:g.1A>T":
                results[h] = _hit("intron_variant")
            elif h == "NC_000001.11:g.5A>T":
                results[h] = _hit("synonymous_variant")
            elif h == "NC_000001.11:g.6A>T":
                results[h] = _hit("missense_variant")
        return results, {}

    def fake_recoder(hgvs_strings, *, api_url, timeout_seconds):
        assert hgvs_strings == ["NM_000000.1:c.1A>T"]
        return ({"NM_000000.1:c.1A>T": ["NC_000001.11:g.5A>T", "NC_000001.11:g.6A>T"]}, "")

    monkeypatch.setattr(mod, "_run_batches_concurrent", fake_run_batches)
    monkeypatch.setattr(mod, "run_variant_recoder", fake_recoder)

    out = mod.get_functional_consequence(
        ["NC_000001.11:g.1A>T", "NM_000000.1:c.1A>T"],
        api_url="https://rest.ensembl.org",
        timeout_seconds=10,
        batch_size=200,
    )

    assert out["NC_000001.11:g.1A>T"][0] == "intron_variant"
    assert out["NM_000000.1:c.1A>T"][0] == "missense_variant"  # most severe of synonymous vs missense


# ---------------------------------------------------------------------------
# main() integration — skip/limit and keep-existing
# ---------------------------------------------------------------------------


def test_main_applies_skip_and_limit(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["variant_urn", "mapped_hgvs_g"],
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v1", "mapped_hgvs_g": "NC_000001.11:g.1A>T"})
        writer.writerow({"variant_urn": "v2", "mapped_hgvs_g": "NC_000001.11:g.2C>G"})
        writer.writerow({"variant_urn": "v3", "mapped_hgvs_g": "NC_000001.11:g.3G>A"})

    def fake_get(hgvs_strings, *, api_url, timeout_seconds, batch_size, max_workers=1):
        return {
            h: _hit("missense_variant" if h.endswith("2C>G") else "synonymous_variant")
            for h in dict.fromkeys(hgvs_strings)
        }

    monkeypatch.setattr(mod, "get_functional_consequence", fake_get)

    mod.main([str(in_path), str(out_path), "--skip", "1", "--limit", "1", "--row-batch-size", "1"])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert [r["variant_urn"] for r in rows] == ["v2"]
    assert rows[0]["vep.most_severe_mutational_consequence"] == "missense_variant"
    assert rows[0]["vep.error"] == ""
    assert rows[0]["vep.access_date"] == date.today().isoformat()


def test_main_keep_existing_preserves_annotated_fills_blank(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    today = date.today().isoformat()

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "variant_urn", "mapped_hgvs_g",
                "vep.most_severe_mutational_consequence", "vep.access_date", "vep.error",
            ],
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        # Already annotated — should be preserved.
        writer.writerow({
            "variant_urn": "v1", "mapped_hgvs_g": "NC_000001.11:g.1A>T",
            "vep.most_severe_mutational_consequence": "synonymous_variant",
            "vep.access_date": "2025-01-01", "vep.error": "",
        })
        # Blank annotation — should be newly annotated.
        writer.writerow({
            "variant_urn": "v2", "mapped_hgvs_g": "NC_000001.11:g.2C>G",
            "vep.most_severe_mutational_consequence": "",
            "vep.access_date": "", "vep.error": "",
        })

    looked_up: list[str] = []

    def fake_get(hgvs_strings, *, api_url, timeout_seconds, batch_size, max_workers=1):
        looked_up.extend(hgvs_strings)
        return {h: _hit("missense_variant") for h in hgvs_strings}

    monkeypatch.setattr(mod, "get_functional_consequence", fake_get)

    mod.main([str(in_path), str(out_path), "--keep-existing"])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert looked_up == ["NC_000001.11:g.2C>G"]

    assert rows[0]["variant_urn"] == "v1"
    assert rows[0]["vep.most_severe_mutational_consequence"] == "synonymous_variant"
    assert rows[0]["vep.access_date"] == "2025-01-01"

    assert rows[1]["variant_urn"] == "v2"
    assert rows[1]["vep.most_severe_mutational_consequence"] == "missense_variant"
    assert rows[1]["vep.access_date"] == today


def test_main_keep_existing_makes_no_api_calls_when_all_annotated(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "variant_urn", "mapped_hgvs_g",
                "vep.most_severe_mutational_consequence", "vep.access_date", "vep.error",
            ],
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({
            "variant_urn": "v1", "mapped_hgvs_g": "NC_000001.11:g.1A>T",
            "vep.most_severe_mutational_consequence": "synonymous_variant",
            "vep.access_date": "2025-01-01", "vep.error": "",
        })

    call_count: list[int] = []

    def fake_get(hgvs_strings, *, api_url, timeout_seconds, batch_size, max_workers=1):
        call_count.append(1)
        return {}

    monkeypatch.setattr(mod, "get_functional_consequence", fake_get)

    mod.main([str(in_path), str(out_path), "--keep-existing"])

    assert call_count == [], "get_functional_consequence should not be called when all rows already annotated"

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert rows[0]["vep.most_severe_mutational_consequence"] == "synonymous_variant"
    assert rows[0]["vep.access_date"] == "2025-01-01"


# ---------------------------------------------------------------------------
# annotate_row — access_date alignment
# ---------------------------------------------------------------------------


def test_access_date_is_pipe_aligned_for_multiple_candidates():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T|NC_000001.11:g.2C>G"}
    cache = {
        "NC_000001.11:g.1A>T": _hit("missense_variant"),
        "NC_000001.11:g.2C>G": _hit("synonymous_variant"),
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_g"], access_date="2026-05-29",
    )
    assert out["vep.access_date"] == "2026-05-29|2026-05-29"


def test_access_date_preserves_empty_candidate_slots():
    row = {"mapped_hgvs_g": "NC_000001.11:g.1A>T||NC_000001.11:g.3G>A"}
    cache = {
        "NC_000001.11:g.1A>T": _hit("missense_variant"),
        "NC_000001.11:g.3G>A": _hit("intron_variant"),
    }
    out = mod.annotate_row(
        row, cache, col_prefix="vep", hgvs_cols=["mapped_hgvs_g"], access_date="2026-05-29",
    )
    assert out["vep.access_date"] == "2026-05-29||2026-05-29"
