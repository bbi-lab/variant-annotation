import csv
import threading
import time

import pytest

from src import add_variant_position_alleles as mod


def _write_tsv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "mapped_hgvs_g", "mapped_hgvs_c", "mapped_hgvs_p"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_annotate_variants_preserves_input_order_with_concurrency(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {"variant_urn": "v1", "mapped_hgvs_g": "g1", "mapped_hgvs_c": "c1", "mapped_hgvs_p": "p1"},
            {"variant_urn": "v2", "mapped_hgvs_g": "g2", "mapped_hgvs_c": "c2", "mapped_hgvs_p": "p2"},
            {"variant_urn": "v3", "mapped_hgvs_g": "g3", "mapped_hgvs_c": "c3", "mapped_hgvs_p": "p3"},
        ],
    )

    active = 0
    max_active = 0
    lock = threading.Lock()

    def fake_parse_hgvs(hgvs_value, resolve_missing_ref_alleles=False):
        nonlocal active, max_active
        with lock:
            active += 1
            max_active = max(max_active, active)

        text = (hgvs_value or "").strip()
        if text.endswith("1"):
            time.sleep(0.08)
        else:
            time.sleep(0.02)

        with lock:
            active -= 1

        # start, stop, ref, alt, touches_intronic_region, spans_intron
        return "1", "1", "A", "T", False, False

    monkeypatch.setattr(mod, "_parse_hgvs", fake_parse_hgvs)

    mod.annotate_variants(str(inp), str(out), max_workers=3)

    rows = _read_tsv(out)
    assert [r["variant_urn"] for r in rows] == ["v1", "v2", "v3"]
    assert max_active > 1


def test_annotate_variants_invalid_max_workers_raises(tmp_path):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [{"variant_urn": "v1", "mapped_hgvs_g": "", "mapped_hgvs_c": "", "mapped_hgvs_p": ""}],
    )

    with pytest.raises(ValueError, match="max_workers must be >= 1"):
        mod.annotate_variants(str(inp), str(out), max_workers=0)
