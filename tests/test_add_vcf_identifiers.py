import csv
import threading
import time

import pytest

from src import add_vcf_identifiers as mod


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

        # start, stop, ref, alt, touches_intronic_region, spans_intron, chromosome
        return "1", "1", "A", "T", False, False, "1"

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


def test_annotate_variants_skip_and_limit(tmp_path):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {"variant_urn": "v1", "mapped_hgvs_g": "", "mapped_hgvs_c": "", "mapped_hgvs_p": ""},
            {"variant_urn": "v2", "mapped_hgvs_g": "", "mapped_hgvs_c": "", "mapped_hgvs_p": ""},
            {"variant_urn": "v3", "mapped_hgvs_g": "", "mapped_hgvs_c": "", "mapped_hgvs_p": ""},
            {"variant_urn": "v4", "mapped_hgvs_g": "", "mapped_hgvs_c": "", "mapped_hgvs_p": ""},
        ],
    )

    mod.annotate_variants(str(inp), str(out), max_workers=1, skip=1, limit=2)

    rows = _read_tsv(out)
    assert [r["variant_urn"] for r in rows] == ["v2", "v3"]


def test_annotate_variants_invalid_skip_or_limit_raises(tmp_path):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [{"variant_urn": "v1", "mapped_hgvs_g": "", "mapped_hgvs_c": "", "mapped_hgvs_p": ""}],
    )

    with pytest.raises(ValueError, match="skip must be >= 0"):
        mod.annotate_variants(str(inp), str(out), skip=-1)

    with pytest.raises(ValueError, match="limit must be >= 0"):
        mod.annotate_variants(str(inp), str(out), limit=-1)


def test_aa_3to1_conversion():
    """Test amino acid 3-letter to 1-letter conversion."""
    assert mod._aa_3to1("Ala") == "A"
    assert mod._aa_3to1("Arg") == "R"
    assert mod._aa_3to1("Leu") == "L"
    assert mod._aa_3to1("Pro") == "P"
    assert mod._aa_3to1("*") == "*"
    assert mod._aa_3to1("Ter") == "*"  # Termination -> *
    assert mod._aa_3to1("A") == "A"  # Already 1-letter
    assert mod._aa_3to1("ala") == "A"  # Case insensitive
    assert mod._aa_3to1("ALA") == "A"  # Case insensitive


def test_extract_chromosome_from_hgvs():
    """Test chromosome extraction from HGVS accession codes."""
    assert mod._extract_chromosome_from_hgvs("NC_000001.11:g.100A>T") == "1"
    assert mod._extract_chromosome_from_hgvs("NC_000022.11:g.1000C>G") == "22"
    assert mod._extract_chromosome_from_hgvs("NC_000023.11:g.1A>G") == "X"  # X chromosome
    assert mod._extract_chromosome_from_hgvs("NC_000024.10:g.1A>G") == "Y"  # Y chromosome
    assert mod._extract_chromosome_from_hgvs("NC_012920.1:m.1A>G") == "M"  # Mitochondrial
    assert mod._extract_chromosome_from_hgvs("NM_000001:c.100A>T") is None  # Transcript, no chromosome
    assert mod._extract_chromosome_from_hgvs(None) is None
    assert mod._extract_chromosome_from_hgvs("") is None


def test_protein_hgvs_with_1letter_codes():
    """Test that protein HGVS parsing uses 1-letter amino acid codes."""
    # Test single amino acid substitution
    start, stop, ref, alt, _, _, _ = mod._parse_hgvs("NP_000001.2:p.Pro656Leu")
    assert ref == "P"
    assert alt == "L"
    assert start == "656"
    
    # Test synonymous variant (ref == alt should repeat)
    start, stop, ref, alt, _, _, _ = mod._parse_hgvs("NP_000001.2:p.Pro656Pro")
    assert ref == "P"
    assert alt == "P"
    
    # Test deletion — alt should be empty string
    start, stop, ref, alt, _, _, _ = mod._parse_hgvs("NP_000001.2:p.Pro656del")
    assert ref == "P"
    assert alt == ""

    # Test termination alt -> *
    start, stop, ref, alt, _, _, _ = mod._parse_hgvs("NP_000001.2:p.Pro656Ter")
    assert ref == "P"
    assert alt == "*"

    # Test range deletion — alt should be empty string
    start, stop, ref, alt, _, _, _ = mod._parse_hgvs("NP_000001.2:p.Pro656_Leu660del")
    assert ref == "P_L"
    assert alt == ""


def test_genomic_chromosome_extraction():
    """Test full parsing with chromosome extraction for genomic variants."""
    start, stop, ref, alt, _, _, chrom = mod._parse_hgvs("NC_000001.11:g.1000A>T")
    assert chrom == "1"
    assert ref == "A"
    assert alt == "T"


def test_pipe_delimited_hgvs_columns_produce_pipe_delimited_output(tmp_path):
    """Pipe-delimited HGVS values (protein reverse translations) should produce
    pipe-delimited output columns with one entry per DNA candidate."""
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"

    with open(inp, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "mapped_hgvs_g", "mapped_hgvs_c", "mapped_hgvs_p"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({
            "variant_urn": "v1",
            "mapped_hgvs_g": "NC_000010.11:g.87864473A>T|NC_000010.11:g.87864474C>G",
            "mapped_hgvs_c": "NM_000314.8:c.4A>T|NM_000314.8:c.5C>G",
            "mapped_hgvs_p": "NP_000305.3:p.Thr2Ser",
        })

    mod.annotate_variants(str(inp), str(out), max_workers=1)

    rows = _read_tsv(out)
    row = rows[0]

    # Genomic column: two candidates, both on chr 10
    assert row["mapped_hgvs_g_chromosome"] == "10|10"
    assert row["mapped_hgvs_g_start"] == "87864473|87864474"
    assert row["mapped_hgvs_g_ref"] == "A|C"
    assert row["mapped_hgvs_g_alt"] == "T|G"

    # Transcript column: two candidates
    assert row["mapped_hgvs_c_start"] == "4|5"
    assert row["mapped_hgvs_c_ref"] == "A|C"

    # Protein column: single value, no pipe
    assert "|" not in row["mapped_hgvs_p_start"]
    assert row["mapped_hgvs_p_ref"] == "T"  # Thr -> T
    assert row["mapped_hgvs_p_alt"] == "S"  # Ser -> S

