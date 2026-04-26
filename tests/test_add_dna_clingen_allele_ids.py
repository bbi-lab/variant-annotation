import csv
import threading
import time

import pytest

from src import add_dna_clingen_allele_ids as mod


def _write_tsv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "variant_urn",
                "raw_hgvs_nt",
                "raw_hgvs_pro",
                "mapped_hgvs_c",
                "mapped_hgvs_g",
                "clingen_allele_id",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_dna_row_copies_existing_clingen_id(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v1",
                "raw_hgvs_nt": "NM_000001.1:c.1A>G",
                "raw_hgvs_pro": "",
                "mapped_hgvs_c": "NM_000001.1:c.1A>G",
                "mapped_hgvs_g": "NC_000001.11:g.100A>G",
                "clingen_allele_id": "CA111",
            }
        ],
    )

    # Should not need network for this case.
    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", lambda *args, **kwargs: None)

    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    assert rows[0]["dna_clingen_allele_id"] == "CA111"


def test_protein_row_pipe_candidates_preserve_positions(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v_pro",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala1Val",
                "mapped_hgvs_c": "NM_1:c.10A>G|NM_1:c.11A>G|NM_1:c.12A>G",
                "mapped_hgvs_g": "NC_1:g.100A>G|NC_1:g.101A>G|NC_1:g.102A>G",
                "clingen_allele_id": "PA123",  # Protein variant uses PA prefix
            }
        ],
    )

    def fake_query(hgvs, max_retries=3):
        if hgvs == "NM_1:c.10A>G":
            return {"@id": "https://reg.genome.network/allele/CA10"}
        if hgvs == "NC_1:g.101A>G":
            return {"@id": "https://reg.genome.network/allele/CA101"}
        return None

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", fake_query)

    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    # candidate1 -> c hit (CA10), candidate2 -> c miss then g hit (CA101), candidate3 -> miss
    assert rows[0]["dna_clingen_allele_id"] == "CA10|CA101|"


def test_single_candidate_without_existing_id_uses_c_then_g(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v2",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "",
                "mapped_hgvs_c": "NM_2:c.20A>G",
                "mapped_hgvs_g": "NC_2:g.200A>G",
                "clingen_allele_id": "",
            }
        ],
    )

    calls = []

    def fake_query(hgvs, max_retries=3):
        calls.append(hgvs)
        if hgvs == "NC_2:g.200A>G":
            return {"@id": "https://reg.genome.network/allele/CA200"}
        return None

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", fake_query)

    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    assert calls == ["NM_2:c.20A>G", "NC_2:g.200A>G"]
    assert rows[0]["dna_clingen_allele_id"] == "CA200"


def test_mismatched_pipe_lengths_are_padded(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v3",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala1Val",
                "mapped_hgvs_c": "NM_3:c.30A>G|NM_3:c.31A>G",
                "mapped_hgvs_g": "NC_3:g.300A>G",
                "clingen_allele_id": "PA789",  # Protein variant uses PA prefix
            }
        ],
    )

    def fake_query(hgvs, max_retries=3):
        if hgvs == "NM_3:c.30A>G":
            return {"@id": "https://reg.genome.network/allele/CA30"}
        return None

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", fake_query)

    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    # second candidate has c only and no hit -> empty placeholder retained
    assert rows[0]["dna_clingen_allele_id"] == "CA30|"


def test_protein_single_candidate_does_not_reuse_existing_id(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v_pro_single",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala1Val",
                "mapped_hgvs_c": "NM_4:c.40A>G",
                "mapped_hgvs_g": "NC_4:g.400A>G",
                "clingen_allele_id": "PA456",  # Protein variant uses PA prefix
            }
        ],
    )

    def fake_query(hgvs, max_retries=3):
        if hgvs == "NM_4:c.40A>G":
            return {"@id": "https://reg.genome.network/allele/CA40"}
        return None

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", fake_query)

    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    assert rows[0]["dna_clingen_allele_id"] == "CA40"


def test_dna_row_with_invalid_ca_prefix_raises_error(tmp_path, monkeypatch):
    """DNA rows must have clingen_allele_id starting with 'CA'."""
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v_invalid",
                "raw_hgvs_nt": "NM_000001.1:c.1A>G",
                "raw_hgvs_pro": "",
                "mapped_hgvs_c": "NM_000001.1:c.1A>G",
                "mapped_hgvs_g": "NC_000001.11:g.100A>G",
                "clingen_allele_id": "PA999",  # Wrong prefix for DNA variant
            }
        ],
    )

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", lambda *args, **kwargs: None)

    import pytest
    with pytest.raises(ValueError, match="Invalid ClinGen allele ID prefix for DNA variant"):
        mod.add_dna_clingen_allele_ids(str(inp), str(out))


def test_protein_row_with_valid_pa_prefix(tmp_path, monkeypatch):
    """Protein rows with valid PA prefix should succeed."""
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v_pro_pa",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala1Val",
                "mapped_hgvs_c": "NM_1:c.10A>G",
                "mapped_hgvs_g": "NC_1:g.100A>G",
                "clingen_allele_id": "PA12345",  # Valid protein prefix
            }
        ],
    )

    def fake_query(hgvs, max_retries=3):
        if hgvs == "NM_1:c.10A>G":
            return {"@id": "https://reg.genome.network/allele/CA10"}
        return None

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", fake_query)

    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    # Protein row should do new lookup, not reuse PA prefix
    assert rows[0]["dna_clingen_allele_id"] == "CA10"


def test_protein_row_with_invalid_ca_prefix_raises_error(tmp_path, monkeypatch):
    """Protein rows must have clingen_allele_id starting with 'PA'."""
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v_pro_invalid",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala1Val",
                "mapped_hgvs_c": "NM_1:c.10A>G",
                "mapped_hgvs_g": "NC_1:g.100A>G",
                "clingen_allele_id": "CA999",  # Wrong prefix for protein variant
            }
        ],
    )

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", lambda *args, **kwargs: None)

    import pytest
    with pytest.raises(ValueError, match="Invalid ClinGen allele ID prefix for protein variant"):
        mod.add_dna_clingen_allele_ids(str(inp), str(out))


def test_empty_clingen_id_does_not_raise_error(tmp_path, monkeypatch):
    """Rows with empty clingen_allele_id should not raise validation errors."""
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v_empty",
                "raw_hgvs_nt": "NM_000001.1:c.1A>G",
                "raw_hgvs_pro": "",
                "mapped_hgvs_c": "NM_000001.1:c.1A>G",
                "mapped_hgvs_g": "NC_000001.11:g.100A>G",
                "clingen_allele_id": "",  # Empty is OK
            }
        ],
    )

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", lambda *args, **kwargs: None)

    # Should not raise
    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    assert rows[0]["dna_clingen_allele_id"] == ""


def test_concurrent_processing_preserves_row_order(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v1",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala1Val",
                "mapped_hgvs_c": "NM_1:c.10A>G",
                "mapped_hgvs_g": "NC_1:g.100A>G",
                "clingen_allele_id": "PA1",
            },
            {
                "variant_urn": "v2",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala2Val",
                "mapped_hgvs_c": "NM_2:c.20A>G",
                "mapped_hgvs_g": "NC_2:g.200A>G",
                "clingen_allele_id": "PA2",
            },
            {
                "variant_urn": "v3",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala3Val",
                "mapped_hgvs_c": "NM_3:c.30A>G",
                "mapped_hgvs_g": "NC_3:g.300A>G",
                "clingen_allele_id": "PA3",
            },
        ],
    )

    active = 0
    max_active = 0
    lock = threading.Lock()

    def fake_query(hgvs, max_retries=3):
        nonlocal active, max_active
        with lock:
            active += 1
            max_active = max(max_active, active)
        # Intentionally make first row slower so completion order differs.
        if "NM_1" in hgvs:
            time.sleep(0.08)
        else:
            time.sleep(0.02)
        with lock:
            active -= 1
        # Resolve via c. lookup
        suffix = hgvs.split(":", 1)[0].split("_")[-1]
        return {"@id": f"https://reg.genome.network/allele/CA{suffix}"}

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", fake_query)

    mod.add_dna_clingen_allele_ids(str(inp), str(out), max_workers=3)

    rows = _read_tsv(out)
    # Output must remain in the original input order.
    assert [r["variant_urn"] for r in rows] == ["v1", "v2", "v3"]
    # Work should have been concurrent (more than one request active at once).
    assert max_active > 1


def test_invalid_max_workers_raises(tmp_path):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v1",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "",
                "mapped_hgvs_c": "",
                "mapped_hgvs_g": "",
                "clingen_allele_id": "",
            }
        ],
    )

    with pytest.raises(ValueError, match="max_workers must be >= 1"):
        mod.add_dna_clingen_allele_ids(str(inp), str(out), max_workers=0)


def test_placeholder_ids_from_clingen_response_are_treated_as_null(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v_placeholder",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala1Val",
                "mapped_hgvs_c": "NM_1:c.10A>G",
                "mapped_hgvs_g": "NC_1:g.100A>G",
                "clingen_allele_id": "PA1",
            }
        ],
    )

    def fake_query(hgvs, max_retries=3):
        if hgvs == "NM_1:c.10A>G":
            return {"@id": "_:CA123"}
        if hgvs == "NC_1:g.100A>G":
            return {"@id": "_:PA999"}
        return None

    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", fake_query)

    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    assert rows[0]["dna_clingen_allele_id"] == ""


def test_existing_placeholder_id_is_treated_as_empty_for_validation_and_reuse(tmp_path, monkeypatch):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    _write_tsv(
        inp,
        [
            {
                "variant_urn": "v_existing_placeholder",
                "raw_hgvs_nt": "NM_000001.1:c.1A>G",
                "raw_hgvs_pro": "",
                "mapped_hgvs_c": "NM_000001.1:c.1A>G",
                "mapped_hgvs_g": "NC_000001.11:g.100A>G",
                "clingen_allele_id": "_:CA123",  # Should be treated as empty
            }
        ],
    )

    # No API hit, so output should remain empty and should not raise prefix error.
    monkeypatch.setattr(mod, "_query_clingen_by_hgvs", lambda *args, **kwargs: None)

    mod.add_dna_clingen_allele_ids(str(inp), str(out))

    rows = _read_tsv(out)
    assert rows[0]["dna_clingen_allele_id"] == ""
