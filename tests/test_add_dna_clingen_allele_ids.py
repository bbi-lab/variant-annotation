import csv

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
