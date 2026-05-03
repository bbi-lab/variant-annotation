import csv

from src.normalize_intra_codon_haplotypes import haplotypes_to_delins, normalize_haplotype_to_delins


def test_normalize_haplotype_to_delins_same_codon_with_gap_uses_target_sequence():
    def fake_mapper(accession, components, target_sequence, row_label, dcd):
        assert accession is None
        assert components == ["1A>G", "3G>T"]
        return ["c.1A>G", "c.3G>T"]

    # Position 2 remains unchanged and should be copied from target_sequence.
    rewritten = normalize_haplotype_to_delins(
        raw_hgvs_nt="c.[1A>G;3G>T]",
        target_sequence="ATGCC",
        row_label="1",
        dcd=None,
        mapper=fake_mapper,
    )

    assert rewritten == "c.1_3delinsGTT"


def test_normalize_haplotype_to_delins_skips_intronic_or_utr_components():
    rewritten = normalize_haplotype_to_delins(
        raw_hgvs_nt="c.[76+1G>A;77A>T]",
        target_sequence="ATGCC",
        row_label="1",
        dcd=None,
    )
    assert rewritten is None


def test_normalize_haplotype_to_delins_requires_single_codon():
    def fake_mapper(accession, components, target_sequence, row_label, dcd):
        return ["c.1A>G", "c.4C>T"]

    rewritten = normalize_haplotype_to_delins(
        raw_hgvs_nt="c.[1A>G;4C>T]",
        target_sequence="ATGCC",
        row_label="1",
        dcd=None,
        mapper=fake_mapper,
    )
    assert rewritten is None


def test_haplotypes_to_delins_adds_orig_column_and_rewrites_rows(tmp_path):
    input_file = tmp_path / "in.tsv"
    output_file = tmp_path / "out.tsv"

    with input_file.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "raw_hgvs_nt", "target_sequence"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "variant_urn": "v1",
                "raw_hgvs_nt": "c.[1A>G;3G>T]",
                "target_sequence": "ATGCC",
            }
        )
        writer.writerow(
            {
                "variant_urn": "v2",
                "raw_hgvs_nt": "c.[1A>G;4C>T]",
                "target_sequence": "ATGCC",
            }
        )

    # Monkeypatch via module-level import to avoid network/dcd calls in tests.
    from src import normalize_intra_codon_haplotypes as mod

    def fake_mapper(accession, components, target_sequence, row_label, dcd):
        if components == ["1A>G", "3G>T"]:
            return ["c.1A>G", "c.3G>T"]
        return ["c.1A>G", "c.4C>T"]

    original_mapper = mod._map_components_to_mapped_c
    mod._map_components_to_mapped_c = fake_mapper
    try:
        stats = haplotypes_to_delins(str(input_file), str(output_file))
    finally:
        mod._map_components_to_mapped_c = original_mapper

    assert stats["rows"] == 2
    assert stats["rewritten"] == 1

    with output_file.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert "orig_raw_hgvs_nt" in rows[0]
    assert rows[0]["orig_raw_hgvs_nt"] == "c.[1A>G;3G>T]"
    assert rows[0]["raw_hgvs_nt"] == "c.1_3delinsGTT"

    assert rows[1]["orig_raw_hgvs_nt"] == "c.[1A>G;4C>T]"
    assert rows[1]["raw_hgvs_nt"] == "c.[1A>G;4C>T]"
