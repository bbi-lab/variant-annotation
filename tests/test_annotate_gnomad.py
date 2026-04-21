from src.annotate_gnomad import GnomadRecord, annotate_row


def test_annotate_row_prefers_first_matching_candidate():
    records = {
        "CA2": GnomadRecord(
            caid="CA2",
            allele_count=10,
            allele_number=100,
            allele_frequency=0.1,
            minor_allele_frequency=0.1,
            faf95_max=0.09,
            faf95_max_ancestry="nfe",
        ),
        "CA3": GnomadRecord(
            caid="CA3",
            allele_count=5,
            allele_number=100,
            allele_frequency=0.05,
            minor_allele_frequency=0.05,
            faf95_max=None,
            faf95_max_ancestry="",
        ),
    }

    row = {"dna_clingen_allele_id": "CA1|CA2|CA3"}
    out = annotate_row(row, records, "gnomad.v4.1", "dna_clingen_allele_id")

    assert out["gnomad.v4.1.minor_allele_frequency"] == "0.1"
    assert out["gnomad.v4.1.allele_frequency"] == "0.1"
    assert out["gnomad.v4.1.allele_count"] == "10"
    assert out["gnomad.v4.1.allele_number"] == "100"
    assert out["gnomad.v4.1.faf95_max"] == "0.09"
    assert out["gnomad.v4.1.faf95_max_ancestry"] == "nfe"


def test_annotate_row_handles_no_match():
    records = {}
    row = {"dna_clingen_allele_id": "CA1|CA2"}
    out = annotate_row(row, records, "gnomad.v4.1", "dna_clingen_allele_id")

    assert out["gnomad.v4.1.minor_allele_frequency"] == ""
    assert out["gnomad.v4.1.allele_frequency"] == ""


def test_annotate_row_custom_dna_column():
    records = {
        "CA7": GnomadRecord(
            caid="CA7",
            allele_count=2,
            allele_number=20,
            allele_frequency=0.1,
            minor_allele_frequency=0.1,
            faf95_max=None,
            faf95_max_ancestry="",
        )
    }
    row = {"my_dna_ids": "CA7"}
    out = annotate_row(row, records, "gnomad.v4.1", "my_dna_ids")

    assert out["gnomad.v4.1.minor_allele_frequency"] == "0.1"
