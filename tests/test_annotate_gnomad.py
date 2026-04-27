import csv

import src.annotate_gnomad as mod

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


def test_main_applies_skip_and_limit(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "dna_clingen_allele_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v1", "dna_clingen_allele_id": "CA1"})
        writer.writerow({"variant_urn": "v2", "dna_clingen_allele_id": "CA2"})
        writer.writerow({"variant_urn": "v3", "dna_clingen_allele_id": "CA3"})
        writer.writerow({"variant_urn": "v4", "dna_clingen_allele_id": "CA4"})

    monkeypatch.setattr(mod, "ensure_local_gnomad_ht", lambda *args, **kwargs: tmp_path / "dummy.ht")
    monkeypatch.setattr(
        mod,
        "load_gnomad_records_for_caids",
        lambda local_ht, caids, cache_dir: {
            caid: GnomadRecord(
                caid=caid,
                allele_count=1,
                allele_number=10,
                allele_frequency=0.1,
                minor_allele_frequency=0.1,
                faf95_max=None,
                faf95_max_ancestry="",
            )
            for caid in caids
        },
    )

    mod.main([
        str(in_path),
        str(out_path),
        "--gnomad-ht-uri",
        "dummy-uri",
        "--skip",
        "1",
        "--limit",
        "2",
    ])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert [r["variant_urn"] for r in rows] == ["v2", "v3"]
    assert [r["gnomad.v4_1.minor_allele_frequency"] for r in rows] == ["0.1", "0.1"]
