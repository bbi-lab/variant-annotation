import csv

from src import reverse_translate_protein_variants as rtpv


class _FakeConnection:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def _write_tsv(path, rows):
    fieldnames = ["variant_urn", "mapped_hgvs_p", "mapped_hgvs_c", "mapped_hgvs_g", "mapping_error", "assayed_variant_level"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_reverse_translate_preserves_row_order_and_updates_targets(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {
            "variant_urn": "v0",
            "mapped_hgvs_p": "NP_000001.1:p.Ala1Val",
            "mapped_hgvs_c": "",
            "mapped_hgvs_g": "",
            "mapping_error": "",
            "assayed_variant_level": "",
        },
        {
            "variant_urn": "v1",
            "mapped_hgvs_p": "NP_000001.1:p.Ala2Val",
            "mapped_hgvs_c": "NM_000001.1:c.2A>G",
            "mapped_hgvs_g": "NC_000001.11:g.2A>G",
            "mapping_error": "",
            "assayed_variant_level": "dna",
        },
        {
            "variant_urn": "v2",
            "mapped_hgvs_p": "NP_000001.1:p.Ala3Val",
            "mapped_hgvs_c": "",
            "mapped_hgvs_g": "",
            "mapping_error": "",
            "assayed_variant_level": "",
        },
    ]
    _write_tsv(input_path, rows)

    monkeypatch.setenv("UTA_DB_URL", "postgresql://user:pass@localhost:5432/uta")
    monkeypatch.setattr(rtpv, "_find_reverse_translate_cli", lambda: "/usr/bin/reverse-translate-variants")
    monkeypatch.setattr(rtpv.psycopg2, "connect", lambda **kwargs: _FakeConnection())
    monkeypatch.setattr(rtpv, "_resolve_transcript_accession", lambda *args, **kwargs: "NM_000001.1")

    def fake_batch(cli_path, rows, **kwargs):
        translated = []
        for row in rows:
            translated.append(
                {
                    "transcript": row["transcript"],
                    "hgvs_p": row["hgvs_p"],
                    "hgvs_c": f"{row['transcript']}:c.123A>G",
                    "hgvs_g": "NC_000001.11:g.123A>G",
                }
            )
        return translated, []

    monkeypatch.setattr(rtpv, "_run_reverse_translate_batch", fake_batch)
    monkeypatch.setattr(rtpv, "_populate_derived_hgvs_columns", lambda *args, **kwargs: None)

    rtpv.reverse_translate_protein_variants(str(input_path), str(output_path))

    out_rows = _read_tsv(output_path)
    assert [r["variant_urn"] for r in out_rows] == ["v0", "v1", "v2"]
    assert out_rows[0]["mapped_hgvs_c"] == "NM_000001.1:c.123A>G"
    assert out_rows[1]["mapped_hgvs_c"] == "NM_000001.1:c.2A>G"
    assert out_rows[2]["mapped_hgvs_c"] == "NM_000001.1:c.123A>G"


def test_reverse_translate_adds_error_when_transcript_unresolved(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {
            "variant_urn": "v0",
            "mapped_hgvs_p": "NP_000001.1:p.Ala1Val",
            "mapped_hgvs_c": "",
            "mapped_hgvs_g": "",
            "mapping_error": "existing",
            "assayed_variant_level": "",
        }
    ]
    _write_tsv(input_path, rows)

    monkeypatch.setenv("UTA_DB_URL", "postgresql://user:pass@localhost:5432/uta")
    monkeypatch.setattr(rtpv, "_find_reverse_translate_cli", lambda: "/usr/bin/reverse-translate-variants")
    monkeypatch.setattr(rtpv.psycopg2, "connect", lambda **kwargs: _FakeConnection())
    monkeypatch.setattr(rtpv, "_resolve_transcript_accession", lambda *args, **kwargs: "")
    monkeypatch.setattr(rtpv, "_run_reverse_translate_batch", lambda *args, **kwargs: ([], []))
    monkeypatch.setattr(rtpv, "_populate_derived_hgvs_columns", lambda *args, **kwargs: None)

    rtpv.reverse_translate_protein_variants(str(input_path), str(output_path))

    out_rows = _read_tsv(output_path)
    assert len(out_rows) == 1
    assert "existing" in out_rows[0]["mapping_error"]
    assert "Unable to resolve transcript accession" in out_rows[0]["mapping_error"]


def test_derive_joined_hgvs_fields_tracks_intronic_and_spans(monkeypatch):
    # Simulate parser behavior for an intron-spanning del and a standard substitution.
    def fake_parse_hgvs(candidate, resolve_missing_ref_alleles=True):
        if "76+1_77-1del" in candidate:
            return ("76+1", "77-1", "AG", "", True, True)
        return ("90", "90", "A", "G", False, False)

    monkeypatch.setattr(rtpv, "_parse_hgvs", fake_parse_hgvs)

    joined = "NM_000001.1:c.76+1_77-1del|NM_000001.1:c.90A>G"
    start, stop, ref, alt, touches_intronic, spans_intron = rtpv._derive_joined_hgvs_fields(
        joined,
        resolve_missing_ref_alleles=True,
    )

    assert start == "76+1|90"
    assert stop == "77-1|90"
    assert ref == "AG|A"
    assert alt == "|G"
    assert touches_intronic == "true|false"
    assert spans_intron == "true|false"


def test_build_parser_defaults_and_flags():
    parser = rtpv._build_parser()
    args = parser.parse_args(["input.tsv", "output.tsv"])

    assert args.input_file == "input.tsv"
    assert args.output_file == "output.tsv"
    assert args.assembly == "GRCh38"
    assert args.strict_ref_aa is True
    assert args.resolve_missing_ref_alleles is True
    assert args.allow_length_changing_stop_candidates is True


def test_resolve_transcript_accession_prefers_transcript_in_protein_hgvs(monkeypatch):
    row = {
        "mapped_hgvs_p": "NM_000001.1:p.Arg17His",
        "fallback_col": "NM_999999.1:c.1A>G",
    }

    # Should not call RefSeq-protein lookup when mapped_hgvs_p already has transcript accession.
    monkeypatch.setattr(
        rtpv,
        "_resolve_transcript_from_refseq_protein_id",
        lambda connection, protein_accession: "NM_SHOULD_NOT_BE_USED.1",
    )

    transcript = rtpv._resolve_transcript_accession(
        row,
        mapped_hgvs_p_col="mapped_hgvs_p",
        transcript_fallback_columns=("fallback_col",),
        uta_connection=object(),
    )
    assert transcript == "NM_000001.1"


def test_resolve_transcript_accession_uses_refseq_mapping_before_fallback(monkeypatch):
    row = {
        "mapped_hgvs_p": "NP_000001.1:p.Arg17His",
        "fallback_col": "NM_999999.1:c.1A>G",
    }

    monkeypatch.setattr(
        rtpv,
        "_resolve_transcript_from_refseq_protein_id",
        lambda connection, protein_accession: "NM_000010.2",
    )

    transcript = rtpv._resolve_transcript_accession(
        row,
        mapped_hgvs_p_col="mapped_hgvs_p",
        transcript_fallback_columns=("fallback_col",),
        uta_connection=object(),
    )
    assert transcript == "NM_000010.2"


def test_resolve_transcript_accession_falls_back_when_refseq_lookup_empty(monkeypatch):
    row = {
        "mapped_hgvs_p": "NP_000001.1:p.Arg17His",
        "fallback_col": "NM_123456.7:c.42G>A",
    }

    monkeypatch.setattr(
        rtpv,
        "_resolve_transcript_from_refseq_protein_id",
        lambda connection, protein_accession: "",
    )

    transcript = rtpv._resolve_transcript_accession(
        row,
        mapped_hgvs_p_col="mapped_hgvs_p",
        transcript_fallback_columns=("fallback_col",),
        uta_connection=object(),
    )
    assert transcript == "NM_123456.7"
