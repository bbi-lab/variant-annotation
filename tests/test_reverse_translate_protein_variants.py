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
    assert "Unable to resolve transcript accession" in out_rows[0]["reverse_translation_error"]


def test_derive_joined_hgvs_fields_tracks_intronic_and_spans(monkeypatch):

    # Simulate parser behavior for an intron-spanning del and a standard substitution.
    def fake_parse_hgvs(candidate, resolve_missing_ref_alleles=True):
        if "76+1_77-1del" in candidate:
            return ("76+1", "77-1", "AG", "", True, True, "1")
        return ("90", "90", "A", "G", False, False, "1")

    monkeypatch.setattr(rtpv, "_parse_hgvs", fake_parse_hgvs)

    joined = "NM_000001.1:c.76+1_77-1del|NM_000001.1:c.90A>G"
    start, stop, ref, alt, touches_intronic, spans_intron, chromosome, warnings = rtpv._derive_joined_hgvs_fields(
        joined,
        resolve_missing_ref_alleles=True,
    )

    assert start == "76+1|90"
    assert stop == "77-1|90"
    assert ref == "AG|A"
    assert alt == "|G"
    assert touches_intronic == "true|false"
    assert spans_intron == "true|false"
    assert chromosome == "1|1"


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


# ---------------------------------------------------------------------------
# WT codon helpers
# ---------------------------------------------------------------------------

import pytest


@pytest.mark.parametrize(
    "hgvs_p, expected",
    [
        ("NP_000001.1:p.Met1Met", ("Met", 1, "Met")),
        ("NM_000001.1:p.Glu2Glu", ("Glu", 2, "Glu")),
        ("p.Trp3Trp", ("Trp", 3, "Trp")),
        # Synonymous shorthand "=" should be expanded to the ref AA
        ("NP_000001.1:p.Ala5=", ("Ala", 5, "Ala")),
        # Non-WT (ref != alt) should still parse correctly
        ("NP_000001.1:p.Ala1Val", ("Ala", 1, "Val")),
        # Cannot parse bare string without p. prefix
        ("Met1Met", None),
        # Empty string
        ("", None),
    ],
)
def test_parse_protein_hgvs_aa_change(hgvs_p, expected):
    assert rtpv._parse_protein_hgvs_aa_change(hgvs_p) == expected


@pytest.mark.parametrize(
    "transcript, aa_position, codon, expected",
    [
        ("NM_000001.1", 1, "ATG", "NM_000001.1:c.1_3delinsATG"),
        ("NM_000001.1", 2, "GAA", "NM_000001.1:c.4_6delinsGAA"),
        ("NM_000002.1", 100, "TGG", "NM_000002.1:c.298_300delinsTGG"),
    ],
)
def test_build_wt_c_hgvs(transcript, aa_position, codon, expected):
    assert rtpv._build_wt_c_hgvs(transcript, aa_position, codon) == expected


@pytest.mark.parametrize(
    "ref_aa3, mode, expected_codon",
    [
        ("Met", "unambiguous", "ATG"),
        ("Trp", "unambiguous", "TGG"),
        ("Glu", "unambiguous", None),  # not unambiguous
        ("Met", "all", "ATG"),
        ("Trp", "all", "TGG"),
        ("Met", "none", None),
    ],
)
def test_get_wt_codon_for_mode_unambiguous_cases(ref_aa3, mode, expected_codon, monkeypatch):
    # "all" mode for Glu would query UTA; patch it to avoid a real DB call.
    monkeypatch.setattr(rtpv, "_lookup_codon_from_uta", lambda t, p, c: "GAA")
    codon = rtpv._get_wt_codon_for_mode(ref_aa3, 1, "NM_000001.1", mode, object())
    assert codon == expected_codon


def test_get_wt_codon_for_mode_all_queries_uta_for_ambiguous(monkeypatch):
    monkeypatch.setattr(rtpv, "_lookup_codon_from_uta", lambda t, p, c: "GAG")
    codon = rtpv._get_wt_codon_for_mode("Glu", 2, "NM_000001.1", "all", object())
    assert codon == "GAG"


def test_wt_codon_mode_requires_include_indels(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(input_path, [])
    monkeypatch.setenv("UTA_DB_URL", "postgresql://user:pass@localhost:5432/uta")
    monkeypatch.setattr(rtpv, "_find_reverse_translate_cli", lambda: "/usr/bin/reverse-translate-variants")

    with pytest.raises(ValueError, match="--include-indels"):
        rtpv.reverse_translate_protein_variants(
            str(input_path),
            str(output_path),
            wt_codon_mode="unambiguous",
            include_indels=False,
        )


def test_wt_codon_mode_invalid_value(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    _write_tsv(input_path, [])
    monkeypatch.setenv("UTA_DB_URL", "postgresql://user:pass@localhost:5432/uta")
    monkeypatch.setattr(rtpv, "_find_reverse_translate_cli", lambda: "/usr/bin/reverse-translate-variants")

    with pytest.raises(ValueError, match="wt_codon_mode must be"):
        rtpv.reverse_translate_protein_variants(
            str(input_path),
            str(output_path),
            wt_codon_mode="invalid",
            include_indels=True,
        )


def test_reverse_translate_wt_codon_mode_unambiguous_appends_met_candidate(tmp_path, monkeypatch):
    """For p.Met1Met with wt_codon_mode='unambiguous', a WT codon delins should be appended."""
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {
            "variant_urn": "wt0",
            "mapped_hgvs_p": "NP_000001.1:p.Met1Met",
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
    monkeypatch.setattr(rtpv, "_populate_derived_hgvs_columns", lambda *args, **kwargs: None)

    # Patch the hgvs library at the module level so the imports inside the function get mocks.
    import hgvs.dataproviders.uta
    import hgvs.assemblymapper
    import hgvs.parser as hgvs_parser_mod
    monkeypatch.setattr(hgvs.dataproviders.uta, "connect", lambda url: object())
    monkeypatch.setattr(hgvs.assemblymapper, "AssemblyMapper", lambda *a, **kw: object())
    monkeypatch.setattr(hgvs_parser_mod, "Parser", lambda: object())

    # CLI returns empty (no synonymous variants for Met)
    monkeypatch.setattr(rtpv, "_run_reverse_translate_batch", lambda *a, **kw: ([{"hgvs_c": "", "hgvs_g": ""}], []))
    # g. mapping returns a fake genomic string
    monkeypatch.setattr(rtpv, "_do_map_c_to_g", lambda c, p, m: "NC_000001.11:g.1_3delinsATG")

    rtpv.reverse_translate_protein_variants(
        str(input_path),
        str(output_path),
        include_indels=True,
        wt_codon_mode="unambiguous",
    )

    out_rows = _read_tsv(output_path)
    assert len(out_rows) == 1
    assert out_rows[0]["mapped_hgvs_c"] == "NM_000001.1:c.1_3delinsATG"
    assert out_rows[0]["mapped_hgvs_g"] == "NC_000001.11:g.1_3delinsATG"


def test_reverse_translate_wt_codon_not_duplicated_when_already_present(tmp_path, monkeypatch):
    """If CLI already returned the WT codon as a synonymous candidate, don't duplicate it."""
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {
            "variant_urn": "syn0",
            "mapped_hgvs_p": "NP_000001.1:p.Met1Met",
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
    monkeypatch.setattr(rtpv, "_populate_derived_hgvs_columns", lambda *args, **kwargs: None)
    monkeypatch.setattr(rtpv, "_do_map_c_to_g", lambda c, p, m: "NC_000001.11:g.1_3delinsATG")

    # CLI already returned the WT delins
    monkeypatch.setattr(
        rtpv,
        "_run_reverse_translate_batch",
        lambda *a, **kw: ([{"hgvs_c": "NM_000001.1:c.1_3delinsATG", "hgvs_g": "NC_000001.11:g.1_3delinsATG"}], []),
    )

    import hgvs.dataproviders.uta
    import hgvs.assemblymapper
    import hgvs.parser as hgvs_parser_mod
    monkeypatch.setattr(hgvs.dataproviders.uta, "connect", lambda url: object())
    monkeypatch.setattr(hgvs.assemblymapper, "AssemblyMapper", lambda *a, **kw: object())
    monkeypatch.setattr(hgvs_parser_mod, "Parser", lambda: object())

    rtpv.reverse_translate_protein_variants(
        str(input_path),
        str(output_path),
        include_indels=True,
        wt_codon_mode="unambiguous",
    )

    out_rows = _read_tsv(output_path)
    # Should NOT be doubled
    assert out_rows[0]["mapped_hgvs_c"] == "NM_000001.1:c.1_3delinsATG"


def test_build_parser_wt_codon_mode_default():
    parser = rtpv._build_parser()
    args = parser.parse_args(["input.tsv", "output.tsv"])
    assert args.wt_codon_mode == "none"


def test_build_parser_wt_codon_mode_choices():
    parser = rtpv._build_parser()
    for choice in ("none", "unambiguous", "all"):
        args = parser.parse_args(["input.tsv", "output.tsv", "--wt-codon-mode", choice])
        assert args.wt_codon_mode == choice

