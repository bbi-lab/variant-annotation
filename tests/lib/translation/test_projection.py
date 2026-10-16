# Copyright (C) 2026 Jeremy Stone
# SPDX-License-Identifier: AGPL-3.0-or-later

import random
from dataclasses import dataclass

import pytest

hgvs_exceptions = pytest.importorskip("hgvs.exceptions")

from variant_annotation.lib.translation.projection import coding_to_genomic, genomic_to_coding  # noqa: E402

# ---------------------------------------------------------------------------
# Rules, against a coordinate translator that answers from a table
# ---------------------------------------------------------------------------


class _TableCoordinates:
    """Answers each projection from ``answers``; anything else is an hgvs rejection of the variant."""

    def __init__(self, answers: dict[tuple[str, str], str]):
        self.answers = answers

    def _answer(self, op: str, hgvs: str) -> str:
        if (op, hgvs) not in self.answers:
            raise hgvs_exceptions.HGVSInvalidVariantError(f"no {op} for {hgvs}")
        return self.answers[(op, hgvs)]

    def c_to_g(self, c):
        return self._answer("c_to_g", c)

    def c_to_g_literal(self, c):
        return self._answer("c_to_g_literal", c)

    def g_to_c(self, g, tx):
        return self._answer("g_to_c", g)

    def g_to_c_literal(self, g, tx):
        return self._answer("g_to_c_literal", g)

    def c_to_n(self, c):
        return self._answer("c_to_n", c)

    def c_to_p(self, c):
        raise NotImplementedError


IN_EXON_C = "NM_1.1:c.151_153delinsGCG"
IN_EXON_G = "NC_000016.10:g.1315654_1315656delinsGCG"


@pytest.mark.unit
def test_a_change_within_one_exon_projects_to_its_normalized_genomic_form():
    coordinates = _TableCoordinates(
        {
            ("c_to_g_literal", IN_EXON_C): IN_EXON_G,
            ("c_to_n", IN_EXON_C): "NM_1.1:n.201_203delinsGCG",
            ("c_to_g", IN_EXON_C): IN_EXON_G,
        }
    )
    assert coding_to_genomic(IN_EXON_C, coordinates) == IN_EXON_G


@pytest.mark.unit
@pytest.mark.parametrize(
    "hgvs_c, literal, n",
    [
        # A codon split by a 117-base intron (TP53 codon 25).
        (
            "NM_000546.6:c.73_75delinsAGG",
            "NC_000017.11:g.7676403_7676522delinsCCT",
            "NM_000546.6:n.215_217delinsAGG",
        ),
        ("NM_000546.6:c.73_75del", "NC_000017.11:g.7676403_7676522del", "NM_000546.6:n.215_217del"),
        # An insertion between two exons' bases.
        ("NM_1.1:c.1903_1904insGAA", "NC_000001.11:g.1000_1500insTTC", "NM_1.1:n.1953_1954insGAA"),
    ],
)
def test_a_change_across_an_exon_junction_has_no_genomic_projection(hgvs_c, literal, n):
    coordinates = _TableCoordinates({("c_to_g_literal", hgvs_c): literal, ("c_to_n", hgvs_c): n})
    assert coding_to_genomic(hgvs_c, coordinates) is None


@pytest.mark.unit
def test_a_single_position_skips_the_span_check():
    coordinates = _TableCoordinates(
        {("c_to_g_literal", "NM_1.1:c.74T>G"): "NC_1.1:g.100A>C", ("c_to_g", "NM_1.1:c.74T>G"): "NC_1.1:g.100A>C"}
    )
    assert coding_to_genomic("NM_1.1:c.74T>G", coordinates) == "NC_1.1:g.100A>C"


@pytest.mark.unit
def test_a_change_with_an_intronic_offset_is_taken_as_written():
    hgvs_c = "NM_1.1:c.123+2_123+5del"
    coordinates = _TableCoordinates(
        {
            ("c_to_g_literal", hgvs_c): "NC_1.1:g.1002_1005del",
            ("c_to_n", hgvs_c): "NM_1.1:n.173+2_173+5del",
            ("c_to_g", hgvs_c): "NC_1.1:g.1002_1005del",
        }
    )
    assert coding_to_genomic(hgvs_c, coordinates) == "NC_1.1:g.1002_1005del"


@pytest.mark.unit
def test_an_inversion_is_written_as_the_delins_on_its_trimmed_span():
    hgvs_c = "NM_1.1:c.151_153delinsCAT"
    coordinates = _TableCoordinates(
        {
            ("c_to_g_literal", hgvs_c): "NC_1.1:g.1000_1002delinsATG",
            ("c_to_n", hgvs_c): "NM_1.1:n.201_203delinsCAT",
            ("c_to_g", hgvs_c): "NC_1.1:g.1001_1002inv",
        }
    )
    assert coding_to_genomic(hgvs_c, coordinates) == "NC_1.1:g.1001_1002delinsTG"


@pytest.mark.unit
def test_an_inversion_that_cannot_be_trimmed_keeps_the_unnormalized_delins():
    hgvs_g = "NC_1.1:g.1000_1002delinsATG"
    coordinates = _TableCoordinates(
        {
            ("g_to_c_literal", hgvs_g): "NM_1.1:c.-3_-1delinsCAT",
            ("c_to_n", "NM_1.1:c.-3_-1delinsCAT"): "NM_1.1:n.48_50delinsCAT",
            ("g_to_c", hgvs_g): "NM_1.1:c.-2_-1inv",
        }
    )
    assert genomic_to_coding(hgvs_g, "NM_1.1", coordinates) == "NM_1.1:c.-3_-1delinsCAT"


@pytest.mark.unit
def test_a_variant_hgvs_rejects_has_no_projection():
    assert coding_to_genomic("NM_1.1:c.99999A>G", _TableCoordinates({})) is None


@pytest.mark.unit
@pytest.mark.parametrize(
    "error",
    [
        ConnectionError("cdot-rest unreachable"),
        hgvs_exceptions.HGVSError("Permanently lost connection to postgresql://uta"),
        # hgvs's sequence fetchers wrap every fetch failure, network errors included, in this.
        hgvs_exceptions.HGVSDataNotAvailableError("Failed to fetch NC_000017.11 (Connection refused)"),
    ],
)
def test_a_failure_to_reach_reference_data_is_raised(error):
    class _Unreachable(_TableCoordinates):
        def c_to_g_literal(self, c):
            raise error

    with pytest.raises(type(error)):
        coding_to_genomic(IN_EXON_C, _Unreachable({}))


@pytest.mark.unit
def test_a_genomic_change_covering_a_whole_intron_has_no_coding_projection():
    hgvs_g = "NC_000017.11:g.7676403_7676522del"
    coordinates = _TableCoordinates(
        {
            ("g_to_c_literal", hgvs_g): "NM_000546.6:c.73_75del",
            ("c_to_n", "NM_000546.6:c.73_75del"): "NM_000546.6:n.215_217del",
        }
    )
    assert genomic_to_coding(hgvs_g, "NM_000546.6", coordinates) is None


# ---------------------------------------------------------------------------
# Property: every codon-level change, through real hgvs mapping on synthetic
# transcripts, either projects to the same genomic sequence change or has no
# projection because it crosses an exon junction.
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_CONTIG = "NC_000001.11"
# Exon lengths in transcript order. With the coding sequence starting at transcript base 11, the three
# junctions fall after codon phases 0, 1 and 2.
_EXON_LENGTHS = (61, 61, 61, 100)
_INTRON_LENGTH = 117
_CDS_START_I = 10
_CDS_END_I = 280
_GENOME_OFFSET = 500


def _revcomp(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass
class _Exon:
    tx_start_i: int
    tx_end_i: int
    alt_start_i: int
    alt_end_i: int


class _SyntheticTranscripts:
    """A minimal hgvs data provider over one synthetic multi-exon transcript per strand."""

    def __init__(self, seed: int = 7):
        rng = random.Random(seed)
        span = sum(_EXON_LENGTHS) + _INTRON_LENGTH * (len(_EXON_LENGTHS) - 1)
        self.genome = "".join(rng.choice("ACGT") for _ in range(2 * _GENOME_OFFSET + span))
        self.exons: dict[str, list[_Exon]] = {}
        self.strands = {"NM_PLUS.1": 1, "NM_MINUS.1": -1}
        self.transcripts: dict[str, str] = {}
        for tx_ac, strand in self.strands.items():
            exons, tx_start, genomic_start = [], 0, _GENOME_OFFSET
            # Genomic order runs with the transcript on the plus strand and against it on the minus strand.
            lengths = _EXON_LENGTHS if strand == 1 else tuple(reversed(_EXON_LENGTHS))
            placed = []
            for length in lengths:
                placed.append((genomic_start, genomic_start + length))
                genomic_start += length + _INTRON_LENGTH
            ordered = placed if strand == 1 else list(reversed(placed))
            for (alt_start, alt_end), length in zip(ordered, _EXON_LENGTHS):
                exons.append(_Exon(tx_start, tx_start + length, alt_start, alt_end))
                tx_start += length
            self.exons[tx_ac] = exons
            spliced = "".join(
                self.genome[e.alt_start_i : e.alt_end_i] for e in sorted(exons, key=lambda e: e.alt_start_i)
            )
            self.transcripts[tx_ac] = spliced if strand == 1 else _revcomp(spliced)

    # hgvs data provider interface ------------------------------------------------------------------

    def data_version(self):
        return "synthetic"

    def schema_version(self):
        return "1.1"

    def get_assembly_map(self, assembly_name):
        return {_CONTIG: "1"}

    def get_seq(self, ac, start_i=None, end_i=None):
        sequence = self.genome if ac == _CONTIG else self.transcripts[ac]
        return sequence[start_i:end_i]

    def get_tx_mapping_options(self, tx_ac):
        return [{"tx_ac": tx_ac, "alt_ac": _CONTIG, "alt_aln_method": "splign"}]

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        return {
            "hgnc": "SYNTH",
            "cds_start_i": _CDS_START_I,
            "cds_end_i": _CDS_END_I,
            "tx_ac": tx_ac,
            "alt_ac": alt_ac,
            "alt_aln_method": alt_aln_method,
        }

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        strand = self.strands[tx_ac]
        return [
            {
                "tx_ac": tx_ac,
                "alt_ac": alt_ac,
                "alt_strand": strand,
                "alt_aln_method": alt_aln_method,
                "ord": order,
                "tx_start_i": exon.tx_start_i,
                "tx_end_i": exon.tx_end_i,
                "alt_start_i": exon.alt_start_i,
                "alt_end_i": exon.alt_end_i,
                "cigar": f"{exon.tx_end_i - exon.tx_start_i}=",
            }
            for order, exon in sorted(enumerate(self.exons[tx_ac]), key=lambda item: item[1].alt_start_i)
        ]

    def get_tx_identity_info(self, tx_ac):
        return {
            "tx_ac": tx_ac,
            "alt_ac": tx_ac,
            "alt_aln_method": "transcript",
            "cds_start_i": _CDS_START_I,
            "cds_end_i": _CDS_END_I,
            "lengths": list(_EXON_LENGTHS),
            "hgnc": "SYNTH",
        }

    def get_pro_ac_for_tx_ac(self, tx_ac):
        return None

    def get_gene_info(self, gene):
        raise NotImplementedError

    def get_acs_for_protein_seq(self, seq):
        raise NotImplementedError

    def get_similar_transcripts(self, tx_ac):
        raise NotImplementedError

    def get_tx_for_gene(self, gene):
        raise NotImplementedError

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        raise NotImplementedError

    # Model -----------------------------------------------------------------------------------------

    def genomic_index(self, tx_ac: str, c: int) -> tuple[int, int]:
        """The genomic index of coding position ``c`` and the transcript-order index of its exon."""
        tx = _CDS_START_I + c - 1
        for order, exon in enumerate(self.exons[tx_ac]):
            if exon.tx_start_i <= tx < exon.tx_end_i:
                offset = tx - exon.tx_start_i
                strand = self.strands[tx_ac]
                return (exon.alt_start_i + offset if strand == 1 else exon.alt_end_i - 1 - offset), order
        raise ValueError(c)

    def expected_genome(self, tx_ac: str, start: int, end: int, alt: str, insertion: bool) -> str | None:
        """The genome after a coding change, or ``None`` when its bases sit in more than one exon."""
        (g_start, start_exon), (g_end, end_exon) = self.genomic_index(tx_ac, start), self.genomic_index(tx_ac, end)
        if start_exon != end_exon:
            return None

        if self.strands[tx_ac] == -1:
            g_start, g_end, alt = g_end, g_start, _revcomp(alt)
        if insertion:
            return self.genome[: g_start + 1] + alt + self.genome[g_start + 1 :]
        return self.genome[:g_start] + alt + self.genome[g_end + 1 :]

    def apply_genomic(self, hgvs_g: str, parser) -> str:
        variant = parser.parse(hgvs_g)
        start, end = variant.posedit.pos.start.base - 1, variant.posedit.pos.end.base - 1
        edit = variant.posedit.edit
        if edit.type == "ins":
            return self.genome[: start + 1] + edit.alt + self.genome[start + 1 :]
        if edit.type == "dup":
            return self.genome[: end + 1] + self.genome[start : end + 1] + self.genome[end + 1 :]
        if edit.type == "identity":
            return self.genome
        assert edit.type in ("sub", "del", "delins"), f"unexpected edit type in {hgvs_g}"
        return self.genome[:start] + (edit.alt or "") + self.genome[end + 1 :]


def _codon_changes(transcript: str):
    """Codon-level changes like those reverse translation emits, with their coding span and inserted bases."""
    for codon in range(1, (_CDS_END_I - _CDS_START_I) // 3 + 1):
        start, end = 3 * codon - 2, 3 * codon
        ref = transcript[_CDS_START_I + start - 1 : _CDS_START_I + end]
        middle_kept = ("T" if ref[0] != "T" else "C") + ref[1] + ("G" if ref[2] != "G" else "A")
        for alt in {ref, _revcomp(ref), middle_kept, "TGA"}:
            yield f"c.{start}_{end}delins{alt}", start, end, alt, False
        yield f"c.{start}_{end}del", start, end, "", False
        if end < _CDS_END_I - _CDS_START_I:
            yield f"c.{end}_{end + 1}insGCA", end, end + 1, "GCA", True


@pytest.fixture(scope="module")
def synthetic():
    hgvs_parser = pytest.importorskip("hgvs.parser")
    from variant_annotation.lib.clients.coordinates import HgvsMapper

    provider = _SyntheticTranscripts()
    return provider, HgvsMapper(provider), hgvs_parser.Parser()


@pytest.mark.slow
@pytest.mark.parametrize("tx_ac", ["NM_PLUS.1", "NM_MINUS.1"])
def test_every_codon_change_projects_to_the_same_genomic_change_or_none(synthetic, tx_ac):
    provider, coordinates, parser = synthetic
    split = 0
    for change, start, end, alt, insertion in _codon_changes(provider.transcripts[tx_ac]):
        hgvs_c = f"{tx_ac}:{change}"
        expected = provider.expected_genome(tx_ac, start, end, alt, insertion)
        projected = coding_to_genomic(hgvs_c, coordinates)

        if expected is None:
            split += 1
            assert projected is None, f"{hgvs_c} crosses an exon junction but projected to {projected}"
            continue

        assert projected is not None, f"{hgvs_c} sits in one exon but has no projection"
        assert "inv" not in projected and ".[" not in projected, projected
        assert provider.apply_genomic(projected, parser) == expected, f"{hgvs_c} -> {projected}"

        # The genomic form projects back to a coding change with the same effect on the genome.
        back = genomic_to_coding(projected, tx_ac, coordinates)
        assert back is not None, f"{projected} has no coding projection"

    # Every junction lands inside a codon or an insertion site, so some changes must be split.
    assert split > 0
