"""Sequence Ontology molecular-consequence vocabulary — the canonical term set and severity order.

Shared vocabulary: pure data and pure functions, no service dependencies. Any ``lib/`` module and any
external consumer (mavedb-api) may import from here.

The term set is Ensembl's published *calculated variant consequences* table, which is the closed set
VEP can emit and the ranking VEP's own ``most_severe_consequence`` field uses. It is fetched from
``https://rest.ensembl.org/info/variation/consequence_types?rank=1`` and reproduced here so ranking is
deterministic and offline — a severity comparison must not depend on a network call.

Terms are SO terms; SO accessions are carried alongside so a consumer can emit a CURIE
(``SO:0001583``) rather than a bare label. Ensembl's ranking is *not* the SO hierarchy — it is
Ensembl's own severity judgement over SO terms — so ``rank`` is documented as "Ensembl severity", not
"SO severity".

Why this lives in the library rather than in each caller: mavedb-api and the offline CLI pipeline must
agree on which of a transcript's consequence terms is "the" consequence. Two hand-maintained copies of
a 41-term ordered list diverge silently — and did (mavedb-api's copy carried 4 duplicates, 18 terms
Ensembl never emits, two non-canonical spellings, and had ``feature_elongation``/``feature_truncation``
30 places below their true rank).

Vocabulary note: these terms describe a variant's predicted effect on a transcript — its **molecular
consequence**. That is deliberately not "functional consequence": in MaveDB's domain a *functional*
effect is what an assay measured (VA-Spec "functional impact"), which is a different fact from a
different source. Callers should carry the molecular-consequence name through to storage and display.

Maintenance: Ensembl adds terms across releases (the three fine-grained splice terms arrived in 105).
When a release adds or reorders terms, refresh :data:`ENSEMBL_CONSEQUENCE_RANKING` from the endpoint
above and bump :data:`RESOLVER_VERSION` in :mod:`variant_annotation.lib.vep` so downstream caches
know their stored answers were computed under a superseded ranking.
"""

from __future__ import annotations

from typing import Iterable, Optional

#: Ensembl's calculated variant consequences, most severe first, as ``(SO term, SO accession)``.
#: Verified against ``/info/variation/consequence_types?rank=1`` on Ensembl release 116 (2026-08-12).
#: Order is load-bearing: index == severity rank - 1.
ENSEMBL_CONSEQUENCE_RANKING: tuple[tuple[str, str], ...] = (
    ("transcript_ablation", "SO:0001893"),
    ("splice_acceptor_variant", "SO:0001574"),
    ("splice_donor_variant", "SO:0001575"),
    ("stop_gained", "SO:0001587"),
    ("frameshift_variant", "SO:0001589"),
    ("stop_lost", "SO:0001578"),
    ("start_lost", "SO:0002012"),
    ("transcript_amplification", "SO:0001889"),
    ("feature_elongation", "SO:0001907"),
    ("feature_truncation", "SO:0001906"),
    ("inframe_insertion", "SO:0001821"),
    ("inframe_deletion", "SO:0001822"),
    ("missense_variant", "SO:0001583"),
    ("protein_altering_variant", "SO:0001818"),
    ("splice_donor_5th_base_variant", "SO:0001787"),
    ("splice_region_variant", "SO:0001630"),
    ("splice_donor_region_variant", "SO:0002170"),
    ("splice_polypyrimidine_tract_variant", "SO:0002169"),
    ("incomplete_terminal_codon_variant", "SO:0001626"),
    ("start_retained_variant", "SO:0002019"),
    ("stop_retained_variant", "SO:0001567"),
    ("synonymous_variant", "SO:0001819"),
    ("coding_sequence_variant", "SO:0001580"),
    ("mature_miRNA_variant", "SO:0001620"),
    ("5_prime_UTR_variant", "SO:0001623"),
    ("3_prime_UTR_variant", "SO:0001624"),
    ("non_coding_transcript_exon_variant", "SO:0001792"),
    ("intron_variant", "SO:0001627"),
    ("NMD_transcript_variant", "SO:0001621"),
    ("non_coding_transcript_variant", "SO:0001619"),
    ("coding_transcript_variant", "SO:0001968"),
    ("upstream_gene_variant", "SO:0001631"),
    ("downstream_gene_variant", "SO:0001632"),
    ("TFBS_ablation", "SO:0001895"),
    ("TFBS_amplification", "SO:0001892"),
    ("TF_binding_site_variant", "SO:0001782"),
    ("regulatory_region_ablation", "SO:0001894"),
    ("regulatory_region_amplification", "SO:0001891"),
    ("regulatory_region_variant", "SO:0001566"),
    ("intergenic_variant", "SO:0001628"),
    ("sequence_variant", "SO:0001060"),
)

#: The SO terms alone, most severe first. The ordered vocabulary most callers want.
CONSEQUENCE_TERMS: tuple[str, ...] = tuple(term for term, _ in ENSEMBL_CONSEQUENCE_RANKING)

#: SO term -> SO accession, for emitting CURIEs alongside labels.
SO_ACCESSIONS: dict[str, str] = dict(ENSEMBL_CONSEQUENCE_RANKING)

#: SO term -> Ensembl severity rank (1 = most severe). Built from the tuple so the two cannot drift.
CONSEQUENCE_RANKS: dict[str, int] = {term: index + 1 for index, term in enumerate(CONSEQUENCE_TERMS)}

#: Rank assigned to a term Ensembl does not publish — sorts after every known term without discarding
#: it. A new Ensembl release emitting an unknown term degrades to "least severe", never to "no answer".
UNRANKED = len(CONSEQUENCE_TERMS) + 1


def is_known_term(term: str) -> bool:
    """Whether ``term`` is in Ensembl's published consequence set."""
    return term in CONSEQUENCE_RANKS


def rank(term: str) -> int:
    """Ensembl severity rank of ``term`` (1 = most severe), or :data:`UNRANKED` if unpublished.

    Total over any string, so it is safe as a ``sort`` key on raw upstream output.
    """
    return CONSEQUENCE_RANKS.get(term, UNRANKED)


def most_severe(terms: Iterable[str]) -> Optional[str]:
    """The most severe of ``terms`` by Ensembl rank, or ``None`` when ``terms`` is empty.

    Unpublished terms are ranked last but never dropped: an unrecognised term is still a real answer
    from VEP, so a set consisting only of unknown terms returns one of them (the first in input order)
    rather than ``None``. ``None`` means "nothing to rank", which callers read as "VEP said nothing" —
    conflating the two would turn a vocabulary gap into a spurious missing annotation.
    """
    ordered = sort_by_severity(terms)
    return ordered[0] if ordered else None


def sort_by_severity(terms: Iterable[str]) -> list[str]:
    """``terms`` ordered most severe first, de-duplicated, unpublished terms last.

    Ties among unpublished terms keep their input order (``sorted`` is stable), so the output is
    deterministic for a given input.
    """
    seen: dict[str, None] = dict.fromkeys(terms)
    return sorted(seen, key=rank)
