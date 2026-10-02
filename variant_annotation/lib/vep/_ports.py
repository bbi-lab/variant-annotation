"""Consumer-owned narrow protocols for VEP resolution's backing services.

Defined here (in the consumer) so that ``clients/`` and mavedb-api can satisfy them structurally
without importing from ``vep/``. The arrows point down: ``_core`` imports these; ``clients/`` does not.

The protocols are deliberately pitched at the *raw response* level rather than returning resolved
consequences. The resolution rule — which ``transcript_consequences`` entry to read, which of its terms
wins — is the part that must be identical in every caller, so it stays in ``_core`` where it is pure and
unit-testable. Only the HTTP round trip is injected.
"""

from __future__ import annotations

from typing import Any, Mapping, Protocol, Sequence


class VepLookup(Protocol):
    """Ensembl VEP HGVS endpoint."""

    def annotate_hgvs(self, hgvs: Sequence[str], *, refseq: bool) -> Sequence[Mapping[str, Any]]:
        """Return VEP's raw response entries for ``hgvs``, unmodified.

        Each entry is expected to carry ``input`` (the submitted string, VEP echoes it back),
        ``most_severe_consequence``, and optionally ``transcript_consequences``. Entries for inputs VEP
        could not resolve may be absent from the response entirely — the caller detects those by
        comparing against what it submitted, so an implementation must not fabricate placeholders.

        ``refseq`` requests RefSeq transcript ids in ``transcript_consequences`` (Ensembl's
        ``refseq=1``). Implementations must pass it through rather than deciding for themselves: which
        transcript set to ask for is a function of the input accessions, which ``_core`` partitions on.

        Raises:
            Exception: any transport or HTTP failure, after whatever retry policy the implementation
                applies. ``_core`` catches broadly, attributes the failure to the batch's inputs, and
                reports them ``ERRORED`` — so a failure must surface as a raise, never as an empty
                response, which is indistinguishable from a genuine set of misses.
        """
        ...


class VariantRecoder(Protocol):
    """Ensembl Variant Recoder endpoint — the fallback for HGVS that VEP cannot parse."""

    def to_genomic(self, hgvs: Sequence[str]) -> Mapping[str, Sequence[str]]:
        """Return a mapping from each input HGVS to its genomic (``NC_``) HGVS equivalents.

        An input with no genomic equivalent is simply absent from the mapping — that is a genuine
        empty, not an error. One input may map to several genomic forms (an ambiguous protein change
        recodes to every codon-level DNA change encoding it), so the value is a sequence.

        Raises:
            Exception: any transport or HTTP failure, after retries. See :meth:`VepLookup.annotate_hgvs`.
        """
        ...


class ReferenceSequence(Protocol):
    """Transcript reference-sequence facts, used to recognise reference-identical input.

    Optional: :func:`~variant_annotation.lib.vep.resolve_consequences` skips reference-identical
    detection when no implementation is supplied.
    """

    def coding_interval_reference(self, transcript: str, start: int, stop: int) -> str | None:
        """Return the reference bases of ``transcript`` over the inclusive 1-based c. interval.

        Returns ``None`` when the interval cannot be resolved (unknown transcript, out of range).
        """
        ...
