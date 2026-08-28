"""Molecular-consequence resolution over Ensembl VEP responses.

Two layers, both public through the package:

- The **kernel** (:func:`resolve_entry`, :func:`requested_transcript`, :func:`partition_batches`,
  :func:`combine_recoded`) is pure — no I/O, no ports. It is the part that must produce identical
  answers in every caller, so it is separable and directly unit-testable against literal VEP payloads.
- The **orchestration** (:func:`resolve_consequences`) drives the kernel over the injected ports,
  handling batching, the Recoder fallback, and failure attribution.

A caller that already owns its own transport (mavedb-api's worker runs batches on an asyncio loop with
its own progress reporting and rate limiting) can import the kernel alone and keep its transport, and
still be guaranteed the same consequence for the same input. A caller with no transport opinion uses
:func:`resolve_consequences`.

Why the transcript match matters
--------------------------------
VEP annotates every transcript overlapping a variant's genomic position, not just the transcript the
input HGVS names. Its top-level ``most_severe_consequence`` is therefore the worst call across *all*
of them — a single BRCA1 coding variant returns 368 ``transcript_consequences`` entries, and the
headline term routinely describes a transcript with a different reading frame or exon structure than
the one the assay was performed against. Resolution against the input's own transcript is the
difference between "this variant is synonymous in the gene we assayed" and "this variant is
splice-disrupting in some other overlapping isoform".
"""

from __future__ import annotations

import logging
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Iterable, Mapping, Optional, Sequence

from ..accessions import (
    extract_accession,
    looks_like_refseq_transcript_accession,
    looks_like_transcript_accession,
    strip_accession_version,
)
from ..hgvs.no_change import parse_coding_delins, states_no_change
from ..sequence_ontology import most_severe, sort_by_severity
from ._ports import ReferenceSequence, VariantRecoder, VepLookup
from .types import (
    NO_CHANGE_TERM,
    ConsequenceOutcome,
    ConsequenceResolution,
    ConsequenceSource,
    VepConfig,
    VepInput,
)

logger = logging.getLogger(__name__)

#: Cap on a stored error description. Bounded so a verbose upstream HTML error page cannot land whole
#: in a database column or a CSV cell.
_MAX_ERROR_LENGTH = 300


def _truncate_error(message: str) -> str:
    """Bound an error description's length, marking it when truncated.

    Transport-level, not format-level: escaping for a particular sink (TSV delimiters, JSON) belongs to
    that sink's adapter, not here.
    """
    text = " ".join((message or "").split())
    if len(text) <= _MAX_ERROR_LENGTH:
        return text
    return text[:_MAX_ERROR_LENGTH] + "..."


def requested_transcript(vep_input: VepInput) -> Optional[str]:
    """The version-stripped transcript ``vep_input``'s consequence should be read from.

    An explicit ``VepInput.transcript`` wins; otherwise the transcript is inferred from the HGVS
    accession, which succeeds for coding and protein input and fails for genomic input (an ``NC_``
    accession names a chromosome). ``None`` means no transcript-specific resolution is possible and the
    caller will fall back to VEP's cross-transcript headline.

    A protein accession is *not* usable here even though it identifies the right gene product: VEP keys
    ``transcript_consequences`` on transcript ids, so ``NP_…`` would never match. Callers holding a
    protein allele should resolve the consequence from its coding representation instead.
    """
    if vep_input.transcript:
        return strip_accession_version(vep_input.transcript.strip())

    accession = strip_accession_version(extract_accession(vep_input.hgvs))
    return accession if accession and looks_like_transcript_accession(accession) else None


def needs_refseq_transcripts(vep_input: VepInput) -> bool:
    """Whether ``vep_input`` must be queried with Ensembl's ``refseq=1`` transcript set.

    Keyed on the *effective* transcript from :func:`requested_transcript`, not on the HGVS prefix. That
    distinction is the whole point: a genomic ``NC_…:g.`` input carrying a RefSeq transcript hint needs
    ``refseq=1`` so ``transcript_consequences`` comes back with ``NM_`` ids the hint can match. Keying
    on the HGVS prefix alone — as a naive implementation does — sends that input to the Ensembl
    transcript set, guaranteeing the match fails and the answer silently degrades to ``most_severe``.
    """
    transcript = requested_transcript(vep_input)
    return bool(transcript) and looks_like_refseq_transcript_accession(transcript or "")


def partition_batches(inputs: Sequence[VepInput], *, batch_size: int) -> list[tuple[tuple[VepInput, ...], bool]]:
    """Split ``inputs`` into ``(batch, refseq)`` request groups.

    Inputs are grouped by which transcript set they need before being chunked, because ``refseq`` is a
    per-request flag: mixing both kinds in one batch would force one group to be queried against the
    wrong transcript set. Order within each group is preserved so a caller can correlate progress
    reporting with its input order.

    Duplicate HGVS strings are *not* collapsed here — two inputs may share an HGVS string but request
    different transcripts, which are genuinely different questions. :func:`resolve_consequences`
    de-duplicates at the wire level instead, by the ``(hgvs, refseq)`` pair actually sent.
    """
    if batch_size < 1:
        raise ValueError(f"batch_size must be at least 1, got {batch_size}")

    refseq_inputs = tuple(i for i in inputs if needs_refseq_transcripts(i))
    ensembl_inputs = tuple(i for i in inputs if not needs_refseq_transcripts(i))

    batches: list[tuple[tuple[VepInput, ...], bool]] = []
    for group, refseq in ((refseq_inputs, True), (ensembl_inputs, False)):
        for start in range(0, len(group), batch_size):
            batches.append((group[start : start + batch_size], refseq))
    return batches


def resolve_entry(vep_input: VepInput, entry: Mapping[str, Any]) -> ConsequenceResolution:
    """Resolve one VEP response ``entry`` into a consequence for ``vep_input``. Pure.

    Resolution order:

    1. When a transcript is requested (:func:`requested_transcript`), scan ``transcript_consequences``
       for the entry whose ``transcript_id`` matches it with versions stripped. On a hit, the answer is
       the most severe of *that entry's* ``consequence_terms`` and the source is ``TRANSCRIPT``.
    2. Otherwise fall back to the top-level ``most_severe_consequence``, source ``MOST_SEVERE``.
    3. With neither, the outcome is ``ABSENT`` — VEP answered and had nothing for this input.

    A matched transcript entry carrying an empty ``consequence_terms`` falls through to step 2 rather
    than resolving to ``ABSENT``: the match proves the transcript is in VEP's set, so an empty term list
    is a gap in that entry, not a statement that the variant has no consequence anywhere.

    Only the *first* matching transcript entry is used. VEP emits one entry per (transcript, variant)
    pair, so a version-stripped transcript id matches at most once in practice; if a response ever
    carried two versions of the same transcript, taking the first is deterministic given VEP's stable
    response order.
    """
    wanted = requested_transcript(vep_input)

    if wanted:
        for candidate in entry.get("transcript_consequences") or ():
            if not isinstance(candidate, Mapping):
                continue
            candidate_id = str(candidate.get("transcript_id") or "")
            if strip_accession_version(candidate_id) != wanted:
                continue

            terms = sort_by_severity(str(t) for t in (candidate.get("consequence_terms") or ()))
            if not terms:
                break

            return ConsequenceResolution(
                input=vep_input,
                outcome=ConsequenceOutcome.RESOLVED,
                consequence_terms=terms,
                most_severe_consequence=terms[0],
                source=ConsequenceSource.TRANSCRIPT,
                matched_transcript=candidate_id or None,
            )

    headline = entry.get("most_severe_consequence")
    if headline:
        return ConsequenceResolution(
            input=vep_input,
            outcome=ConsequenceOutcome.RESOLVED,
            consequence_terms=[str(headline)],
            most_severe_consequence=str(headline),
            source=ConsequenceSource.MOST_SEVERE,
        )

    return ConsequenceResolution(input=vep_input, outcome=ConsequenceOutcome.ABSENT)


def combine_recoded(vep_input: VepInput, resolutions: Iterable[ConsequenceResolution]) -> ConsequenceResolution:
    """Collapse the resolutions of ``vep_input``'s recoded genomic forms into one answer. Pure.

    Variant Recoder can map one input to several genomic expressions (a protein change recodes to every
    codon-level DNA change encoding it). Each is queried independently, and the combined answer is the
    most severe consequence across those that resolved — the conservative reading, since any of the
    genomic forms is a valid spelling of the input.

    The source is always ``MOST_SEVERE``: a recoded genomic query carries no transcript, so no
    transcript-specific claim can be made about the result even when the original input named one.

    With nothing resolved, an ``ERRORED`` member makes the whole answer ``ERRORED`` (unknown, retry)
    and only an all-``ABSENT`` set makes it ``ABSENT`` (a confirmed empty). A partial failure must not
    be recorded as a negative.
    """
    resolved_terms: list[str] = []
    first_error: Optional[str] = None

    for resolution in resolutions:
        if resolution.outcome is ConsequenceOutcome.RESOLVED and resolution.most_severe_consequence:
            resolved_terms.append(resolution.most_severe_consequence)
        elif resolution.outcome is ConsequenceOutcome.ERRORED and first_error is None:
            first_error = resolution.error

    if resolved_terms:
        winner = most_severe(resolved_terms)
        return ConsequenceResolution(
            input=vep_input,
            outcome=ConsequenceOutcome.RESOLVED,
            consequence_terms=[winner] if winner else [],
            most_severe_consequence=winner,
            source=ConsequenceSource.MOST_SEVERE,
        )

    if first_error is not None:
        return ConsequenceResolution(input=vep_input, outcome=ConsequenceOutcome.ERRORED, error=first_error)

    return ConsequenceResolution(input=vep_input, outcome=ConsequenceOutcome.ABSENT)


def reference_identical_resolution(vep_input: VepInput) -> ConsequenceResolution:
    """The resolution for an input that describes no sequence change. Pure.

    Carries :data:`~variant_annotation.lib.vep.types.NO_CHANGE_TERM` and source
    ``REFERENCE_IDENTICAL`` — an explicit label, so downstream can tell "wild type by construction"
    apart from "annotation has not run" and from ``synonymous_variant``.
    """
    return ConsequenceResolution(
        input=vep_input,
        outcome=ConsequenceOutcome.RESOLVED,
        consequence_terms=[NO_CHANGE_TERM],
        most_severe_consequence=NO_CHANGE_TERM,
        source=ConsequenceSource.REFERENCE_IDENTICAL,
    )


def _dispatch_batches(
    batches: Sequence[tuple[tuple[VepInput, ...], bool]],
    *,
    vep: VepLookup,
    max_workers: int,
) -> tuple[dict[tuple[str, bool], Mapping[str, Any]], dict[tuple[str, bool], str]]:
    """Run every batch and return ``(entries, errors)`` keyed by ``(hgvs, refseq)``.

    A batch-level failure is attributed to that batch's inputs only, so one bad batch degrades its own
    members to ``ERRORED`` while every other batch's results stand. Aborting the whole run on a single
    transport blip would discard work already paid for against a rate-limited service.
    """
    entries: dict[tuple[str, bool], Mapping[str, Any]] = {}
    errors: dict[tuple[str, bool], str] = {}

    def run_one(
        batch: tuple[VepInput, ...], refseq: bool
    ) -> tuple[tuple[VepInput, ...], bool, Sequence[Mapping[str, Any]], Optional[str]]:
        hgvs = list(dict.fromkeys(i.hgvs for i in batch))
        try:
            return batch, refseq, vep.annotate_hgvs(hgvs, refseq=refseq), None
        except Exception as exc:  # noqa: BLE001 - any transport failure is an unknown, not a negative
            logger.warning("VEP request failed for a batch of %d HGVS: %s", len(hgvs), exc)
            return batch, refseq, (), f"VEP request failed: {exc}"

    if max_workers > 1 and len(batches) > 1:
        with ThreadPoolExecutor(max_workers=min(max_workers, len(batches))) as pool:
            completed = list(pool.map(lambda args: run_one(*args), batches))
    else:
        completed = [run_one(batch, refseq) for batch, refseq in batches]

    for batch, refseq, response, error in completed:
        for entry in response:
            submitted = entry.get("input")
            if submitted:
                entries[(str(submitted), refseq)] = entry
        if error is not None:
            for vep_input in batch:
                key = (vep_input.hgvs, refseq)
                if key not in entries:
                    errors[key] = _truncate_error(error)

    return entries, errors


def resolve_consequences(
    inputs: Sequence[VepInput],
    *,
    vep: VepLookup,
    recoder: Optional[VariantRecoder] = None,
    reference: Optional[ReferenceSequence] = None,
    config: Optional[VepConfig] = None,
) -> list[ConsequenceResolution]:
    """Resolve a molecular consequence for every input, one resolution per input in input order.

    Three passes:

    1. Query VEP, partitioned by transcript set, and resolve each answered input against its own
       transcript (:func:`resolve_entry`).
    2. For inputs VEP neither answered nor errored on — protein HGVS it cannot parse, exotic coding
       forms it rejects — recode to genomic via ``recoder`` and query again, combining the recovered
       forms (:func:`combine_recoded`). Skipped when ``recoder`` is ``None`` or
       ``config.recode_misses`` is false; those inputs stay ``ABSENT``.
    3. For inputs still unresolved, check whether they describe no sequence change and label them
       ``REFERENCE_IDENTICAL`` (requires ``reference``).

    An input whose request failed is reported ``ERRORED`` and is never recoded: the answer is unknown,
    so spending a Recoder round trip on it would be attributing a fallback to a question that was never
    asked. Callers must retry those rather than storing the outcome.

    ``reference`` is consulted only for inputs that reached step 3, so the no-change check costs nothing
    for the overwhelming majority of inputs that VEP resolves normally.
    """
    settings = config or VepConfig()
    if not inputs:
        return []

    entries, errors = _dispatch_batches(
        partition_batches(inputs, batch_size=settings.batch_size),
        vep=vep,
        max_workers=settings.max_workers,
    )

    resolutions: dict[int, ConsequenceResolution] = {}
    unanswered: list[tuple[int, VepInput]] = []

    for index, vep_input in enumerate(inputs):
        key = (vep_input.hgvs, needs_refseq_transcripts(vep_input))
        entry = entries.get(key)
        if entry is not None:
            resolution = resolve_entry(vep_input, entry)
            if resolution.outcome is ConsequenceOutcome.RESOLVED:
                resolutions[index] = resolution
                continue
            # VEP answered but classified nothing: still a miss worth recoding.
            unanswered.append((index, vep_input))
            continue

        error = errors.get(key)
        if error is not None:
            resolutions[index] = ConsequenceResolution(input=vep_input, outcome=ConsequenceOutcome.ERRORED, error=error)
            continue

        unanswered.append((index, vep_input))

    if unanswered and recoder is not None and settings.recode_misses:
        resolutions.update(_resolve_via_recoder(unanswered, vep=vep, recoder=recoder, settings=settings))

    for index, vep_input in unanswered:
        if index in resolutions:
            continue
        if reference is not None and _is_reference_identical(vep_input, reference):
            resolutions[index] = reference_identical_resolution(vep_input)
        else:
            resolutions[index] = ConsequenceResolution(input=vep_input, outcome=ConsequenceOutcome.ABSENT)

    return [resolutions[index] for index in range(len(inputs))]


def _resolve_via_recoder(
    unanswered: Sequence[tuple[int, VepInput]],
    *,
    vep: VepLookup,
    recoder: VariantRecoder,
    settings: VepConfig,
) -> dict[int, ConsequenceResolution]:
    """Recode the unanswered inputs to genomic HGVS and re-query VEP with them.

    Returns resolutions only for inputs the fallback settled. An input the Recoder had no genomic form
    for is left out, so the caller's reference-identical check and ``ABSENT`` default still apply. An
    input whose Recoder request *failed* is returned ``ERRORED`` — unknown, not a confirmed empty.
    """
    to_recode = list(dict.fromkeys(vep_input.hgvs for _, vep_input in unanswered))

    # Batch the Recoder POST the way VEP is batched: Ensembl bounds a single POST, so one oversized
    # request would 400 and strand every miss as ERRORED. A failed chunk marks only its own HGVS unknown
    # (as _dispatch_batches does for VEP), so one bad batch does not sink the others.
    recoded: dict[str, Sequence[str]] = {}
    recoder_errors: dict[str, str] = {}
    for start in range(0, len(to_recode), settings.batch_size):
        chunk = to_recode[start : start + settings.batch_size]
        try:
            recoded.update(recoder.to_genomic(chunk))
        except Exception as exc:  # noqa: BLE001 - a Recoder failure leaves those answers unknown
            logger.warning("Variant Recoder request failed for %d HGVS: %s", len(chunk), exc)
            error = _truncate_error(f"Variant Recoder request failed: {exc}")
            for hgvs in chunk:
                recoder_errors[hgvs] = error

    # Genomic forms carry no transcript, so every recoded query is a plain most-severe lookup.
    genomic_inputs = [VepInput(hgvs=genomic) for forms in recoded.values() for genomic in dict.fromkeys(forms)]

    genomic_entries: dict[tuple[str, bool], Mapping[str, Any]] = {}
    genomic_errors: dict[tuple[str, bool], str] = {}
    if genomic_inputs:
        genomic_entries, genomic_errors = _dispatch_batches(
            partition_batches(genomic_inputs, batch_size=settings.batch_size),
            vep=vep,
            max_workers=settings.max_workers,
        )

    per_genomic: dict[str, ConsequenceResolution] = {}
    for genomic_input in genomic_inputs:
        key = (genomic_input.hgvs, False)
        entry = genomic_entries.get(key)
        if entry is not None:
            per_genomic[genomic_input.hgvs] = resolve_entry(genomic_input, entry)
        elif key in genomic_errors:
            per_genomic[genomic_input.hgvs] = ConsequenceResolution(
                input=genomic_input,
                outcome=ConsequenceOutcome.ERRORED,
                error=genomic_errors[key],
            )

    by_hgvs: dict[str, list[tuple[int, VepInput]]] = defaultdict(list)
    for index, vep_input in unanswered:
        by_hgvs[vep_input.hgvs].append((index, vep_input))

    settled: dict[int, ConsequenceResolution] = {}
    for hgvs, members in by_hgvs.items():
        forms = [f for f in dict.fromkeys(recoded.get(hgvs) or ()) if f in per_genomic]
        if forms:
            for index, vep_input in members:
                settled[index] = combine_recoded(vep_input, (per_genomic[f] for f in forms))
        elif hgvs in recoder_errors:
            # The Recoder request covering this HGVS failed: unknown, not a confirmed empty. An HGVS the
            # Recoder simply had no genomic form for is left out instead, to fall through to ABSENT.
            for index, vep_input in members:
                settled[index] = ConsequenceResolution(
                    input=vep_input, outcome=ConsequenceOutcome.ERRORED, error=recoder_errors[hgvs]
                )

    return settled


def _is_reference_identical(vep_input: VepInput, reference: ReferenceSequence) -> bool:
    """Whether ``vep_input`` describes no sequence change, per the transcript reference.

    Recognises the coding ``delins`` form the mapping pipeline produces for an unchanged interval
    (``NM_000049.4:c.12_14delinsGCT`` where the transcript already reads ``GCT`` there). Explicit HGVS
    no-change forms (``c.123=``) are recognised without consulting the reference — the notation itself
    is the claim.
    """
    if states_no_change(vep_input.hgvs):
        return True

    delins = parse_coding_delins(vep_input.hgvs)
    if delins is None:
        return False

    transcript, start, stop, alt = delins
    if stop - start + 1 != len(alt):
        return False

    current = reference.coding_interval_reference(transcript, start, stop)
    return current is not None and current.upper() == alt.upper()
