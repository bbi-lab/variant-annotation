"""Public types for the VEP consequence-resolution API."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

#: Label for a reference-identical ("no change") input — an HGVS expression describing no sequence
#: change, which VEP rejects because it is not a sequence variant. Deliberately **not** an SO term:
#: ``synonymous_variant`` means "a change that leaves the amino acid intact", which is a different
#: claim from "no change was described". Reference-identical rows are real and load-bearing in
#: functional-assay data (a wild-type control well), so they get an explicit label rather than a null
#: that cannot be told apart from "not yet annotated".
NO_CHANGE_TERM = "no_change"


class ConsequenceSource(str, Enum):
    """Which resolution path produced a consequence — provenance, not severity.

    Recorded on every resolution so a stored consequence can be audited without re-querying, and so a
    consumer can tell a transcript-specific answer from a cross-transcript one.

    - ``TRANSCRIPT`` — read from the ``transcript_consequences`` entry matching the input's transcript.
      The trustworthy case: the consequence applies to the transcript the allele actually lives on.
    - ``MOST_SEVERE`` — VEP's top-level ``most_severe_consequence``, which is the worst call across
      *every* transcript overlapping the position. Used only when no transcript could be matched
      (genomic input with no transcript hint, or a transcript absent from VEP's set). Treat as a
      lower-confidence answer: it regularly reports a consequence that does not apply to the
      transcript of interest.
    - ``REFERENCE_IDENTICAL`` — the input describes no sequence change; see :data:`NO_CHANGE_TERM`.
      No VEP call was made.
    """

    TRANSCRIPT = "transcript"
    MOST_SEVERE = "most_severe"
    REFERENCE_IDENTICAL = "reference_identical"


class ConsequenceOutcome(str, Enum):
    """Whether an input got an answer, a confirmed non-answer, or no answer at all.

    The three states are kept distinct because they demand different handling downstream, and
    collapsing any two of them loses information that cannot be recovered later:

    - ``RESOLVED`` — a consequence was determined. ``most_severe_consequence`` is set.
    - ``ABSENT`` — VEP was queried successfully and returned no consequence for this input. A settled
      negative; re-querying under the same source version will return the same nothing.
    - ``ERRORED`` — the VEP or Recoder *request* failed (transport error, non-200 after retries). The
      answer is **unknown**, not negative. Callers must retry rather than storing a null, and must not
      let this overwrite a previously held consequence.
    """

    RESOLVED = "resolved"
    ABSENT = "absent"
    ERRORED = "errored"


@dataclass(frozen=True)
class VepInput:
    """One HGVS expression to resolve a molecular consequence for.

    ``hgvs`` may be genomic (``NC_…:g.``), coding (``NM_…:c.`` / ``ENST…:c.``), or protein
    (``NP_…:p.``).

    ``transcript`` names the transcript the consequence should be read from. When omitted it is
    inferred from ``hgvs`` itself, which works for coding and protein input but **not** for genomic
    input — an ``NC_`` accession is a chromosome, not a transcript, so a genomic variant with no hint
    can only fall back to ``most_severe`` across every overlapping transcript. Pass the transcript
    explicitly for genomic input whenever the caller knows it: in mavedb-api that is the coding allele
    paired to the genomic one through ``projection_group``, which is exactly the transcript the assay
    was performed against. Versions are ignored when matching, so ``NM_000049`` and ``NM_000049.4``
    are equivalent here.
    """

    hgvs: str
    transcript: Optional[str] = None


@dataclass(frozen=True)
class ConsequenceResolution:
    """The resolved molecular consequence for one :class:`VepInput`.

    ``consequence_terms`` is every term VEP reported for the matched transcript, severity-ordered.
    It is a single-element list when the answer came from ``most_severe`` (the top-level field carries
    one term, not a set) and empty when the outcome is not ``RESOLVED``.

    ``most_severe_consequence`` is the single headline term — ``consequence_terms[0]``, or the
    top-level VEP field, or :data:`NO_CHANGE_TERM`. ``None`` unless the outcome is ``RESOLVED``.

    ``matched_transcript`` is the versioned transcript id as VEP reported it, set only when
    ``source`` is ``TRANSCRIPT``. Retained rather than assumed equal to ``input.transcript``: the match
    is version-insensitive, so recording what VEP actually used shows which transcript *version*
    produced the call.

    ``error`` carries a sanitized failure description when the outcome is ``ERRORED``, and is ``None``
    otherwise. It is never populated for ``ABSENT`` — a confirmed empty is not a failure.
    """

    input: VepInput
    outcome: ConsequenceOutcome
    consequence_terms: list[str] = field(default_factory=list)
    most_severe_consequence: Optional[str] = None
    source: Optional[ConsequenceSource] = None
    matched_transcript: Optional[str] = None
    error: Optional[str] = None


@dataclass(frozen=True)
class VepConfig:
    """Resolution behaviour knobs — everything except the port implementations.

    ``batch_size`` is capped by Ensembl at 200 HGVS notations per ``/vep/human/hgvs`` POST; the
    default sits at that ceiling because every request costs a round trip against a rate-limited
    public service.

    ``max_workers`` bounds concurrent in-flight batches. Ensembl rate-limits per client, so raising
    this trades throughput against 429s; the conservative default keeps a single-threaded caller's
    behaviour unchanged.

    ``recode_misses`` controls the Variant Recoder fallback. VEP cannot parse protein HGVS and rejects
    some exotic coding forms, so an input it could not resolve is recoded to genomic and re-queried.
    The recovered consequence is necessarily ``MOST_SEVERE`` — a genomic query carries no transcript —
    so a caller that would rather record ``ABSENT`` than a cross-transcript guess can turn it off.
    """

    batch_size: int = 200
    max_workers: int = 1
    recode_misses: bool = True
