"""variant_annotation.vep — public API for molecular-consequence resolution.

Composition roots (``pipeline/``, mavedb-api) import from here.

Two entry points, both guaranteeing the same answer for the same input:

- :func:`resolve_consequences` — full orchestration over the injected ports. Use it when you have no
  transport opinion.
- :func:`resolve_entry` plus :func:`partition_batches` / :func:`requested_transcript` /
  :func:`combine_recoded` — the pure kernel. Use it when you own your transport (an asyncio loop, your
  own rate limiter, your own progress reporting) and only need the resolution rule.
"""

from ._core import (
    combine_recoded,
    needs_refseq_transcripts,
    partition_batches,
    reference_identical_resolution,
    requested_transcript,
    resolve_consequences,
    resolve_entry,
)
from ._ports import ReferenceSequence, VariantRecoder, VepLookup
from .types import (
    NO_CHANGE_TERM,
    ConsequenceOutcome,
    ConsequenceResolution,
    ConsequenceSource,
    VepConfig,
    VepInput,
)

#: Version of the resolution *rule* — the logic in ``_core`` plus the ranking in
#: ``lib.sequence_ontology``. Bump it whenever a change would make this library return a different
#: consequence for an unchanged input and an unchanged Ensembl release: a new or reordered term in
#: ``ENSEMBL_CONSEQUENCE_RANKING``, a change to the transcript-matching rule, a change to how the
#: Recoder fallback combines forms.
#:
#: Consumers that cache or persist resolutions should key their staleness check on
#: ``(ensembl_release, RESOLVER_VERSION)``, not on the Ensembl release alone. The release captures
#: upstream change; this captures ours. Keying on the release alone means a fix shipped here is
#: invisible to every already-annotated record — they look current and are silently never recomputed.
#:
#: Not tied to the package version: a packaging or CLI change must not invalidate stored answers, and a
#: resolution-rule change must invalidate them even in a patch release.
RESOLVER_VERSION = "1"

__all__ = [
    "NO_CHANGE_TERM",
    "RESOLVER_VERSION",
    "ConsequenceOutcome",
    "ConsequenceResolution",
    "ConsequenceSource",
    "ReferenceSequence",
    "VariantRecoder",
    "VepConfig",
    "VepInput",
    "VepLookup",
    "combine_recoded",
    "needs_refseq_transcripts",
    "partition_batches",
    "reference_identical_resolution",
    "requested_transcript",
    "resolve_consequences",
    "resolve_entry",
]
