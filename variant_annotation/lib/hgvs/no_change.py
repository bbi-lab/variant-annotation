"""Recognising HGVS expressions that describe no sequence change. Pure — no I/O, no service calls.

Reference-identical rows are ordinary in functional-assay data (a wild-type control), but they are not
sequence variants, so Ensembl VEP rejects them. Telling them apart from genuine variants is what lets a
consumer label them explicitly instead of leaving a null that cannot be distinguished from "annotation
has not run yet".

Deliberately free of any dependency on ``hgvs/parse.py``, which still requires ``src/`` on the path (see
the transitional-state section in ``docs/architecture.md``). The mavedb-api-facing modules import from
here, so it must stay importable from the installed library alone.
"""

from __future__ import annotations

import re
from typing import Optional

#: Explicit HGVS "no change" notation: a trailing ``=`` asserting the sequence is unchanged
#: (``NM_000049.4:c.123=``, ``NC_000003.12:g.10141916C=``, ``NP_000537.3:p.Met4=``).
_NO_CHANGE_RE = re.compile(r"^[^:]+:[cgnpmr]\.[^:]*=\s*$")

#: Coding ``delins`` over an explicit interval, e.g. ``NM_000049.4:c.256_257delinsAA``. Whether such an
#: expression is reference-identical cannot be decided from the notation alone — it depends on what the
#: transcript actually reads over that interval — so this only extracts the parts a caller needs to ask.
_CODING_DELINS_RE = re.compile(
    r"^(?P<transcript>[^:]+):c\.(?P<start>\d+)_(?P<stop>\d+)delins(?P<alt>[ACGTNacgtn]+)\s*$"
)


def states_no_change(hgvs: str) -> bool:
    """Whether ``hgvs`` explicitly asserts no sequence change via HGVS ``=`` notation.

    The notation is itself the claim, so no reference sequence is needed to decide it. Covers all
    sequence types, including protein ``p.Met4=`` forms — which are the cases that cannot be recoded to
    an alternate synonymous codon at all when the residue has no codon redundancy (Met, Trp).
    """
    return bool(_NO_CHANGE_RE.match((hgvs or "").strip()))


def parse_coding_delins(hgvs: str) -> Optional[tuple[str, int, int, str]]:
    """Parse a coding interval ``delins`` into ``(transcript, start, stop, alt)``.

    Returns ``None`` when ``hgvs`` is not of that form. The interval is 1-based and inclusive, as
    written in the expression.

    A ``delins`` is the form the mapping pipeline emits for an unchanged interval — the replacement
    bases can equal the reference — so a caller that can read the transcript's reference bases uses this
    to decide whether the expression is reference-identical.
    """
    match = _CODING_DELINS_RE.match((hgvs or "").strip())
    if match is None:
        return None
    return (
        match.group("transcript"),
        int(match.group("start")),
        int(match.group("stop")),
        match.group("alt"),
    )
