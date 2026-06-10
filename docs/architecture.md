# Library architecture

`variant_annotation/lib/` is the installable library core. Scripts in `src/` are
thin composition roots and CLI wiring. The separation is strict: library code never
imports from `src/`; scripts are the only place that wire everything together.

## Why this library exists

Variant annotation needs to produce consistent results whether it is being run
through the MaveDB API or being applied manually to data in the lab. Having the
logic live in one place — as an installable library — means both callers use
exactly the same implementation. The CLI tools in `src/` make the same library
callable from the command line without duplicating any logic.

The split between this library and `dcd_mapping` reflects that division of
responsibility: `dcd_mapping` produces variant mappings at the assay level;
`variant_annotation` handles the annotation and equivalence-class work that
follows, and is the shared layer both the API and manual pipelines call.

## Layers

```
variant_annotation/lib/          ← installed library
  translation/                   ← feature module (ports-and-adapters pattern)
    _ports.py                    ← narrow protocols
    _core.py                     ← business logic; imports ports, not clients
    types.py                     ← public input/output types
    __init__.py                  ← re-exports public API
  clients/                       ← protocol implementations
    uta.py                       ← UTA-backed TranscriptSource
    coordinates.py               ← HGVS-mapper-backed CoordinateTranslator
  pipeline/                      ← CSV composition boundary (CLI-only)
    reverse_translate_step.py    ← ColumnConfig + process_rows
  hgvs/                          ← shared vocabulary module
  accessions.py                  ← shared vocabulary module
  consequence.py                 ← shared vocabulary module

src/                             ← composition roots + CLI entry points
  reverse_translate_protein_variants.py   ← wires clients into process_rows
  add_vcf_identifiers.py         ← (not yet refactored to import downward)
  ...
```

The dependency arrows are strictly downward:

```
src/  →  pipeline/  →  translation/  →  _ports.py
                     →  accessions.py        ↑
                     →  consequence.py       │
                                       shared vocabulary
src/  →  clients/  →  accessions.py
                    →  hgvs/
```

`_core.py` imports from `_ports.py` and shared vocabulary only — never from
`clients/`. `clients/` never imports from `translation/`. `pipeline/` is the
composition boundary between the CSV layer and the library.

There are four layers in the system:

- **Feature modules** (`translation/`) — own a domain operation, have external
  service dependencies, follow the `_core/_ports/types/__init__` layout.
- **Shared vocabulary modules** (`hgvs/`, `accessions.py`, `consequence.py`) —
  pure functions and types with no service dependencies; importable from anywhere
  in the library without causing cycles.
- **Protocol implementations** (`clients/`) — concrete implementations of
  protocols defined by feature modules.
- **Composition roots** (`src/`) — wire clients into library calls, own all I/O
  and CLI concerns, hold no business logic.

## Feature module structure

Feature modules own a domain operation and depend on external services. They
follow the same four-file layout. `translation/` is the reference:

**`types.py`** — public input and output types as plain dataclasses. No imports
from `_ports.py` or `clients/`; no I/O. Both the core logic and callers import
from here, so keeping it dependency-free prevents cycles and makes the types
independently importable for tests or type annotations.

**`_ports.py`** — narrow `Protocol` definitions for any external services the
module needs. Private (underscore prefix) because callers never import from it
directly — implementations satisfy it structurally. Owned by the consumer module,
not by `clients/`, so the dependency arrow always points toward the consumer.

**`_core.py`** — business logic. Imports from `_ports.py` and `types.py` and
shared library vocabulary (`accessions`, `consequence`, `hgvs/`). Never imports
from `clients/` — that would collapse the layer boundary. Private because the
public surface is `__init__.py`; internal organisation can change without
breaking callers.

**`__init__.py`** — re-exports everything callers need and declares `__all__`.
This is the stable public surface. Callers always import from the package, not
from private submodules. When internal files are reorganised, only this file
needs updating to keep the public API intact.

Modules without external service dependencies (e.g. a pure transformation step)
can skip `_ports.py` and `_core.py` and expose their functions directly from
`__init__.py`. The four-file layout is for modules that need dependency
injection, not a universal requirement.

## Shared vocabulary modules

`hgvs/`, `accessions.py`, and `consequence.py` are shared vocabulary: pure
functions and types that multiple modules in the library need. They have no
external service dependencies and no ports. Any `lib/` module can import from them
without risk of creating a dependency cycle.

`hgvs/` graduated from a single utility function into a directory because HGVS
parsing logic was needed by `clients/`, `pipeline/`, and the CLI scripts — enough
distinct consumers that a shared home earned its place. `accessions.py` and
`consequence.py` are still single files because their scope hasn't grown past that.

**Recognising when a new shared vocabulary module is warranted:** when the same
parsing, classification, or type logic appears in more than one `lib/` module, or
when a `_core.py` is accumulating utility functions that aren't specific to its
domain operation. Extract the shared concept into its own module rather than
making one feature module depend on another. The name should describe the
vocabulary, not the feature that first needed it.

The same principle applies at different scales: a single `.py` file is the right
shape for a focused concept (`accessions.py`); a directory is warranted when the
concept has distinct sub-concerns worth separating (`hgvs/parse.py`,
`hgvs/fields.py`).

## Public API shape

Each feature module exposes a public API through its `__init__.py`, which
re-exports the types and functions callers need and lists them in `__all__`.
Private submodules use a leading underscore. Callers always import from the
package, not from private submodules — internal layout can change without breaking
callers.

Shared vocabulary modules may skip the `__init__.py` indirection when the module
is a single file and callers importing directly from it is fine. When the module
grows into a directory, add an `__init__.py` to restore a stable import surface.

External consumers (such as mavedb-api) import only from feature modules and
`clients/`. They never import from `pipeline/` or `src/`.

## Ports and adapters

`_ports.py` in each `lib/` module defines narrow `Protocol` interfaces **owned by
the consumer**. Implementations live in `clients/` or in the external consumer
(mavedb-api), and satisfy the protocols structurally without importing from
`_ports.py`.

**Why extract UTA and the HGVS mapper behind protocols rather than calling them
directly?** Three reasons:

1. **Testability.** UTA is a PostgreSQL service; the HGVS mapper depends on
   seqrepo. Both require real infrastructure. Extracting them behind protocols
   means the core logic can be tested with stub implementations that return
   controlled data — no database required. The real clients only touch tests that
   specifically exercise the service integration.

2. **The API may eventually manage its own UTA connection.** Currently mavedb-api
   uses the library's `clients/uta.py` directly. If and when the API moves to its
   own UTA connection management (e.g. to share a connection across multiple jobs
   or align with its existing `AlleleTranslator` lifecycle), the protocol makes
   that swap a drop-in replacement with no changes to the core logic.

3. **The UTA connection-slot constraint is real.** Opening a new UTA connection
   per variant (or per call) exhausts the server's reserved connection budget in
   production. Protocols make it explicit that the caller controls connection
   lifetime — the library never opens its own UTA connection; that decision
   belongs to the composition root.

Consumer ownership prevents circular imports: `translation/` defines the protocols
it needs; `clients/` implements them; neither knows about the other. The
composition root (`src/` or mavedb-api) is the only place that sees both sides.

## Composition roots (`src/`)

`src/` is the fourth layer: the only place in the system that sees all the other
layers at once. A composition root instantiates concrete clients, injects them
into library calls, and handles everything I/O-facing — file handles, env vars,
argparse or click wiring, CSV streaming, progress logging. It owns no business
logic. If logic in `src/` looks reusable or testable in isolation, it belongs in
`lib/` instead.

Composition roots are the glue that makes the library usable as a CLI tool. They
are not a migration destination — a script does not graduate into `lib/` over
time. The layer boundary is permanent: `lib/` is the library, `src/` is how a
human operator invokes it.

**What a composition root does:**
- Reads env vars and validates preconditions (`UTA_DB_URL`, etc.)
- Instantiates clients (`UtaClient`, `HgvsMapper`) with the right connection
  lifetime for the job
- Opens input/output files and manages CSV streaming
- Calls into `lib/` (directly, or via a `pipeline/` adapter for complex CSV work)
- Logs progress and handles top-level errors

**What it does not do:**
- Touch row dicts with column name strings — that is `ColumnConfig`'s job
- Contain conditional logic about variant types or annotation rules
- Import from another `src/` script

`src/reverse_translate_protein_variants.py` is the reference. It instantiates
`UtaClient` and `HgvsMapper`, opens the input and output files, and delegates all
logic to `process_rows`. The only decisions it makes are operational: batch size,
progress interval, when to flush.

**`pipeline/` as a CSV sub-layer.** When a script's row-level CSV logic is complex
enough to test independently — column name mappings, row classification, merging
results back — it can be extracted into a `pipeline/` module. `pipeline/` sits
between `src/` and `lib/`: it imports from the library and is imported by `src/`,
but it is still CLI-only and mavedb-api never touches it. Not every script needs
one; a script that reads one column and writes one column does not.

| Concern | Lives in |
|---|---|
| Argparse / click wiring | `src/` |
| Env var loading, file handles, streaming | `src/` |
| Client instantiation | `src/` |
| Complex CSV row schema and batching | `pipeline/` |
| Business logic | `variant_annotation/lib/` |

## ColumnConfig — the CSV schema adapter pattern

When a `pipeline/` module handles CSV rows, all column name knowledge belongs in
a single dataclass: the schema adapter. `ColumnConfig` in
`pipeline/reverse_translate_step.py` is the reference implementation.

**Why centralise column names.** CSV column names are user-facing — they default
to the pipeline's own conventions but are overridable via CLI flags. Without a
central adapter, column name strings scatter across read calls, write calls, and
field derivation logic. Renaming a column or making it configurable then requires
hunting down every reference. The adapter makes each column name a single
configurable field, and every access to that column goes through a method.

**Structure.** A schema adapter is a dataclass whose fields are column name
strings. It exposes read methods that extract typed values from row dicts, and
write methods that merge results back. No code outside the adapter touches a row
dict with a raw column name string.

```python
@dataclass
class ColumnConfig:
    # column name fields — one per logical column, defaulting to pipeline conventions
    input_hgvs: str = "mapped_hgvs_p"
    result_column: str = "result"
    error_column: str = "error"

    # read methods — extract typed values from a row
    def read_hgvs(self, row: dict) -> str | None:
        return (row.get(self.input_hgvs) or "").strip() or None

    # write methods — merge results back into a row
    def write_result(self, row: dict, value: str) -> None:
        row[self.result_column] = value

    def write_error(self, row: dict, message: str) -> None:
        row[self.error_column] = message
```

The column name fields are passed directly to the constructor from parsed CLI
args, which makes every column name a first-class CLI option with no additional
wiring.

**The pipe-join boundary.** When candidates are stored as pipe-delimited strings
in a column (multiple values per cell), the join and split happen only inside the
adapter — `write_result` joins, a `_read_candidates` method splits. Everything
between those two points is `list[str]`. Splitting or joining anywhere else —
in `_core.py`, in the CLI, anywhere outside the adapter — is a violation.

**When to use this pattern.** Any `pipeline/` module that reads more than a couple
of named columns benefits from it. The signal that you need it: you find yourself
passing column name strings as function arguments, or writing `row["some_column"]`
in more than one place. Extract those into a `ColumnConfig` dataclass, give each
column a field with the right default, and route all row access through its
methods.

## `hgvs/parse.py` transitional state

`variant_annotation/lib/hgvs/parse.py` currently re-exports `parse_hgvs` and
`apply_vcf_anchor` from `src/add_vcf_identifiers.py`. The functions live in a
script that has not yet been refactored to import downward from the library.

This means `hgvs/parse.py`, `hgvs/fields.py`, and `pipeline/` all require `src/`
to be on the Python path — satisfied in tests via `pythonpath = [".", "src"]` and
at runtime via Docker. The mavedb-api-facing surface (`translation/` and
`clients/`) does not import from `hgvs/parse.py` and is unaffected.

When `add_vcf_identifiers.py` is refactored to import downward, the function
bodies move into `hgvs/parse.py` and become canonical. New HGVS parsing logic
belongs in `hgvs/parse.py` now — do not add it to `add_vcf_identifiers.py`.

## Test structure

Tests mirror the layers. `tests/lib/` mirrors `variant_annotation/lib/`
directory-for-directory (`tests/lib/translation/`, `tests/lib/pipeline/`, …);
`tests/scripts/` holds the tests for the composition roots in `src/`. Library
unit tests stub services through the protocol fixtures in `tests/lib/conftest.py`
(`stub_transcripts`, `stub_coordinates`) rather than patching import paths —
the same dependency-injection seam the ports exist for.

Three markers, declared in `[tool.pytest.ini_options]` in `pyproject.toml`,
classify every test:

- **`unit`** — exercises a single lib function or class in isolation with all
  dependencies stubbed. No real services, no subprocess, no CLI, no file pipeline.
- **`integration`** — exercises multiple components together: a CLI script run
  end-to-end (`main()`, `CliRunner`, or the script's entry function with real
  file I/O), or a `pipeline/` step wired to real collaborators.
- **`slow`** — additive, never standalone. A test is `slow` *and* `unit`/
  `integration` when it pays real `time.sleep`, subprocess startup, or
  large-fixture cost. Lets `-m "not slow"` trim the CI inner loop.

The markers are the contract for *which* suite runs when: `pytest -m unit` for the
fast inner loop, `pytest -m integration` before merging, `pytest -m "not slow"` to
skip the expensive concurrency tests. That contract only holds if every test is
labelled by **what it actually exercises** — not by which directory or file it
lives in.

**Markers are applied per test, not per file.** A module-level
`pytestmark = pytest.mark.<kind>` is correct *only* when every test in the file is
the same category. That is true throughout `tests/lib/` — those files are
uniformly `unit` — so they carry a single module-level marker.

It is usually **not** true in `tests/scripts/`. A single CLI script's test file
typically has a mixed mandate: pure-helper unit tests (a parser, a formatter, a
row classifier called directly with literal inputs and stubbed collaborators) sit
alongside end-to-end integration tests (`mod.main()` / `CliRunner` over real temp
files). In a mixed file, a blanket module-level marker mislabels half the suite —
so each test function (or test class) carries its own `@pytest.mark.unit` /
`@pytest.mark.integration` decorator, and there is no module-level marker. Class
groupings get one decorator on the class; standalone functions get one each.

Two invariants make the partition checkable:

- **Exhaustive** — every collected test carries `unit` or `integration`.
  `pytest -m "not unit and not integration" --collect-only` must collect nothing.
- **Exclusive** — no test carries both. `pytest -m "unit and integration"
  --collect-only` must collect nothing.

When a script is migrated into the library, its tests move and are re-marked to
match their new home — see the test-migration section in
[migration-guide.md](migration-guide.md).

## Invariants not enforced by tooling

- `_core.py` must never import from `clients/`. No import check enforces this.
- `ColumnConfig` is the only place that reads or writes named CSV columns.
- `pipeline/` is a CLI-only module; mavedb-api must not import from it.
- UTA-backed providers must be constructed once per job and reused across variants
  (see the UTA connection-slot exhaustion incident documented in lore).
- Every test carries exactly one of `unit` / `integration`; mixed `tests/scripts/`
  files mark per test, never with a blanket module-level marker.
