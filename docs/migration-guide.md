# Migrating a script to the library format

This guide covers moving scripts in `src/` that mix business logic, I/O, and CLI
wiring toward the four-layer structure described in [architecture.md](architecture.md).
Not every script needs a full migration — read the decision section before starting.

## When migration is worth doing

- Another caller (mavedb-api, another script, tests) needs to reuse the logic
- The script is hard to test without spinning up real services
- The same parsing or classification logic has appeared in more than one script
- The core logic is deeply tangled with CSV or file concerns

A script that is purely mechanical (filter columns, reorder rows, join TSVs) has
no business logic worth isolating and can stay as-is.

## Deciding which kind of extraction

Before writing any new files, identify what you actually have:

**The logic needs to be called from mavedb-api or reused across multiple callers**
→ Extract a feature module. Full path below.

**The same utility function (parsing, classification, type) appears in more than
one module**
→ Extract shared vocabulary. Lighter path below.

**The logic is self-contained to this script but is tangled with I/O**
→ No `lib/` changes needed. Slim the script into a proper composition root: pull
the logic into a plain function, push all I/O and env concerns to `main()`. The
split question is: *does this know how data is stored or transported?* If yes, it
stays in `src/`. If no, it can be a plain function within the same file.

---

## Path 1: Extracting a feature module

Use this when the logic needs to be callable from mavedb-api or from outside the
CLI. `variant_annotation/lib/translation/` is the reference for every step.

### 1. Sort the script into buckets

- **Pure logic** — transforms data without I/O or service calls
- **Service calls** — talks to UTA, HGVS mapper, ClinGen, etc.
- **I/O + wiring** — file handles, CSV streaming, argparse/click, env vars

Pure logic goes in `_core.py`. Service calls get abstracted behind a protocol in
`_ports.py` and implemented in `clients/`. I/O stays in `src/`.

### 2. Define input/output types in `types.py`

```python
# variant_annotation/lib/my_step/types.py
from dataclasses import dataclass

@dataclass
class MyInput:
    hgvs: str

@dataclass
class MyResult:
    input: MyInput
    ...

@dataclass
class MyError:
    input: MyInput
    error: str
```

No service imports, no I/O. Both `_core.py` and callers import from here —
keeping it dependency-free prevents cycles.

### 3. Define service protocols in `_ports.py`

Only needed if the core logic calls an external service. Define the narrowest
interface possible — only the methods `_core.py` actually calls.

```python
# variant_annotation/lib/my_step/_ports.py
from typing import Protocol

class MyServiceSource(Protocol):
    def lookup_thing(self, key: str) -> str | None: ...
```

The protocol is owned by the consumer (`my_step/`), not by `clients/`. The
implementation satisfies it structurally without importing from `_ports.py`.

### 4. Write the core logic in `_core.py`

```python
# variant_annotation/lib/my_step/_core.py
from ._ports import MyServiceSource
from .types import MyInput, MyResult, MyError

def process(
    inputs: list[MyInput],
    *,
    service: MyServiceSource,
) -> tuple[list[MyResult], list[MyError]]:
    ...
```

Imports: `_ports.py`, `types.py`, shared vocabulary (`accessions`, `hgvs/`, etc.).
Never imports from `clients/` — that collapses the layer boundary.

### 5. Add a client implementation in `clients/`

```python
# variant_annotation/lib/clients/my_service.py
class MyServiceClient:
    def lookup_thing(self, key: str) -> str | None:
        # call the real service
        ...
```

Satisfies `MyServiceSource` structurally. Does not import from `_ports.py`.

### 6. Expose the public API in `__init__.py`

```python
# variant_annotation/lib/my_step/__init__.py
from ._core import process
from ._ports import MyServiceSource
from .types import MyInput, MyResult, MyError

__all__ = ["process", "MyServiceSource", "MyInput", "MyResult", "MyError"]
```

Callers import from `variant_annotation.lib.my_step`, never from private
submodules.

### 7. Slim the script into a composition root

```python
# src/my_script.py
from variant_annotation.lib.clients.my_service import MyServiceClient
from variant_annotation.lib.my_step import MyInput, process

def main():
    # parse args / load env
    service = MyServiceClient(...)
    inputs = [MyInput(...) for row in read_csv(...)]
    results, errors = process(inputs, service=service)
    write_csv(results)
```

The script instantiates clients, handles I/O, and delegates everything else.

### 8. Test each layer independently

Library tests stub the protocol directly — no patching of import paths:

```python
class _FakeService:
    def lookup_thing(self, key: str) -> str | None:
        return "stub"

def test_process():
    results, errors = process([MyInput(hgvs="NP_000001.1:p.Met1Val")], service=_FakeService())
    assert len(results) == 1
```

Composition root tests use `monkeypatch` to swap client constructors, as in
`test_reverse_translate_protein_variants.py`.

---

## Path 2: Extracting shared vocabulary

Use this when the same function or type is needed in more than one `lib/` module
and has no external service dependencies.

The extraction is simpler than a feature module — no ports, no injection, no
`__init__.py` required for a single-file module.

1. Create a new file in `variant_annotation/lib/` named after the concept, not
   the feature that first needed it (e.g. `accessions.py`, not
   `translation_utils.py`).
2. Move the shared functions/types there.
3. Update all callers to import from the new module.
4. If the concept grows to have distinct sub-concerns, promote it to a directory
   with its own `__init__.py` (as `hgvs/` did).

The only rule: shared vocabulary modules must have no service dependencies and no
imports from feature modules or `clients/`. They sit at the base of the dependency
graph and everything else imports from them.

---

## Migrating a script's tests

Tests follow the code. When logic moves out of `src/` and into `lib/`, the tests
that cover that logic move with it — and get re-marked to match their new home.
The marker taxonomy and the per-test discipline are described in the test-structure
section of [architecture.md](architecture.md); this is how a migration touches them.

**1. Move the unit tests into `tests/lib/`.** The tests that called the script's
helper functions directly — parsers, classifiers, formatters — now belong beside
the module they test, under `tests/lib/<module>/`, mirroring the new library
layout. They stay `unit`.

**2. Switch stubbing from monkeypatch to protocol stubs.** A composition-root test
swaps client *constructors* with `monkeypatch`. A library unit test stubs the
*protocol* directly — pass a fake satisfying `MyServiceSource` into the core call,
or reuse a fixture in `tests/lib/conftest.py`. The migration is the moment to drop
the import-path patching; the ports exist precisely so you don't need it.

**3. Leave the integration tests in `tests/scripts/`.** Once the script is a slim
composition root, its remaining tests exercise the end-to-end wiring — `main()`
over real temp files, CLI argument handling, streaming. Those stay in
`tests/scripts/` and are marked `integration`.

**4. Re-mark both files precisely.** A migration usually changes a file's mandate:

- The `tests/scripts/` file may go from **mixed** to **uniform** once its unit
  tests leave — collapse the per-test decorators into a single module-level
  `pytestmark = pytest.mark.integration`.
- The new `tests/lib/` file is uniformly `unit` — give it one module-level
  `pytestmark = pytest.mark.unit`, not per-test decorators.
- If the script file is still mixed afterward, keep decorating per test.

**5. Verify the partition still holds.** After moving tests, confirm nothing slipped
through the cracks:

```
pytest -m "not unit and not integration" --collect-only   # must collect nothing
pytest -m "unit and integration" --collect-only            # must collect nothing
```

Re-marking is part of the migration changeset, not a follow-up. A moved test with a
stale marker is the same defect as a moved function with a stale import.

---

## What not to do

**Don't create protocols for stdlib or subprocess.** File I/O, `subprocess`,
`csv`, `pathlib` — these are not services and don't need protocols.

**Don't put column-name logic in `_core.py`.** If the extraction involves CSV
rows, column-name mapping belongs in a `ColumnConfig`-style adapter in
`pipeline/`, not in the core module. The core function takes typed inputs and
returns typed outputs; the adapter translates between those types and raw row
dicts. See the ColumnConfig section in [architecture.md](architecture.md) for the
full pattern.

**Don't name a shared vocabulary module after the feature that first needed it.**
`hgvs.py` is right; `translation_hgvs_utils.py` is wrong. The name should
describe the concept so other modules can discover and reuse it.

---

## Reference

- [architecture.md](architecture.md) — the four layers and their invariants
- `variant_annotation/lib/translation/` — reference feature module
- `variant_annotation/lib/hgvs/` — reference shared vocabulary module
- `variant_annotation/lib/pipeline/reverse_translate_step.py` — reference
  `pipeline/` CSV adapter
- `src/reverse_translate_protein_variants.py` — reference composition root
