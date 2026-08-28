# Review 2 of `plan_disulfide_cysteines.md`

**Lens:** software correctness, API surface, regression risk.
**Reviewed at:** `6-disulfide-cysteines`, working tree at `b5e8409` plus untracked docs. Date 2026-08-28.
**Method:** every claim below was checked against the code, most of them by building a scratch copy of the repository at `/tmp/.../scratchpad/proto`, implementing the plan as literally as the text allows (new `calculations/disulfides.py`, parser hook, `disulfide_partner`, `side_chain_pka`, key-based `Atom.charge`, `Structure.disulfides`, `run_metadata` key) and running the suite against it. Nothing in `/mnt/shared/prodes` was modified.

Baseline: `python -m pytest -q` -> **257 passed in 100.88s** **[verified 2026-08-28]**.
Prototype: `PYTHONPATH=proto/src:proto python -m pytest -q` -> **4 failed, 253 passed in 141.90s** **[verified 2026-08-28]**.

---

## Verdict

The plan is unusually well grounded: I re-checked every number in it and did not find one that is wrong. The physics, the choke point (`Atom.charge`), the decision not to add a feature column, the MAJOR bump and the propagation story are all correct, and a faithful implementation reproduces the headline results exactly (1GDW formal charge at pH 8.5 goes from -1.0 to +7.0, pI 8.899 to 10.342, ARH96693 output bit-identical to `main`). What it gets wrong is smaller and entirely in the maintenance corner: the terminus half of the `Atom.charge` rewrite converts a today-wrong-number into a crash on realistic pKa-file input; the test-impact section asserts that `tests/test_residue.py` is unaffected when two of its assertions break; the "the `_pka` cache cannot be stale" invariant holds only for the one call site the plan happens to use, not for the public function it is exposing; and two of the Phase-4 test cases cannot be written the way they are described because `Residue.pkas` dereferences `self.chain`.

**Do not merge as written.** Fix findings 1-5 (all small, all localised) and it is good to implement.

---

## Findings

### 1. BLOCKING - reading the terminus pKa by key turns a wrong number into a crash

**Plan says** (Phase 1, `Atom.charge`): "terminus: look up the `"N+"` or `"C-"` entry by key rather than taking `pkas[-1]`." It specifies a `None` fallback for the side chain ("if that is `None`, the branch is skipped") but says nothing about the terminus.

**Actually true.** `Structure.redo_pkas` (`src/prodes/core/structure.py:203-206`) overwrites `residue._pka` wholesale with whatever list the pKa file supplies. If that list has no `N+`/`C-` key and the residue is terminal, a key lookup returns `None`, and the very next line is `pka > ph`.

**Evidence** (prototype, 1GDW + `tests/data/1GDW.pka` with the `N+` entry for residue 1 removed) **[verified 2026-08-28]**:

```
prototype: res1 pkas [{'LYS': 11.34}] terminus N -> TypeError: '>' not supported between instances of 'NoneType' and 'int'
main:      res1 pkas [{'LYS': 11.34}] terminus N -> charge(7) = 2.0
```

This is not a contrived input. `redo_pkas` keys on residue *number* alone (`structure.py:204`), a defect the plan itself lists under Risks. In a two-chain structure, chain B residue 1 is terminal and inherits chain A residue 1's entry list; if chain A residue 1 is a non-terminal titratable residue, that list has no terminus key and the whole run dies. The workflow the issue names is antibodies, which are always multi-chain. `tests/data/1GDW_hpp_pka.json` already shows the shape in the wild: residue 5 carries `[{"N+": 7.541}]` and nothing else **[verified 2026-08-28]**.

**Fix.** In the terminus branch, fall back to `data.n_term_pka()` / `data.c_term_pka()` when the keyed entry is absent, and log once. Never index a `None`. Say so in the plan, because "look it up by key" reads as complete and is not.

---

### 2. SERIOUS - the test-impact section is wrong about `tests/test_residue.py`

**Plan says** (Phase 4): "`tests/test_data.py`, `tests/test_residue.py` — check whether the `CYS` entry assertions still hold; `general_data.json` is not being edited, so they should."

**Actually true.** The reasoning is the wrong one: those assertions do not read `general_data.json`, they read a parsed residue. `tests/test_residue.py:6` parses `1GDW_h.pdb.zip` and `residues["CYS"]` resolves to **CYS 6**, which is bonded to CYS 128 **[verified 2026-08-28]**. Two assertions break:

* `tests/test_residue.py:141` `residues["CYS"].pkas == [{"CYS": 8.33}]` -> `None`
* `tests/test_residue.py:197` `round(residues["CYS"].charge(14), 3) == -1` -> `0`

`tests/test_data.py` is genuinely unaffected: it asserts on the JSON, not on a residue.

**Fix.** List both `tests/test_residue.py` assertions explicitly in the plan's edit list, with the same "these are the numbers the issue exists to move" comment the plan already promises for `test_structure.py`. Better: point `residues["CYS"]` at a free thiol so the test keeps testing the default pKa, and add the bonded case as a separate assertion.

---

### 3. SERIOUS - `assign_disulfides` does not invalidate `_pka`, so the invariant holds only by accident

**Plan says** (Phase 1, `Residue`): "The existing `_pka` cache makes ordering load-bearing: `pkas` is computed once, on first access, and cached. Detection runs inside the parser before any caller can touch `pkas`, so the cache cannot be stale."

**Actually true.** That is a property of the *one* call site, not of the function. `assign_disulfides(structure, ssbond_pairs=())` is described as "the entry point" of a public module; nothing stops a caller running it on a structure whose `pkas` have already been read, and the plan's own Phase-2 acceptance and Phase-3 count reporting invite exactly that.

**Evidence** (prototype) **[verified 2026-08-28]**: read `.pkas` on a cysteine first (caching `[{'CYS': 8.33}]`), then call `assign_disulfides(structure)` -> `disulfide_partner` is set, `pkas` is still `[{'CYS': 8.33}]`, and the SG stays titratable.

Two further inaccuracies in the same paragraph:

* "`pkas` is computed once ... and cached" is false for the case this change creates. `residue.py:42` caches only when the list is non-empty, so a bonded mid-chain cysteine recomputes on every access and returns `None` every time. Harmless, but the docstring the plan wants to write would be wrong.
* "Phase 2's change to `redo_pkas` is the only other writer of `_pka`" - `_pka` is also a constructor argument (`residue.py:6`), so a caller can pass one in.

**Fix.** `assign_disulfides` sets `residue._pka = None` on both members of each pair it marks. One line, and it makes the ordering claim true rather than circumstantial.

---

### 4. SERIOUS - two Phase-4 test cases cannot be written the way they are described

**Plan says** (Phase 1 acceptance, Phase 4): "A hand-built structure with one disulfide and one free cysteine"; "A contested sulfur: three `SG` atoms where two candidate pairs share one atom".

**Actually true.** `Residue.pkas` does `self == self.chain.residues[0]` (`residue.py:38`). A `Residue` built without a chain raises on the first `.pkas` access.

**Evidence** (prototype) **[verified 2026-08-28]**: `Residue("CYS", Structure("hand"), 1).pkas` -> `AttributeError: 'NoneType' object has no attribute 'residues'`.

The geometry-only cases (`geometric_disulfides`, contested sulfur, 2.4 A pair) are fine on bare `Atom` objects, because detection never touches `pkas`. Anything asserting on charge needs a real chain.

**Fix.** Write the mixed-cysteine and terminal-cysteine cases as small PDB text fixtures parsed by `PDBparser`, which is what the plan already proposes for the two `SSBOND` fixtures. Drop "hand-built structure" from the Phase-1 acceptance list, or state that "hand-built" means "hand-written PDB text".

---

### 5. SERIOUS - `read_ssbond_line` robustness is unspecified, and the plan's own fixtures will trip it

**Plan says**: "lines whose record name is `SSBOND` are parsed by column into `(chain1, number1, chain2, number2, sym1, sym2)`", and "a record carrying a symmetry operator other than `1555` in either symmetry field ... is skipped".

**Actually true.** The PDB columns are chainID1 at col 16 (`line[15]`), seqNum1 17-21, chainID2 col 30, seqNum2 32-35, sym1 60-65, sym2 67-72. A hand-written fixture line that stops after the second residue is 35 characters long; `line[59:65]` is then `""`, which is "other than 1555", so under the rule as written **every bond in that fixture is silently skipped**. A line shorter than 16 characters raises `IndexError` rather than anything a user can act on. My prototype only avoided this because I added `or "1555"` defensively **[verified 2026-08-28]** - the plan does not ask for it.

**Fix.** Specify: an empty symmetry field means `1555`; a line too short to hold the residue fields raises a `ValueError` naming the file and line. Then the unit test the plan wants ("`read_ssbond_line` on a single real `SSBOND` line, including one with a symmetry operator") should gain a truncated line as a third case.

---

### 6. SERIOUS - records-first suppresses geometry wholesale, including for bonds the records omit

**Plan says**: "`assign_disulfides(structure, ssbond_pairs=())` ... Uses the `SSBOND` records when the file has any, otherwise geometry", and Phase 4 pins the behaviour: "a geometric pair the record set does not name ... is *not* called because records win when present."

**Actually true** as a design choice, but it is all-or-nothing on the *file*, not per residue. One surviving `SSBOND` line in a header that has been trimmed (a very common state for files that have been through a modelling pipeline) suppresses detection of every other bridge in the structure. That fails in the direction the issue is about: spurious negative charge.

**Fix.** Either let records win only for the residues they name and let geometry fill the rest - which is strictly safer and costs nothing, since the two agree wherever both fire - or keep the plan's rule and add a warning when geometry finds pairs the records do not name, so the suppression is visible. The plan already logs the inverse disagreement; this one is the more dangerous direction and is currently silent.

---

### 7. MINOR - `run_metadata` should take the count as a keyword argument

**Plan says**: "`run_metadata` gains `"disulfides": n`".

**Actually true and safe.** Adding a top-level metadata key does **not** break `tests/test_output_bundle.py`: `test_the_run_record_says_what_was_run` compares only `record["settings"]` exactly (`test_output_bundle.py:103-112`), and `tests/test_screening.py:337` is a substring check. Verified: prototype with `"disulfides"` added to `run_metadata` -> `pytest tests/test_output_bundle.py` **32 passed** **[verified 2026-08-28]**.

The only risk is the signature. `run_metadata` is a public function of a public module, and the repo already has an out-of-tree consumer (`tests/test_pka.py:14-16` records that biochai imports `convert_propka`). `run.py:553` is the only in-repo caller **[verified 2026-08-28: grep over `src/`, `tests/`, `scripts/`]**.

**Fix.** `def run_metadata(pdb_file, settings, n_points, ep, disulfides=0)`. Say so in the plan.

---

### 8. MINOR - the Phase-3 print cannot go where the plan puts it

**Plan says**: "`prepare_structure` prints the disulfide count alongside the existing `calculating {name}` line."

**Actually true.** `print(f"calculating {structure.name}")` is in `calculate()` at `run.py:529`, after `construct_surface_grid`; `prepare_structure` is `run.py:453-461`. A print inside `prepare_structure` lands *before* the "calculating" line, with a whole SASA pass between them, not alongside it. No test asserts on stdout (`grep -n "capsys\|capfd\|stdout" tests/*.py` -> only `test_parser.py:111`) **[verified 2026-08-28]**, so this is cosmetic - but pick one place and say which.

---

### 9. MINOR - `write_pdb` drops `SSBOND`, and nothing in the plan or the tests notices

**Actually true.** `write_pdb` (`parser.py:232-360`) emits `ATOM`/`TER`/`END` only. A structure written out and read back loses its authoritative records and falls back to geometry. On 1GDW that gives the same four pairs, so the loss is invisible there; on a structure where records and geometry disagree it is not. Note that the *bundle* is unaffected: `output.bundle_files` copies the input file's text via `read_pdb_text` (`output.py:174`), so `SSBOND` records travel with the bundle intact.

`write_pdb` has no in-repo caller outside `tests/test_parser.py` **[verified 2026-08-28]**.

**Fix.** One sentence in the `write_pdb` docstring, or write the records back out. Not worth blocking on, but it belongs in the plan's Risks so the next reader does not rediscover it.

---

### 10. MINOR - the version argument is right; the version *stamp* will lag

The MAJOR bump is correct by the repo's own rule: the reported quantity changes for any disulfide-containing structure above pH ~8. But `prodes_version()` (`output.py:113-119`) reads `importlib.metadata.version("prodes")`, which comes from installed distribution metadata rather than from `pyproject.toml` at runtime. In this environment the editable install reports **5.0.0** while `prodes.__file__` points at the working tree **[verified 2026-08-28: `pip show prodes` -> Version: 5.0.0]**. So every bundle produced in a dev environment that has not been reinstalled after the bump carries the old version. `tests/test_dependency_versions.py` checks dependency pins and the Python version, not the project version, so nothing catches it. Worth one line in Phase 3: bump *and* reinstall.

---

### 11. MINOR - insertion codes are silently dropped, on both sides

`_read_line` reads `residue_number = int(line[22:26])` (`parser.py:127`) and ignores the iCode column, and `_add_atom` merges residues by number, so `100` and `100A` collapse into one residue. `SSBOND` carries icode1/icode2 at columns 22 and 36, which `record_disulfides` keyed on `(chain, number)` cannot represent. Pre-existing, not caused by this change - but the issue names antibodies, and Kabat/Chothia numbering is built on insertion codes. The plan already carries a "one quirk worth knowing" note about `altLoc`; this deserves the same treatment.

---

### 12. MINOR - there is no multi-chain fixture in the repository at all

All five shipped structures are single chain: `1GDW`, `1GDW_h`, `1GPB`, `ARH96693`, `ARH98503` all report exactly one chain `A` **[verified 2026-08-28]**. So inter-chain disulfides - the entire antibody case the issue is about, and the one place `record_disulfides`' `(chain, number)` keying differs from `redo_pkas`' number-only keying - would ship with no coverage.

**Fix.** Make one of the two hand-written fixtures two-chain with an inter-chain bridge. It costs nothing extra and is the only way the plan's claim that "the disulfide code ... resolves residues by chain and number and does not inherit the defect" gets tested.

---

### 13. MINOR - module placement is acceptable but introduces a new dependency edge

`src/prodes/io/parser.py` today imports only `prodes.core.*`; `prodes.calculations.*` is reached from `core` through function-local imports to dodge cycles (`atom.py:70`, `structure.py:82`). A top-level `from prodes.calculations.disulfides import ...` in `parser.py` is a new `io -> calculations` edge. It does **not** create a cycle - the prototype imports cleanly **[verified 2026-08-28]** - and `calculations/` is a defensible home for something that is fundamentally a distance calculation. Two conditions: `disulfides.py` must never import `prodes.io`, and it should carry a "why this exists" module docstring in the style of `output.py` and `feature_dictionary.py`, which the plan's sketch does not have (it has one line of prose in the plan and no docstring in the design). The plan is compliant with the house rules on `__init__.py` (untouched) and leading underscores (all names public).

---

### 14. MINOR - hard-coding PROPKA's 99.99 sentinel

Phase 2 suppresses the debug line "unless the value is PROPKA's 99.99 sentinel". Comparing a parsed float against a literal is fragile and predictor-specific. Suppress instead for any value that cannot titrate over the reachable pH range (say, outside 0-14), which covers `convert_hpp`'s `14.0`/`0` clamps and `convert_pypka`'s "Not" handling as well.

---

## Tests that break

Confirmed by running the full suite against a faithful prototype **[verified 2026-08-28: 4 failed, 253 passed]**.

| Test | Why it breaks | In the plan's list? |
|---|---|---|
| `tests/test_residue.py::test_pkas` (line 141) | `residues["CYS"]` is CYS 6 of `1GDW_h`, bonded to CYS 128. `.pkas` becomes `None`; the assertion expects `[{"CYS": 8.33}]`. | **No** - the plan says it should hold |
| `tests/test_residue.py::test_charge` (line 197) | `round(residues["CYS"].charge(14), 3)` goes `-1` -> `0`, same residue. | **No** |
| `tests/test_structure.py::test_isoelectric_point` (line 39) | `8.899` -> `10.342` on `1GDW_h`. | Yes |
| `tests/test_structure.py::test_charge` (line 62) | pH 11 assertion `-13` -> `-5.0`. pH 7 assertion is unaffected. | Yes |

Nothing else moves. In particular these all pass untouched, which is the regression guard working: `test_sasa.py` and `test_feature_dictionary.py` (ARH96693 has no disulfides), `test_pka.py` (all eight 1GDW cysteines come back as PROPKA's 99.99, so file-driven and structure-driven agree), `test_output_bundle.py`, `test_screening.py` (its 1GDW run is at pH 7, where a free thiol is formally neutral anyway), `test_parallel.py`, `test_parser.py`, `test_atom.py`, `test_data.py`, `test_grid_wizard.py`, `test_geometry.py`, `test_read.py`, `test_cli.py`, `test_run_prodes.py`.

`tests/test_structure.py::test_dipole` survives because it is evaluated at pH 7 with formal charges, where a cysteine contributes zero either way.

---

## What the plan got right

Everything below I re-derived independently and it holds.

* **`Atom.charge` really does have to change, and the plan's diagnosis of why is exact.** `atom.py:81/87/94/97` index `self.residue.pkas[0]`; with the side-chain entry gone that is `None[0]` for a mid-chain cysteine and the C-terminal 2.34 for a terminal one.
* **The latent positional bug is real.** With `_pka = [{"N+": 8.35}, {"LYS": 11.34}]` on 1GDW residue 1, `NZ.charge(10)` is `0` and the backbone `N.charge(10)` is `1`; with the normal order they are `1.0` and `0` **[verified 2026-08-28]**. Fixing it as part of this change is the right call, and `docs/visualisation/plan_apbs_vs_prodes_comparison_REVIEW1.md:293` independently found the same defect.
* **Every fixture number.** SG counts 8/8/9/3/4, detected pairs 4/4/0/0/0, zero `SSBOND` records in any shipped file, the four 1GDW distances (2.022, 2.030, 2.031, 2.050) and the 6.541 A nearest non-bonded pair **[verified 2026-08-28]**. Also worth recording: ARH96693's closest SG-SG pair is 3.68 A and ARH98503's is 5.15 A, so both AlphaFold models are comfortably clear of the window.
* **The 1GDW charge table and the pI.** Formal charge at pH 8.5 `-1.0` -> `7.0`, pI `8.899` -> `10.342` **[verified 2026-08-28, prototype]**.
* **Every PROPKA claim.** `tests/data/1GDW.pka` yields 45 groups, all eight cysteines at `99.99`, and **zero** titratable-or-terminal residues omitted **[verified 2026-08-28]**. The issue's premise really is wrong, and the plan is right to say so and to say it in the README.
* **`Atom.charge` is the single choke point, and the plan's propagation claim is complete.** `Residue.pkas` (`residue.py:34`) is the only reader of `["pka"]` from `general_data.json` **[verified 2026-08-28: grep]**. `Residue.charge`, `Structure.charge`, `Structure.dipole`, `Structure.isoelectric_point`, `run.charged_atom_arrays`, the two surface-grid charge filters (`run.py:202`, `run.py:295`) and `geometry.map_ep_to_plane` (`geometry.py:283`) all route through it. `Residue.charged_atoms` dispatches on the entry key and needs no change, as the plan says.
* **`viz.py` needs no change.** `charged_residue_underlay` hard-codes `ASP+GLU` and `LYS+ARG+HIS` and never renders cysteine, so the picture cannot disagree with the numbers.
* **Multiprocessing is a non-issue, for the reason the plan implies rather than states.** No structure object is ever pickled: `run_tasks` maps a module-level function over integer indices and the children inherit state through fork. The three charge filters run in the *parent*, so the workers receive arrays from which bonded SG atoms are already absent. `disulfide_partner` creates a `Residue` <-> `Residue` reference cycle, which is shared copy-on-write and never copied. Worth one sentence in the plan so a future reader does not have to re-derive it.
* **Not adding a feature column.** Correct: `tests/test_feature_dictionary.py` pins the column list from the two YAML files, so a new column is a much larger contract change than the issue asks for.
* **"ARH96693 features bit-identical to `main`" is achievable and achieved.** Full-feature CSVs from `main` and the prototype diff clean **[verified 2026-08-28]**.
