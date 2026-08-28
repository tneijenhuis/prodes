# Implementation plan: stop charging disulfide-bonded cysteines

**Status:** plan only, nothing implemented.
**Branch:** `6-disulfide-cysteines`, cut from `main` at `b5e8409`.
**Issue:** [#6](https://github.com/datacatalysis/prodes/issues/6), "Cystines are charged as free thiols, distorting charge features above pH 8".
**Date:** 2026-08-28.
**Provenance:** first draft, written against the code at `b5e8409`. Every number below marked **[verified 2026-08-28]** was recomputed in the `prodes` conda environment rather than read out of the issue.

---

## Goal

A cysteine whose `SG` is bonded to another cysteine's `SG` has no thiol proton and cannot titrate. Prodes currently gives every `CYS` a pKa of 8.33 regardless of structure, so both halves of every disulfide are titrated as free thiols. Above about pH 8 this puts a large spurious negative charge on the protein.

After this change, a cysteine that Prodes judges to be in a disulfide carries no side-chain pKa and therefore no side-chain charge, at any pH, by any route.

---

## Background

### The data model, and where charge comes from

`PDBparser.parse` builds `Structure -> Chain -> Residue -> Atom`. Only lines whose record name equals the parser's `identifier` are read, and that defaults to `"ATOM"`; every other record type in the file, `SSBOND` included, is discarded without being looked at. **[verified 2026-08-28]**

Two objects decide charge, and they consult different sources:

| Caller | What it reads | Keyed by |
|---|---|---|
| `Residue.pkas` | `data.residue_data(self.name)["pka"]`, plus the terminus defaults | residue *type*, from `general_data.json` |
| `Residue.charged_atoms` | `self.pkas`, dispatching on the key of each entry | the pKa entry's key |
| `Atom.charge` | `data.residue_data(self.residue_name)["potential_charge"]` and `["charged_atoms"]` to decide *whether* the atom is charged, then `self.residue.pkas[0]` (side chain) or `self.residue.pkas[-1]` (terminus) for the *value* | residue type for the decision, list **position** for the value |

`Residue.pkas` returns a list of single-entry dicts such as `[{"LYS": 10.5}, {"N+": 9.69}]`, and returns `None` for a residue that is neither titratable nor terminal. The side-chain entry is appended before the terminus entry, which is the only reason `pkas[0]` and `pkas[-1]` currently pick the right ones.

`Structure.charge` sums `Residue.charge`, which sums `Atom.charge` over heavy atoms. `Structure.dipole` and every surface and shell electrostatic feature go through `Atom.charge` as well, so a fix at the pKa level propagates to all of them without touching any feature code.

### Why `Atom.charge` has to change too

Setting a bonded cysteine's pKa to `None`, as the issue proposes, is not sufficient on its own, and on its own it crashes.

`Atom.charge` asks `residue_data(self.residue_name)` whether the atom is a charged atom of its residue type. For `SG` in a `CYS` that answer is always yes, whatever this particular residue's pKa list says. It then indexes `self.residue.pkas[0]` for the value. With the side-chain entry removed, `pkas` is either `None` (`TypeError: 'NoneType' object is not subscriptable`) or, for a terminal cysteine, `[{"C-": 2.34}]`, in which case `pkas[0]` hands the C-terminal pKa of 2.34 to the `SG` atom and gives it a full negative charge at every pH above 2.3. Verified by emptying a residue's `_pka` and calling `charge`: `IndexError` / wrong value, never a silent zero. **[verified 2026-08-28]**

So `Atom.charge` must stop reading the pKa by list position and start reading it by group name. That is the same defect as an existing latent bug: with the entries in the other order, `[{"N+": 8.35}, {"LYS": 11.34}]`, the `NZ` atom of Lys 1 reads the N-terminal pKa and comes out at charge 0 at pH 10 instead of 1.0. **[verified 2026-08-28]** Both predictors we ship converters for happen to emit the side chain first, so this has never bitten, but it is the mechanism this change depends on and it should be made correct rather than relied upon.

### How bad it is now

`tests/data/1GDW.pdb.zip` is lysozyme: 130 residues, 8 cysteines, all 8 in disulfides, no free thiol. Its `SG`-`SG` distances **[verified 2026-08-28]**:

| Pair | Distance (A) |
|---|---|
| 77-95 | 2.022 |
| 6-128 | 2.030 |
| 30-116 | 2.031 |
| 65-81 | 2.050 |
| next nearest non-bonded pair (65-95) | 6.541 |

Bonded and non-bonded are separated by more than 4 A, so the cutoff is not a delicate choice on this structure.

Net charge of 1GDW as shipped, and with the 8 cysteines silenced **[verified 2026-08-28]**:

| pH | Formal now | Formal fixed | Average now | Average fixed |
|---|---|---|---|---|
| 7.0 | 7.000 | 7.000 | 6.735 | 7.092 |
| 7.4 | 7.000 | 7.000 | 6.178 | 7.019 |
| 8.5 | **-1.000** | **7.000** | 1.936 | 6.709 |
| 9.0 | **-1.000** | **7.000** | -0.458 | 6.133 |

Isoelectric point moves from **8.899 to 10.342**. At pH 8.5 the formal charge is wrong by 8 e on a protein whose real formal charge is +7, i.e. the reported number has the wrong sign. This reproduces the issue's arithmetic on a real structure.

### The issue's PROPKA premise is wrong, and this matters

The issue says PROPKA "generally omits bonded cysteines from its titratable list", and asks for that to be confirmed first. It was, and it is not what PROPKA does.

`tests/data/1GDW.pka`, a genuine propka3.4.0 run committed with the repository, lists **all eight** bonded cysteines in its summary table, each with a pKa of `99.99` **[verified 2026-08-28]**. That is PROPKA's sentinel for "this group does not titrate", not a prediction. `convert_propka` reads it as the float 99.99 and `redo_pkas` applies it, and because `neg_charge(99.99, pH)` underflows to zero and `99.99 < pH` is false, the residue ends up neutral at every reachable pH. The current code therefore already gets a PROPKA-supplied cystine right, by accident, through a magic number.

Consequences for this plan:

1. **The bug is a no-pKa-file bug.** With `--pka` from PROPKA the charge features are already close to correct. Without it they are badly wrong. That is the opposite of the issue's framing, and it should be stated in the README, because "use PROPKA" is the recommended default workflow and users following it are less affected than users who are not.
2. **Point 4 of the issue** ("make `redo_pkas` not silently fall back for residues the pKa file omits") is still worth doing, but it is not the safety net for cystines that the issue thought it was. On 1GDW, PROPKA omits **zero** titratable residues **[verified 2026-08-28]**, so the warning is silent on our own fixture, which is the right behaviour for a warning.
3. The fix must not fight the sentinel. A structure-detected disulfide and a 99.99 from the file must agree, and they do: both mean no charge.

### What the shipped fixtures can and cannot test

**[verified 2026-08-28]**

| Fixture | ATOM lines | `SSBOND` records | `SG` atoms | Detected pairs |
|---|---|---|---|---|
| `1GDW.pdb.zip` | 1022 | 0 | 8 | 4 |
| `1GDW_h.pdb.zip` | 2003 | 0 | 8 | 4 |
| `1GPB.pdb.zip` | 6699 | 0 | 9 | 0 |
| `ARH96693.pdb.zip` | 479 | 0 | 3 | 0 |
| `ARH98503.pdb.zip` | 3106 | 0 | 4 | 0 |

So: 1GDW covers all-bonded, 1GPB covers many-free-thiols-no-bonds, and the two AlphaFold models cover the no-disulfide case that must not change. **Nothing shipped has an `SSBOND` record and nothing shipped mixes bonded and free cysteines in one structure**, so both of those need a small hand-written fixture. `1GDW.pdb.zip` also has no alternate locations (`altLoc` is blank on every atom), so the altLoc question below cannot be answered from the fixtures either.

### One quirk worth knowing before writing the detection

The parser takes the atom name from columns 13-17, which includes the `altLoc` column: `name = line[12:17].strip()`. An `SG` atom with `altLoc` `A` therefore parses as the atom name `"SG A"`, not `"SG"`. Such an atom already matches nothing in `charged_atoms` and already carries no charge today, so it is not a regression, but the detection must be written knowing that a disordered cysteine will simply not be seen. Out of scope to fix; worth a comment so the next reader does not think it was overlooked.

---

## Phase 1: detect the disulfides and remove the side-chain pKa

### New module, `src/prodes/calculations/disulfides.py`

Pure functions over parsed objects, no I/O.

```
SG_SG_DISTANCE_ANGSTROM = 2.05      # mean Cys SG-Cys SG bond length
SG_SG_TOLERANCE_ANGSTROM = 0.3      # accepted window is 1.75 to 2.35 A
```

* `sulfur_atoms(structure)` — the `SG` atoms of every `CYS` residue, in parse order.
* `geometric_disulfides(structure)` — returns a list of `(residue, residue)` pairs. All `SG`-`SG` distances inside the window are collected, sorted shortest first, and assigned greedily so that each `SG` takes at most one partner and the shortest candidate wins a contested sulfur. The count of cysteines is small (single digits to low tens), so a plain pairwise loop over `numpy` coordinates is used; no spatial index is warranted.
* `record_disulfides(structure, ssbond_pairs)` — resolves `(chain, number)` pairs read from `SSBOND` records to `Residue` objects, ignoring any pair naming a residue that is not in the coordinates.
* `assign_disulfides(structure, ssbond_pairs=())` — the entry point. Uses the `SSBOND` records when the file has any, otherwise geometry. Sets `residue.disulfide_partner` on both members of each pair, stores the list on `structure.disulfides`, and returns it.

**Records beat geometry, but disagreement is reported.** `SSBOND` is authoritative when present, as the issue asks. When a record names a pair whose `SG`-`SG` distance is outside the window, or names a residue absent from the coordinates, the discrepancy is logged at warning level and the record is still honoured for the pairs that resolve. A record carrying a symmetry operator other than `1555` in either symmetry field names a crystallographic copy that is not in this file's coordinates and is skipped.

Rationale for records-first: an `SSBOND` record encodes the depositor's chemistry, and a structure can legitimately hold a strained or poorly refined bond that sits outside a 0.3 A window. Rationale for logging the disagreement: a stale record on a reduced or mutated model is the one case where trusting it is wrong, and silence there would be indistinguishable from success.

### Parser hook

`PDBparser._read_pdb` gains a second thing to collect while it walks the file: lines whose record name is `SSBOND` are parsed by column into `(chain1, number1, chain2, number2, sym1, sym2)` in a small function `read_ssbond_line`, testable on a single string. After the atom loop and the existing C-terminus assignment, `_read_pdb` calls `assign_disulfides`.

Detection therefore happens for **every** parse, including `prodes.read()` and every test that builds a structure, not only inside the pipeline's `prepare_structure`. A `Structure` that has been parsed from a file is then always self-consistent about its cysteines, which is the property that keeps this from being re-broken by a future caller that builds a structure its own way.

The parser does not read `SSBOND` when the caller passes a non-default `identifier` (`parse(file, identifier="HETATM")`), because in that mode it is not reading protein atoms at all. Geometry still runs and finds nothing, since there are no `CYS` residues.

### `Residue`

* New attribute `disulfide_partner = None`, set in `__init__`.
* `pkas` skips the side-chain entry when `disulfide_partner is not None`. Terminus entries are unaffected: a bonded cysteine at a chain terminus still has its titratable backbone group.
* New property `side_chain_pka` — the value of the entry in `pkas` whose key is this residue's name, or `None` if there is none. This is the accessor `Atom.charge` will use.

The existing `_pka` cache makes ordering load-bearing: `pkas` is computed once, on first access, and cached. Detection runs inside the parser before any caller can touch `pkas`, so the cache cannot be stale. That ordering is stated in the docstring, and Phase 2's change to `redo_pkas` is the only other writer of `_pka`.

### `Atom.charge`

Replace both positional reads:

* side chain: `self.residue.side_chain_pka`, and if that is `None`, the atom carries no side-chain charge and the branch is skipped. This is what makes a bonded cysteine's `SG` neutral, and it is the single point through which the whole fix propagates to the surface and shell features.
* terminus: look up the `"N+"` or `"C-"` entry by key rather than taking `pkas[-1]`.

`Residue.charged_atoms` already dispatches on the entry key and needs no change; with no `CYS` entry present it simply never appends the `SG`.

### Acceptance for phase 1

* 1GDW reports 4 disulfides; each of the 8 cysteines has a partner; the partner relation is symmetric.
* 1GPB reports 0 despite having 9 cysteines.
* 1GDW formal charge at pH 8.5 is +7.0 rather than -1.0, and its isoelectric point is 10.342 rather than 8.899.
* A hand-built structure with one disulfide and one free cysteine titrates the free one and not the bonded one.
* `ARH96693` features are bit-identical to `main`.

---

## Phase 2: make the pKa file interaction honest

### `Structure.redo_pkas`

Two changes.

**A bonded cysteine keeps no side-chain pKa, whatever the file says.** When the supplied dict gives a `CYS` entry for a residue that detection marked as bonded, that entry is dropped and the drop is logged at debug level unless the value is PROPKA's 99.99 sentinel, which agrees with us and is not worth a line. Structure beats file here because the file's value is a prediction about a group that does not exist. If the remaining entries for that residue are empty, `_pka` is left at `None` so that the property recomputes the terminus defaults rather than being pinned to an empty list, which `Atom.charge` cannot read.

**Residues the file omits are reported, not silently defaulted.** After applying the dict, count the residues that are titratable by default (or terminal) and absent from it, and emit one warning naming the count and the first few. Bonded cysteines are excluded from that count, since their absence is correct. The message says plainly that those residues fall back to the textbook pKa, which is the fact the current silence hides.

This is one log line per run at most, and on `1GDW.pka` it is silent, because PROPKA covers all 45 titratable groups. **[verified 2026-08-28]**

### Acceptance for phase 2

* A pKa file that gives a bonded cysteine a real pKa of 9.0 does not make it charged.
* A pKa file with a titratable residue removed produces a warning naming that residue.
* `1GDW.pka` produces no warning.
* The end-to-end pKa test in `tests/test_pka.py` still passes.

---

## Phase 3: make the count visible, and say so in the docs

* `prepare_structure` prints the disulfide count alongside the existing `calculating {name}` line, in the same plain style the pipeline already uses.
* `run_metadata` gains `"disulfides": n`, so every bundle's `prodes_run.json` records what detection found. This is a metadata field, not a feature column: the feature set is pinned by `features_reduced.yaml` and `features_full_only.yaml` and asserted column by column in `tests/test_feature_dictionary.py`, so adding a column is a much larger contract change than this issue calls for.
* README: a subsection under **pKa values and protonation states** stating what is detected, the cutoff, that `SSBOND` wins when present, that AlphaFold models have no `SSBOND` and are handled geometrically, and that PROPKA already reports cystines as non-titratable so a run with `--pka` changes far less than one without. The "No structure preparation" bullet, which currently lists "form disulfide bonds" among the things Prodes does not do, is now half wrong and needs correcting: Prodes does not *form* them, but it does *recognise* them.
* `pyproject.toml` version 5.0.0 -> **6.0.0**. The reported quantity changes: formal charge, isoelectric point, average charge and every electrostatic feature move for any structure with a disulfide, at any pH where a thiol would have titrated. Per the versioning rule that is MAJOR, and the version is stamped into every bundle's `prodes_run.json` by `prodes_version`, so the drift would be visible to users. Also add the change to the README's improvements list.

---

## Phase 4: tests, `tests/test_disulfides.py` plus edits

New module for detection and charge consequences:

* `read_ssbond_line` on a single real `SSBOND` line, including one with a symmetry operator.
* Geometry on 1GDW: 4 pairs, the exact partner numbers, symmetry of the relation.
* Geometry on 1GPB: 9 cysteines, 0 pairs.
* A contested sulfur: three `SG` atoms where two candidate pairs share one atom, and the shorter wins.
* A pair just outside the window (2.4 A) is not called.
* `SSBOND` beats geometry: a hand-written fixture whose record names a pair the geometry would not find, and the inverse, a geometric pair the record set does not name, which is *not* called because records win when present.
* An `SSBOND` naming a residue that is not in the coordinates is skipped with a warning rather than raising.
* Charge consequences on 1GDW at pH 7.0, 8.5 and 9.0, against the table above.
* A mixed structure: bonded cysteine neutral, free cysteine titrating normally at pH 9.
* A bonded cysteine that is also the C-terminal residue keeps its C-terminal charge and gains no side-chain charge, which is the case that crashes without the `Atom.charge` change.

Fixtures: two small hand-written PDB files under `tests/data/`, one with `SSBOND` records and one mixing a bonded pair with a free thiol. Written as plain text files, small enough to read in a diff.

Edits to existing tests:

* `tests/test_structure.py` — `test_charge` at pH 11 and `test_isoelectric_point` change on 1GDW. Update to the new values with a comment saying why, since these are exactly the numbers the issue exists to move.
* `tests/test_pka.py` — add the two `redo_pkas` cases from phase 2.
* `tests/test_data.py`, `tests/test_residue.py` — check whether the `CYS` entry assertions still hold; `general_data.json` is not being edited, so they should.

Regression guard: `tests/data/ARH96693_prodes_orig_output.csv` has no disulfides, so `tests/test_sasa.py` must still pass untouched. If it does not, the change has leaked somewhere it should not have.

---

## Files touched

| File | Change |
|---|---|
| `src/prodes/calculations/disulfides.py` | new: detection, pairing, assignment |
| `src/prodes/io/parser.py` | collect `SSBOND` lines, call `assign_disulfides` |
| `src/prodes/core/residue.py` | `disulfide_partner`, `side_chain_pka`, skip bonded `CYS` in `pkas` |
| `src/prodes/core/atom.py` | read pKas by group name rather than list position |
| `src/prodes/core/structure.py` | `disulfides` attribute; `redo_pkas` hardening and warning |
| `src/prodes/run.py` | report the count; pass it to `run_metadata` |
| `src/prodes/output.py` | `disulfides` in the run metadata |
| `pyproject.toml` | 5.0.0 -> 6.0.0 |
| `README.rst` | disulfide subsection, correct the "no structure preparation" bullet, improvements list |
| `tests/test_disulfides.py` | new |
| `tests/data/*.pdb` | two small hand-written fixtures |
| `tests/test_structure.py`, `tests/test_pka.py` | updated expectations, new pKa-file cases |

## Risks

* **A 2.05 +/- 0.3 A window is a chemistry choice, not a fact.** On 1GDW the nearest non-bonded pair is at 6.5 A, so nothing is marginal there, but a low-resolution or coarse model could put a real bond at 2.4 A and be missed. Missing a bond fails safe in the sense that it reproduces today's behaviour, and the logged count is how a user notices.
* **Records-first can be wrong.** A reduced structure that kept its `SSBOND` records would be treated as oxidised. The distance check that logs a warning is the mitigation; making it an error would break files that are merely strained.
* **This changes shipped numbers.** Any model already fitted on Prodes features from disulfide-containing structures at pH above 8 will need refitting. That is the point of the change, but it is a MAJOR version bump and should be called out in the release notes, not just the diff.
* **Adjacent bug, deliberately not fixed here:** `redo_pkas` keys by residue number alone, so in a multi-chain structure a pKa predicted for chain A residue 30 is also applied to chain B residue 30. For an antibody, which is the workflow the issue names, that is wrong. It is independent of disulfides and deserves its own issue rather than being folded into this one; the disulfide code itself resolves residues by chain and number and does not inherit the defect.
