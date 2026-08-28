# Implementation review 1: disulfide-bonded cysteines

**Reviewed:** working tree on `6-disulfide-cysteines`, uncommitted, against `main` at `b5e8409`.
**Lens:** correctness, adversarially. Assume there is a bug.
**Date:** 2026-08-28. Everything tagged **[verified 2026-08-28]** was run, not reasoned from the diff.
**Environment:** `conda activate prodes`, run from `/mnt/shared/prodes`; experiments in a scratch copy, source under review untouched.

---

## Verdict

**Ship.** No blocking defect was found. The chemistry, the numbers and the blast radius are all right: I recomputed every asserted charge, the isoelectric point and the whole 1AO6 fixture from the raw PDB text with an independent script and all of them agree exactly; the fix reaches every charge path; a 400-structure fuzz of `assign_disulfides` produced no crash and no invariant violation; the pKa-file rewrite is bit-identical to `main` on the shipped PROPKA fixture and the ARH96693 regression is bit-identical on all 106 features. The findings below are one cosmetic bug in a log message, four places where a guard exists but no test would notice if it were deleted, and a handful of plan divergences that are improvements or immaterial. None of them changes a reported number.

---

## Findings

### 1. Minor — the "two SG atoms" warning prints the residue name twice

`disulfides.py:cysteine_sulfurs` formats `"CYS %s carries %d SG atoms..."` with `residue_label(residue)`, and `residue_label` already begins with the residue name.

Observed output **[verified 2026-08-28]**, from a synthetic file with two `SG` atoms in one residue:

```
WARNING CYS CYS A1 carries 2 SG atoms, so two residues have been read as one; ...
```

**Fix:** drop the literal `CYS ` from the format string, or pass `residue.number` instead of `residue_label(residue)`.

### 2. Serious (tests) — four guards survive deletion with the whole suite green

Mutation runs, each against `tests/test_disulfides.py tests/test_residue.py tests/test_structure.py tests/test_pka.py tests/test_output_bundle.py` (95 tests, baseline green) **[verified 2026-08-28]**:

| Mutation | Result |
|---|---|
| `record_disulfides`: symmetry-operator filter replaced by `if False:` | **95 passed** |
| `record_disulfides`: self-pair filter (`(chain1, number1) == (chain2, number2)`) replaced by `if False:` | **95 passed** |
| `record_disulfides`: duplicate-claim guard (`id(first) in claimed or ...`) replaced by `if False:` | **95 passed** |
| `cysteine_sulfurs`: `sulfurs.extend(found[:1])` → `sulfurs.extend(found)` | **95 passed** |

The first two survive for the same reason: the only test that exercises either, `test_a_symmetry_record_does_not_suppress_the_real_bonds`, uses a record that is *both* a self-pair *and* symmetry `2655`, so each filter covers for the other. That is the 1ERT record verbatim, but it means the two rules of plan §3.2 are tested as one. A record with a non-identity operator naming **two different** residues — the case that matters, because honouring it would bond two cysteines 30 Å apart — is currently caught only by code no test defends. It does work today **[verified 2026-08-28]**: a `2655` record naming `A1` and `A2` at 30 Å yields `structure.disulfides == []`.

The duplicate-claim guard matters because without it `len(structure.disulfides)` double-counts, and that number is printed to the user and written into `prodes_run.json`. It also works today: two identical `SSBOND A1-A2` records give one bond plus `WARNING two SSBOND records claim the same cysteine, CYS A1 or CYS A2; keeping the first` **[verified 2026-08-28]**.

`found[:1]` matters because without it a residue that merged two cysteines under one number can be paired **with itself** by geometry, which would silence a genuinely titratable thiol.

**Fix:** four small tests — a `2655` record over two distinct residues, a `1555` self-pair, a duplicated record asserting `len(structure.disulfides) == 1`, and a residue with two `SG` atoms 2.05 Å apart asserting no self-bond. All four are two-line additions to the existing synthetic-file helpers in `tests/test_disulfides.py`.

### 3. Minor — a record dropped for its symmetry operator is dropped silently

`record_disulfides` `continue`s on a non-identity operator and on a self-pair with no log at all. Plan §3.2 filters 1 and 2 say "skipped" and do not require a warning, so this is plan-conformant, but it interacts with `read_ssbond_line`: a line truncated *inside* the symmetry field, rather than before it, yields a partial string that is not `"1555"` and the record vanishes without a trace.

**[verified 2026-08-28]** `read_ssbond_line` on a line ending at column 69 returns `('A', 1, 'A', 2, '1555', '555')`, and `'555'` is then treated as a crystallographic operator.

Impact is low: a cysteine no record claimed falls through to geometry, so a wrongly dropped record usually still produces the bond. It only loses a bond that geometry cannot see.

**Fix:** treat a symmetry field shorter than four characters as the identity, and log the symmetry skip at debug.

### 4. Minor — `assign_disulfides` is not idempotent: stale partners survive re-detection

The function clears each newly bonded residue's `_pka` cache but never clears `disulfide_partner` on residues that are *no longer* bonded. Calling it a second time with different records leaves the first call's partners in place. `tests/test_disulfides.py::test_assign_disulfides_can_be_called_after_the_pkas_have_been_read` works around this by resetting `disulfide_partner = None` on every residue by hand before the second call, so the test name overstates what is guaranteed.

Not reachable in the shipped pipeline — the parser is the only caller and calls it once — but the docstring says "Safe to call on a structure whose pKas have already been read", which reads as a promise about arbitrary call sites.

**Fix:** one loop at the top of `assign_disulfides` setting `residue.disulfide_partner = None` and `residue._pka = None` for every `CYS`, or narrow the docstring.

### 5. Minor — `mark_disulfide_partner` throws away the whole pKa list, not just the side chain

`self._pka = None` discards any values `redo_pkas` had already merged onto that residue, including its terminus. Only bites when `assign_disulfides` runs *after* `redo_pkas`, which `prepare_structure` never does (parse assigns, then redo merges), so no shipped path loses data. It becomes reachable the moment anyone re-detects on a structure that has had a pKa file applied.

**Fix:** drop only the entry keyed by the residue's own name instead of clearing the list, or accept it and say so in the docstring.

### 6. Minor — plan divergences, all benign, one of them a wrong claim in the plan

| Plan | Implementation | Verdict |
|---|---|---|
| §5.1 "Logged at debug level" for a real pKa on a bridged cysteine | logged at `logger.info` | Undocumented, immaterial |
| §4.3 lists `disulfide_partner`, `side_chain_pka`, `group_pka` | also adds `set_group_pka`, `mark_disulfide_partner` on `Residue` and `titratable_groups` on `Structure` | Improvement; both are needed by §5.1, but they are new public API the plan did not budget for |
| §4.1 lists four functions | also `residue_label`, `sulfur_distance`, `candidate_pairs`, `residues_by_chain_and_number`, `report_record_distance`; all public, no leading underscores | Improvement, and `candidate_pairs` is what makes the greedy assignment directly testable |
| §7.2 "one at 1.5 A is not and is logged" | the test uses 0.5 Å | Immaterial |
| §5.3 "A pKa file with a titratable group removed warns, **naming it**" | `test_omitted_groups_are_reported` asserts only the phrase `"are not in the pKa file"`, not that the group is named | Minor test gap; the code does name the first five |
| §4.2 "Under `parse(file, identifier="HETATM")` ... geometry finds nothing because there are no CYS residues" | **false.** A file whose CYS residues are `HETATM` records parses into `CYS` residues and geometry bonds them: `parse(..., identifier="HETATM")` on two `HETATM` cysteines 2.05 Å apart returns **1 disulfide** **[verified 2026-08-28]** | The *behaviour* is defensible (arguably correct); the plan's stated reasoning is wrong. `SSBOND` records are still ignored on that path, so records-vs-geometry precedence silently differs by identifier |

### 7. Minor — `tests/data/test.csv` still records the old 1GDW numbers

Rows carry `8.899` for the isoelectric point and `6.734965...` for the average charge. Nothing references the file **[verified 2026-08-28]**: `grep -rn "test\.csv" --include=*.py --include=*.rst --include=*.md .` returns nothing. Pre-existing dead fixture, but it is now a committed record of the behaviour this change corrects.

**Fix:** delete it, in this change or another.

### 8. Minor (environment, not code) — the installed metadata still says 5.0.0

`pyproject.toml` is bumped to 6.0.0 but the editable install in the `prodes` env is `__editable__.prodes-5.0.0.pth`, so `prodes_version()` returns `5.0.0` **[verified 2026-08-28]** and every bundle written from this working tree stamps the old version into `prodes_run.json`. Nothing in the suite catches it. Re-install before generating anything users will see.

---

## Adversarial results: nothing broke

### Charge paths

Every path that produces a charge was traced and exercised.

* `Atom.charge` side-chain branch, formal and non-formal — reads `residue.side_chain_pka`, `None` short-circuits to 0. This is the single branch the whole fix rests on.
* `Atom.charge` terminus branch, formal and non-formal — reads `group_pka("N+"|"C-")`, `None` branch present.
* `Residue.charge`, `Residue.charged_atoms`, `Structure.charge`, `Structure.dipole`, `Structure.isoelectric_point`, `run.charged_atom_arrays`, the surface-grid phase and the shell phase all reach charge through those two branches and needed no change.

`Residue.charged_atoms` correctly stops listing a bridged `SG`, because it dispatches on the entry key and there is no `CYS` entry to dispatch on.

Two guards in `Atom.charge` are provably dead but harmless: `self.name in ("N", "C")` (terminus is only ever `"N"` or `"C"`), and reading the terminus by name rather than `pkas[-1]` (with the §5.1 merge in place, the terminus is always last). Both mutations pass the suite **[verified 2026-08-28]**; keep them, they are cheap.

### Hostile inputs to detection

**[verified 2026-08-28]**, synthetic files through `PDBparser().parse`:

| Input | Result |
|---|---|
| no cysteines | `[]` |
| one cysteine | `[]` |
| cysteine with no `SG` | `[]`, no crash |
| duplicate `SSBOND` records | 1 bond, warned |
| contradictory records `A1-A2` and `A1-A3` | `A1-A2` kept, warned, `A3` left to geometry |
| record naming a non-`CYS` residue | skipped with warning; the other cysteines still pair geometrically |
| blank chain IDs, in `ATOM` and `SSBOND` alike | resolves correctly to chain `""` |
| negative residue numbers (`-5`, `-4`) | parsed and bonded |
| three sulfurs at an exact distance tie | deterministic; `sorted(key=...[0])` never compares `Atom` objects, so no `TypeError` |
| coincident sulfurs (0.00 Å) | not bonded, warned |
| residue carrying two `SG` atoms | warned, first `SG` used, no self-bond |
| empty file / `SSBOND` records but no `ATOM` records | `IndexError` at `residues[-1].terminus = "C"` — **pre-existing**, that line is unchanged by this diff and fires before `assign_disulfides` |
| `read_ssbond_line` on `""`, `"SSBOND"`, non-numeric fields | `None` plus a warning, never raises |

**Fuzz:** 400 random structures (1–10 residues, 2 chains, repeated and negative numbers, mixed `CYS`/`ALA`, 0–4 records with random symmetry operators and dangling residue references). No exception, and on every one: no residue partnered with itself, the partner relation symmetric, only `CYS` partnered, no residue in two bonds, `len(partnered) == 2 * len(disulfides)`, and every bridged `SG` neutral at pH 14 both formal and non-formal. **[verified 2026-08-28]**

### `redo_pkas`

**[verified 2026-08-28]**, on 1GDW and 1AO6:

* Bridged cysteine given `{"CYS": 9.0}` — dropped, `pkas is None`, `charge(11) == 0`, logged at info.
* Called twice — second call merges over the first, no accumulation.
* Residue number not in the structure — no-op, no crash.
* Empty dict — no-op.
* Non-numeric value — propagates to a `TypeError` in `Atom.charge`. Unreachable through `read_pka`, which casts with `float()`.
* Group the residue does not have (`{"SER": 10.0}` on a `PRO`) — dropped with a warning; previously this raised `TypeError` in `charged_atoms`.
* Duplicate residue numbers across chains (1AO6 `A34`/`B34`) — both receive the value. Pre-existing chain-blind keying, plan §10.1, unchanged by this diff and not made worse.

**`id(residue)` bookkeeping is sound.** Every residue whose id enters `applied` is held alive for the whole call by `self.residues` (a NumPy object array holds strong references), and `titratable_groups()` reads from the same array, so no id can be recycled between the two passes. The same argument covers `claimed` in `record_disulfides` and `paired` in `geometric_disulfides`.

**It preserves the old behaviour where it should**, and the strongest evidence is end-to-end: 1GDW run with the shipped PROPKA file, `--full-features`, `--ionic-strength 0`, `main` versus branch — **0 of 106 features differ** **[verified 2026-08-28]**. Merging instead of replacing, and dropping the 99.99 entries instead of applying them, changes nothing a user sees.

### Numbers in the tests, recomputed from scratch

An independent script read `1GDW.pdb.zip` directly, applied Henderson–Hasselbalch to the residue table in `general_data.json`, detected disulfides by an independent 1.6–2.5 Å greedy pairing, and bisected for the pI. It shares no code with `prodes`. **[verified 2026-08-28]**

| Quantity | Recomputed | Asserted |
|---|---|---|
| disulfides | 4: 77–95 (2.022 Å), 6–128 (2.030), 30–116 (2.031), 65–81 (2.050) | same, `test_lysozyme_has_four_disulfides` |
| charge(7.0, formal) / (7.0, average) | 7.0 / 7.092 | 7.0 / 7.092 |
| charge(8.5, formal) / (8.5, average) | 7.0 / 6.709 | 7.0 / 6.709 |
| charge(9.0, formal) / (9.0, average) | 7.0 / 6.133 | 7.0 / 6.133 |
| charge(11, formal) | −5.0 | −5, `test_structure.py::test_charge` |
| isoelectric point | 10.342 (8.899 with cystines titrating) | 10.342 |

The 1AO6 fixture was validated the same way, straight from the zip **[verified 2026-08-28]**: 1156 residues, 2 chains, 70 cysteines, 34 `SSBOND` records all `1555/1555`, 34 geometric bonds with the longest at 2.048 Å and the nearest non-bonded pair at 5.951 Å, free thiols exactly `A34` and `B34`, and the record set equal to the geometric set. No insertion codes, no altLocs, no `CYS` missing an `SG`. 187 KB, 0.41 s to parse. Every claim in plan §7.1 holds.

README numbers restated in §6 of the plan also check out: on the reduced 54-feature set, 1GDW with versus without the PROPKA file now differs in **19** features, with the pI moving 10.342 → 10.565 and formal charge +7 → +8 **[verified 2026-08-28]** — exactly what the README now says.

### Multiprocessing

`parallel.py` forks and passes only an integer index; workers read `SURFACE_GRID_STATE` / `AVERAGE_CHARGE_GRID_STATE` / `SHELL_STATE`, which the parent fills before the pool exists. `disulfide_partner` is set at parse time, long before any fork, and the charged-atom arrays are built in the parent, so the new per-residue state is inherited through copy-on-write and never crosses a pickle boundary. Empirically, 1GDW with `PRODES_N_WORKERS=1` and `PRODES_N_WORKERS=6` produce **identical values on all 106 features** **[verified 2026-08-28]**.

### Regression guard

ARH96693, `--full-features --ionic-strength 0`, `main` versus branch: **0 of 106 features differ** **[verified 2026-08-28]** — exact equality, stronger than the tolerance-based comparison in `tests/test_sasa.py`. `test_output_bundle.py` additionally pins `record["disulfides"] == 0` for it.

Full suite: **297 passed in 119 s**. `ruff check src tests scripts`, `black --check src tests scripts` and `mypy src` all clean **[verified 2026-08-28]**.

---

## Acceptance criteria

| # | Criterion | Status | Evidence |
|---|---|---|---|
| 1 | A test structure with known disulfides reports the correct number detected | **Met** | 1GDW 4 bonds with exact partner numbers; 1AO6 34 bonds; both recomputed independently from the raw PDB **[verified 2026-08-28]**. `test_lysozyme_has_four_disulfides`, `test_albumin_is_detected_from_its_ssbond_records`, `test_geometry_alone_reproduces_the_albumin_records` |
| 2 | Formal charge and isoelectric point for a cystine-rich structure change in the expected direction and magnitude at pH 8.5 | **Met** | 1GDW formal charge at pH 8.5 goes −1.0 → **+7.0**, pI 8.899 → **10.342**; both recomputed from scratch and matching **[verified 2026-08-28]**. `test_lysozyme_charge_is_no_longer_dragged_down_by_its_cystines[8.5-True-7.0]`, `test_lysozyme_isoelectric_point`, `test_structure.py::test_isoelectric_point` |
| 3 | A free-thiol Cys in the same structure is still titrated normally | **Met** | 1AO6 `A34` keeps `[{"CYS": 8.33}]` and titrates (0 at pH 7, −1 at pH 9 and 14) while `A53` is bonded to `A62` and stays 0. `test_residue.py::test_a_free_cysteine_keeps_the_thiol_pka`, `test_a_free_cysteine_titrates_and_a_bonded_one_does_not` |
| 4 | Behaviour is unchanged for structures with no disulfides | **Met** | ARH96693 bit-identical on all 106 features against `main` **[verified 2026-08-28]**; 1GPB reports 9 cysteines and 0 bonds (`test_free_cysteines_are_not_paired`); `tests/test_sasa.py` untouched and green |
| 5 | The pKa-file interaction from point 4 is covered by a test | **Met** | `test_omitted_groups_are_reported`, `test_a_complete_pka_file_is_reported_silently`, `test_a_group_the_file_omits_keeps_its_default`, `test_a_pka_for_a_group_the_residue_does_not_have_is_dropped`, `test_a_predicted_pka_cannot_titrate_a_bonded_cysteine`. Note the criterion changed shape: PROPKA does **not** omit bridged cysteines, it writes 99.99, which plan §2.4 established and `test_propka_reports_bonded_cysteines_as_not_titratable` now pins. The group counts behind it check out: `1GDW.pka` covers 46 groups over 45 residues, and the structure has 38 titratable groups once the 8 bridged cysteines are excluded, so the omission warning is correctly silent **[verified 2026-08-28]** |

---

## What is solid — do not churn it

* **The detection criterion and the pairing.** 2.5 Å upper bound, 1.6 Å lower bound, greedy shortest-first over a `numpy` pairwise matrix, one partner per sulfur. Correct on both real fixtures, correct under 400 fuzzed structures, and `candidate_pairs` sorts by distance alone so it never tries to order `Atom` objects.
* **Per-cysteine precedence.** Records first, geometry over what is left. `test_an_incomplete_record_set_keeps_the_bonds_it_forgot` and `test_a_record_wins_the_cysteines_it_names` both fail under the obvious mutations, and `test_geometry_alone_reproduces_the_albumin_records` is a genuine cross-check on 34 real bonds.
* **Where the fix is applied.** Removing the side-chain entry in `Residue.pkas` rather than special-casing `charge()` is the right level: it reaches `Atom.charge`, `Residue.charged_atoms`, `Structure.charge`, `dipole`, `isoelectric_point`, the grid and the shell with no feature code touched. Mutating that one condition fails 15 tests across four modules.
* **Keying pKa lookups by group name.** Necessary, correct, and defensively written. Reverting the side-chain half to `pkas[0]` fails `test_a_bonded_cysteine_at_a_chain_terminus_keeps_its_terminal_charge`, which is exactly the case the plan predicted would crash.
* **The `redo_pkas` merge.** Reverting it to list replacement fails three tests, and it is bit-identical to `main` on the shipped PROPKA fixture. Both properties at once is what you want from this rewrite.
* **The 1AO6 fixture.** Every property claimed for it is true, it is the only thing in the repository that can test a free thiol and a cystine side by side, and 187 KB / 0.41 s is a fair price.
* **Docs and version.** The README numbers were actually recomputed and are right; 6.0.0 is the correct bump for a changed reported quantity; the dated caveat on `identification_of_propka_dependent_features.md` is the honest way to handle an analysis that cannot be re-run.
