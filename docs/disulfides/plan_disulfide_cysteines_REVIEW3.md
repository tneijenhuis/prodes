# Adversarial review 3 — `plan_disulfide_cysteines.md`

**Reviewer lens:** refute. **Date:** 2026-08-28. **Code reviewed at:** `b5e8409`, branch `6-disulfide-cysteines` (verified: `git branch --show-current` = `6-disulfide-cysteines`, `git log -1` = `b5e8409`, working tree clean apart from untracked docs). **[verified 2026-08-28]**

## Verdict

I set out to find a wrong number and could not find one in the plan's quantitative core. Every SG-SG distance, every cell of the charge table, both isoelectric points, the whole fixture inventory, the PROPKA 99.99 claim, the zero-omissions claim and both crash modes reproduce exactly. That part of the plan is unusually well checked. What it does contain is a wrong *inference*: the plan asserts that `tests/test_residue.py` should survive untouched because `general_data.json` is not being edited — and that reasoning is a non-sequitur, because `test_residue.py` draws its `residues["CYS"]` from 1GDW_h, whose every cysteine is disulfide-bonded. Two assertions there fail under the plan, and the file is missing from the Files-touched table. Beyond that, the plan's own `Atom.charge` rewrite opens a new, reachable failure mode that its Phase-2 warning is structurally incapable of catching, it leaves three stale documentation claims uncorrected (including one the plan itself relies on), and it presents an alternative-free "has to change" where an alternative exists and is already in use in the code. **Do not merge as written.** All findings are edits to the plan, not a redesign; the science and the scope are otherwise sound.

## Re-verification table

| # | Claim in plan | Plan's value | Recomputed | Verdict |
|---|---|---|---|---|
| 1 | Parser reads only `identifier` lines; `SSBOND` discarded | — | `_read_pdb` filters `line[0:6].strip() == self.identifier`; nothing else is inspected | correct |
| 2 | SG-SG 77-95 | 2.022 | 2.022 | correct |
| 3 | SG-SG 6-128 | 2.030 | 2.030 | correct |
| 4 | SG-SG 30-116 | 2.031 | 2.031 | correct |
| 5 | SG-SG 65-81 | 2.050 | 2.050 | correct |
| 6 | next non-bonded pair 65-95 | 6.541 | 6.541 | correct |
| 7 | "separated by more than 4 A" | >4 | 6.541 − 2.050 = 4.491 | correct |
| 8 | Window 2.05 ± 0.3 | 1.75–2.35 | 1.75–2.35 | correct |
| 9 | Formal charge now @ 7.0 / 7.4 / 8.5 / 9.0 | 7 / 7 / −1 / −1 | 7.000 / 7.000 / −1.000 / −1.000 | correct |
| 10 | Formal charge fixed @ same | 7 / 7 / 7 / 7 | 7.000 / 7.000 / 7.000 / 7.000 | correct |
| 11 | Average charge now | 6.735 / 6.178 / 1.936 / −0.458 | identical to 3 dp | correct |
| 12 | Average charge fixed | 7.092 / 7.019 / 6.709 / 6.133 | identical to 3 dp | correct |
| 13 | pI moves 8.899 -> 10.342 | | 8.899 -> 10.342 | correct |
| 14 | 1GDW: 130 residues, 8 Cys, all bonded, no free thiol | | 130 residues, 8 CYS, 4 pairs, 0 free | correct |
| 15 | `1GDW.pka` lists all 8 Cys at 99.99 | 8 | 8 (lines 208-215 of the file) | correct |
| 16 | PROPKA omits **zero** titratable residues on 1GDW | 0 | 0 (45 default-titratable residue numbers, 45 propka keys, symmetric difference empty) | correct |
| 17 | "PROPKA covers all **45 titratable groups**" | 45 groups | **46 groups over 45 residues** (LYS 1 and N+ 1 share a residue) | **wrong** — see F6 |
| 18 | Fixture ATOM lines 1022 / 2003 / 6699 / 479 / 3106 | | 1022 / 2003 / 6699 / 479 / 3106 | correct |
| 19 | Fixture SSBOND records all 0 | | 0 everywhere; record types present are only ATOM/TER/CRYST1/END | correct |
| 20 | Fixture SG counts 8 / 8 / 9 / 3 / 4 | | 8 / 8 / 9 / 3 / 4 | correct |
| 21 | Detected pairs 4 / 4 / 0 / 0 / 0 | | 4 / 4 / 0 / 0 / 0 under the plan's own window (1GPB min SG-SG 8.309; ARH96693 3.676; ARH98503 5.148) | correct |
| 22 | No altLoc in any fixture | | `line[16]` is blank on every ATOM line of all five | correct |
| 23 | `SG` with altLoc A parses as `"SG A"` | | `_read_line` on a synthetic altLoc line returns name `'SG A'` | correct |
| 24 | Removing the side-chain entry gives `TypeError` when `pkas` is `None` | TypeError | `TypeError: 'NoneType' object is not subscriptable` | correct |
| 25 | "emptying a residue's `_pka`" gives `IndexError` | IndexError | `IndexError: list index out of range` — but **unreachable by the plan's own design** | correct-but-irrelevant, see F9 |
| 26 | Terminal Cys with `[{"C-": 2.34}]` gives SG a full negative charge above 2.3 | −1 | `SG.charge(9, formal=True)` = −1.0; `charge(3)` = −1.0; average −0.99999978 | correct |
| 27 | Latent ordering bug: reversed order makes Lys 1 `NZ` read 0 at pH 10 instead of 1.0 | 0 vs 1.0 | NZ 1.0 -> 0 | correct but incomplete — the `N` atom simultaneously goes 0 -> 1 (F10) |
| 28 | "Both predictors we ship converters for" emit side chain first | 2 predictors | **3 converters shipped** (`convert_hpp`, `convert_propka`, `convert_pypka`); README line 393 says "three predictors" | **wrong count**, F7 |
| 29 | `redo_pkas` keys by residue number alone (multi-chain defect) | | confirmed in `Structure.redo_pkas`; `convert_propka` also drops the chain ID | correct |
| 30 | `prepare_structure` prints the count "alongside the existing `calculating {name}` line" | | that `print` is in `calculate()` (`run.py:526`), two calls *after* `prepare_structure` | **wrong**, F8 |
| 31 | `tests/test_residue.py` CYS assertions "should" still hold | hold | **two fail** (F1) | **wrong** |
| 32 | ARH96693 / ARH98503 are "AlphaFold models" | | files carry no header, no `REMARK`, only `CRYST1 1.000 1.000 1.000 P 1`; consistent with a predicted model but not stated anywhere in the file | unverified inference stated as fact (minor) |

Additional numbers I computed that the plan does not state, and should:

| Quantity (1GDW, default settings, pH 7) | Value |
|---|---|
| Reduced features (54) changed by the fix | **1** — `Isoelectric point` only |
| Full features (105) changed by the fix | **20** — `Isoelectric point`, `Average charge`, and the 18 `SurfEp*Average` columns |
| Features changed by the fix when `--pka` is supplied | **0** |
| PROPKA-vs-default feature delta, today, default ionic strength | 19 of 54 (README says 18; 18 is the value at `--ionic-strength 0`) |
| pI with PROPKA, before and after the fix | 10.565 both |

---

## Findings

### F1 — BLOCKING. `tests/test_residue.py` breaks, and the plan's reason for saying it will not is a non-sequitur.

**Plan says** (Phase 4): "`tests/test_data.py`, `tests/test_residue.py` — check whether the `CYS` entry assertions still hold; `general_data.json` is not being edited, so they should."

**Actually true.** `tests/test_residue.py:6` builds its `residues` dict from `tests/data/1GDW_h.pdb.zip`, taking the first residue of each type. The first CYS encountered is **CYS 6**, which is bonded to CYS 128 at 2.030 A. **[verified 2026-08-28]** Under the plan, `residues["CYS"].pkas` becomes `None` and its side-chain charge becomes 0, so:

* `test_pkas`, line 141: `assert residues["CYS"].pkas == [{"CYS": 8.33}]` -> fails (`None`).
* `test_charge`, line 195: `assert round(residues["CYS"].charge(14), 3) == -1` -> fails (0).
* line 170 (`residues["CYS"].charge(8) == 0`) survives.

The stated reason is wrong on its own terms: the change does not come from `general_data.json`, it comes from `Residue.pkas` gaining a structural precondition, so whether the JSON is edited is irrelevant to whether these assertions hold. The **Files touched** table omits `tests/test_residue.py` entirely while Phase 4 mentions it — an internal contradiction.

**Fix.** Move `tests/test_residue.py` into the Files-touched table as *updated*. Either (a) pick a free-thiol cysteine for the CYS row by sourcing it from 1GPB, or (b) keep 1GDW_h and assert the new truth (`pkas is None`, `charge(14) == 0`) plus a second assertion on a free thiol from 1GPB, so the file tests both states rather than losing coverage of the titrating one.

### F2 — BLOCKING. The `Atom.charge` terminus rewrite has an unhandled, reachable case, and the Phase-2 warning cannot catch it.

**Plan says:** "terminus: look up the `"N+"` or `"C-"` entry by key rather than taking `pkas[-1]`." It never says what happens when there is no such entry.

**Actually true.** That case is reachable today, without any hand-crafted input, through `redo_pkas`. `redo_pkas` replaces a residue's entire `_pka` list with whatever the file holds for that residue number. A pKa file that names the side chain of the N-terminal or C-terminal residue but not its terminus therefore deletes the terminus pKa. Demonstrated on 1GDW **[verified 2026-08-28]**:

```
before                      [{'LYS': 10.5}, {'N+': 9.69}]
redo_pkas({1: [{"LYS": 11.34}]})
after                       [{'LYS': 11.34}]
N atom .charge(7, formal=True)  ->  1     # from pkas[-1] == LYS 11.34, not from N+
```

Today this is silently wrong (right at pH 7, wrong at pH 10: it returns 1 where N+ 9.69 gives 0). After the rewrite it is either a `KeyError`, a `None` dereference, or a silent zero, depending on which the implementer picks — the plan does not say, so two people build two different things.

Worse, the Phase-2 warning is specified as "count the **residues** that are titratable by default (or terminal) and absent from it". Residue 1 *is* in the dict, so no warning fires. The warning is at residue granularity while the defect is at group granularity.

**Fix.** (a) State the behaviour explicitly: a terminal atom whose residue has no `"N+"`/`"C-"` entry falls back to `data.n_term_pka()` / `data.c_term_pka()` and logs, or raises — pick one and write it down. (b) Re-specify the Phase-2 warning to count **groups** (side chain, N+, C-) rather than residues, so a file that supplies a side chain and drops the terminus is reported. This is the same root cause as F6.

### F3 — SERIOUS. Three documentation claims go stale and none is in the plan's README list.

**Plan says** (Phase 3) the README needs a new disulfide subsection, a correction to the "No structure preparation" bullet, and an entry in the improvements list. That is the whole list.

**Actually true.** Three further statements become false:

1. **`README.rst:383`**: "On 1GDW, feeding in PROPKA values changes **18 of the 54 features**, moving the isoelectric point from **8.9** to 10.6 and the formal charge from +7 to +8." After the fix the no-PROPKA pI is 10.342, so this becomes "from 10.3 to 10.6". Separately, the "18" is already stale at the shipped default: I measure **19** at `--ionic-strength 0.15` and 18 at `--ionic-strength 0`. **[verified 2026-08-28]** This sentence is the README's headline argument for using PROPKA and the fix cuts its force roughly in half; leaving it is the least honest option available.
2. **`docs/identification_of_propka_dependent_features.md`**, cited from that same README line: 819 AlphaFold structures, 22 of 56 features PROPKA-dependent. AlphaFold models place bonded SG atoms at bonding distance (that is the plan's own premise for geometric detection), so the no-PROPKA arm of that comparison moves and the report's headline numbers are no longer reproducible from the shipped code.
3. **README "Output compatibility"**: "The **values** are identical only with ``--ionic-strength 0``." After this change they are not identical to the original Prodes even at zero ionic strength, for any structure with a disulfide. This is an explicit compatibility contract and it needs a second clause.

**Fix.** Add all three to the Phase-3 README bullet. For (2), either regenerate the report or add a dated note at its head saying which prodes version produced it.

### F4 — SERIOUS. "Records win when present" is under-specified and its failure mode is not the one the plan mitigates.

**Plan says:** "Uses the `SSBOND` records when the file has any, otherwise geometry", and a Phase-4 test pins "a geometric pair the record set does not name, which is *not* called because records win when present." The Risks section covers only the inverse case (a reduced structure that kept its records).

**Two problems.**

* **Ambiguity.** "When the file has any" is evaluated *before* or *after* the plan's own filters (symmetry operator not `1555`; residue absent from the coordinates)? A file with two `SSBOND` records, both crystallographic copies, has records but yields zero pairs. Under one reading detection returns 0 and geometry never runs; under the other geometry finds the real bonds. The plan does not say which. Two implementers, two behaviours, and the difference is silent.
* **The unmitigated failure.** Incomplete `SSBOND` sets are common (truncated files, chain extractions, hand-edited models, output from tools that write only intra-chain records). A file with 1 of 4 records loses 3 real disulfides, and the plan's design guarantees geometry will not rescue them. The plan's warning fires on records that *disagree* with geometry, not on geometry that finds bonds the records omit — so this case is silent by construction. The Risks section names the opposite hazard and not this one.

**Fix.** Either (a) state the fallback is per-file and evaluated on the raw record count, and add a warning when geometry finds a pair the records omit; or (b) take the union — records add pairs, geometry adds pairs, records win only where the two disagree about a *contested* sulfur. (b) is safer and costs one paragraph. Decide, do not list.

### F5 — SERIOUS (scope). The `Atom.charge` rewrite is presented as forced; the terminus half of it is not needed at all, and the side-chain half has an alternative the code already uses.

**Plan says:** "Why `Atom.charge` has to change too ... Setting a bonded cysteine's pKa to `None`, as the issue proposes, is not sufficient on its own, and on its own it crashes."

**Actually true.**

* The **side-chain** half is genuinely required *if* you remove the entry — confirmed, `TypeError`. But the plan documents on the previous page that PROPKA already gets cystines right through a 99.99 sentinel, and the current `Atom.charge` handles that correctly with no change at all. I reproduced the entire fix by giving bonded cysteines `_pka = [{"CYS": 99.99}]` and touching nothing else: formal charge at pH 8.5 = +7.0, pI = 10.342, 1 of 54 reduced and 20 of 105 full features changed at pH 7 — the same numbers the plan predicts. **[verified 2026-08-28]** So the rewrite is a *choice* between a clean model and a magic number, not a necessity. That is a defensible choice; presenting it as forced when the plan itself documented the alternative two pages earlier is the one place the plan overstates.
* The **terminus** half (`pkas[-1]` -> key lookup) is not required by this issue at all. With `side_chain_pka` returning `None`, a bonded C-terminal cysteine skips the side-chain branch and the terminus branch reads `pkas[-1] == {"C-": 2.34}`, which is correct. The plan admits the ordering bug "has never bitten". It is a real bug worth fixing, but it is an independent bug in a `MAJOR` release fixing a different thing — exactly the treatment the plan gives the multi-chain `redo_pkas` defect, which it correctly defers.

**Fix.** Rewrite the section as "why we prefer removing the entry over a sentinel", listing both options and the reason for the choice. Either split the terminus fix into its own commit with its own test and note in the plan that it is a separate defect ridden along, or move it to its own issue. F2 raises the cost of keeping it, since it is the terminus lookup that introduces the new unhandled case.

### F6 — MINOR. "45 titratable groups" is 46 groups over 45 residues.

**Plan says:** "on `1GDW.pka` it is silent, because PROPKA covers all 45 titratable groups."

**Actually true.** 45 residue numbers, **46** groups: 8 ASP + 3 GLU + 1 C- + 1 HIS + 8 CYS + 6 TYR + 5 LYS + 13 ARG + 1 N+. Residue 1 carries both `LYS` and `N+`. Prodes' own defaults also yield 46 groups over 45 residues. **[verified 2026-08-28]** `tests/test_pka.py` asserts `len(pkas) == 45`, which is the count of *keys*, and the plan appears to have read that number as a group count.

Not just a wording slip: the residue-vs-group distinction is the direct cause of F2's blind spot.

**Fix.** "all 46 titratable groups on 45 residues", and specify the Phase-2 warning at group granularity.

### F7 — MINOR. "Both predictors we ship converters for" — there are three.

`src/prodes/io/pka_converter.py` ships `convert_hpp`, `convert_propka` and `convert_pypka`, and `README.rst:393` says "Prodes ships converters for three predictors". **[verified 2026-08-28]** The side-chain-first ordering claim is verifiable for propka (`1GDW.pka`: `LYS 1` at line 222, `N+ 1` at line 240) and for the shipped H++ fixture (`1GDW_hpp_pka.json` key `686` = `[{"LYS": 10.197}, {"C-": 0.0}]`). It is **not** verifiable for pypka — no pypka fixture is shipped, and pypka orders by its own table.

**Fix.** "All three converters we ship happen to emit the side chain first on the files we have; there is no pypka fixture, so that one is inferred." That is also the honest framing for F5's argument, which leans on this claim.

### F8 — MINOR. `prepare_structure` cannot print "alongside the existing `calculating {name}` line".

That `print` lives in `calculate()` at `run.py:526`, two calls after `prepare_structure` and after `construct_surface_grid` — which is the slowest phase. **[verified 2026-08-28]** A line emitted from `prepare_structure` appears seconds *before* `calculating 1GDW`, unlabelled, which is worse than useless in a batch log.

**Fix.** Return the count from `prepare_structure` and print it from `calculate()` on the same line: `calculating 1GDW (4 disulfides)`. It is needed in `calculate()` anyway for `run_metadata`.

### F9 — MINOR. The verification sentence for the crash claim describes an experiment that matches neither stated failure mode.

**Plan says:** "Verified by emptying a residue's `_pka` and calling `charge`: `IndexError` / wrong value, never a silent zero."

**Actually true.** Emptying `_pka` to `[]` does give `IndexError` **[verified 2026-08-28]** — but `[]` is not a state the plan's design produces. `Residue.pkas` only caches when `len(pkas) > 0`, so a bonded non-terminal cysteine yields `None`, and the reachable error is `TypeError`. The `IndexError` is an artefact of the experiment. Small, but it is the one place the plan cites a verification that does not test the thing being claimed.

**Fix.** "Verified by making `pkas` return `None`: `TypeError`. Verified for a terminal cysteine by setting `_pka = [{'C-': 2.34}]`: SG comes out at −1.0 at every pH above 2.34."

### F10 — MINOR. The latent ordering bug is described half-way.

With `[{"N+": 8.35}, {"LYS": 11.34}]` at pH 10, `NZ` goes 1.0 -> 0 (as the plan says) **and** the backbone `N` goes 0 -> 1, because `pkas[-1]` then holds the LYS value. **[verified 2026-08-28]** Net residue charge is accidentally unchanged; the *positions* of the charge move, so `Dipole` and every surface and shell potential feature shift while `Formal charge` does not. Worth a sentence, because "charge 0 instead of 1.0" understates what a reader would look for if this ever fired.

### F11 — MINOR. No rule against two SG atoms in one residue.

`geometric_disulfides` is specified as a pairwise loop over `sulfur_atoms(structure)` with a greedy assignment. Nothing forbids the two atoms of a candidate pair belonging to the same `Residue`. That happens whenever a residue ends up with two SG atoms — a multi-`MODEL` file, since `PDBparser._add_atom` keys residues by `(chain, number)` and appends the second model's atoms to the first model's `Residue` objects rather than starting new ones. The 1.75 A floor makes a same-residue pair unlikely, not impossible, and a residue paired with itself would set `disulfide_partner = self`.

**Fix.** One line in the plan: candidate pairs must come from different `Residue` objects. Costs nothing, removes a class of nonsense.

### F12 — MINOR. `side_chain_pka` and the cache need explicit `None` handling.

`Residue.pkas` returns `None`, not `[]`, for a residue with no entries, and for a bonded non-terminal cysteine it will return `None` on **every** access, because `_pka` is only assigned when the list is non-empty — so the cache the plan calls "computed once, on first access, and cached" is not a cache for exactly the residues this change creates. Functionally harmless (`Atom.charge` is called a handful of times per atom; the EP phases build charge arrays once), but the plan asserts a caching property that does not hold for the new case, and `side_chain_pka` must guard `pkas is None` or it raises on the first bonded cysteine it sees.

**Fix.** State that `side_chain_pka` returns `None` when `pkas` is `None`, and drop or qualify the "computed once and cached" sentence.

### F13 — MINOR. The PROPKA sentinel claim is inferred from one 3.4.0 file and stated as a fact about PROPKA.

`tests/data/1GDW.pka` is a propka3.4.0 run and does show 99.99 for all eight cystines. **[verified 2026-08-28]** But `environment.yml:15` pins `propka=3.5.1`, and propka is **not installed in the `prodes` conda env** (`pip show propka` -> not found), so nothing in this repository or this environment can confirm that 3.5.1 still writes 99.99, or that it does so for every disulfide in every structure. The plan writes it as settled ("That is PROPKA's sentinel for 'this group does not titrate'") and then builds a README recommendation on it.

**Fix.** Attribute it: "propka 3.4.0 writes 99.99 for both partners of every disulfide in our one fixture; we have not run 3.5.1." Or install propka and check — it is one command and the plan's Phase-2 debug-log special case for 99.99 depends on the answer.

### F14 — MINOR. `run_metadata`'s signature change is an unremarked API break.

`run_metadata(pdb_file, settings, n_points, ep)` gains a parameter. The plan notes that `biochai` is an outside caller of `convert_propka`, so the package is treated as having external consumers; it does not check whether anything imports `run_metadata`. The `MAJOR` bump covers it, but it belongs in the risk list next to "this changes shipped numbers". Adding the key at the top level (not inside `settings`) does *not* break `tests/test_output_bundle.py::test_the_run_record_says_what_was_run`, which asserts equality on `record["settings"]` only. **[verified 2026-08-28]**

### F15 — MINOR (honesty, in the plan's favour). The impact at the default pH is much smaller than the plan's framing implies.

The plan's headline table is at pH 8.5 and 9.0. At the shipped default of pH 7, on 1GDW, the fix changes **1 of the 54** reduced features (`Isoelectric point`) and **20 of the 105** full features (`Isoelectric point`, `Average charge`, and the 18 `SurfEp*Average` columns). With `--pka`, it changes **0**. **[verified 2026-08-28]** The plan's Phase-3 justification is correctly qualified ("at any pH where a thiol would have titrated"), so this is not a wrong claim — but a reader deciding whether to refit a model will want the pH-7 numbers, and they are not in the plan. The `MAJOR` bump remains justified: `Isoelectric point` is a reported quantity, it moves by 1.4 pH units, and it moves at *every* pH because the pI search sweeps 0-14.

### F16 — MINOR. Insertion codes are not mentioned anywhere.

`SSBOND` carries `iCode` in columns 22 and 36; the plan's `read_ssbond_line` returns `(chain1, number1, chain2, number2, sym1, sym2)` and drops them. That is *consistent* with `PDBparser._read_line`, which reads the residue number from `line[22:26]` and drops the iCode too — so residues 100 and 100A merge into one `Residue` and an `SSBOND` naming 100A resolves to that merged object. Pre-existing, no worse than today, but the plan flags the altLoc quirk and not this one, which is the closer analogue and bites `record_disulfides` directly.

**Fix.** One sentence next to the altLoc paragraph.

---

## Issue acceptance criteria

| # | Criterion | Met? | Notes |
|---|---|---|---|
| 1 | A test structure with known disulfides reports the correct number detected | **Yes** | Phase 1 acceptance and Phase 4 test: 1GDW -> 4 pairs, exact partner numbers, symmetry. Recomputed: 4 pairs (6-128, 30-116, 65-81, 77-95). Phase 3 also exposes the count in `prodes_run.json`. |
| 2 | `Formal charge` and `Isoelectric point` for a cystine-rich structure change in the expected direction and magnitude at pH 8.5 | **Yes** | Table verified cell for cell; Phase 4 pins pH 7.0/8.5/9.0. One gap: the *feature* `Formal charge` is computed at the run pH, and no test runs the pipeline at pH 8.5 — the assertions are on `Structure.charge` directly. Adequate, but say so. |
| 3 | A free-thiol Cys in the same structure is still titrated normally | **Partly** | No shipped fixture mixes bonded and free cysteines — the plan says so plainly, which is to its credit — and covers it with a hand-written fixture. So the criterion is met on a synthetic structure, not on "the same structure" as criterion 2. 1GPB (9 free thiols, 0 pairs) covers the free-thiol half separately. Acceptable; the deviation should be stated rather than left for a reader to notice. |
| 4 | Behaviour is unchanged for structures with no disulfides | **Yes, but under-tested as specified** | ARH96693 bit-identical is the right guard and 0 pairs is confirmed. But the `Atom.charge` rewrite (F5) also touches every *terminus* charge in the package, for every structure, disulfide or not, and the plan's only regression guard is a structure with no pKa file. Add a no-disulfide structure **with** `--pka` to the guard — 1GDW with `1GDW.pka` is the obvious candidate and, per my measurement, must come out bit-identical (0 of 54 features change). |
| 5 | The pKa-file interaction from point 4 is covered by a test | **Yes as scoped, no as needed** | Phase 2 gives two tests (bonded Cys with a real pKa in the file; a titratable residue removed from the file). Both are right. But per F2 and F6 the warning is specified at residue granularity, so the case that actually loses data today — a file naming a residue's side chain but not its terminus — is neither warned about nor tested. Extend the criterion's coverage before calling it met. |

Issue point 3 ("set the pKa to `None` ... Prefer this over special-casing in `charge()`") is **deviated from**: the plan special-cases in `Atom.charge` via `side_chain_pka`. The deviation is argued and is correct on the merits — but see F5 for why the argument as written is stronger than the facts support.

Issue point 5 ("expose the count") is **met** (stdout line plus `prodes_run.json`), and the plan's reasoning for keeping it out of the feature columns — `features_reduced.yaml` / `features_full_only.yaml` are asserted column by column in `tests/test_feature_dictionary.py` — is correct and is the right call.

---

## What would break for an existing user

Handled by the plan: changed feature values (`MAJOR` bump, release-note callout), the need to refit models.

Not handled: (a) the three stale documentation claims of F3, one of which is the package's own headline argument for the recommended workflow; (b) `run_metadata`'s signature (F14); (c) a user whose pKa file omits terminus groups now hits the new, unspecified code path of F2. Nothing here is expensive to fix, but (a) and (c) should not ship silently.
