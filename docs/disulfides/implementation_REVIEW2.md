# Implementation review 2: disulfide-bonded cysteines (issue #6)

**Reviewer lens:** maintainer and scientist. **Date:** 2026-08-28. **Branch:** `6-disulfide-cysteines`, uncommitted working tree.
**Environment:** `conda activate prodes`, run from `/mnt/shared/prodes` (a copy at a scratch path; the installed editable package resolves to `/mnt/shared/prodes/src`, so every number below was produced by the code under review).

## Verdict

The code is right and the science is defensible. The fix sits at the pKa level exactly as the issue asked, so it propagates to every consumer without feature-code changes; the positional-to-keyed pKa lookup in `Atom.charge` repairs a real latent bug and is properly tested; `redo_pkas`'s replacement-to-merge change fixes a silent data loss that predates this issue. The 2.5 Å cutoff is not a guess — I confirmed it independently against both tools the README names. Every headline number in the plan reproduces exactly, 297 tests pass, and pre-commit is clean including the untracked files. What is not right is the documentation: two sentences in `README.rst` are wrong as written, one of them in the release note that tells an existing user whether their pipeline is affected, and it understates the change at the default pH 7. **Ship the code. Do not ship `README.rst` as written** — findings 1 and 2 are one-line edits and must land first. Nothing else blocks.

---

## Findings

### 1. Serious — the release note understates the blast radius, and it is the sentence users will read

`README.rst:688`: "corrected the charge-derived features of every structure with a disulfide **above about pH 8**".

At the **default** pH of 7, 20 of the 105 full features move on 1GDW **[verified 2026-08-28]**. I ran `calculate` twice on `tests/data/1GDW.pdb.zip` with `full_features=True`, once as shipped and once with `prodes.io.parser.assign_disulfides` stubbed to a no-op:

```
Isoelectric point:      8.899  ->  10.342
Average charge:         6.7350 ->   7.0924
SurfEpSumAverage:    1574.219  -> 1645.483
NSurfPosEpAverage:       5168  ->    5235
NSurfNegEpAverage:       1866  ->    1803
... 15 more SurfEp*Average columns
```

`Formal charge` is unchanged at pH 7, which is what "above about pH 8" is true of. But `Isoelectric point` is not a pH-dependent quantity at all and moves by 1.443 units at every pH, and the whole `*Average` electrostatic block moves at pH 7. A user reading only the improvements bullet concludes a pH-7 pipeline is unaffected. It is not.

The same framing weakens `README.rst:448`, "Below about pH 8 the error is small" — true of the formal charge, false of the isoelectric point, which is the feature that paragraph itself names two clauses earlier.

**Fix.** `README.rst:688`: "…which changes the charge-derived features of every structure with a disulfide, at every pH. Above about pH 8 the formal charge changes too; below it the isoelectric point and the fractional-charge features still move." `README.rst:448`: qualify "the error is small" as "the *formal*-charge error is small".

### 2. Serious — "the two cannot disagree" is false

`README.rst:453`: "That is the same cutoff PROPKA and PDB2PQR use, so the two cannot disagree about which cysteines are bridged."

The cutoff claim is correct **[verified 2026-08-28]**: `pip download propka==3.5.1 --no-binary :all:` → `propka/bonds.py:24  DISULFIDE_DISTANCE = 2.5`, applied in `_find_bonds_for_atoms` and setting `atom.cysteine_bridge`, with `propka/group.py:503-504` writing `99.99` for a bridged atom; `pip download pdb2pqr` → `pdb2pqr-3.7.1/pdb2pqr/config.py:162  BONDED_SS_LIMIT = 2.5`, used at `biomolecule.py:697`.

The *conclusion* does not follow. Same threshold, different rule, and there are three ways they diverge:

- **The record route.** PROPKA never reads `SSBOND`. `record_disulfides` honours a record whose sulfurs are 15 Å apart — the branch the implementation's own `test_a_record_is_honoured_even_when_the_sulfurs_are_far_apart` pins. Prodes calls that cysteine bridged; PROPKA gives it a real predicted pKa.
- **Greedy one-partner assignment.** PROPKA's `_find_bonds_for_atoms` runs over every pair and sets `cysteine_bridge` on both members of *every* S-S pair within 2.5 Å. Three sulfurs mutually inside the window: PROPKA bridges all three, prodes bridges only the shortest pair and titrates the third at 8.33.
- **The 1.6 Å floor.** Prodes has one, PROPKA has none. Duplicated coincident sulfurs are bridged by PROPKA and not by prodes — again pinned by the implementation's own `test_coincident_sulfurs_are_reported_rather_than_bonded`.

**Fix.** End the sentence at "…the same cutoff PROPKA and PDB2PQR use." If you want the reassurance, say what is actually true: "so the geometric route agrees with PROPKA on which pairs are close enough. It can still differ where an `SSBOND` record overrules the distance, or where three cysteines cluster inside the cutoff and Prodes takes only the shortest pair."

### 3. Minor — the `docs/identification_of_propka_dependent_features.md` note claims the wrong kind of bound

The new note says "the dependence reported below is an upper bound". That document's headline metric is a **count** ("PROPKA-dependent: **22**"), and the count does not shrink. Measured on 1GDW, reduced feature set **[verified 2026-08-28]**:

| | PROPKA changes | pI without pKa file → with |
|---|---|---|
| detection disabled (old behaviour) | 19 of 54 | 8.899 → 10.565 |
| detection on (this branch) | 19 of 54 | 10.342 → 10.565 |

Same 19 features respond either way. What shrinks is the *magnitude* — the pI gap falls from 1.666 to 0.223 units, which is what would raise the R² figures the document tabulates.

**Fix.** "…so the per-feature R² values below understate the agreement, and the magnitudes of the differences are an upper bound. The set of PROPKA-dependent features is unchanged: on 1GDW it is 19 of 54 either way."

### 4. Minor — stale number left in a file this change edited

`tests/test_pka.py:11` still says "On 1GDW a real PROPKA prediction moves **18** of the 54 reduced features". The README was corrected to 19 in the same commit; 19 is the measured value, and it was 19 before this change too (table in finding 3), so the old 18 was already wrong. **Fix:** change 18 to 19 in the module docstring.

### 5. Minor — the bundle manifest no longer describes the bundle

`README.rst:457` tells the user the disulfide count is "recorded in ``prodes_run.json``". Both places that list what `prodes_run.json` contains still say "version, settings and time of the run": `src/prodes/output.py:88` (the `README.txt` written into every bundle) and `README.rst:189`. **Fix:** "version, settings, disulfide count and time of the run" in both.

### 6. Minor — `assign_disulfides` is not idempotent, and the count it reports can contradict the structure

It never clears prior assignments. On 1AO6, parsed (34 bonds), then `assign_disulfides(structure, ())` called again with `geometric_disulfides` stubbed to return nothing **[verified 2026-08-28]**:

```
second pass bonds: 0   structure.disulfides: 0
residues still carrying a partner after a pass that found none: 68
their pkas: None
```

`structure.disulfides` reports 0 while 68 cysteines are still silently non-titratable — and `structure.disulfides` is precisely the number `run.py:529` prints and the README (line 457) tells the user to sanity-check. Unreachable through the shipped pipeline (the parser calls it once), but the function's docstring explicitly advertises being called on an already-parsed structure. **Fix:** three lines at the top of `assign_disulfides` — for every `CYS` residue, `residue.disulfide_partner = None; residue._pka = None`.

### 7. Minor — the runtime version still says 5.0.0

`pyproject.toml` is 6.0.0, but `prodes_version()` reads `importlib.metadata.version("prodes")`, which comes from the installed distribution. In this environment **[verified 2026-08-28]**: `pip show prodes` → `Version: 5.0.0`, and `python -c "from prodes.output import prodes_version; print(prodes_version())"` → `5.0.0`. Every bundle produced from this branch right now stamps `"prodes_version": "5.0.0"` next to a README that says 6.0. `tests/test_dependency_versions.py` compares files against each other, not against the environment, so nothing catches it. Plan review 2 raised this; the refined plan dropped it. **Fix:** `pip install -e .` before generating any reference output or tagging, and consider a test asserting `prodes_version()` matches `pyproject.toml`.

### 8. Minor — `redo_pkas` silently changed which duplicate wins, and the new fixture makes it reachable

`convert_propka` keys by residue number only and *appends* on collision (`pka_converter.py:100-103`), so a two-chain PROPKA run gives `pkas[34] = [{"CYS": a}, {"CYS": b}]`. Old code stored the list and `Atom.charge` read `pkas[0]` → chain A's value for both chains. New code applies both entries to the same group in order **[verified 2026-08-28]**:

```
before: [{'ASP': 3.86}]
redo_pkas({18: [{"ASP": 3.11}, {"ASP": 7.77}]})
after : [{'ASP': 7.77}]   side_chain_pka = 7.77
```

First-wins became last-wins. Both are wrong, and neither warns — the duplicate counts as "applied", so the new omission warning stays silent. This is adjacent bug §10.1 of the plan, correctly out of scope, but 1AO6 is now a committed two-chain fixture and the same commit ships chain-aware disulfide resolution next to chain-blind pKa resolution. **Fix:** open the issue before merging, and say in the pKa section of the README that predicted values are applied by residue number without regard to chain.

### 9. Minor — the "reported before it is discarded" promise is invisible

`structure.py:256` logs at `INFO` when a pKa file gives a bridged cysteine a real (non-99.99) value. Prodes ships no logging configuration, so the default `logging.lastResort` handler emits `WARNING` and above only — I confirmed a `logger.warning` from the disulfide module does reach stderr with no configuration **[verified 2026-08-28]**, which means the `INFO` line does not. The `redo_pkas` docstring says "a real value from another predictor is reported before it is discarded". In practice it is discarded silently. This is the one place where Prodes overrules a third-party predictor's number. **Fix:** promote it to `WARNING`, or change the docstring to say "logged at INFO, which the default handler does not show".

### 10. Minor — the limitation list is honest but not complete

Two reachable cases are missing from `README.rst:463-469`:

- **Alternate location indicators.** `parser._read_line` uses `name = line[12:17].strip()`, which swallows the altLoc column, so an `SG` with altLoc `A` parses under the atom name `"SG  A"`. `cysteine_sulfurs` matches `atom.name == "SG"` and never sees it. Verified on a synthetic two-cysteine file with SG atoms 2.05 Å apart and altLoc `A` **[verified 2026-08-28]**: `disulfides = 0`, atom names `['C   A', 'CA  A', 'CB  A', 'N   A', 'O   A', 'SG  A']`, and the first cysteine keeps `[{'CYS': 8.33}, {'N+': 9.69}]`. The charge is unaffected (a mis-named `SG` is not in `charged_atoms`, so it was already neutral), but the **reported count** is wrong, and the count is the diagnostic the README asks users to trust. This is the plan's adjacent bug §10.4; the README does not mention it. Note also that `cysteine_sulfurs`'s "more than one SG" warning cannot fire in this case, because the two names differ.
- **Hydrophobicity is untouched.** A cystine still carries the free-thiol `CYS` value from `hydrophobicity.json`, so `CYSSurfFrac` and the whole MHP block are unchanged by this fix. The issue called that out of scope; the README should say so, since the section reads as "here is everything Prodes now knows about a cystine".

### 11. Minor — a second behaviour change ships in 6.0.0 with no user-facing note

`redo_pkas` went from replace to merge, gained two drop rules and an omission warning. That changes results for any file that names a residue's side chain but not its terminus, and it is a public method with an external caller (`biochai`, per `tests/test_pka.py`). The improvements list has one bullet, about disulfides. **Fix:** a second bullet — "pKa files merged per titratable group — a predicted side-chain value no longer deletes the residue's terminal pKa, entries for groups a residue does not have are dropped with a warning, and groups the file omits are reported."

### 12. Minor — the one doctored committed fixture has no provenance

`tests/data/1AO6.pdb.zip` (187 KB) holds `ATOM: 9198, SSBOND: 34, TER: 2, CRYST1: 1, END: 1` **[verified 2026-08-28]** — waters and every header record stripped. `tests/test_disulfides.py`'s docstring says the synthetic cases are written inline "so that the doctoring a test depends on is visible in the test itself", and then depends on a doctored committed file whose doctoring is recorded only in §7.1 of the plan. A maintainer diffing against RCSB 1AO6 will not find that. **Fix:** four lines in the `test_disulfides.py` docstring, or a `tests/data/README.md`: PDB ID, source, download date, and the exact record filter applied. (1GDW and 1GPB have the same gap; this is the moment to close it.)

*The fixture itself is worth its weight.* It is the only structure that gives a free thiol and a bonded cysteine in one file with `SSBOND` records that agree with the geometry, it is the repository's first multi-chain fixture, and it costs 0.58 s to parse **[verified 2026-08-28]**; the whole of `test_disulfides.py` plus `test_residue.py` runs in 3.5 s. Keep it.

### 13. Minor — grammar and unsupported claims

- `run.py:529` prints "1 disulfide bonds detected" for a single bond.
- `README.rst:452`: a stretched record "usually means a strained bond and occasionally means a reduced structure carrying stale records" — an unsupported frequency claim, and the code comment in `report_record_distance` leans the other way. Drop "usually / occasionally".

### 14. Minor — test gaps

The tests pin behaviour rather than shape, and the fixture choices are right. What is missing:

- No end-to-end assertion that a real disulfide-containing structure reports a **nonzero** count in a bundle. `test_the_run_record_says_what_was_run` asserts 0 on ARH96693 and `test_the_run_record_carries_the_disulfide_count` passes a fabricated 17. Nothing connects `calculate` on a structure with bonds to the number in `prodes_run.json`.
- `Structure.titratable_groups` is new public API with no direct test.
- The "two SSBOND records claim the same cysteine, keeping the first" branch (`disulfides.py:186`) is unreached by any test.
- `read_ssbond_line` with a blank chain ID (single-chain files with an empty column 16) is untested, and it is the case where a record must still match a residue the parser stored under chain name `""`.
- The exact boundaries are untested: 2.50 Å is inclusive and 1.60 Å is inclusive in the code; the tests probe 2.45/2.55 and 0.5.
- `read_ssbond_line` lives in `io/parser.py` but is tested only from `tests/test_disulfides.py`. Defensible, but a reader of `tests/test_parser.py` will not know it is covered.

---

## Documentation: every touched claim, recomputed

All measured on this branch **[verified 2026-08-28]**.

| Claim | Location | Recomputed | Verdict |
|---|---|---|---|
| PROPKA changes "19 of the 54 features" on 1GDW | `README.rst:383` | 19 of 54 numeric columns | **correct** |
| "isoelectric point from 10.34 to 10.57" | `README.rst:383` | 10.342 → 10.565 | **correct** |
| "formal charge from +7 to +8" | `README.rst:383` | 7.0 → 8.0 | **correct** |
| formal charge at pH 8.5 "-1 when it is +7" | `README.rst:448` | old -1.0, new +7.0 | **correct** |
| "isoelectric point at 8.9 when it is 10.3" | `README.rst:448` | 8.899 → 10.342 | **correct** |
| "Below about pH 8 the error is small" | `README.rst:448` | true of formal charge; pI is off by 1.443 at every pH, `Average charge` by 0.357 e at pH 7 | **misleading** (finding 1) |
| cutoff "2.5 Å … the same cutoff PROPKA and PDB2PQR use" | `README.rst:453` | propka 3.5.1 `bonds.py:24 = 2.5`; pdb2pqr 3.7.1 `config.py:162 = 2.5` | **correct** |
| "so the two cannot disagree about which cysteines are bridged" | `README.rst:453` | three named divergences | **wrong** (finding 2) |
| PROPKA "reports a bridged cysteine as 99.99 … Prodes has always passed that through" | `README.rst:461` | `propka/group.py:503-504`; all 8 bridged cysteines are 99.99 in `tests/data/1GDW.pka`; `neg_charge(99.99, ph) ≈ 0` | **correct** |
| "recorded in ``prodes_run.json``" | `README.rst:457` | present as a top-level key; but both manifests still omit it | **correct but incomplete** (finding 5) |
| improvements: "above about pH 8" | `README.rst:688` | 20 of 105 features move at pH 7 | **wrong** (finding 1) |
| output compatibility: "identical only with ``--ionic-strength 0``, and then only for a structure with no disulfide bonds" | `README.rst:697` | ARH96693 has 0 disulfides and `tests/test_sasa.py` passes untouched | **correct** |
| "Almost no structure preparation … it *recognises* existing disulfide bonds" | `README.rst:675` | accurate | **correct** |
| "the dependence reported below is an upper bound" | `docs/identification_of_propka_dependent_features.md:7` | count unchanged (19 of 54 before and after); magnitude shrinks | **wrong as phrased** (finding 3) |
| "moves 18 of the 54 reduced features" | `tests/test_pka.py:11` | 19 | **wrong** (finding 4) |
| `disulfides.py` docstring: "moved the formal charge at pH 8.5 from +7 to -1 and the isoelectric point by 1.4 units" | module docstring | 8.899 vs 10.342 = 1.443 | **correct** |
| `atom.py` docstring: a residue whose terminus entry comes first "hands the terminal pKa to the side chain" | `charge()` docstring | a bridged terminal CYS has `pkas == [{"C-": 2.34}]`; old `pkas[0]` gave the SG 2.34 | **correct**, and tested |
| `MAX_SG_SG_DISTANCE_ANGSTROM` comment: 6VXX at 2.361 Å, metallothionein-2 at 2.050 Å | `disulfides.py` | not re-fetched; both were verified in plan review 1 and the refined plan | **not independently re-verified** |
| plan §7.1: 1AO6 has 34 bonds, longest 2.048 Å, nearest non-bonded 5.951 Å, free thiols A34/B34 | fixture | 34 bonds, 1.979–2.048 Å, nearest non-bonded 5.951 Å, free `[('A', 34), ('B', 34)]` | **correct** |

---

## What is solid

- **The fix is in the right place.** A bridged cysteine loses its side-chain pKa entry in `Residue.pkas`, so `charged_atoms`, `Atom.charge`, `Structure.charge`, `dipole`, `isoelectric_point`, the surface grid and the shell phase all follow with no feature-code change. That is what the issue asked for and it is why the change is 451 lines rather than a thousand.
- **`Atom.charge` positional → keyed is a genuine repair, not collateral.** The old `pkas[0]` / `pkas[-1]` reads were a latent bug independent of disulfides, and `test_a_bonded_cysteine_at_a_chain_terminus_keeps_its_terminal_charge` pins the exact failure.
- **`redo_pkas` merge.** Fixes real silent data loss (a side-chain prediction deleting the residue's terminal pKa) and finally makes the README's long-standing promise true.
- **The cutoff is sourced, not invented,** and both sources check out. The reasoning in the `MAX_SG_SG_DISTANCE_ANGSTROM` comment — why not tighter, why not wider — is the best kind of constant documentation in this repo.
- **Per-cysteine precedence** is the right rule and the 1ERT justification is real; the synthetic test for it is clearer than a committed fixture would have been.
- **All plan numbers reproduce exactly**: 1GDW 4 bonds `{(6,128),(30,116),(65,81),(77,95)}`, formal/average charge table at 7.0/7.4/8.5/9.0, pI 8.899 → 10.342; 1GPB 9 cysteines and 0 bonds; 1AO6 70 cysteines, 34 bonds, free A34/B34.
- **297 tests pass in 117 s**; `pre-commit run --all-files` passes with the untracked files staged (black, ruff, mypy).
- **The regression guard held.** `tests/test_sasa.py` against `ARH96693_prodes_orig_output.csv` is untouched and passes, which is the evidence that nothing leaked into disulfide-free structures.
- **Fits the repo.** No leading-underscore helpers in the new module, every function has a docstring in the house voice, constants carry a why-comment, `print` for pipeline progress and `logging` for anomalies matches existing practice, and warnings do reach a bare CLI user through `logging.lastResort` (confirmed).
- **6.0.0 is the correct bump.** The reported quantity changes and the version is stamped into every bundle.

## What breaks for an existing user

1. **Feature values** for any structure with a disulfide, at every pH — not only above 8. Communicated in "Output compatibility" and by the major bump; **understated** in the improvements bullet (finding 1).
2. **Any model fitted on Prodes ≤5 features** from disulfide-containing structures needs refitting. Stated in the plan's risk section; not in the README. Worth one sentence.
3. **`Structure.redo_pkas` semantics** (replace → merge, entries dropped, omission warning). Not communicated anywhere user-facing (finding 11), and the duplicate-key winner silently flipped (finding 8).
4. **New stderr warnings** on any pKa file that does not cover every titratable group. Desirable, but new and unannounced.
5. **New stdout line** from `calculate`. Trivial unless someone parses stdout.
6. `run_metadata` gained a keyword argument with a default and `prodes_run.json` gained a key — both additive.
