# Refined plan: stop charging disulfide-bonded cysteines

**Status:** plan, ready to implement. Supersedes `plan_disulfide_cysteines.md`, which must not be edited: it and the three reviews are the audit trail.
**Branch:** `6-disulfide-cysteines`, cut from `main` at `b5e8409`.
**Issue:** [#6](https://github.com/datacatalysis/prodes/issues/6), "Cystines are charged as free thiols, distorting charge features above pH 8".
**Date:** 2026-08-28.
**Provenance:** original plan, then three independent reviews (`_REVIEW1` structural chemistry, `_REVIEW2` software correctness, `_REVIEW3` adversarial). All three returned *do not merge as written*. No headline number in the original plan was wrong; every defect was in the criterion, the precedence rule, or the blast radius. This revision changes the detection cutoff, changes the precedence rule from per-file to per-cysteine, replaces `redo_pkas`'s list-replacement with a per-group merge, and adds a real fixture that the original plan said did not exist.

Everything marked **[verified 2026-08-28]** was recomputed here, after the reviews, not carried over.

---

## 1. Goal

A cysteine whose `SG` is bonded to another cysteine's `SG` has no thiol proton and cannot titrate. Prodes gives every `CYS` a pKa of 8.33 regardless of structure, so both halves of every disulfide titrate as free thiols. Above about pH 8 this puts a large spurious negative charge on the protein.

After this change, a cysteine Prodes judges to be in a disulfide carries no side-chain pKa and therefore no side-chain charge, at any pH, by any route.

---

## 2. Background

### 2.1 The object graph and the two charge paths

`PDBparser.parse` builds `Structure -> Chain -> Residue -> Atom`. Only lines whose record name equals the parser's `identifier` are read, and that defaults to `"ATOM"`; every other record type, `SSBOND` included, is discarded unread. **[verified 2026-08-28]**

Charge is decided in two places that consult different sources:

| Caller | What it reads | Keyed by |
|---|---|---|
| `Residue.pkas` | `data.residue_data(self.name)["pka"]` plus terminus defaults | residue *type*, from `general_data.json` |
| `Residue.charged_atoms` | `self.pkas`, dispatching on each entry's key | the entry's key |
| `Atom.charge` | `residue_data(self.residue_name)["potential_charge"]` and `["charged_atoms"]` decide *whether* the atom is charged; then `self.residue.pkas[0]` (side chain) or `self.residue.pkas[-1]` (terminus) for the *value* | residue type for the decision, list **position** for the value |

`Residue.pkas` returns a list of single-entry dicts, `[{"LYS": 10.5}, {"N+": 9.69}]`, or `None` for a residue that is neither titratable nor terminal. The side-chain entry is appended before the terminus entry, and that ordering is the only reason `pkas[0]` and `pkas[-1]` currently pick the right entries.

`Structure.charge` sums `Residue.charge`, which sums `Atom.charge` over heavy atoms. `Structure.dipole`, `Structure.isoelectric_point`, `run.charged_atom_arrays`, the surface-grid phase and the shell phase all go through `Atom.charge`. A fix at the pKa level therefore reaches every electrostatic feature without touching feature code.

### 2.2 Why `Atom.charge` must change, and why the alternative was rejected

Setting a bonded cysteine's pKa to `None`, as the issue asks, does not work on its own: it crashes.

`Atom.charge` asks the residue *type* table whether `SG` is a charged atom of a `CYS`. That answer is always yes, whatever this residue's pKa list holds. It then indexes `pkas[0]`. With the side-chain entry removed, `pkas` is either `None` (`TypeError: 'NoneType' object is not subscriptable`) or, for a terminal cysteine, `[{"C-": 2.34}]`, in which case `pkas[0]` hands the C-terminal pKa to the `SG` atom and gives it a full negative charge above pH 2.3. **[verified 2026-08-28]**

Review 3 argued the rewrite is avoidable: give a bonded cysteine a sentinel pKa of 99.99, the mechanism the code already relies on when PROPKA supplies one, and nothing else needs to change. That is true numerically and it is rejected on design grounds:

* The issue asks for the pKa to be `None`, "matching how `ALA` and other non-titratable residues are already handled", precisely so that the fix propagates rather than being a special case.
* A sentinel leaves `Residue.pkas` reporting a titratable group that does not chemically exist, and leaves `Residue.charged_atoms` listing the `SG` as a charged atom carrying `-0.0`. Every future consumer of those two APIs inherits the lie.
* 99.99 is a magic number belonging to a third-party file format. Adopting it as internal state means prodes cannot later distinguish "PROPKA said this does not titrate" from "prodes decided this does not titrate".

The related change to the *terminus* lookup is kept for a different reason, given in 2.3.

### 2.3 `redo_pkas` replaces rather than merges, and that loses data today

`Structure.redo_pkas` assigns `residue._pka = pka_dict[number]`, replacing the whole list. A file that names a residue's side chain but not its terminus therefore deletes the terminus pKa. On 1GDW, applying `{1: [{"LYS": 11.34}]}` leaves residue 1 with `pkas == [{"LYS": 11.34}]`, and the backbone `N` atom, reading `pkas[-1]`, then takes its charge from the lysine value. **[verified 2026-08-28]** The N-terminal pKa of 9.69 is gone and nothing reports it.

This is why the terminus half of the lookup is in scope: keying by name without fixing the replacement would turn a silent wrong number into a crash (review 2's blocking finding), and fixing the replacement without keying by name would leave `pkas[-1]` reading whichever entry happens to be last. The two changes are one change.

The README already promises the merge semantics that the code does not implement: "Residues that appear in the file get the predicted value; every other residue keeps its default." This makes that true one level down, per titratable group rather than per residue.

### 2.4 What PROPKA actually does, verified against the pinned version

The issue says PROPKA "generally omits bonded cysteines from its titratable list". It does not.

PROPKA 3.5.1, the version pinned in `environment.yml`, was installed in a throwaway virtualenv and run on both shipped structures. **[verified 2026-08-28]**

* On **1GDW**, all eight bonded cysteines appear in the summary with a pKa of `99.99`, reproducing the committed propka3.4.0 file exactly.
* On **1GPB**, whose nine cysteines are all free thiols, every one gets a real predicted pKa between 8.32 and 13.79, and none is 99.99.
* In the propka source, `99.99` is assigned in exactly one place, `propka/group.py`, guarded by the atom's `cysteine_bridge` flag, and the disulfide criterion is `DISULFIDE_DISTANCE = 2.5` in `propka/bonds.py`. PROPKA never reads `SSBOND` records; it detects geometrically. **[verified 2026-08-28]**

So 99.99 is PROPKA's marker for "this cysteine is bridged", nothing else, and PROPKA's own detection agrees with the geometric result on both structures: 8 of 8 bridged on 1GDW, 0 of 9 on 1GPB. That is a free external cross-check of our detector.

Three consequences:

1. **This is a no-pKa-file bug.** With `--pka` from PROPKA the charge features are already close to right, by way of a magic number. Without it they are badly wrong. That inverts the issue's framing and belongs in the README, since "use PROPKA" is the recommended workflow.
2. Issue point 4, warning about residues a pKa file omits, is still worth doing but is not the cystine safety net the issue thought. On 1GDW, PROPKA omits zero titratable groups. **[verified 2026-08-28]**
3. A detected disulfide and a 99.99 from the file agree, so the two mechanisms never fight.

### 2.5 How bad it is now

`tests/data/1GDW.pdb.zip` is mutant human lysozyme R21G: 130 residues, 8 cysteines, all 8 bridged, no free thiol. Its `SG`-`SG` distances are 2.022, 2.030, 2.031 and 2.050 A, and the nearest non-bonded pair is at 6.541 A. **[verified 2026-08-28]**

Net charge as shipped, against the same structure with the eight cysteines silenced **[verified 2026-08-28]**:

| pH | Formal now | Formal fixed | Average now | Average fixed |
|---|---|---|---|---|
| 7.0 | 7.000 | 7.000 | 6.735 | 7.092 |
| 7.4 | 7.000 | 7.000 | 6.178 | 7.019 |
| 8.5 | **-1.000** | **7.000** | 1.936 | 6.709 |
| 9.0 | **-1.000** | **7.000** | -0.458 | 6.133 |

Isoelectric point moves from **8.899 to 10.342**. At pH 8.5 the reported formal charge has the wrong sign on a protein whose real formal charge is +7.

---

## 3. The detection criterion

### 3.1 Cutoff: one upper bound at 2.5 A

The original plan proposed 2.05 +/- 0.3 A, an accepted window of 1.75 to 2.35 A. **That upper bound is too tight and would miss real bonds.** Review 1 measured 142 annotated disulfides across 12 entries and found the bonded range runs to 2.362 A. Re-verified here on 6VXX, the SARS-CoV-2 spike trimer, a 2.8 A cryo-EM structure: `Cys391`-`Cys525` is modelled at **2.361, 2.361 and 2.362 A** in the three protomers, while every other bonded pair in the entry sits between 2.018 and 2.058 A and the nearest non-bonded pair is above 3.6 A. **[verified 2026-08-28]** A 2.35 A ceiling silently drops three real bonds from that one structure.

The cutoff adopted is **2.5 A**, which is what PROPKA 3.5.1 and PDB2PQR both use, verified in the PROPKA source above. Using the same number as the tool the README tells users to run means the two cannot disagree about which cysteines are bridged.

**It must not go higher.** The AlphaFold model of human metallothionein-2 (`AF-P02795-F1`), a protein with twenty metal-binding cysteines and zero disulfides, places `Cys15`-`Cys29` at **2.050 A**, dead centre of the real bond distribution, with the rest of the metal cluster at 3.32 A and above. **[verified 2026-08-28]** No distance cutoff can exclude that false positive, and raising the cutoff toward the cluster distances would convert a single false positive into a dozen. 2.5 A keeps a 0.8 A margin below the cluster cloud.

A lower bound of **1.6 A** is also applied. No real S-S bond is that short; a pair below it means coincident or duplicated atoms, which happens in modelling output. Such a pair is logged and not paired.

**Rejected: a chi-3 dihedral filter.** Biotite screens on `CB-SG-SG-CB` at 90 +/- 10 degrees. Review 1 measured that dihedral across the same 142 real bonds and found it spans 37 to 174 degrees with 54 per cent outside 80 to 100 degrees. The filter would reject more real bonds than it saves.

### 3.2 Precedence: per cysteine, not per file

The original plan said records win when the file has any, and geometry runs only when it has none. That is wrong at file granularity, and 1ERT shows why. **[verified 2026-08-28]**

* **1ERT** (reduced human thioredoxin) has exactly one `SSBOND` record: `CYS A 73` to `CYS A 73`, symmetry operators `2655` and `1555`. It is a self-pair across a crystallographic symmetry axis, so after filtering it names nothing in this file's coordinates. Under file-granularity precedence, the presence of that one record suppresses geometry entirely.
* **1ERU** (oxidised human thioredoxin) has the same symmetry record plus the real `CYS A 32` to `CYS A 35` at 2.024 A. Its true answer is one disulfide; 1ERT's true answer is zero, its closest `SG`-`SG` pair being 3.921 A.

The rule adopted: an `SSBOND` record is authoritative **for the two cysteines it names**. Records are applied first; geometry then runs over the cysteines that no surviving record claimed. A file whose record set is incomplete keeps its remaining bonds, and a file whose only records are symmetry pairs falls through to geometry entirely.

Record filtering, in order:

1. Either symmetry field present and not `1555` — skipped, the partner is a crystallographic copy not in these coordinates. A **blank or truncated** symmetry field is read as `1555`, because hand-written and trimmed PDB files routinely stop before column 72.
2. Both sides naming the same chain and residue number — skipped, necessarily a symmetry self-pair.
3. Either side naming a residue not present in the coordinates, or not a `CYS` — skipped with a warning.
4. A surviving record whose measured `SG`-`SG` distance exceeds the geometric cutoff — **honoured**, with a warning naming the pair and the distance. The depositor's chemistry wins over our window, but a stale record on a reduced or mutated model is exactly the case where trusting it silently would be wrong.

---

## 4. Phase 1: detect the bonds, and remove the side-chain pKa

### 4.1 New module, `src/prodes/calculations/disulfides.py`

Pure functions over parsed objects, no I/O, no leading-underscore helpers.

```
MAX_SG_SG_DISTANCE_ANGSTROM = 2.5   # PROPKA 3.5.1 and PDB2PQR use the same value
MIN_SG_SG_DISTANCE_ANGSTROM = 1.6   # below this the atoms are coincident, not bonded
```

* `cysteine_sulfurs(structure)` — the `SG` atoms of every `CYS` residue, in parse order. Warns when one residue carries more than one `SG`, which is the symptom of merged alternate locations or insertion codes (see 7.2).
* `geometric_disulfides(sulfurs)` — every `SG`-`SG` pair inside the window, sorted shortest first, assigned greedily so each sulfur takes at most one partner and the shortest candidate wins a contested sulfur. Cysteine counts are small, tens at most, so a pairwise `numpy` distance matrix is used with no spatial index.
* `record_disulfides(structure, records)` — resolves parsed `SSBOND` records to `Residue` pairs, applying the four filters of 3.2.
* `assign_disulfides(structure, records=())` — the entry point. Applies records first, runs geometry over the cysteines left unclaimed, sets `residue.disulfide_partner` on both members of every pair, clears each affected residue's `_pka` cache, stores the list on `structure.disulfides` and returns it.

The `_pka` clear matters: review 2 demonstrated a stale `[{"CYS": 8.33}]` surviving re-detection when `assign_disulfides` is called on a structure whose `pkas` had already been read. It is one line and it makes the function correct from any call site, not only from inside the parser.

### 4.2 Parser

`read_ssbond_line(line)` parses one `SSBOND` record by column into `(chain1, number1, chain2, number2, sym1, sym2)`, tolerating a line truncated before the symmetry fields, and returns `None` for a line it cannot read rather than raising. It is testable on a single string.

`PDBparser._read_pdb` collects `SSBOND` lines during the pass it already makes over the file, and calls `assign_disulfides` after the atom loop and the existing C-terminus assignment. Detection therefore happens on **every** parse, including `prodes.read()`, not only inside the pipeline's `prepare_structure`, so a parsed `Structure` is always self-consistent about its cysteines.

`SSBOND` is only collected when the parser is reading `ATOM` records, its default. Under `parse(file, identifier="HETATM")` the records are ignored, because they refer to residues that parse is not building. Geometry still runs, and on a file whose cysteines are `HETATM` records it will find their bonds; that is harmless and arguably right, and it is noted here only because an earlier draft of this plan claimed it finds nothing, which is wrong.

### 4.3 `Residue`

* New attribute `disulfide_partner = None`, set in `__init__`.
* `pkas` omits the side-chain entry when `disulfide_partner is not None`. Terminus entries are untouched: a bridged cysteine at a chain end still has its titratable backbone group.
* New property `side_chain_pka` — the value of the entry in `pkas` whose key equals this residue's name, else `None`.
* New method `group_pka(key)` — the value of the entry under `key`, else `None`. Used for `"N+"` and `"C-"`.

### 4.4 `Atom.charge`

Both positional reads are replaced by keyed ones, each with an explicit `None` branch:

* side chain: `self.residue.side_chain_pka`. `None` means the group does not exist and the atom carries no side-chain charge. This single branch is what makes a bridged cysteine's `SG` neutral everywhere.
* terminus: `self.residue.group_pka("N+")` or `group_pka("C-")`. `None` means no terminus charge rather than a crash. With the merge of 5.1 in place this should be unreachable, and it is written defensively because review 2 demonstrated the crash.

`Residue.charged_atoms` already dispatches on the entry key and needs no change; with no `CYS` entry present it never appends the `SG`.

### 4.5 Acceptance

* 1GDW: 4 disulfides, all 8 cysteines partnered, relation symmetric.
* 1GPB: 0 disulfides despite 9 cysteines.
* 1AO6: 34 disulfides, and exactly `A34` and `B34` unpartnered.
* 1GDW formal charge at pH 8.5 is +7.0 not -1.0; isoelectric point 10.342 not 8.899.
* ARH96693 features bit-identical to `main`.

---

## 5. Phase 2: make the pKa file interaction honest

### 5.1 `Structure.redo_pkas` merges per group

Replacement becomes a merge. For each residue named in the dict, each supplied entry is applied to the matching group of that residue's current pKa list; groups the file does not name keep their default. Then:

* An entry naming the residue's own name where the residue is a **bridged cysteine** is dropped: structure beats file, because the file's value is a prediction about a group that does not exist. Logged at debug level, except for the 99.99 sentinel, which agrees with us and is not worth a line.
* An entry naming a group the residue does not have at all, for example a `SER` pKa from H++, is dropped with a warning. Today such an entry reaches `Residue.charged_atoms` and raises `TypeError` on `atom.name in None`.

### 5.2 Omitted groups are reported

After applying the dict, count the titratable groups, side chains and termini, that the file did not name, and emit one warning giving the count and the first few. Counted per **group**, not per residue: `tests/data/1GDW.pka` covers 45 residues but **46 groups**, residue 1 carrying both `LYS` and `N+`. **[verified 2026-08-28]** A residue-level count cannot see a missing terminus, which is the case that loses data (2.3).

Bridged cysteines are excluded from the count; their absence is correct. On `1GDW.pka` the warning is silent, because PROPKA covers all 46 groups.

### 5.3 Acceptance

* A pKa file giving a bridged cysteine a real pKa of 9.0 does not make it charged.
* A pKa file naming a side chain but not the terminus leaves the terminus pKa at its default, and the backbone atom's charge is unchanged.
* A pKa file with a titratable group removed warns, naming it.
* `1GDW.pka` warns about nothing, and `tests/test_pka.py` passes unchanged.

---

## 6. Phase 3: visibility, docs, version

* `calculate` prints the disulfide count immediately after `prepare_structure`, before the slow phases, in the plain style the pipeline already uses. Not in `prepare_structure` itself and not in the parser, so a library call to `prodes.read` stays quiet.
* `run_metadata` gains a `disulfides` count, passed as a keyword argument. It is a top-level metadata key, not inside `settings`, so `tests/test_output_bundle.py`, which compares `record["settings"]` exactly, is unaffected. **[verified 2026-08-28]** It is deliberately not a feature column: the feature set is pinned by `features_reduced.yaml`, `features_full_only.yaml` and asserted column by column in `tests/test_feature_dictionary.py`, and adding a column is a far larger contract change than this issue calls for.
* `README.rst`, under **pKa values and protonation states**: a subsection saying what is detected, the 2.5 A cutoff and that it matches PROPKA, that `SSBOND` wins per cysteine, that AlphaFold models carry no `SSBOND` and are handled geometrically, that PROPKA already reports bridged cysteines as non-titratable so a `--pka` run changes far less than a default one, and the two honest limitations: metal-bound cysteines are still titrated, and a model can place two cysteines at bonding distance without a bond.
* `README.rst`, the "No structure preparation" bullet currently lists "form disulfide bonds" among what Prodes does not do. Prodes still does not *form* them but now *recognises* them; correct the bullet.
* `README.rst`, the PROPKA impact paragraph quotes "18 of the 54 features" and "isoelectric point from 8.9 to 10.6" against the default run. The default run is what changes here, so both numbers are recomputed after implementation and restated.
* `README.rst`, "Output compatibility" promises values identical to the original Prodes at `--ionic-strength 0`. That promise now holds only for structures without disulfides; state the exception.
* `docs/identification_of_propka_dependent_features.md` reports an 819-structure analysis run against the old behaviour. The analysis cannot be re-run here, so it gets a dated note saying which behaviour it was measured against.
* `pyproject.toml`: 5.0.0 -> **6.0.0**. The reported quantity changes: formal charge, isoelectric point, average charge and every electrostatic feature move for any structure with a disulfide at any pH where a thiol would have titrated, and the isoelectric point moves at every pH. The version is stamped into every bundle's `prodes_run.json`, so the drift would be visible. Add the change to the README improvements list.

---

## 7. Phase 4: tests

### 7.1 New fixture

`tests/data/1AO6.pdb.zip`, human serum albumin, trimmed to `CRYST1`, `SSBOND`, `ATOM`, `TER` and `END` records. 187 KB zipped, the same order as the existing `1GPB.pdb.zip` at 120 KB. **[verified 2026-08-28]** It is the only fixture that carries all four properties this change needs at once:

| Property | Value |
|---|---|
| Cysteines | 70, over chains A and B |
| Disulfides by geometry | 34, longest 2.048 A, nearest non-bonded pair 5.951 A |
| Free thiols | exactly `A34` and `B34`, the textbook free cysteine of serum albumin |
| `SSBOND` records | 34, all `1555/1555`, all intra-chain, agreeing exactly with geometry |

It supplies the issue's third acceptance criterion, a free thiol titrating normally in the same structure as the bonded ones, on a real structure rather than a synthetic one. No shipped fixture could do that: 1GDW is all-bonded and the other three have no disulfides at all. It is also the repository's first multi-chain fixture.

The synthetic cases are built inline with `tmp_path` rather than committed, so the doctoring is visible in the test that depends on it.

### 7.2 `tests/test_disulfides.py`

* `read_ssbond_line` on a real record, on one with a non-`1555` symmetry operator, on one truncated before the symmetry fields, and on a malformed line.
* 1GDW: 4 pairs, exact partner numbers, symmetry of the relation.
* 1GPB: 9 cysteines, 0 pairs.
* 1AO6: 34 pairs from records; `A34` and `B34` unpartnered; the same 34 pairs found by geometry when the records are stripped, which is the cross-check that the two routes agree on real data.
* Contested sulfur: three `SG` atoms, two candidate pairs sharing one, shorter wins.
* Boundary: a pair at 2.45 A is called, one at 2.55 A is not, one at 1.5 A is not and is logged.
* Per-cysteine precedence: a fixture whose only record is a symmetry self-pair still gets its geometric bonds, the 1ERT case.
* An `SSBOND` naming a residue absent from the coordinates is skipped with a warning, not an exception.
* An `SSBOND` naming a pair 15 A apart is honoured, with a warning, and the cysteines it consumed are not re-paired by geometry.
* Charge consequences on 1GDW at pH 7.0, 8.5 and 9.0 against the table in 2.5.
* On 1AO6: a bridged cysteine is neutral at pH 9 and `A34` titrates normally.
* A bridged cysteine that is also the C-terminal residue keeps its C-terminal charge and gains no side-chain charge. This is the case that crashes without 4.4.

### 7.3 Existing tests

Confirmed to break, all on 1GDW or 1GDW_h, all expected **[verified by review 2's prototype run: 253 passed, 4 failed]**:

| Test | Why |
|---|---|
| `tests/test_residue.py::test_pkas` | `residues["CYS"]` is Cys 6 of 1GDW_h, bridged; asserts `[{"CYS": 8.33}]` |
| `tests/test_residue.py::test_charge` | same residue; asserts `charge(14) == -1` |
| `tests/test_structure.py::test_charge` | `charge(11)` on all-bridged lysozyme |
| `tests/test_structure.py::test_isoelectric_point` | 8.899 -> 10.342 |

The original plan claimed the `test_residue.py` assertions "should hold because `general_data.json` is not being edited". That was a non-sequitur and both reviews 2 and 3 caught it. Each is updated to the new value with a comment naming this issue, and `tests/test_residue.py` is added to the files-touched list. The two `test_residue` cases are additionally *extended* rather than merely retargeted: the free-thiol assertion moves to a free cysteine so that the textbook pKa stays covered.

Regression guard: `tests/data/ARH96693_prodes_orig_output.csv` holds no disulfides, so `tests/test_sasa.py` must pass untouched. If it does not, the change has leaked.

---

## 8. Files touched

| File | Change |
|---|---|
| `src/prodes/calculations/disulfides.py` | new: detection, pairing, record filtering, assignment |
| `src/prodes/io/parser.py` | `read_ssbond_line`; collect records; call `assign_disulfides` |
| `src/prodes/core/residue.py` | `disulfide_partner`, `side_chain_pka`, `group_pka`, skip bridged `CYS` in `pkas` |
| `src/prodes/core/atom.py` | read pKas by group name with `None` branches, not by list position |
| `src/prodes/core/structure.py` | `disulfides` attribute; `redo_pkas` merges per group, drops inapplicable entries, warns on omissions |
| `src/prodes/run.py` | report the count; pass it to `run_metadata` |
| `src/prodes/output.py` | `disulfides` in the run record |
| `pyproject.toml` | 5.0.0 -> 6.0.0 |
| `README.rst` | disulfide subsection; correct the structure-preparation bullet; restate the PROPKA impact numbers; note the compatibility exception; improvements list |
| `docs/identification_of_propka_dependent_features.md` | dated note on which behaviour it measured |
| `tests/data/1AO6.pdb.zip` | new fixture |
| `tests/test_disulfides.py` | new |
| `tests/test_residue.py`, `tests/test_structure.py` | updated expectations, free-thiol coverage moved to a free cysteine |
| `tests/test_pka.py` | merge, drop and omission-warning cases |

---

## 9. Known limitations, stated rather than hidden

These are real and are not fixed here. Each is written into the README or a code comment so a user meets it before it bites.

* **Metal-bound cysteines are still titrated as free thiols.** A zinc-finger or metallothionein cysteine is deprotonated and coordinated, not protonated, and no `SG`-`SG` criterion can see that. Metallothionein-2's cluster sits at 3.3 A and up, well outside any bonding cutoff. Prodes will charge those cysteines. Detecting them needs the metal, which lives in `HETATM` records the parser discards.
* **Thioether and other non-disulfide cysteine links** are not detected either: cytochrome c's `Cys14` and `Cys17` bond to the heme, 8.7 A apart, unreachable by any `SG`-`SG` cutoff.
* **A model can place two cysteines at bonding distance without a bond.** The metallothionein AlphaFold model does exactly that at 2.050 A. The reported count is how a user notices.
* **Low-confidence AlphaFold regions may model a real bond too long to detect.** Verified for a preproinsulin model with two of three bonds at 3.47 and 3.53 A by review 1. Missing a bond fails safe: it reproduces today's behaviour.
* **A monomer model cannot show an inter-chain bond.** An AlphaFold IgG1 heavy chain leaves the cysteines that would bond to the light chain and to the partner heavy chain looking free, and Prodes will titrate them.
* **Selenocysteine (`SEC`) raises `KeyError`** in `residue_data` and stops the run. Pre-existing, unrelated to this change.

## 10. Adjacent bugs found while reviewing, deliberately not fixed here

Each deserves its own issue. None is a prerequisite for this change.

1. **`redo_pkas` keys on residue number alone**, so in a multi-chain structure a pKa predicted for chain A residue 30 is applied to chain B residue 30 as well. For an antibody, the workflow issue #6 names, that is wrong. The disulfide code resolves residues by chain and number and does not inherit the defect.
2. **The parser ignores `MODEL`/`ENDMDL`.** A 20-model NMR file, 1PIT, parses to 58 residues with 16 901 atoms piled into one of them and a formal charge of **-1362**. **[verified 2026-08-28]** Any NMR ensemble currently produces nonsense silently.
3. **Insertion codes are dropped.** The parser never reads column 27, so `CYS A 100` and `ALA A 100A` merge into a single `Residue`. Antibody Kabat numbering makes this reachable for the named workflow. Detection here warns when a residue carries more than one `SG`, which is the visible symptom, but does not fix residue identity.
4. **Alternate locations are dropped by accident.** `name = line[12:17].strip()` includes the `altLoc` column, so an `SG` with `altLoc` `A` parses under the atom name `"SG A"`, matches nothing in `charged_atoms`, and is silently uncharged today. Detection will not see it either.
5. **`propka` is pinned in `environment.yml` at 3.5.1 but is not installed** in the local `prodes` environment. **[verified 2026-08-28]** `tests/test_dependency_versions.py` compares files against each other, not against the environment, so nothing catches it.

## 11. Risks

* **The cutoff is a chemistry choice.** 2.5 A matches PROPKA and PDB2PQR and is the most defensible single number available, but it is still a threshold, and section 9 lists what falls either side of it.
* **Records-first can be wrong.** A reduced structure that kept its `SSBOND` records is treated as oxidised. The distance warning is the mitigation; making it an error would break files that are merely strained.
* **This changes shipped numbers.** Any model fitted on Prodes features from disulfide-containing structures above pH 8 needs refitting. That is the point of the change, and it is why the version goes to 6.0.0.
