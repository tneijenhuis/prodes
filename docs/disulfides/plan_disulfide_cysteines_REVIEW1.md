# Review 1 of `plan_disulfide_cysteines.md` — structural biology and physical chemistry

**Reviewer lens:** protein chemistry and structural correctness. Does the proposed criterion find real disulfides on real files, and are the resulting numbers defensible?
**Date:** 2026-08-28. Reviewed at `b5e8409` on branch `6-disulfide-cysteines`. No source, test or plan file was modified.
**Method:** every number in the plan was recomputed in the `prodes` conda environment. External evidence comes from the propka 3.5.1 and pdb2pqr 3.7.1 sources and from 12 PDB entries and 10 AlphaFold DB models downloaded for this review.

---

## Verdict

The plan's arithmetic is impeccable — every recomputable figure in it reproduces exactly, its diagnosis that `Atom.charge` must stop indexing `pkas` by list position is correct and is the load-bearing part of the change, and its correction of the issue's PROPKA premise is right for the right reason. The chemistry model, "a bonded cysteine has no titratable group", is the correct one and matches what PDB2PQR does. Where it falls down is the detection criterion and its precedence rule. The 2.05 ± 0.3 Å window is tighter than every production tool, tighter than the two tools Prodes already interoperates with, and demonstrably misses three real annotated disulfides in a 2.8 Å cryo-EM structure; the file-level "records beat geometry" gate has a constructible failure on a real PDB entry; and the plan is silent on insertion codes and multi-MODEL files, both of which break the `(chain, number)` key it proposes to resolve `SSBOND` records with. **Do not merge as written.** Findings 1, 3, 4, 5 and 6 are cheap to fix at the plan stage and expensive to fix after the code is written.

---

## Findings

### 1. Serious — the 2.35 Å upper bound is too tight. Use 2.5 Å.

**Plan says:** window `2.05 ± 0.3` Å, i.e. 1.75–2.35 Å. The risk section concedes "a low-resolution or coarse model could put a real bond at 2.4 Å and be missed", and keeps the window anyway on the grounds that missing a bond fails safe.

**Actually true:** it is not hypothetical, and the failure is already in the PDB. Across 142 `SSBOND`-annotated `1555/1555` disulfides in 12 entries **[verified 2026-08-28]**:

| Entry | resolution | bonds | min SG–SG | max SG–SG | nearest **non**-bonded SG–SG |
|---|---|---|---|---|---|
| 2VB1 | 0.65 Å | 4 | 2.026 | 2.038 | 6.373 |
| 1AKI | 1.50 Å | 4 | 1.970 | 2.018 | 6.650 |
| 1GDW | 1.80 Å | 4 | 2.022 | 2.050 | 6.541 |
| 1CRN | atomic | 3 | 2.004 | 2.051 | 7.566 |
| 1PIT | NMR | 3 | 2.071 | 2.086 | 6.724 |
| 3HHR | 2.80 Å | 8 | 1.996 | 2.033 | 6.446 |
| 1IGT | 2.80 Å | 17 | 2.013 | 2.041 | 5.251 |
| 1HZH | 2.70 Å | 15 | 2.020 | 2.046 | 11.556 |
| 6VXX | 2.80 Å cryo-EM | 36 | 2.018 | **2.362** | 7.945 |
| 7CN4 | 2.93 Å | 45 | 2.019 | 2.037 | 9.364 |
| 1ERU / 2TRX | X-ray | 3 | 2.024 | 2.090 | 9.68 / 28.27 |

Overall: **max bonded 2.362 Å, min non-bonded 5.251 Å.** 6VXX (SARS-CoV-2 spike ectodomain) has Cys391–Cys525 at 2.361, 2.361 and 2.362 Å in protomers A, B and C. All three are annotated disulfides. All three fall **outside** the plan's window.

What the field uses **[verified 2026-08-28, from source]**:

| Tool | criterion |
|---|---|
| PROPKA 3.5.1 | `DISULFIDE_DISTANCE = 2.5` (`propka/bonds.py:24`), any S–S pair by element |
| PDB2PQR 3.7.1 | `BONDED_SS_LIMIT = 2.5` (`pdb2pqr/config.py:162`), SG–SG |
| Biotite example | 2.05 ± 0.05 Å **and** χ3 = 90 ± 10° |
| Biopython | no built-in; the community threshold is 3.0 Å |
| PyMOL | `distance <= connect_cutoff (0.35) + 0.2 + (vdw1+vdw2)/2` — effectively ~2.2–2.4 Å depending on PyMOL's S radius (*unverified*, formula from the PyMOL wiki, radius not checked) |
| PDB intermolecular-disulfide survey (PMC9323673) | 2.5 Å, chosen "to account for inaccuracies in experimental data" |

**Fix:** replace the symmetric window with a single upper cutoff of **2.5 Å**, matching PROPKA and PDB2PQR — the two tools Prodes already reads pKa files from, so agreement between Prodes' own detection and a supplied `.pka` becomes exact rather than approximate. Keep a nominal 1.5 Å floor as a corruption check only; there is no physical reason for a 1.75 Å lower bound. The measured 2.36 → 5.25 Å gap means 2.5 Å costs nothing in false positives on experimental structures. Do not go above 2.5 (see finding 2).

---

### 2. Serious — metal-ligating cysteines are the real false-positive risk, and they bound the cutoff from above.

**Plan says:** nothing. Metal sites are not mentioned.

**Actually true:** a bare SG–SG distance rule cannot tell a disulfide from a metal-bridged thiolate pair, and on apo models it gets it wrong **[verified 2026-08-28]**:

- **AF-P02795, human metallothionein-2** — 20 cysteines, *zero* disulfides in reality (every one is a Zn/Cd thiolate ligand). AlphaFold, which places no metals, models **Cys15–Cys29 at 2.050 Å** — a hard false positive dead centre in the plan's window. Seven further pairs sit at 3.32–3.50 Å.
- **AF-P10109, human ferredoxin-1** — the [2Fe-2S] ligand pairs at 3.318 and 3.364 Å (pLDDT 95–98).
- **AF-P08046, EGR1 zinc fingers** — no SG pair under 3.5 Å.

So the metal-cluster S···S floor sits around 3.3 Å. Anything at or below 2.5 Å is safe against it; 3.0 Å is marginal; the 3.0 Å threshold Biopython users commonly quote is not.

The false positive on an apo metallothionein model is arguably the *less* wrong answer — a Zn-ligating cysteine is a thiolate whose charge is compensated by a metal that Prodes does not model at all (HETATM is not read), so charging it as a free thiol with pKa 8.33 is already wrong. But this must be a stated design position, not something a user discovers.

**Fix:** cap the cutoff at 2.5 Å, and say in the README subsection Phase 3 already plans that metal-coordinating cysteines are not recognised as such and that a metal-site model may occasionally be called as a disulfide.

---

### 3. Serious — "records beat geometry" must be per-cysteine, not per-file.

**Plan says:** `assign_disulfides` "uses the `SSBOND` records when the file has any, otherwise geometry."

**Actually true:** a file can have `SSBOND` records that resolve to nothing, and the gate then switches geometry off entirely. This is not a thought experiment **[verified 2026-08-28]**:

```
1ERT (reduced human thioredoxin, 1.70 A) — its ONLY SSBOND record:
SSBOND   1 CYS A   73    CYS A   73                          2655   1555  2.05
1ERU (oxidised human thioredoxin) — two records:
SSBOND   1 CYS A   32    CYS A   35                          1555   1555  2.02
SSBOND   2 CYS A   73    CYS A   73                          1555   2655  1.99
```

Both entries carry a symmetry-copy record, which the plan correctly skips. But strip the 32–35 line from 1ERU and the plan's gate sees "the file has records", disables geometry, skips the one record it has for symmetry, and reports **0 disulfides** where the true answer is 1 — the Cys32–Cys35 bond is right there at 2.024 Å.

**Fix:** honour records for the cysteines a resolvable `1555/1555` record names, and run geometry over the cysteines that no such record names. Identical behaviour on well-formed files, no cliff on partial ones.

Evidence that this costs nothing: across 1GDW, 2VB1, 1AKI, 1CRN, 1PIT, 3HHR, 1IGT, 1HZH, 6VXX, 7CN4, 1ERU and 2TRX there is not one geometric pair under 3.0 Å without a matching record, and not one `1555/1555` record without matching geometry **[verified 2026-08-28]**. Modern wwPDB `SSBOND` lists are complete; the per-cysteine rule only ever helps on damaged or edited input.

Related, smaller: the plan proposes skipping any record with "a symmetry operator other than 1555 in either symmetry field". All 144 records I looked at have both fields populated **[verified 2026-08-28]**, but a line truncated by a third-party writer loses them. Treat blank/absent symmetry fields as `1555`, not as "skip".

---

### 4. Serious — a record can name the same residue on both sides. Reject that explicitly.

**Plan says:** `record_disulfides` resolves `(chain, number)` pairs and ignores pairs naming a residue absent from the coordinates. Nothing about a record whose two ends are the same residue.

**Actually true:** 1ERT's only record is `CYS A 73 — CYS A 73` (see above). The plan's symmetry check happens to catch this one, but nothing in the design forbids self-pairing, and a self-pair would set `residue.disulfide_partner = residue` and silence a residue that is a genuinely free, reactive surface thiol — Cys73 of human thioredoxin is the glutathionylation site.

**Fix:** an explicit guard that the two ends must be different `Residue` objects, skipped with a warning. Also add the insertion-code fields (columns 22 and 36) to the `read_ssbond_line` return value the plan specifies; they are needed for finding 5 and are free to parse.

---

### 5. Serious — insertion codes silently merge residues, so `(chain, number)` is not a key.

**Plan says:** `record_disulfides` "resolves `(chain, number)` pairs read from `SSBOND` records to `Residue` objects".

**Actually true:** `PDBparser._read_line` reads `residue_number = int(line[22:26].strip())` and never touches column 27 (index 26), the insertion code. A two-residue test file with `CYS A 100` followed by `ALA A 100A` parses as **one** `Residue`, named `CYS`, holding all ten atoms **[verified 2026-08-28]**:

```
CYS 100 terminus: N atoms: ['N','CA','C','O','CB','SG','N','CA','C','O']
```

For the antibody workflow the issue is written around this is not an edge case: Kabat and Chothia numbering routinely use H100A…H100Z, and a `SSBOND` naming `CYS H 100` cannot be distinguished from one naming `CYS H 100A`.

**Fix:** at minimum, state the limitation and have the parser warn when it sees a duplicate `(chain, number)`. The cleaner route, which also sidesteps the merge entirely: resolve each `SSBOND` record by finding the SG pair whose distance is closest to the record's declared `Length` field (columns 74–78, populated in all 144 records I checked), and fall back to number matching only when that is ambiguous.

---

### 6. Serious — multi-MODEL files. Specify that detection iterates residues, and say the format is unsupported.

**Plan says:** nothing about `MODEL`. `sulfur_atoms(structure)` is described as "the `SG` atoms of every `CYS` residue, in parse order".

**Actually true:** `_read_pdb` has no `MODEL` handling and the result is already catastrophic, before any disulfide work. On 1PIT (BPTI, 20 NMR models) Prodes parses 17 780 atoms into 58 residues and puts **16 901 of them into the last residue**, `ALA 58`, because `_make_new_residue` fires only for unseen residue numbers and `self.current_residue` is never reset at a model boundary. Formal charge at pH 7 comes out as **-1362** **[verified 2026-08-28]**.

Why it matters here: if `sulfur_atoms` iterates `structure.residues` it happens to see only model 1's eight SGs and gets the right three pairs (2.071, 2.080, 2.086 Å). If it iterates `structure.atoms` filtering on `residue_name == "CYS" and name == "SG"` — the natural NumPy-friendly implementation, and the one `run.py` already uses for charged atoms — it sees all 120 SGs and produces cross-model pairs, including pairs inside a single `Residue`. The plan does not say which, and picks the wrong default by omission.

**Fix:** state `structure.residues` explicitly in the plan, and add the "an SG may not bond to a residue that is itself" guard from finding 4 as a second line of defence. Separately, raise the multi-MODEL parser defect as its own issue, next to the multi-chain `redo_pkas` keying bug the plan already flags — a `-1362` charge on BPTI is worse than anything in issue #6.

---

### 7. Serious — the AlphaFold claim is over-general, and no shipped fixture tests it.

**Plan says:** "AlphaFold models contain no `SSBOND` records but do place the `SG` atoms at bonding distance, so geometric detection works and is the only viable route." The plan's own fixture table then shows that both AlphaFold fixtures have zero disulfides, so nothing tests the claim.

**Actually true** for high-confidence regions, and false for low-confidence ones **[verified 2026-08-28, AlphaFold DB v6]**:

| Model | disulfides found (Å) |
|---|---|
| AF-P61626 human lysozyme | 2.018, 2.032, 2.039, 2.054 (all pLDDT 99) |
| AF-P00698 hen lysozyme | 2.016–2.047 (pLDDT 99) |
| AF-P00760 bovine trypsin | 2.024–2.035, all six |
| AF-P02769 BSA | 2.011–2.066, all seventeen |
| AF-P01344 IGF-2 | 2.002, 2.030, 2.051 — all three, even at pLDDT 60 |
| **AF-P01308 human preproinsulin** | B7–A7 at **2.086**; B19–A20 at **3.467**; A6–A11 at **3.527** (pLDDT 49–56) |

Two of insulin's three disulfides are modelled at 3.5 Å and are unrecoverable by any cutoff that stays safe against metal sites (finding 2).

Second limitation, specific to the workflow the issue names: an AlphaFold model of a *single* antibody chain physically cannot form interchain disulfides. **AF-P01857 (IgG1 heavy chain constant region)** has 11 SG atoms; geometry finds the 3 intradomain bonds and leaves 5 free — the light-chain and hinge cysteines, which are all bonded in the real IgG1 **[verified 2026-08-28]**. Monomer models will still over-charge them.

**Fix:** narrow the README wording to high-confidence AlphaFold regions; state the low-pLDDT and monomer-model limits explicitly; make the logged disulfide count the user's sanity check, which Phase 3 already provides. If a positive AlphaFold fixture is wanted, AF-P61626 is human lysozyme and would pair naturally with 1GDW.

---

### 8. Minor, and a "do not churn" — do **not** add a CB–SG–SG–CB dihedral filter.

Biotite's example gates on χ3 = 90 ± 10°, so a reviewer or refiner may push for it. Measured on the same 142 annotated disulfides **[verified 2026-08-28]**: |χ3| ranges **37.1° to 174.1°**, median 89.1°, and **77 of 142 (54%) fall outside 80–100°**. The extremes are ordinary bonds — 1IGT 214–128 at 37.1°, 7CN4 480–488 at 174.1°, 6VXX 391–525 at 45.2°. A dihedral filter would reject over half of real disulfides.

Distance alone is unambiguous because the bonded/non-bonded gap is ~2.9 Å wide (finding 1). The plan is right to use distance only, and the plan should say *why*, so this does not get relitigated.

---

### 9. Minor — the plan's PROPKA analysis is correct; two details to tighten.

**Verified in full [verified 2026-08-28]:**
- `tests/data/1GDW.pka` lists all eight cysteines with `99.99` in both the detailed table (lines 94–108) and the summary (lines 208–215), with model pKa 9.00.
- `convert_propka` returns `{6: [{'CYS': 99.99}], ...}` and `redo_pkas` applies it. `neg_charge(99.99, ph)` does not overflow and `99.99 < ph` is false at every reachable pH, so the residue is neutral by both routes.
- **99.99 has exactly one meaning in PROPKA output.** In propka 3.5.1 it is assigned in a single place: `Group.calculate_total_pka` returns early with `self.pka_value = 99.99` if and only if `self.atom.cysteine_bridge` (`propka/group.py:502–505`). There is no second meaning.
- PROPKA is **geometry-only**: the string `SSBOND` does not appear anywhere in the package. `_find_bonds_for_atoms` sets `cysteine_bridge` on any pair of `S` **elements** within 2.5 Å — keyed on element, not atom name, so in principle a MET SD could trip it.
- With `1GDW.pka` applied, 1GDW gives formal +8 / +7 / +6 / +6 at pH 7.0 / 7.4 / 8.5 / 9.0 and pI **10.565**, against −1 and pI 8.899 without it. The plan's headline conclusion — this is a *no-pKa-file* bug — is right.

**Two corrections.** The plan says PROPKA "covers all 45 titratable groups". It is 45 *residues* carrying 46 *groups*; residue 1 carries both `LYS 11.34` and `N+ 8.35`. Zero omissions either way, so the conclusion stands. Second, note that PROPKA's summary emits `N+` **last**, after `LYS 1`, which is why `convert_propka` happens to produce side-chain-first order; that ordering is an accident of PROPKA's output layout, not a guarantee, which strengthens the plan's own case for key-based lookup.

---

### 10. Minor — the plan's latent ordering bug is wrong at both ends, not one.

**Plan says:** with `[{"N+": 8.35}, {"LYS": 11.34}]`, "the `NZ` atom of Lys 1 reads the N-terminal pKa and comes out at charge 0 at pH 10 instead of 1.0."

**Also true, and unstated:** the backbone `N` reads `pkas[-1]` = 11.34 and returns **1** where the correct answer with 8.35 is **0** **[verified 2026-08-28]** on 1GDW residue 1. Both reads are positional and both are wrong in the reversed case. The proposed fix (key lookup for both the side chain and the terminus) already covers it; just describe it correctly so the test is written for both.

---

### 11. Minor — `pkas` stops caching for bonded cysteines, in a hot loop.

With the side-chain entry skipped and no terminus entry, `len(pkas) == 0`, so `_pka` stays `None` and the property recomputes on every access. `Atom.charge` is called once per atom per surface point inside `Point.electrostatic_potential` (`src/prodes/core/point.py:110` and `:116`). `prodes.data.data` is loaded once at import, so each recompute is a dict lookup plus two NumPy element comparisons — small, but free to avoid.

**Fix:** cache the empty result (store `[]` and have `pkas` return `self._pka or None`), or memoise `side_chain_pka`. Worth one sentence in the plan since the plan already makes the `_pka` cache ordering load-bearing.

---

### 12. Minor — the fixture biology labels.

**1GDW is not "lysozyme" generically.** It is **mutant human lysozyme**, R21G (`SEQADV 1GDW GLY A 21 UNP P61626 ARG 39 ENGINEERED MUTATION`), 1.80 Å, with one Na⁺ ion, `TITLE: CRYSTAL STRUCTURE OF MUTANT HUMAN LYSOZYME SUBSTITUTED AT LEFT-HANDED HELICAL POSITIONS` **[verified 2026-08-28, from the deposited file]**. The 130-residue length and the 6–128 / 30–116 / 65–81 / 77–95 pairing are the human pattern; hen egg-white lysozyme is 129 residues with 6–127 / 30–115 / 64–80 / 76–94. The plan's numbers are right, but the bare label invites the next reader to check them against HEWL and conclude they are wrong. Say "human lysozyme (R21G mutant)".

**The deposited 1GDW does carry its four `SSBOND` records** — lengths 2.03 / 2.03 / 2.05 / 2.02, matching my measured 2.030 / 2.031 / 2.050 / 2.022. Only the shipped fixture is stripped to `CRYST1` / `ATOM` / `TER` / `END`. The plan proposes hand-writing a synthetic `SSBOND` fixture; re-zipping the real 1GDW with its header would be a better one, since it lets the records-versus-geometry agreement test run on real data instead of on coordinates someone typed.

**1GPB checks out:** glycogen phosphorylase b, 1.90 Å, **zero** `SSBOND` records, 9 SG atoms, PLP and water as the only HETATM **[verified 2026-08-28]**. Nine genuinely free thiols in a cytosolic enzyme. Good negative control, as the plan says. (Unrelated: its PLP is a Schiff base on Lys680, which Prodes titrates as a free lysine. Out of scope, not a disulfide problem.)

---

### 13. Minor — quantify the dipole change, since it moves and no acceptance criterion covers it.

The plan's acceptance criteria stop at formal charge, average charge and pI. `Structure.dipole` also moves, because it iterates `residue.charged_atoms(ph)` and a bonded Cys's SG stops appearing there **[verified 2026-08-28]**, 1GDW:

| pH | dipole formal, now → fixed | dipole partial, now → fixed |
|---|---|---|
| 7.0 | 162.053 → 162.053 | 156.020 → 156.201 |
| 8.5 | 169.046 → **162.053** | 165.512 → 164.528 |
| 9.0 | 169.046 → **162.053** | 173.390 → 170.064 |
| 11.0 | 294.782 → **280.684** | 270.642 → 260.978 |

Small in relative terms, but `Dipole` is a shipped feature that changes, so it belongs in the release note and ideally in one assertion.

---

## What the plan got right

Do not churn any of this.

1. **Every recomputable number.** SG–SG distances 2.022 / 2.030 / 2.031 / 2.050 with the next nearest pair at 6.541 Å; the charge table at pH 7.0 / 7.4 / 8.5 / 9.0 (formal 7 / 7 / −1 / −1 now, 7 / 7 / 7 / 7 fixed; average 6.735 / 6.178 / 1.936 / −0.458 now, 7.092 / 7.019 / 6.709 / 6.133 fixed); pI 8.899 → 10.342; the fixture table (ATOM counts 1022 / 2003 / 6699 / 479 / 3106, zero `SSBOND` in all five, `altLoc` blank on every atom of all five). All **[verified 2026-08-28]**.
2. **The diagnosis that setting the pKa to `None` alone is not enough, and crashes.** Verified precisely: a **C-terminal** cysteine with the side-chain entry dropped has `pkas == [{'C-': 2.34}]`, and `Atom.charge` hands 2.34 to the SG and returns a **full −1.0 at pH 7**. A non-terminal one whose `pkas` returns `None` raises `TypeError: 'NoneType' object is not subscriptable`. Moving both reads to key lookup is the correct fix and is the actual substance of the change; the disulfide detection is the easy half.
3. **Correcting the issue's PROPKA premise** and drawing the right consequence from it: the bug bites hardest on the *unprepared* workflow, and that belongs in the README because it inverts the issue's framing.
4. **"No titratable group" over "very high pKa".** This is the chemically correct model: a cystine sulfur has no S–H, so it is not a weak acid with an inaccessible pKa, it is not an acid. PDB2PQR takes the same view — its `CYX` patch removes `HG` outright (`pdb2pqr/aa.py:412–420`). PROPKA's 99.99 exists only because its output is a fixed-width table with one row per group; it is a formatting artefact, not a model, and the plan is right to treat structural detection as beating the file.
5. **Distance-only criterion, no dihedral** (finding 8).
6. **Detection in the parser, not in `prepare_structure`**, so `prodes.read()` (`src/prodes/__init__.py:54`) and hand-built test structures are consistent with pipeline runs. Correct call.
7. **The altLoc note.** Verified: an SG with `altLoc` `A` parses as the atom name `"SG A"`, so a disordered cysteine is invisible to both `charged_atoms` and detection. Correctly identified as pre-existing and out of scope.
8. **Deferring the multi-chain `redo_pkas` keying bug** to its own issue, and saying that the disulfide code resolves by chain and number so it does not inherit the defect.
9. **5.0.0 → 6.0.0.** Correct under the versioning rule: the reported quantity changes, and `prodes_version` is stamped into every bundle.
10. **Not adding a feature column.** Metadata-only is the right call given `tests/test_feature_dictionary.py` asserts the column set.

The existing suite passes at `b5e8409` (`tests/test_structure.py test_pka.py test_residue.py test_atom.py test_parser.py`, 44 passed) **[verified 2026-08-28]**, and `test_structure.py:39` (`pI == 8.899`) and `:62` (`charge(11) == -13`) are the two assertions the plan correctly identifies as needing new values. For the record, the new pH-11 formal charge is **−5.0** **[verified 2026-08-28]**, which the plan does not state.

---

## What the plan is silent about and should cover

**Cysteines that are neither free thiols nor cystines** and will still carry pKa 8.33 after this change. The README subsection Phase 3 already plans is the right place for one paragraph:

- **Metal ligands** — zinc fingers, Fe–S clusters, metallothionein (finding 2).
- **Thioether links to heme.** Cytochrome *c* Cys14 and Cys17 are covalently attached to the haem vinyls and are not titratable. Their SG–SG distance is 8.712 Å in 1HRC and 8.459 Å in 5CYT **[verified 2026-08-28]**, so no distance criterion will ever pair them. This is a permanent residual error, small in magnitude, but it should be named rather than left for a user to find.
- **Trisulfide bridges**, a documented IgG2 product-quality attribute. The extra sulfur pushes SG–SG to roughly 4 Å, so a trisulfide reads as two free thiols.
- **Glutathionylated, palmitoylated and prenylated cysteines** — same story.
- **HETATM cysteines.** `PDBparser.parse` defaults to `identifier="ATOM"`, so modified cysteines deposited as HETATM (`CSO`, `CSD`, `CME`, `OCS`, `SNC`) are not read at all. No regression, but the reported disulfide count will not match the file's `SSBOND`/`LINK` count for such entries, and the log line from Phase 3 is where a user would notice the discrepancy.
- **Selenocysteine is a non-issue** and should be dismissed in one sentence so nobody spends effort on it: `data.residue_data("SEC")` raises `KeyError`, so Prodes cannot compute charge or mass for *any* structure containing `SEC` **[verified 2026-08-28]**. The parse succeeds; the first charge call dies.

**Isoelectric point range.** `Structure.isoelectric_point` bisects a hard-coded 0–14 window and returns the last midpoint with no warning if the charge never crosses zero. Silencing cysteines moves every pI upward — 1GDW goes 8.899 → 10.342 — so more structures will approach the ceiling and a very basic, cysteine-rich protein could saturate silently. Either warn on saturation or say in the README that a reported pI near 14 means "did not cross zero".

**Hydrophobicity.** The issue explicitly raised a second `CYS` entry and deferred it; the plan does not mention it at all, which will get it rediscovered. One sentence recording it as deferred is enough. Worth noting when it is picked up: the Miyazawa–Jernigan value Prodes ships for `CYS` (7.93, near Ile at 8.83) comes from a PDB-derived contact statistic dominated by *buried, disulfide-bonded* cysteines, so the single entry is already biased toward the cystine form. That is an argument for leaving `hydrophobicity.json` alone in this change, not for splitting it.

**A mixed real acceptance fixture.** Every acceptance criterion in the plan is either an all-bonded structure (1GDW), an all-free one (1GPB, the AlphaFold pair), or a hand-written toy. **1HZH** (IgG1 b12, 2.70 Å) is 32 SG atoms, 15 disulfides and 2 free cysteines — a real structure that exercises the mixed case, the interchain case and the antibody workflow the issue is written around, in one file. If its size is acceptable zipped, it is a better acceptance fixture than anything typed by hand.
