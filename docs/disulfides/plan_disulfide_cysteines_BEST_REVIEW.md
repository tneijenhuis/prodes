# Which review of the disulfide plan did the better work

Judgement on the three independent reviews of `plan_disulfide_cysteines.md`, written after re-verifying their contested claims. Date 2026-08-28.

## Ranking

1. **Review 1**, structural chemistry.
2. **Review 2**, software correctness.
3. **Review 3**, adversarial.

**This ranking is not a filter.** It says which reviewer did the better work, not which findings survive. Every one of the three found at least one thing the other two missed, and the refined plan takes material from all three. Dropping review 3's findings because it ranks last would lose the terminus-deletion data loss and the 45-versus-46 group count, both of which changed the design.

All three reached the same verdict: do not merge as written. None of them found a wrong number in the original plan's headline tables. The ranking therefore turns on what each found *beyond* the document, not on who caught the author out.

## Why review 1 wins

It found the one defect that would have shipped a quietly wrong feature: **the 2.35 A upper bound misses real disulfides.** Nothing inside this repository could have revealed that. Every shipped fixture has its bonds between 2.02 and 2.05 A, comfortably inside the original window, so any amount of careful work on the repo alone would have concluded the cutoff was fine. Review 1 went and measured 142 annotated disulfides across twelve entries, found the real range runs to 2.362 A, and produced a decisive counter-example.

Re-verified here: 6VXX, the SARS-CoV-2 spike trimer at 2.8 A, models `Cys391`-`Cys525` at 2.361, 2.361 and 2.362 A in the three protomers, while all its other bonds sit between 2.018 and 2.058 A. The original window drops three real bonds from that entry alone.

It then did the harder half: it established the ceiling as well as the floor. The AlphaFold model of metallothionein-2, a protein with zero disulfides, puts two metal-binding cysteines at 2.050 A, so the cutoff cannot be widened toward the metal-cluster distances without turning one false positive into a dozen. That converted "pick a bigger number" into "pick 2.5, the number PROPKA and PDB2PQR already use", which is both defensible and consistent with the tool the README tells users to run.

Its second decisive finding, the 1ERT case, killed the precedence rule. A file whose only `SSBOND` record is a symmetry self-pair would, under the original file-granularity rule, have had geometry suppressed entirely. That is a real deposited structure, not a hypothetical.

It also produced the negative result that saved work: the chi-3 dihedral filter that biotite uses would reject 54 per cent of real bonds, so the plan was right not to have one and now says why.

## What review 2 uniquely contributed, and must be kept

Review 2 used the best *method* of the three: it built a scratch prototype implementing the plan literally and ran the suite against it, 253 passed and 4 failed. That is how it found the only true blocking implementation defect, and it is why its test-breakage list is exact rather than predicted.

Kept in the refined plan:

* **The terminus lookup has no `None` fallback.** Demonstrated `TypeError: '>' not supported between NoneType and int` where `main` returns 2.0. The refined plan writes the branch defensively even though the merge of section 5.1 should make it unreachable.
* **The `_pka` cache goes stale for any caller other than the parser.** One line to fix, and without it `assign_disulfides` is only correct from the single call site the author had in mind.
* **Two planned tests could not be written as described**, because `Residue.pkas` dereferences `self.chain`, so a hand-built residue raises `AttributeError`. The refined plan builds its synthetic cases as real parsed files under `tmp_path` instead.
* **A truncated `SSBOND` line yields an empty symmetry field**, which under the original rule silently drops every bond. The refined plan reads a blank symmetry field as `1555`.
* Confirmed safe, so nobody has to re-check it: adding a `disulfides` key to `run_metadata` does not break `test_output_bundle.py`, which compares only `record["settings"]` exactly.

## What review 3 uniquely contributed, and must be kept

* **`redo_pkas` deletes a terminus pKa** when a file names a residue's side chain but not its terminus, and the residue-level omission warning the plan proposed could never catch it, because the residue *is* in the dict. This is the finding that turned a keyed-lookup tidy-up into the per-group merge of section 5.1, which is now the largest single change in the plan.
* **45 residues, 46 groups.** The original plan read `len(pkas) == 45` off a test and reported it as a group count. Small in itself, and it is the same confusion that produced the residue-level warning above, so it was worth catching twice.
* **Stale documentation the plan did not list**: the README's PROPKA impact numbers, the output-compatibility promise, and the 819-structure analysis doc, all measured against behaviour this change alters.
* **The scope challenge**, which is the finding I disagreed with and had to settle rather than accept. See below.
* The observation that at pH 7 the fix moves few features, which the refined plan states rather than letting the reader assume the change is dramatic at every pH.

## The one contested finding, and how it was settled

Review 3 argued that the `Atom.charge` rewrite is opportunistic: the whole fix can be had by giving a bridged cysteine a sentinel pKa of 99.99, the mechanism the code already relies on when PROPKA supplies one, with no change to `Atom.charge` at all. It reproduced this and it works numerically.

Settled against it, on grounds that are stated in the refined plan rather than assumed. The issue asks for the pKa to be `None` specifically so the fix propagates instead of being special-cased; a sentinel leaves `Residue.pkas` reporting a group that does not chemically exist and leaves `Residue.charged_atoms` listing the `SG` as charged carrying `-0.0`; and adopting a third-party format's magic number as internal state means prodes can never distinguish "PROPKA said so" from "prodes decided so".

The *terminus* half of the rewrite, which review 3 called an unrelated bug bundled in, is kept for a reason review 3 itself supplied: its own terminus-deletion finding means keying by name without the merge would crash, and the merge without keying by name would leave `pkas[-1]` reading whichever entry happened to be last. The two are one change, not two.

## Re-verification table

Every contested or design-changing claim, re-checked here rather than adjudicated from the documents.

| Claim | Raised by | Re-check result |
|---|---|---|
| PROPKA 3.5.1 uses a 2.5 A disulfide cutoff | R1 | **Confirmed.** `DISULFIDE_DISTANCE = 2.5` in `propka/bonds.py` of the pinned version |
| 99.99 is set in exactly one place, iff `cysteine_bridge` | R1 | **Confirmed.** One assignment, in `propka/group.py` |
| 6VXX has real bonds at 2.361-2.362 A | R1 | **Confirmed.** Three protomers, `Cys391`-`Cys525`, at 2.361, 2.361, 2.362 A; all other bonds 2.018-2.058 A |
| AlphaFold metallothionein-2 has a 2.050 A non-bond | R1 | **Confirmed.** `Cys15`-`Cys29` at 2.050 A; rest of the cluster 3.32 A and above; the protein has no disulfides |
| 1ERT's only `SSBOND` is a symmetry self-pair | R1 | **Confirmed.** `CYS A 73 - CYS A 73`, operators `2655`/`1555`; closest real `SG`-`SG` pair in the file is 3.921 A, so its true count is 0 |
| 1ERU carries the real `Cys32`-`Cys35` bond | R1 | **Confirmed.** 2.024 A, plus the same symmetry record |
| The deposited 1GDW carries its four `SSBOND` records | R1 | **Confirmed.** Four records, all `1555/1555`, lengths 2.03, 2.03, 2.05, 2.02, matching the measured distances exactly |
| PROPKA gives free cysteines real pKa values | mine | **Confirmed.** On 1GPB all nine come back between 8.32 and 13.79, none 99.99, so the sentinel is disulfide-specific |
| A 20-model NMR file parses catastrophically | R1 | **Confirmed.** 1PIT: 58 residues, 16 901 atoms in one of them, formal charge -1362 |
| `test_residue.py::test_pkas` and `::test_charge` break | R2, R3 | **Confirmed.** `residues["CYS"]` is Cys 6 of 1GDW_h, bridged; the plan's claim that they hold was wrong |
| A pKa file naming a side chain but not `N+` deletes the terminus pKa | R3 | **Confirmed.** `{1: [{"LYS": 11.34}]}` leaves `pkas == [{"LYS": 11.34}]`; the 9.69 is gone and the backbone `N` reads the lysine value |
| `1GDW.pka` covers 45 residues but 46 groups | R3 | **Confirmed.** 45 keys, 46 entries; residue 1 carries both `LYS` and `N+` |
| `propka` is pinned but not installed locally | R3 | **Confirmed.** `ModuleNotFoundError`, and `propka3` is not on the path |
| The sentinel approach reproduces the fix numerically | R3 | **Confirmed as fact, rejected as design.** See above |
| All four SG-SG distances, both pI values, the 16 charge-table cells, the fixture inventory, the 99.99 claims, the zero-omissions claim, the `pkas[0]` ordering bug | R2 and R3 independently | **All confirmed exact.** No wrong number in the original plan's tables |

The usual expectation for this exercise is at least one wrong number in the plan. There was not one, and both reviewers who went looking said so plainly rather than manufacturing an objection. The errors were all in judgement: a cutoff too tight, a precedence rule at the wrong granularity, and a blast radius under-counted.

## Scorecard

| Criterion | R1 chemistry | R2 correctness | R3 adversarial |
|---|---|---|---|
| Verified truth | Excellent, all claims held, most required outside data | Excellent, prototype-driven | Excellent, found two real counting errors |
| Blocking value | Highest: the cutoff would have shipped wrong | High: the only true implementation crash | High: the terminus data loss |
| Unique coverage | Cutoff, precedence, chi-3, the whole limitations section | Stale cache, untestable tests, truncated records, metadata safety | Terminus deletion, group count, stale docs, scope challenge |
| Actionability | High, gave the number to use and the source for it | Highest, gave the failing traceback and the one-line fix | High, though the scope challenge needed adjudication |

## What to do with the result

Review 1 forms the spine of section 3, the detection criterion, and section 9, the limitations. Review 2's findings are grafted into sections 4.1, 4.4 and 7.3. Review 3's terminus finding drives section 5.1, which is the largest structural change from the original plan, and its documentation findings drive section 6.

Order of implementation is the order of the refined plan: detection and the pKa removal first, because it has the sharpest test, then the pKa-file merge, then visibility and documentation.

## Meta lesson

The two highest-value findings came from stepping outside the repository in opposite directions.

Review 1 found the cutoff defect only because it downloaded structures the repository does not contain. Every shipped fixture has its disulfides between 2.02 and 2.05 A, so the original window looks perfectly calibrated from inside this repo, and a reviewer who stayed inside it would have confirmed the wrong number with real evidence. When a plan proposes a numeric threshold, the fixtures at hand are the one dataset guaranteed not to test it.

Review 2 found its blocking defect only because it stopped reading and built the thing. The crash is invisible in the plan's prose, which describes the change correctly; it appears only when the described change meets a code path the plan did not think to trace.
