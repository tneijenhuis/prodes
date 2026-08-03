# Redundant Feature Analysis

## Method

R² (Pearson correlation squared) was computed between all pairs of the 105 numeric
features produced by the original Prodes algorithm across **820 proteins**
(1638 rows total, one per protein chain). The full R² matrix is in
`original_feature_set_r_squared_matrix.csv`.

Pair counts above threshold:

| R² threshold | Pairs |
|---|---|
| ≥ 0.99 | 59 |
| ≥ 0.95 | 93 |
| ≥ 0.90 | 166 |

## Findings

### 1. Formal vs Average surface electrostatic potential (R² ≈ 0.997)

Every `SurfEp*Formal` feature is ~0.997 correlated with its `SurfEp*Average`
counterpart. The formal/average charge distinction — which is an implementation
detail of how atom charges are assigned (formal charge states vs pH-averaged
protonation) — makes almost no difference to the resulting surface potentials.

### 2. Mean / Trimean / Median within each stat family (R² ≈ 0.99)

Within each electrostatic-potential family (e.g. `SurfEpMean` vs `SurfEpTrimean`
vs `SurfEpMedian`), the three central-tendency statistics are ~0.99 correlated.
Trimean and median add almost no information beyond the mean. This pattern
repeats for the positive-only and negative-only subsets (`SurfEpPos*`,
`SurfEpNeg*`) and for the hydrophobicity families (`SurfMhp*`, `SurfPosMhp*`,
`SurfNegMhp*`).

### 3. Sum = Mean × N (R² = 1.0)

`ShellEpSumFormal` vs `ShellEpMeanFormal` is **exactly 1.0000**. The sum of
values is perfectly redundant with the mean when the count of values
(`NShellPosEpFormal` / `NShellNegEpFormal`) is also available.

### 4. NShellPosEpFormal = NShellNegEpFormal (R² = 1.0)

The number of positive and negative shell directions are **perfectly
correlated** — they always sum to 120 (the number of sunflower-sphere
directions), so one determines the other exactly.

### 5. Area ≈ NSurfPoints (R² = 0.9998)

Surface area is essentially a linear function of the number of surface grid
points. Keeping both provides almost no additional information.

### 6. Formal charge ≈ Average charge (R² = 0.9976)

`Formal charge` and `Average charge` are nearly identical across these 820
proteins. As with finding 1, this is an implementation detail of the charge
model rather than a biologically meaningful difference.

### 7. Hydrophobicity central-tendency stats (R² ≈ 0.97–0.98)

`SurfMhpMean` / `Trimean` / `Median`, `SurfPosMhpMean` / `Trimean` / `Median`,
and `SurfNegMhpMean` / `Trimean` / `Median` — same pattern as the EP families:
the three central-tendency statistics are near-identical.

### 8. NSurfNegMhp = NSurfPoints − NSurfPosMhp (exact)

Pairwise, the count of negative hydrophobicity surface points is R² = 0.9875 with
`Area` and 0.9865 with `NSurfPoints` — high, but below the ~0.997 seen in findings
1, 3, 5 and 6.

The pairwise figure understates it. The positive and negative MHP counts are
strict complements over the same point set, so

    NSurfPosMhp + NSurfNegMhp = NSurfPoints

holds exactly by construction, and `NSurfPoints` is in turn R² = 0.9998 with
`Area` (finding 5). `NSurfNegMhp` is therefore a near-exact **linear combination**
of `Area` and `NSurfPosMhp`, which no pairwise comparison can reveal.

The two counts are not interchangeable: `NSurfNegMhp` is a protein-size proxy
(R² 0.9875 with `Area`, 0.8579 with `Molecular weight`) because most surface
points have negative MHP, whereas `NSurfPosMhp` is not (0.6976 and 0.4476) and
carries the hydrophobic-patch signal. The negative count is the one to drop.

The `SurfNegMhp` distribution statistics are unaffected — `SurfNegMhpMean`
(R² 0.11 with `Area`) and `SurfNegMhpStd` (0.05) describe the shape of the
negative-MHP distribution rather than how much surface there is.

## Implications

The original 105-feature set contains substantial redundancy. A reduced feature
set could retain essentially the same information content with roughly half the
columns by:

- Keeping only `Mean` (or `Trimean`) from each stat family, dropping `Median`
  and `Sum`.
- Dropping one of `Formal` / `Average` for surface EP (they are ~0.997
  correlated).
- Dropping one of `NShellPosEpFormal` / `NShellNegEpFormal` (perfectly
  redundant).
- Dropping `Sum` features everywhere (recoverable from `Mean × N`).
- Keeping only one of `Area` / `NSurfPoints`.

## Decision: the default feature set

The findings above were acted on by making a reduced 54-feature set the default
output. The original 105 features are still available via `--full-features`; see
the README for how to enable them and when that is advisable.

### Which 51 features are dropped

**1. The average-charge surface electrostatics block — 19 features.** Every `*Average` column: `SurfEpMaxAverage`, `SurfEpMinAverage`, `SurfEpMeanAverage`, `SurfEpTrimeanAverage`, `SurfEpMedianAverage`, `SurfEpSumAverage`, `SurfEpStdAverage`, `NSurfPosEpAverage`, `SurfEpPosMeanAverage`, `SurfEpPosTrimeanAverage`, `SurfEpPosMedianAverage`, `SurfEpPosSumAverage`, `SurfEpPosStdAverage`, `NSurfNegEpAverage`, `SurfEpNegMeanAverage`, `SurfEpNegTrimeanAverage`, `SurfEpNegMedianAverage`, `SurfEpNegSumAverage`, `SurfEpNegStdAverage`. Each is ~0.997 correlated (R²) with its `*Formal` counterpart, which is kept. Formal vs average is a choice of charge model, not a different property of the protein. Skipping this block also skips a whole calculation phase, which is where the reduced set's speed advantage comes from.

**2. Trimean and Median from every statistics family — 18 features.** `SurfEpTrimeanFormal`, `SurfEpMedianFormal`, `SurfEpPosTrimeanFormal`, `SurfEpPosMedianFormal`, `SurfEpNegTrimeanFormal`, `SurfEpNegMedianFormal`, `SurfMhpTrimean`, `SurfMhpMedian`, `SurfPosMhpTrimean`, `SurfPosMhpMedian`, `SurfNegMhpTrimean`, `SurfNegMhpMedian`, `ShellEpTrimeanFormal`, `ShellEpMedianFormal`, `ShellEpPosTrimeanFormal`, `ShellEpPosMedianFormal`, `ShellEpNegTrimeanFormal`, `ShellEpNegMedianFormal`. Within each family the three central-tendency statistics are 0.97-0.99 correlated with the `Mean`, which is kept.

**3. Sum from every statistics family — 9 features.** `SurfEpSumFormal`, `SurfEpPosSumFormal`, `SurfEpNegSumFormal`, `SurfMhpSum`, `SurfPosMhpSum`, `SurfNegMhpSum`, `ShellEpSumFormal`, `ShellEpPosSumFormal`, `ShellEpNegSumFormal`. A sum is the mean times the count, and both the mean and the count are kept, so these are exactly redundant (R² = 1.0000 for the shell pair, where both terms appear in the output).

**4. Average charge — 1 feature.** R² = 0.998 with `Formal charge`, which is kept.

**5. NSurfPoints — 1 feature.** R² = 0.9998 with `Area`, which is kept and is the physically interpretable one of the pair.

**6. NShellNegEpFormal — 1 feature.** The positive and negative shell direction counts always sum to 120, the number of sunflower-sphere directions, so `NShellPosEpFormal` determines this one exactly.

**7. NSurfNegMhp — 1 feature.** See finding 8. `NSurfPosMhp + NSurfNegMhp = NSurfPoints` exactly by construction, and `NSurfPoints` is itself R² = 0.9998 with `Area`, so this count is a linear combination of two features that are kept. Its pairwise R² against `Area` is only 0.9875, below the threshold that governed the other drops, which is why it was retained in the first pass — see the note on the pairwise blind spot below.

**8. NSurfNegEpFormal — 1 feature.** See finding 10. The same complementary-count pattern: the positive and negative surface-EP counts partition the same point set, so this one is R² = 0.99977 predicted by `Area` and `NSurfPosEpFormal` together while its best single partner reaches only 0.8655. It is also the size proxy of the pair (R² 0.8655 with `Area`, against 0.0265 for the positive count), so `NSurfPosEpFormal` is the member that carries signal.

### Which 54 features are kept

- **Whole-molecule** (4) — `Molecular weight`, `Isoelectric point`, `Dipole`, `Formal charge`
- **Surface geometry** (3) — `Area`, `Shape max`, `Shape min`
- **Per-residue surface fractions** (20) — `ALASurfFrac` through `VALSurfFrac`, one per standard amino acid, all unchanged
- **Surface electrostatic potential, formal charges** (9) — `SurfEpMaxFormal`, `SurfEpMinFormal`, `SurfEpMeanFormal`, `SurfEpStdFormal`, `NSurfPosEpFormal`, `SurfEpPosMeanFormal`, `SurfEpPosStdFormal`, `SurfEpNegMeanFormal`, `SurfEpNegStdFormal`
- **Surface hydrophobic potential** (9) — `SurfMhpMax`, `SurfMhpMin`, `SurfMhpMean`, `SurfMhpStd`, `NSurfPosMhp`, `SurfPosMhpMean`, `SurfPosMhpStd`, `SurfNegMhpMean`, `SurfNegMhpStd`
- **Far-field shell electrostatics** (9) — `ShellEpMaxFormal`, `ShellEpminFormal`, `ShellEpMeanFormal`, `ShellEpStdFormal`, `NShellPosEpFormal`, `ShellEpPosMeanFormal`, `ShellEpPosStdFormal`, `ShellEpNegMeanFormal`, `ShellEpNegStdFormal`

Every distinct physical property the original calculates is still represented. What is gone is the second, third and fourth summary statistic of the same distribution, and the second charge model applied to the same surface.

### Why the reduced set is the default

The original package calculates all 105 features and leaves feature selection to the user. Tim Neijenhuis, the original author, makes the fair point that anyone building a model can simply pick the subset they need later, and for an experienced modeller that is exactly right.

Our experience is that this is not what usually happens in practice. A 105-column CSV tends to go into a model as it is. With the dataset sizes typical of this field, often tens to a few hundred proteins rather than tens of thousands, a design matrix in which dozens of columns are ~0.99 correlated with another column does two things: it invites overfitting, and it splits one underlying signal across several collinear columns, so feature importances and regression coefficients become hard to interpret and unstable across resampling. Feature reduction tends to be something people learn about after the first model disappoints.

So the default here is the subset we would select anyway. It is a default, not a restriction: the README lists the conditions under which you should override it with `--full-features`. Every choice above is backed by the findings in this document and by the complete R² matrix in [original_feature_set_r_squared_matrix.csv](original_feature_set_r_squared_matrix.csv), so you can disagree with any individual call and see the numbers it was based on.

## Second pass: the pairwise blind spot

Findings 1-8 came from comparing features two at a time. That screen cannot
detect a feature that is redundant only as a **combination** of two or more
others. `NSurfNegMhp` is the worked example: its best single partner explains
0.9875 of it, which looked survivable, but `NSurfPoints - NSurfPosMhp` gives it
exactly.

A second pass therefore regressed every feature on every *pair* of other
features, using the closed form for a two-regressor OLS fit, which depends only
on the three pairwise correlations:

    R² = (r_y1² + r_y2² - 2·r_y1·r_y2·r_12) / (1 - r_12²)

Every feature was also regressed on all others simultaneously. Same 1638 rows.
The analysis was a one-off; this document is the record of it.

### 9. Size proxies are now fully removed

Only three of the original 105 features are R² ≥ 0.95 against `Area`:
`NSurfPoints` (0.9998), `NSurfNegMhp` (0.9875) and `SurfNegMhpSum` (0.9703).
All three are dropped from the default set, so **no feature that survives is a
rescaled surface area**.

The strongest size relationship left among the kept features is
`Molecular weight` at 0.8411, which is a genuine physical property rather than a
restatement of the surface area, followed by `NSurfPosMhp` at 0.6976 and
`Dipole` at 0.4263.

### 10. NSurfPosEpFormal + NSurfNegEpFormal ≈ NSurfPoints ≈ Area

The same complementary-count pattern as finding 8, one level less exact. The
positive and negative surface-EP counts partition the same point set, so they
sum to `NSurfPoints` apart from points whose formal-charge potential is exactly
zero (on average 32 of ~7000, so the identity is near- rather than fully exact).

Before this pass it was the only remaining collinearity above R² = 0.995 inside
the kept set, and it is invisible pairwise:

| target | predicted by | R² on the pair | R² on best single feature |
|---|---|---|---|
| `NSurfNegEpFormal` | `Area` + `NSurfPosEpFormal` | 0.99977 | 0.8655 (`Area`) |
| `Area` | `NSurfPosEpFormal` + `NSurfNegEpFormal` | 0.99976 | 0.8655 (`NSurfNegEpFormal`) |
| `NSurfPosEpFormal` | `Area` + `NSurfNegEpFormal` | 0.99830 | 0.6059 (`NShellPosEpFormal`) |

These three features carry two dimensions of information, not three. As in
finding 8, the negative count is the size proxy (R² 0.8655 with `Area`) while
the positive count is almost size-independent (0.0265), so `NSurfNegEpFormal`
is the redundant member of the trio. It is dropped from the default set, which
leaves no pair of kept features predicting a third above R² = 0.995.

### 11. The 20 residue surface fractions sum to 1

Each `*SurfFrac` is a residue's area divided by the total surface area, so the
twenty of them sum to exactly 1 by construction (measured mean 1.0000). Any one
of them is therefore determined by the other nineteen, and in exact arithmetic
the twenty fractions plus an intercept are rank-deficient.

In the data they are saved rounded to three decimals, which perturbs the sum
just enough (std 0.0013) that the matrix tests as technically full rank — but
the direction is still severely ill-conditioned, and the redundancy is real
rather than an artefact of this dataset.

This is a compositional-data constraint, not a modelling error, and it is not
resolved by dropping features here: the usual treatments are to omit one residue
as a reference category, apply a log-ratio transform, or rely on regularisation.
Users fitting an unregularised linear model on all twenty fractions plus an
intercept should be aware of it.

### A caveat on reading the pair scan

The scan reports each target's best-scoring pair. Where a single feature already
explains a target almost perfectly — `ShellEpSumFormal` for `ShellEpMeanFormal`,
say — the second slot is filled by an arbitrary passenger such as `ARGSurfFrac`,
which contributes nothing. Compare the pair R² against the best-single R² in the
table above before reading a pair as a genuine two-feature dependency.
