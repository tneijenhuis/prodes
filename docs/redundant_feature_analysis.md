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

### 8. NSurfNegMhp ≈ Area / NSurfPoints (R² ≈ 0.99)

The count of negative hydrophobicity surface points scales almost perfectly
with total surface area.

## Implications

The original 105-feature set contains substantial redundancy. A reduced feature
set could retain the same information content with roughly 40 features by:

- Keeping only `Mean` (or `Trimean`) from each stat family, dropping `Median`
  and `Sum`.
- Dropping one of `Formal` / `Average` for surface EP (they are ~0.997
  correlated).
- Dropping one of `NShellPosEpFormal` / `NShellNegEpFormal` (perfectly
  redundant).
- Dropping `Sum` features everywhere (recoverable from `Mean × N`).
- Keeping only one of `Area` / `NSurfPoints`.
