# Validating the screened electrostatic potential against APBS

Measured 2026-08-26 on the branch that adds ionic screening, using the shipped code path (`prepare_structure`, `construct_surface_grid`, `calculate_surface_grid_features`) rather than a standalone reimplementation, so what is validated is what users get.

Reference: an APBS Poisson-Boltzmann solution at pH 7, 150 mM, protein dielectric 4, sampled at the exact coordinates of the Prodes surface points. Thirteen proteins, 390,427 surface points. The measurement and the reference data are on the `analysis/apbs-electrostatics-comparison` branch.

## Result

| | shipped 4.0.0 | screened, 150 mM |
| --- | --- | --- |
| median overlap with the APBS positive surface (Jaccard) | 0.058 | **0.658** |
| `NSurfPosEpFormal` against APBS, across proteins (Spearman) | 0.641 | **0.995** |
| positive surface, mean across 13 proteins | 6.0 % | **32.4 %** |
| proteins reporting no positive surface at all | 4 | **0** |

The APBS reference calls 33.9 % of the surface positive on average, so the screened value of 32.4 % is close, where the unscreened 6.0 % was not.

A Jaccard of 0.058 means the set of points Prodes called positive and the set APBS called positive were very nearly disjoint. They now overlap substantially.

## Per protein

| protein | Jaccard, shipped | Jaccard, screened | positive %, shipped | positive %, screened | positive %, APBS |
| --- | --- | --- | --- | --- | --- |
| 1AO6 | 0.058 | 0.713 | 2.6 | 37.7 | 45.2 |
| 1BSQ | 0.297 | 0.767 | 10.0 | 34.1 | 33.6 |
| 1CF3 | 0.000 | 0.592 | 0.0 | 26.0 | 26.1 |
| 1F6R | 0.000 | 0.602 | 0.0 | 30.9 | 26.4 |
| 1F8N | 0.371 | 0.699 | 17.3 | 47.2 | 44.4 |
| 1OVT | 0.317 | 0.663 | 16.0 | 45.2 | 48.9 |
| 1TRH | 0.107 | 0.622 | 4.2 | 37.0 | 38.9 |
| 3LA4 | 0.042 | 0.629 | 1.6 | 37.8 | 36.9 |
| 4F5S | 0.002 | 0.680 | 0.1 | 33.1 | 40.8 |
| 4PEP | 0.000 | 0.301 | 0.0 | 1.6 | 3.5 |
| 6FRV | 0.000 | 0.424 | 0.0 | 9.0 | 12.9 |
| 6PO0 | 0.431 | 0.658 | 20.3 | 45.2 | 42.3 |
| AF-P01070 | 0.128 | 0.688 | 5.3 | 36.6 | 40.9 |

Glucose oxidase (1CF3) is the clearest case. Its published Prodes figure is a uniformly red blob and it reported no positive surface anywhere, while APBS calls 26.1 % of that same surface positive. Screened, Prodes reports 26.0 %.

4PEP is the weakest, at a Jaccard of 0.301. It is also the protein whose phosphoserine both tools discard, so some of its real charge is missing from both models.

## What is still not fixed

Screening narrows the gap; it does not close it. The screened potential is still a Coulomb sum at a uniform relative permittivity of 4 with no dielectric boundary, and the screening length used is the one for bulk water rather than the value self consistent with that permittivity. It should not be described as a Poisson-Boltzmann calculation or compared to one in absolute terms.

## The shell features are deliberately unchanged, and this was measured

Screening reaches the `SurfEp` family: 9 of the 54 default features, and 38 of the 105 under `--full-features`, the extra 19 being the `SurfEp*Average` columns computed from partial charges. The nine `ShellEp` features are computed by `geometry.map_ep_to_plane_batch`, which divides the path at the molecular surface and weights the solvent leg with a relative permittivity of 80 against 4 for the protein leg. A distant charge is therefore already damped about twentyfold relative to the surface calculation, and the monopole never dominates the way it did there.

That was checked rather than assumed. A shell value is not a field sample: for each of 120 planes it sums, over every charged atom, the potential from that atom to its own orthogonal projection onto the plane. The APBS equivalent is therefore the APBS field sampled at those same projected points and summed the same way. On 1BSQ:

| | Spearman against APBS | Pearson | positive planes |
| --- | --- | --- | --- |
| shell as shipped, unscreened | **0.877** | 0.898 | 22 of 120 |
| shell with the solvent leg screened | 0.855 | 0.878 | 31 of 120 |
| APBS equivalent | | | 28 of 120 |

The shell was already in good agreement, unlike the surface. Screening it moved the rank correlation slightly the wrong way, improved the positive plane count slightly (22 to 31 against 28) and sign agreement from 92 to 94 per cent. Nothing in that justifies changing a released feature, so the shell is left alone.

This is the opposite of what was expected. The surface and the shell carry symmetric names and look like the same quantity measured at two distances, but they are computed by different models, and only one of them had the defect.

A caveat on scope: this was measured on one structure, because only one APBS grid was retained. It is enough to refute the assumption that the shell has the same defect; it is not enough to prove the shell is well behaved everywhere.

## Reproducing

The comparison script and the sampled reference data live on `analysis/apbs-electrostatics-comparison`. With that branch checked out and the environment from `environment_apbs.yml` available, `python scripts/apbs_comparison.py` regenerates the reference; the numbers above then follow from sampling the shipped potential at the same coordinates.
