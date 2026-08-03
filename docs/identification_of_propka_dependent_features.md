# PROPKA dependency of Prodes surface features

AlphaFold structures compared (present in both runs): **819**
Features analysed (shared by both runs): **56**
PROPKA-dependent: **22**   PROPKA-independent: **34**

## PROPKA-independent features (unchanged for every structure)

- `Molecular weight`
- `Area`
- `ALASurfFrac`
- `ARGSurfFrac`
- `ASNSurfFrac`
- `ASPSurfFrac`
- `CYSSurfFrac`
- `GLNSurfFrac`
- `GLUSurfFrac`
- `GLYSurfFrac`
- `HISSurfFrac`
- `ILESurfFrac`
- `LEUSurfFrac`
- `LYSSurfFrac`
- `METSurfFrac`
- `PHESurfFrac`
- `PROSurfFrac`
- `SERSurfFrac`
- `THRSurfFrac`
- `TRPSurfFrac`
- `TYRSurfFrac`
- `VALSurfFrac`
- `Shape max`
- `Shape min`
- `SurfMhpMax`
- `SurfMhpMin`
- `SurfMhpMean`
- `SurfMhpStd`
- `NSurfPosMhp`
- `SurfPosMhpMean`
- `SurfPosMhpStd`
- `NSurfNegMhp`
- `SurfNegMhpMean`
- `SurfNegMhpStd`

## PROPKA-dependent features

- `Isoelectric point`
- `Dipole`
- `Formal charge`
- `SurfEpMaxFormal`
- `SurfEpMinFormal`
- `SurfEpMeanFormal`
- `SurfEpStdFormal`
- `NSurfPosEpFormal`
- `SurfEpPosMeanFormal`
- `SurfEpPosStdFormal`
- `NSurfNegEpFormal`
- `SurfEpNegMeanFormal`
- `SurfEpNegStdFormal`
- `ShellEpMaxFormal`
- `ShellEpminFormal`
- `ShellEpMeanFormal`
- `ShellEpStdFormal`
- `NShellPosEpFormal`
- `ShellEpPosMeanFormal`
- `ShellEpPosStdFormal`
- `ShellEpNegMeanFormal`
- `ShellEpNegStdFormal`

## R^2 between the two runs

Sorted most PROPKA-sensitive first. Low R^2 means PROPKA changes the feature a lot;
R^2 near 1 with a non-zero `% differ` means it shifts the values but preserves the ranking.

| Feature | R^2 | slope | % of structures differing | mean abs diff | max abs diff | median % diff | n | strength |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `SurfEpMinFormal` | 0.8177 | 0.801 | 72.4 | 0.3351 | 4.2100 | 4.74 | 819 | strong |
| `ShellEpPosStdFormal` | 0.9056 | 0.927 | 60.4 | 0.0053 | 0.0530 | 4.76 | 819 | moderate |
| `ShellEpNegMeanFormal` | 0.9227 | 0.892 | 74.4 | 0.0214 | 0.1680 | 7.59 | 819 | moderate |
| `ShellEpNegStdFormal` | 0.9243 | 0.930 | 72.3 | 0.0072 | 0.0640 | 4.55 | 819 | moderate |
| `ShellEpPosMeanFormal` | 0.9249 | 0.945 | 65.4 | 0.0081 | 0.0840 | 7.62 | 819 | moderate |
| `ShellEpminFormal` | 0.9254 | 0.916 | 73.7 | 0.0391 | 0.2670 | 5.09 | 819 | moderate |
| `NShellPosEpFormal` | 0.9329 | 0.989 | 61.2 | 4.1416 | 48.0000 | 11.76 | 819 | moderate |
| `SurfEpNegStdFormal` | 0.9365 | 0.943 | 74.6 | 0.0353 | 0.3350 | 3.54 | 819 | moderate |
| `ShellEpMeanFormal` | 0.9383 | 0.917 | 74.5 | 0.0272 | 0.1740 | 12.58 | 819 | moderate |
| `SurfEpNegMeanFormal` | 0.9429 | 0.915 | 75.2 | 0.1538 | 1.0040 | 8.64 | 819 | moderate |
| `Isoelectric point` | 0.9448 | 0.983 | 99.9 | 0.2368 | 1.3940 | 3.46 | 819 | moderate |
| `ShellEpMaxFormal` | 0.9465 | 0.938 | 72.8 | 0.0238 | 0.1460 | 9.87 | 819 | moderate |
| `NSurfPosEpFormal` | 0.9522 | 1.056 | 67.8 | 520.3639 | 8279.0000 | 25.45 | 819 | moderate |
| `ShellEpStdFormal` | 0.9522 | 0.948 | 70.0 | 0.0048 | 0.0530 | 2.38 | 819 | moderate |
| `SurfEpMeanFormal` | 0.9529 | 0.939 | 75.2 | 0.1855 | 1.0040 | 11.48 | 819 | moderate |
| `SurfEpStdFormal` | 0.9570 | 0.950 | 74.4 | 0.0252 | 0.2600 | 2.10 | 819 | moderate |
| `SurfEpMaxFormal` | 0.9581 | 0.948 | 73.5 | 0.1558 | 2.0200 | 9.68 | 819 | moderate |
| `SurfEpPosStdFormal` | 0.9586 | 0.965 | 66.8 | 0.0260 | 0.2560 | 6.78 | 819 | moderate |
| `Formal charge` | 0.9701 | 0.946 | 65.7 | 1.3578 | 10.0000 | 11.11 | 819 | moderate |
| `SurfEpPosMeanFormal` | 0.9707 | 0.977 | 67.9 | 0.0357 | 0.2560 | 6.76 | 819 | moderate |
| `NSurfNegEpFormal` | 0.9917 | 0.993 | 68.0 | 527.5201 | 8296.0000 | 0.85 | 819 | weak |
| `Dipole` | 0.9919 | 0.984 | 75.6 | 49.6986 | 555.9310 | 3.51 | 819 | weak |

## Caveats

- Only a reduced feature set of 56 features was tested.
- Prodes rounds its output (mostly 3 dp), so a sub-0.001 PROPKA sensitivity reads as independent.
