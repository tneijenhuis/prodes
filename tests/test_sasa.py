"""Full-pipeline regression tests for the Prodes calculate() function.

How this file works:
  1. Run calculate() on a small test PDB (ARH96693) twice — once with the full
     105-feature set and once with the reduced non-redundant set.
  2. Load a saved reference CSV (ARH96693_prodes_orig_output.csv) that was
     produced by the original, pre-refactoring code.
  3. Compare the freshly calculated output against the reference to catch
     regressions in any calculation (SASA, electrostatic potential,
     lipophilicity, shell features, etc.).

The expensive calculate() calls happen once at module load time (not per test).
Each test function then just reads from the pre-computed DataFrames below.
"""

import tempfile
from pathlib import Path

import pandas as pd

from prodes.run import calculate

PDB_PATH = "tests/data/ARH96693.pdb"
ORIG_OUTPUT = "tests/data/ARH96693_prodes_orig_output.csv"

# Columns that hold integer counts (not floats) — given a looser tolerance.
COUNT_COLUMNS = {
    "NSurfPoints",
    "NSurfPosEpFormal",
    "NSurfNegEpFormal",
    "NSurfPosMhp",
    "NSurfNegMhp",
    "NSurfPosEpAverage",
    "NSurfNegEpAverage",
    "NShellPosEpFormal",
    "NShellNegEpFormal",
}

# --- Run the pipeline once at import time ---
# These are plain variables, not pytest fixtures. calculate() is expensive
# (SASA + 120 shell planes + electrostatic/lipophilic projections), so we
# run it once here and reuse the DataFrames in every test below.

_tmpdir = Path(tempfile.mkdtemp(prefix="prodes_test_"))

_full_out = str(_tmpdir / "full_output.csv")
calculate(PDB_PATH, _full_out, full_features=True)
FULL_OUTPUT = pd.read_csv(_full_out)

_reduced_out = str(_tmpdir / "reduced_output.csv")
calculate(PDB_PATH, _reduced_out, full_features=False)
REDUCED_OUTPUT = pd.read_csv(_reduced_out)

ORIGINAL_OUTPUT = pd.read_csv(ORIG_OUTPUT)


def test_all_columns_present():
    """Every column in the original output should appear in the calculated output."""
    calc_cols = set(FULL_OUTPUT.columns)
    orig_cols = set(ORIGINAL_OUTPUT.columns)
    assert calc_cols == orig_cols, f"Column mismatch: missing={orig_cols - calc_cols}, extra={calc_cols - orig_cols}"


def test_feature_values_match_original():
    """All feature values should match the original output within tolerance."""
    calc = FULL_OUTPUT.iloc[0].drop("ID")
    orig = ORIGINAL_OUTPUT.iloc[0].drop("ID")

    mismatches = []
    for col in orig.index:
        calc_val = float(calc[col])
        orig_val = float(orig[col])

        if col in COUNT_COLUMNS:
            tolerance = 5
        elif "Shell" in col:
            tolerance = 0.05
        else:
            tolerance = 0.01

        if abs(calc_val - orig_val) > tolerance:
            mismatches.append(f"{col}: orig={orig_val}, calc={calc_val}, " f"diff={abs(calc_val - orig_val):.6f}, tol={tolerance}")

    assert len(mismatches) == 0, f"Value mismatches for {len(mismatches)} columns:\n" + "\n".join(mismatches[:20])


# Features that should be present in both full and reduced sets
REDUCED_PRESENT = {
    "Molecular weight",
    "Isoelectric point",
    "Dipole",
    "Formal charge",
    "Area",
    "Shape max",
    "Shape min",
    "SurfEpMaxFormal",
    "SurfEpMinFormal",
    "SurfEpMeanFormal",
    "SurfEpStdFormal",
    "NSurfPosEpFormal",
    "NSurfNegEpFormal",
    "SurfEpPosMeanFormal",
    "SurfEpPosStdFormal",
    "SurfEpNegMeanFormal",
    "SurfEpNegStdFormal",
    "SurfMhpMax",
    "SurfMhpMin",
    "SurfMhpMean",
    "SurfMhpStd",
    "NSurfPosMhp",
    "NSurfNegMhp",
    "SurfPosMhpMean",
    "SurfPosMhpStd",
    "SurfNegMhpMean",
    "SurfNegMhpStd",
    "ShellEpMaxFormal",
    "ShellEpminFormal",
    "ShellEpMeanFormal",
    "ShellEpStdFormal",
    "NShellPosEpFormal",
    "ShellEpPosMeanFormal",
    "ShellEpPosStdFormal",
    "ShellEpNegMeanFormal",
    "ShellEpNegStdFormal",
}

# Features dropped in reduced mode (redundant with others at R² >= 0.95)
REDUCED_ABSENT = {
    "Average charge",
    "NSurfPoints",
    "SurfEpTrimeanFormal",
    "SurfEpMedianFormal",
    "SurfEpSumFormal",
    "SurfEpPosTrimeanFormal",
    "SurfEpPosMedianFormal",
    "SurfEpPosSumFormal",
    "SurfEpNegTrimeanFormal",
    "SurfEpNegMedianFormal",
    "SurfEpNegSumFormal",
    "SurfMhpTrimean",
    "SurfMhpMedian",
    "SurfMhpSum",
    "SurfPosMhpTrimean",
    "SurfPosMhpMedian",
    "SurfPosMhpSum",
    "SurfNegMhpTrimean",
    "SurfNegMhpMedian",
    "SurfNegMhpSum",
    "ShellEpTrimeanFormal",
    "ShellEpMedianFormal",
    "ShellEpSumFormal",
    "ShellEpPosTrimeanFormal",
    "ShellEpPosMedianFormal",
    "ShellEpPosSumFormal",
    "ShellEpNegTrimeanFormal",
    "ShellEpNegMedianFormal",
    "ShellEpNegSumFormal",
    "NShellNegEpFormal",
    # Entire Average-charge surface EP block is skipped
    "SurfEpMaxAverage",
    "SurfEpMinAverage",
    "SurfEpMeanAverage",
    "SurfEpTrimeanAverage",
    "SurfEpMedianAverage",
    "SurfEpSumAverage",
    "SurfEpStdAverage",
    "NSurfPosEpAverage",
    "SurfEpPosMeanAverage",
    "SurfEpPosTrimeanAverage",
    "SurfEpPosMedianAverage",
    "SurfEpPosSumAverage",
    "SurfEpPosStdAverage",
    "NSurfNegEpAverage",
    "SurfEpNegMeanAverage",
    "SurfEpNegTrimeanAverage",
    "SurfEpNegMedianAverage",
    "SurfEpNegSumAverage",
    "SurfEpNegStdAverage",
}


def test_reduced_has_fewer_columns():
    """The reduced feature set should have substantially fewer columns than the full set."""
    assert len(REDUCED_OUTPUT.columns) < len(FULL_OUTPUT.columns)


def test_reduced_contains_core_features():
    """Key non-redundant features should still be present in the reduced set."""
    cols = set(REDUCED_OUTPUT.columns)
    missing = REDUCED_PRESENT - cols
    assert missing == set(), f"Missing core features in reduced set: {missing}"


def test_reduced_drops_redundant_features():
    """Redundant features should be absent from the reduced set."""
    cols = set(REDUCED_OUTPUT.columns)
    present = REDUCED_ABSENT & cols
    assert present == set(), f"Redundant features still present in reduced set: {present}"
