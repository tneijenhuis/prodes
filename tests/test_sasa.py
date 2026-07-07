import pandas as pd
import pytest

from prodes.run import calculate

PDB_PATH = "tests/data/ARH96693.pdb"
ORIG_OUTPUT = "tests/data/ARH96693_prodes_orig_output.csv"

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


@pytest.fixture(scope="module")
def calculated_output(tmp_path_factory):
    """Run the full calculate() pipeline with full_features=True and return the output DataFrame."""
    out_file = str(tmp_path_factory.mktemp("prodes") / "output.csv")
    calculate(PDB_PATH, out_file, full_features=True)
    return pd.read_csv(out_file)


@pytest.fixture(scope="module")
def reduced_output(tmp_path_factory):
    """Run the calculate() pipeline with the default reduced feature set."""
    out_file = str(tmp_path_factory.mktemp("prodes") / "reduced_output.csv")
    calculate(PDB_PATH, out_file, full_features=False)
    return pd.read_csv(out_file)


@pytest.fixture(scope="module")
def original_output():
    """Load the reference output CSV from the original unrefactored code."""
    return pd.read_csv(ORIG_OUTPUT)


def test_all_columns_present(calculated_output, original_output):
    """Every column in the original output should appear in the calculated output."""
    calc_cols = set(calculated_output.columns)
    orig_cols = set(original_output.columns)
    assert calc_cols == orig_cols, (
        f"Column mismatch: missing={orig_cols - calc_cols}, extra={calc_cols - orig_cols}"
    )


def test_feature_values_match_original(calculated_output, original_output):
    """All feature values should match the original output within tolerance."""
    calc = calculated_output.iloc[0].drop("ID")
    orig = original_output.iloc[0].drop("ID")

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
            mismatches.append(
                f"{col}: orig={orig_val}, calc={calc_val}, "
                f"diff={abs(calc_val - orig_val):.6f}, tol={tolerance}"
            )

    assert len(mismatches) == 0, (
        f"Value mismatches for {len(mismatches)} columns:\n"
        + "\n".join(mismatches[:20])
    )


# Features that should be present in both full and reduced sets
REDUCED_PRESENT = {
    "Molecular weight", "Isoelectric point", "Dipole", "Formal charge",
    "Area", "Shape max", "Shape min",
    "SurfEpMaxFormal", "SurfEpMinFormal",
    "SurfEpMeanFormal", "SurfEpStdFormal",
    "NSurfPosEpFormal", "NSurfNegEpFormal",
    "SurfEpPosMeanFormal", "SurfEpPosStdFormal",
    "SurfEpNegMeanFormal", "SurfEpNegStdFormal",
    "SurfMhpMax", "SurfMhpMin", "SurfMhpMean", "SurfMhpStd",
    "NSurfPosMhp", "NSurfNegMhp",
    "SurfPosMhpMean", "SurfPosMhpStd",
    "SurfNegMhpMean", "SurfNegMhpStd",
    "ShellEpMaxFormal", "ShellEpminFormal",
    "ShellEpMeanFormal", "ShellEpStdFormal",
    "NShellPosEpFormal",
    "ShellEpPosMeanFormal", "ShellEpPosStdFormal",
    "ShellEpNegMeanFormal", "ShellEpNegStdFormal",
}

# Features dropped in reduced mode (redundant with others at R² >= 0.95)
REDUCED_ABSENT = {
    "Average charge", "NSurfPoints",
    "SurfEpTrimeanFormal", "SurfEpMedianFormal", "SurfEpSumFormal",
    "SurfEpPosTrimeanFormal", "SurfEpPosMedianFormal", "SurfEpPosSumFormal",
    "SurfEpNegTrimeanFormal", "SurfEpNegMedianFormal", "SurfEpNegSumFormal",
    "SurfMhpTrimean", "SurfMhpMedian", "SurfMhpSum",
    "SurfPosMhpTrimean", "SurfPosMhpMedian", "SurfPosMhpSum",
    "SurfNegMhpTrimean", "SurfNegMhpMedian", "SurfNegMhpSum",
    "ShellEpTrimeanFormal", "ShellEpMedianFormal", "ShellEpSumFormal",
    "ShellEpPosTrimeanFormal", "ShellEpPosMedianFormal", "ShellEpPosSumFormal",
    "ShellEpNegTrimeanFormal", "ShellEpNegMedianFormal", "ShellEpNegSumFormal",
    "NShellNegEpFormal",
    # Entire Average-charge surface EP block is skipped
    "SurfEpMaxAverage", "SurfEpMinAverage",
    "SurfEpMeanAverage", "SurfEpTrimeanAverage", "SurfEpMedianAverage",
    "SurfEpSumAverage", "SurfEpStdAverage",
    "NSurfPosEpAverage", "SurfEpPosMeanAverage", "SurfEpPosTrimeanAverage",
    "SurfEpPosMedianAverage", "SurfEpPosSumAverage", "SurfEpPosStdAverage",
    "NSurfNegEpAverage", "SurfEpNegMeanAverage", "SurfEpNegTrimeanAverage",
    "SurfEpNegMedianAverage", "SurfEpNegSumAverage", "SurfEpNegStdAverage",
}


def test_reduced_has_fewer_columns(reduced_output, calculated_output):
    """The reduced feature set should have substantially fewer columns than the full set."""
    assert len(reduced_output.columns) < len(calculated_output.columns)


def test_reduced_contains_core_features(reduced_output):
    """Key non-redundant features should still be present in the reduced set."""
    cols = set(reduced_output.columns)
    missing = REDUCED_PRESENT - cols
    assert missing == set(), f"Missing core features in reduced set: {missing}"


def test_reduced_drops_redundant_features(reduced_output):
    """Redundant features should be absent from the reduced set."""
    cols = set(reduced_output.columns)
    present = REDUCED_ABSENT & cols
    assert present == set(), f"Redundant features still present in reduced set: {present}"
