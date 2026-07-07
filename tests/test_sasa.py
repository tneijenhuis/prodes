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
    """Run the full calculate() pipeline and return the output DataFrame."""
    out_file = str(tmp_path_factory.mktemp("prodes") / "output.csv")
    calculate(PDB_PATH, out_file)
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
