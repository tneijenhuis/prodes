"""One cached pipeline run per structure and feature set, shared across test modules.

calculate() on the regression structure takes several seconds, and two modules
need both a full and a reduced run of the same file. Running them per module
meant four runs at collection time for two distinct results, and each module
left a temporary directory behind.
"""

import functools
import tempfile
from pathlib import Path

import pandas as pd

from prodes.output import read_features
from prodes.run import DEFAULT_IONIC_STRENGTH_MOLAR, calculate


@functools.cache
def _cached_output(pdb_path: str, full_features: bool, pkas_file: str | None, ionic_strength_molar: float) -> pd.DataFrame:
    """Runs the pipeline once and keeps the result for the rest of the session."""

    with tempfile.TemporaryDirectory(prefix="prodes_tests_") as directory:
        out_file = str(Path(directory) / "bundle.zip")
        calculate(pdb_path, out_file, pkas_file=pkas_file, full_features=full_features, ionic_strength_molar=ionic_strength_molar)

        return read_features(out_file)


def pipeline_output(
    pdb_path: str, full_features: bool, pkas_file: str | None = None, ionic_strength_molar: float = DEFAULT_IONIC_STRENGTH_MOLAR
) -> pd.DataFrame:
    """Returns what calculate() writes for this structure and feature set.

    The first call runs the pipeline; every later call with the same arguments
    gets a copy of that result, so no module can disturb another module's frame.

    Args:
        pdb_path: structure to run, as passed to calculate().
        full_features: True for the full legacy set, False for the reduced one.
        pkas_file: path to a pKa JSON to override the built-in values, or None
            to use them. Cached separately, so a run with pKas and a run
            without are two results rather than one overwriting the other.
        ionic_strength_molar: buffer ionic strength in mol/L. Defaults to what
            the package ships, so tests exercise shipped behaviour. Pass 0 to
            get the unscreened potential of versions before 5.0, which is what
            the regression fixtures were generated with.
    """

    return _cached_output(pdb_path, full_features, pkas_file, ionic_strength_molar).copy()
