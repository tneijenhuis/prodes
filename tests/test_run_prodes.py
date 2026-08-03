"""Tests for the prodes.run_prodes convenience wrapper.

run_prodes is a thin wrapper around prodes.run.calculate. Its failure mode is
silent: an option added to calculate but not to the wrapper is simply dropped,
and the caller gets default behaviour with no error. That is what happened to
full_features and mem_limit_mb. These tests pin the forwarding so it cannot
happen again.
"""

import inspect

import pandas as pd

import prodes
from prodes import run as run_module
from prodes.feature_dictionary import FULL_FEATURE_CODES, REDUCED_FEATURE_CODES

PDB_PATH = "tests/data/ARH96693.pdb.zip"


def test_signature_covers_every_calculate_option():
    """Every parameter calculate accepts must also be accepted by run_prodes."""
    wrapper = set(inspect.signature(prodes.run_prodes).parameters)
    target = set(inspect.signature(run_module.calculate).parameters)

    assert target - wrapper == set(), f"run_prodes is missing calculate options: {sorted(target - wrapper)}"


def test_defaults_match_calculate():
    """A default-valued option must not mean something different via the wrapper."""
    wrapper = inspect.signature(prodes.run_prodes).parameters
    target = inspect.signature(run_module.calculate).parameters

    differing = {
        name: (wrapper[name].default, target[name].default)
        for name in target
        if name in wrapper and wrapper[name].default is not inspect.Parameter.empty and wrapper[name].default != target[name].default
    }

    assert differing == {}, f"defaults disagree between run_prodes and calculate: {differing}"


def test_arguments_are_forwarded(monkeypatch):
    """Values reach calculate unchanged, whatever their position in the signature."""
    captured = {}

    def fake_calculate(pdb_file, out_file, **kwargs):
        captured["pdb_file"] = pdb_file
        captured["out_file"] = out_file
        captured.update(kwargs)

    monkeypatch.setattr(prodes, "calculate", fake_calculate)

    prodes.run_prodes(
        "in.pdb",
        "out.csv",
        pkas_file="x.pka",
        ph=5.5,
        r_probe=1.6,
        hydro_scale="kd",
        full_features=True,
        mem_limit_mb=512,
    )

    assert captured == {
        "pdb_file": "in.pdb",
        "out_file": "out.csv",
        "pkas_file": "x.pka",
        "ph": 5.5,
        "r_probe": 1.6,
        "hydro_scale": "kd",
        "full_features": True,
        "mem_limit_mb": 512,
    }


def test_full_features_reaches_the_output(tmp_path):
    """End to end: the wrapper can actually produce the full feature set."""
    tmpdir = tmp_path

    reduced_out = str(tmpdir / "reduced.csv")
    prodes.run_prodes(PDB_PATH, reduced_out)

    full_out = str(tmpdir / "full.csv")
    prodes.run_prodes(PDB_PATH, full_out, full_features=True)

    reduced_columns = [c for c in pd.read_csv(reduced_out).columns if c != "ID"]
    full_columns = [c for c in pd.read_csv(full_out).columns if c != "ID"]

    assert reduced_columns == list(REDUCED_FEATURE_CODES)
    # As a set: FULL_FEATURE_CODES is reduced-then-dropped, not CSV column order.
    assert set(full_columns) == set(FULL_FEATURE_CODES)
    assert len(full_columns) == len(FULL_FEATURE_CODES)
