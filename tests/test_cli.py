"""Tests for the command line entry point of prodes.run.

parse_arguments is what the CLI user actually meets. It had no tests, and two of
the settings it now covers, the worker count and the chunk size, are read from
the environment deep inside the first calculation phase, so a typo used to
surface a phase late rather than at startup.
"""

import inspect
import os

import pytest

from prodes.run import DEFAULT_IONIC_STRENGTH_MOLAR, calculate, main, parse_arguments


def parse(monkeypatch, *arguments):
    """Runs parse_arguments over a fabricated command line."""

    monkeypatch.setattr("sys.argv", ["prodes", "in.pdb", "out.zip", *arguments])

    return parse_arguments()


def test_defaults_are_the_documented_ones(monkeypatch):
    """A bare invocation gives the reduced feature set and no explicit memory budget."""

    monkeypatch.delenv("PRODES_FULL_FEATURES", raising=False)
    pdb_file, out_file, pkas_file, ph, r_probe, hydro_scale, full_features, mem_limit_mb, ionic_strength_molar = parse(monkeypatch)

    assert (pdb_file, out_file, pkas_file) == ("in.pdb", "out.zip", None)
    assert (ph, r_probe, hydro_scale) == (7, 1.4, "mj_scaled")
    assert full_features is False
    assert mem_limit_mb is None
    assert ionic_strength_molar == DEFAULT_IONIC_STRENGTH_MOLAR


def test_n_workers_flag_sets_the_environment_variable(monkeypatch):
    """--n-workers reaches the workers, which read the setting from the environment."""

    monkeypatch.delenv("PRODES_N_WORKERS", raising=False)
    parse(monkeypatch, "--n-workers", "3")

    assert os.environ["PRODES_N_WORKERS"] == "3"


def test_chunksize_flag_sets_the_environment_variable(monkeypatch):
    """--chunksize likewise, since it is read at dispatch time."""

    monkeypatch.delenv("PRODES_CHUNKSIZE", raising=False)
    parse(monkeypatch, "--chunksize", "5")

    assert os.environ["PRODES_CHUNKSIZE"] == "5"


def test_mem_limit_flag_is_returned_rather_than_exported(monkeypatch):
    """--mem-limit is passed to calculate directly, so it does not need the environment."""

    *_, mem_limit_mb, _ionic_strength = parse(monkeypatch, "--mem-limit", "256")

    assert mem_limit_mb == 256


@pytest.mark.parametrize(
    ("variable", "value"),
    [("PRODES_N_WORKERS", "abc"), ("PRODES_CHUNKSIZE", "0"), ("PRODES_MEM_LIMIT_MB", "-1")],
)
def test_a_bad_setting_fails_at_parse_time(monkeypatch, variable, value):
    """The settings are validated while parsing, before any structure is read."""

    monkeypatch.setenv(variable, value)

    with pytest.raises(ValueError, match=variable):
        parse(monkeypatch)


def test_a_bad_mem_limit_flag_fails_at_parse_time(monkeypatch):
    """A budget given on the command line is checked the same way as the variable."""

    with pytest.raises(ValueError, match="PRODES_MEM_LIMIT_MB must be"):
        parse(monkeypatch, "--mem-limit", "0")


def test_main_forwards_every_parsed_argument(monkeypatch, tmp_path):
    """main passes what it parsed to calculate, in the order calculate expects."""

    captured = {}

    def fake_calculate(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs

    monkeypatch.setattr("prodes.run.calculate", fake_calculate)
    monkeypatch.setattr("sys.argv", ["prodes", str(tmp_path / "in.pdb"), str(tmp_path / "out.zip"), "--ph", "5.5", "--full-features"])

    main()

    # Every position is asserted, not a sample of three. main splats this tuple
    # into calculate positionally, so a parameter inserted in the middle would
    # silently shift every later argument into the wrong slot.
    assert captured["args"] == (
        str(tmp_path / "in.pdb"),
        str(tmp_path / "out.zip"),
        None,
        5.5,
        1.4,
        "mj_scaled",
        True,
        None,
        DEFAULT_IONIC_STRENGTH_MOLAR,
    )
    assert captured["kwargs"] == {}


def test_calculate_signature_matches_the_parsed_order():
    """main calls calculate positionally, so the two orders must not drift apart."""

    assert list(inspect.signature(calculate).parameters) == [
        "pdb_file",
        "out_file",
        "pkas_file",
        "ph",
        "r_probe",
        "hydro_scale",
        "full_features",
        "mem_limit_mb",
        "ionic_strength_molar",
    ]
