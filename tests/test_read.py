"""Tests for the prodes.read entry point.

read() dispatches on the file extension and is the documented way to load a
structure or a pKa file. The pKa branch had no test and no caller, and passed no
argument to read_pka, so every .pka file raised a TypeError instead of loading.
"""

import json

import pytest

import prodes

PDB_PATH = "tests/data/1GDW.pdb.zip"


def test_read_loads_a_zipped_structure():
    """A .pdb.zip is parsed straight through, without unpacking it first."""

    structure = prodes.read(PDB_PATH)

    assert structure.name == "1GDW"
    assert len(structure.atoms) == 1022


def test_read_loads_a_pka_file(tmp_path):
    """A .pka file returns the residue to pKa mapping it holds."""

    pka_file = tmp_path / "structure.pka"
    pka_file.write_text(json.dumps({"5": [{"N+": 7.541}], "8": [{"ARG": 14.0}]}))

    assert prodes.read(str(pka_file)) == {5: [{"N+": 7.541}], 8: [{"ARG": 14.0}]}


def test_read_rejects_an_unsupported_extension(tmp_path):
    """Anything else is refused by name rather than guessed at."""

    with pytest.raises(ValueError, match="File extention not recognized"):
        prodes.read(str(tmp_path / "structure.cif"))
