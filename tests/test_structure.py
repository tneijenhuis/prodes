from collections import Counter

import pytest

from prodes.core.structure import Structure
from prodes.io.parser import PDBparser

structure = PDBparser().parse("tests/data/1GDW_h.pdb.zip")


def test_empty_structure():
    """tests a empty structure"""
    empty_struct = Structure()
    assert empty_struct.name is None


def test_structure_coordinates():
    """tests if coordinates are computed correctly"""

    assert structure.x == 13.243374755381597
    assert structure.y == 14.928170254403136
    assert structure.z == 27.924622309197638


def test_residues():
    assert len(structure.residues) == 130


def test_disulfides_are_found_on_a_parsed_structure():
    """Every structure the parser returns already knows which cysteines are bonded."""

    assert len(structure.disulfides) == 4
    assert {(first.number, second.number) for first, second in structure.disulfides} == {(77, 95), (6, 128), (30, 116), (65, 81)}


def test_atoms():
    assert len(structure.atoms) == 2003


def test_mass():
    """tests the molecular weight of a protein"""

    assert structure.mw == 14601.53


def test_isoelectric_point():
    """tests the isoelectic point of a protein

    1GDW has eight cysteines and all of them are in disulfides. Before those
    were recognised, each was titrated as a free thiol at pKa 8.33 and dragged
    the isoelectric point down by 1.4 units, from 10.342 to 8.899.
    """
    assert structure.isoelectric_point() == 10.342


def test_heavy_atoms():
    """tests if the correct"""
    assert len(structure.heavy_atoms) == 1022


def test_furthest():
    """Tests if the correct atom is returned as furthest"""

    assert structure.furthest_heavy_atom == structure.heavy_atoms[395]


def test_surface_area():
    """Tests if surface area is correct"""
    with pytest.raises(NameError):
        structure.surface_area()


def test_charge():
    """tests the net formal charge at a pH either side of the cysteine pKa

    At pH 7 nothing turns on the cysteines: a free thiol at pKa 8.33 is
    protonated there either way. At pH 11 all eight would have been counted as
    deprotonated, which is where the old -13 came from. All eight are in
    disulfides, so the correct answer is 8 e less negative.
    """

    assert round(structure.charge(7), 3) == 7.0
    assert round(structure.charge(11), 3) == -5


def test_dipole():

    assert round(structure.dipole(7), 3) == 162.053


def test_titratable_groups_excludes_bonded_cysteines():
    """What a pKa predictor could legitimately have an opinion about.

    1GDW has eight cysteines and every one is in a disulfide, so none of them is
    a titratable group. This is what the missing-value warning in redo_pkas
    counts against, so a bonded cysteine appearing here would make a complete
    pKa file look incomplete.
    """

    groups = structure.titratable_groups()
    counts = Counter(key for _, key in groups)

    assert len(groups) == 38
    assert counts == {"ARG": 13, "ASP": 8, "TYR": 6, "LYS": 5, "GLU": 3, "HIS": 1, "N+": 1, "C-": 1}
    assert "CYS" not in counts
