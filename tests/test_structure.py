import pytest

from prodes.core.structure import Structure
from prodes.io.parser import PDBparser


structure = PDBparser().parse("tests/data/1GDW_h.pdb")

def test_empty_structure():
    """tests a empty structure"""
    empty_struct = Structure()
    assert empty_struct.name == None


def test_structure_coordinates():
    """tests if coordinates are computed correctly"""

    assert structure.x == 13.243374755381597
    assert structure.y == 14.928170254403136
    assert structure.z == 27.924622309197638


def test_residues():
    assert len(structure.residues) == 130


def test_atoms():
    assert len(structure.atoms) == 2003


def test_mass():
    """tests the molecular weight of a protein"""

    assert structure.mw == 14601.53


def test_isoelectric_point():
    """tests the isoelectic point of a protein"""
    assert structure.isoelectric_point() == 8.899


def test_heavy_atoms():
    """tests if the correct """
    assert len(structure.heavy_atoms) == 1022


def test_furthest():
    """Tests if the correct atom is returned as furthest"""

    assert structure.furthest_heavy_atom == structure.heavy_atoms[395]


def test_surface_area():
    """Tests if surface area is correct"""
    with pytest.raises(NameError):
        structure.surface_area()

def test_charge():

    assert round(structure.charge(7), 3) == 7.0
    assert round(structure.charge(11), 3) == -13


def test_dipole():

    assert round(structure.dipole(7), 3) == 162.053
