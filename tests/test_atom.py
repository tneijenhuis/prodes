import pytest
from prodes.core.atom import Atom


def test_empty_atom():
    """tests if a empty atom fails"""
    with pytest.raises(TypeError):
        _ = Atom()


def test_atom():
    """tests if a atom is loaded correctly"""

    atom = Atom("ATOM", "NZ", "LYS", "A", 1, 1, 1, 1)
    assert atom.identifier == "ATOM"
    assert atom.name == "NZ"
    assert atom.residue_name == "LYS"
    assert atom.chain_name == "A"
    assert atom.residue_number == 1
    assert atom.x == 1
    assert atom.y == 1
    assert atom.z == 1

    global c_atom, cr_atom, n_atom, o_atom, s_atom

    c_atom = Atom("ATOM", "CA", "LYS", "A", 1, 1, 1, 1)
    cr_atom = Atom("ATOM", "CG", "PHE", "A", 1, 1, 1, 1)
    n_atom = Atom("ATOM", "NZ", "LYS", "A", 1, 1, 1, 1)
    o_atom = Atom("ATOM", "O", "LYS", "A", 1, 1, 1, 1)
    s_atom = Atom("ATOM", "S", "SER", "A", 1, 1, 1, 1)


def test_atom_radius():
    "tests if the correct vdw radius is fetched"

    assert c_atom.radius == 2.0
    assert cr_atom.radius == 1.85
    assert n_atom.radius == 1.5
    assert o_atom.radius == 1.4
    assert s_atom.radius == 1.85


def test_max_area():
    """tests if the correct maximum area is calculated"""

    assert c_atom.max_area() == 145.267244301992
    assert cr_atom.max_area() == 132.73228961416876
    assert n_atom.max_area() == 105.68317686676065
    assert o_atom.max_area() == 98.5203456165759
    assert s_atom.max_area() == 132.73228961416876
# Try a different probe radius
    assert c_atom.max_area(2) == 201.06192982974676


def test_max_n_points():
    """tests if the correct maximum area is calculated"""

    assert c_atom.max_n_points() == 290
    assert cr_atom.max_n_points() == 265
    assert n_atom.max_n_points() == 211
    assert o_atom.max_n_points() == 197
    assert s_atom.max_n_points() == 265
# Try a different probe radius
    assert c_atom.max_n_points(2) == 402


def test_surface_area():
    """tests if surface area returns an error if no surface has been calculated"""

    from prodes.core.structure import Structure
    c_atom.structure = Structure()
    with pytest.raises(NameError):
        c_atom.surface_area()
        cr_atom.surface_area()
        
        