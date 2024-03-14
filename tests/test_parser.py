import pytest
from prodes.io import parser


file_path = "tests/data/1GDW.pdb"
pdb_parser = parser.PDBparser()
pdb_parser.identifier = "ATOM"


def test_main_fail():
    """tests if a incorect file will give an error"""

    with pytest.raises(ValueError):
        pdb_parser.parse(file_path[:-2])


def test_main_accept():
    """tests if a file is parsed correctly"""

    global structure
    structure = pdb_parser.parse(file_path)

    assert len(structure.atoms) == 1022
    assert len(structure.chains) == 1
    assert len(structure.residues) == 130


def test_read():
    """tests if a pdb file read"""

    with open(file_path) as file:
        pdb_parser._read_pdb(file, "test")


def test_parsed_atom():
    """tests parsed atom"""

    atom = structure.atoms[0]

    assert atom.name == "N"
    assert atom.identifier == "ATOM"
    assert atom.chain_name == "A"
    assert atom.residue_number == 1
    assert atom.residue_name == "LYS"
    assert atom.x == 1.134
    assert atom.y == 19.824
    assert atom.z == 22.575
    assert atom.structure == structure


def test_parsed_chain():
    """tests parsed chain"""

    chain = structure.chains[0]
    assert chain.name == "A"


def test_parsed_residue():
    """tests parsed residue"""

    residue = structure.residues[0]

    assert residue.name == "LYS"
    assert residue.number == 1
    assert len(residue.atoms) == 9
