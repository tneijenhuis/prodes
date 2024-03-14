import pytest
from prodes.core.residue import Residue
from prodes.io import parser as ps

structure = ps.PDBparser().parse("tests/data/1GDW_h.pdb")
residues = {}
for residue in structure.residues[1:]:
    if residue.name not in residues:
        residues[residue.name] = residue
for name, res in residues.items():
    print(name, [atom.name for atom in res.atoms])

n_term = structure.residues[0]
c_term = structure.residues[-1]

def test_empty_residue():
    """tests if a empty residue fails"""
    with pytest.raises(TypeError):
        _ = Residue()


def test_attributes():
    "Tests if atributes are correct"

    res = residues["MET"]

    assert res.name == "MET"
    assert res.number == 17
    assert res.structure == structure
    assert res.terminus is None

    n_term = structure.residues[0]
    c_term = structure.residues[-1]

    assert n_term.terminus == "N"
    assert c_term.terminus == "C"


def test_heavy_atoms():
    """tests if correct heavy atoms are called"""
    
    ala_heavies = [
        'N', 'CA', 'C', 'O', 'CB'
    ]
    arg_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG',
        'CD', 'NE', 'CZ', 'NH1', 'NH2'
    ]
    asn_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND2', 'OD1'
    ]
    asp_heavies = [
       'N', 'CA', 'C', 'O', 'CB',
       'CG', 'OD1', 'OD2'
    ]
    cys_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'SG'
    ]
    gln_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE2', 'OE1'
    ]
    glu_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'
    ]
    gly_heavies = [
        'N', 'CA', 'C', 'O'
    ]
    his_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD2', 'ND1',
        'CE1', 'NE2'
    ]
    ile_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'
    ]
    leu_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'
    ]
    lys_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'
    ]
    met_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'
    ]
    phe_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2',
        'CE1', 'CE2', 'CZ'
    ]
    pro_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'
    ]
    ser_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'OG'
    ]
    thr_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG2', 'OG1'
    ]
    trp_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2',
        'CE2', 'CE3', 'NE1', 'CZ2', 'CZ3', 'CH2'
    ]
    tyr_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2',
        'CE1', 'CE2', 'CZ', 'OH'
    ]
    val_heavies = [
        'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'
    ]

    assert [atom.name for atom in residues["ALA"].heavy_atoms] == ala_heavies
    assert [atom.name for atom in residues["ARG"].heavy_atoms] == arg_heavies
    assert [atom.name for atom in residues["ASN"].heavy_atoms] == asn_heavies
    assert [atom.name for atom in residues["ASP"].heavy_atoms] == asp_heavies
    assert [atom.name for atom in residues["CYS"].heavy_atoms] == cys_heavies
    assert [atom.name for atom in residues["GLN"].heavy_atoms] == gln_heavies
    assert [atom.name for atom in residues["GLY"].heavy_atoms] == gly_heavies
    assert [atom.name for atom in residues["GLU"].heavy_atoms] == glu_heavies
    assert [atom.name for atom in residues["HIS"].heavy_atoms] == his_heavies
    assert [atom.name for atom in residues["ILE"].heavy_atoms] == ile_heavies
    assert [atom.name for atom in residues["LEU"].heavy_atoms] == leu_heavies
    assert [atom.name for atom in residues["LYS"].heavy_atoms] == lys_heavies
    assert [atom.name for atom in residues["MET"].heavy_atoms] == met_heavies
    assert [atom.name for atom in residues["PHE"].heavy_atoms] == phe_heavies
    assert [atom.name for atom in residues["PRO"].heavy_atoms] == pro_heavies
    assert [atom.name for atom in residues["SER"].heavy_atoms] == ser_heavies
    assert [atom.name for atom in residues["THR"].heavy_atoms] == thr_heavies
    assert [atom.name for atom in residues["TRP"].heavy_atoms] == trp_heavies
    assert [atom.name for atom in residues["TYR"].heavy_atoms] == tyr_heavies
    assert [atom.name for atom in residues["VAL"].heavy_atoms] == val_heavies

def test_protons():
    """Tests if correct protons are called"""
    
    ala_protons = [
        'H01', 'H02', 'H03', 'H04', 'H12'
    ]
    arg_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H09', 'H10', 'H11', 'H12', 'H13'
    ]
    asn_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06'
    ]
    asp_protons = [
       'H01', 'H02', 'H03', 'H10'
    ]
    cys_protons = [
        'H01', 'H02', 'H03', 'H14'
    ]
    gln_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07', 'H09'
    ]
    glu_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H09'
    ]
    gly_protons = [
        'H01', 'H02', 'H11'
    ]
    his_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07'
    ]
    ile_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H09', 'H10', 'H11'
    ]
    leu_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H09', 'H10', 'H11'
    ]
    lys_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H09', 'H10', 'H11', 'H12', 'H13'
    ]
    met_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H09'
    ]
    phe_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H09'
    ]
    pro_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07'
    ]
    ser_protons = [
        'H01', 'H02', 'H03', 'H04', 'H12'
    ]
    thr_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H14'
    ]
    trp_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H09', 'H10'
    ]
    tyr_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H09'
    ]
    val_protons = [
        'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07',
        'H08', 'H15'
    ]

    assert [atom.name for atom in residues["ALA"].protons] == ala_protons
    assert [atom.name for atom in residues["ARG"].protons] == arg_protons
    assert [atom.name for atom in residues["ASN"].protons] == asn_protons
    assert [atom.name for atom in residues["ASP"].protons] == asp_protons
    assert [atom.name for atom in residues["CYS"].protons] == cys_protons
    assert [atom.name for atom in residues["GLN"].protons] == gln_protons
    assert [atom.name for atom in residues["GLU"].protons] == glu_protons
    assert [atom.name for atom in residues["GLY"].protons] == gly_protons
    assert [atom.name for atom in residues["HIS"].protons] == his_protons
    assert [atom.name for atom in residues["ILE"].protons] == ile_protons
    assert [atom.name for atom in residues["LEU"].protons] == leu_protons
    assert [atom.name for atom in residues["LYS"].protons] == lys_protons
    assert [atom.name for atom in residues["MET"].protons] == met_protons
    assert [atom.name for atom in residues["PHE"].protons] == phe_protons
    assert [atom.name for atom in residues["PRO"].protons] == pro_protons
    assert [atom.name for atom in residues["SER"].protons] == ser_protons
    assert [atom.name for atom in residues["THR"].protons] == thr_protons
    assert [atom.name for atom in residues["TRP"].protons] == trp_protons
    assert [atom.name for atom in residues["TYR"].protons] == tyr_protons
    assert [atom.name for atom in residues["VAL"].protons] == val_protons


def test_pkas():
    """tests weither the correct pka is assigned"""

# testing side chain PKAs
    assert residues['ALA'].pkas is None
    assert residues['ARG'].pkas == [{"ARG": 13.8}]
    assert residues['ASN'].pkas is None 
    assert residues['ASP'].pkas == [{"ASP": 3.86}]
    assert residues['CYS'].pkas == [{"CYS": 8.33}]
    assert residues['GLN'].pkas is None
    assert residues['GLU'].pkas == [{"GLU": 4.25}]
    assert residues['GLY'].pkas is None
    assert residues['HIS'].pkas == [{"HIS": 6.00}]
    assert residues['ILE'].pkas is None
    assert residues['LYS'].pkas == [{"LYS": 10.5}]
    assert residues['MET'].pkas is None
    assert residues['PHE'].pkas is None
    assert residues['PRO'].pkas is None
    assert residues['SER'].pkas is None
    assert residues['THR'].pkas is None
    assert residues['TRP'].pkas is None
    assert residues['TYR'].pkas == [{"TYR": 10.0}]
    assert residues['VAL'].pkas is None

# testing the N and C terminus pka
    assert n_term.pkas == [{"LYS": 10.5}, {"N+": 9.69}]
    assert c_term.pkas == [{"C-": 2.34}]


def test_charge():
    """tests if the correct charge is returned"""

# no charge
    assert residues['ALA'].charge(7) == 0
    assert round(residues['ARG'].charge(14), 3) == 0
    assert residues['ASN'].charge(7) == 0 
    assert round(residues['ASP'].charge(3), 3) == 0
    assert round(residues['CYS'].charge(8), 3) == 0
    assert residues['GLN'].charge(7) == 0
    assert round(residues['GLU'].charge(4), 3) == 0
    assert residues['GLY'].charge(7) == 0
    assert round(residues['HIS'].charge(7), 3) == 0
    assert residues['ILE'].charge(7) == 0
    assert round(residues['LYS'].charge(11), 3) == 0
    assert residues['MET'].charge(7) == 0
    assert residues['PHE'].charge(7) == 0
    assert residues['PRO'].charge(7) == 0
    assert residues['SER'].charge(7) == 0
    assert residues['THR'].charge(7) == 0
    assert residues['TRP'].charge(7) == 0
    assert round(residues['TYR'].charge(7), 3) == 0
    assert residues['VAL'].charge(7) == 0
    assert round(n_term.charge(11), 3) == 0
    assert round(c_term.charge(1), 3)== -0
# Positive charge
    assert round(residues["ARG"].charge(12), 3) == 1
    assert round(residues["LYS"].charge(10), 3) ==1
    assert round(n_term.charge(7), 3) == 2
    assert round(n_term.charge(10), 3) == 1
    assert round(residues["HIS"].charge(5), 3) == 1
# Negative charge
    assert round(residues["GLU"].charge(5), 3) == -1
    assert round(residues["ASP"].charge(14), 3) == -1
    assert round(residues["TYR"].charge(14), 3) == -1
    assert round(residues["CYS"].charge(14), 3) == -1
    assert round(c_term.charge(14), 3) == -1

def test_mass():
    """Tests amino acid masses"""

    assert residues['ALA'].mass == 71.0788
    assert residues['ARG'].mass == 156.1875
    assert residues['ASN'].mass == 114.1038
    assert residues['ASP'].mass == 115.0886
    assert residues['CYS'].mass == 103.1388
    assert residues['GLN'].mass == 128.1307
    assert residues['GLU'].mass == 129.1155
    assert residues['GLY'].mass == 57.0519
    assert residues['HIS'].mass == 137.1411
    assert residues['ILE'].mass == 113.1594
    assert residues['LEU'].mass == 113.1594
    assert residues['LYS'].mass == 128.1741
    assert residues['MET'].mass == 131.1926
    assert residues['PHE'].mass == 147.1766
    assert residues['PRO'].mass == 97.1167
    assert residues['SER'].mass == 87.0782
    assert residues['THR'].mass == 101.1051
    assert residues['TRP'].mass == 186.2132
    assert residues['TYR'].mass == 163.1760
    assert residues['VAL'].mass == 99.1326
