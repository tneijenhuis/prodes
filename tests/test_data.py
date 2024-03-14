import pytest
from prodes import data

def test_hydrophobic_scale():
    assert data.hydrophobic_scale("cw") == {
        'ALA': 0.42, 'ARG': -1.56, 'ASN': -1.03, 'ASP': -0.51, 'CYS': 0.84, 'GLN': -0.96, 'GLU': -0.37,
        'GLY': 0, 'HIS': -2.28, 'ILE': 1.81, 'LEU': 1.8,'LYS': -2.03, 'MET': 1.18, 'PHE': 1.74,
        'PRO': 0.86, 'SER': -0.64, 'THR': -0.26, 'TRP': 1.46, 'TYR': 0.51, 'VAL': 1.34
    }
    with pytest.raises(KeyError):
        data.hydrophobic_scale("K")


def test_residue_data():
    assert data.residue_data("ARG") == {
        'letter': 'R', 'mass': 156.1875, 'potential_charge': 1,
        'pka': 13.8, 'charged_atoms': ['NE', 'NH1', 'NH2'], 'aromatic_carbons': [], 'gly_x_gly': 274.0
    }

    with pytest.raises(KeyError):
        data.hydrophobic_scale("R")
