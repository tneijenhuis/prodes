import importlib.resources
import json

def hydrophobic_scale(scale):
    """Extracts the hydrophobicity scale"""
    
    with importlib.resources.open_text("prodes.data", "hydrophobicity.json") as f:
        data = json.load(f)
    
    try:
        return data[scale]   

    except KeyError:
        raise KeyError(f"Scale {scale} is not found in the list of hydrophobicity scales")


def residue_data(residue):
    """loads the data relating to a specific amino acid residue"""
    aas = all_residues()

    try:
        return aas[residue]
    except KeyError:
        raise KeyError(
            f"{residue} is not found in the list with residues "
            "Please make sure you provide the three letter code all caps"
        )
    

def all_residues():
    """loads all residue data"""
       
    return data["residues"]



def vdw_radius(element):
    """returns the van der waals radius of a particular element"""

    try:
        return data["vdw_radius"][element]

    except KeyError:
        raise KeyError(f"{element} is not found in the list of elements with")


def heavy_atom_mass(atom):
    """returns the mass of a heavy atom"""
    return data["heavy_atom_masses"][atom]


def c_term_pka():
    """returns the standard C terminus pka"""

    return data["C_term_pka"]


def n_term_pka():
    """returns the standard N terminus pka"""

    return data["N_term_pka"]


def _open_general():
    """opens the general data"""
    with importlib.resources.open_text("prodes.data", "general_data.json") as f:
        data = json.load(f)

    return data

data = _open_general()