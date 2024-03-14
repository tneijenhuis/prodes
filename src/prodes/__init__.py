from prodes.io import parser
from prodes.run import calculate

def run_prodes(pdb_file, out_file, pkas_file=None, ph=7, r_probe=1.4, hydro_scale="mj_scales"):
    """Runs Prodes similar as the commantline tool"""

    calculate(pdb_file, out_file, pkas_file, ph, r_probe, hydro_scale)
    

def read(file):
    """Loads files
    List of supported extentions:
    .pdb
    .pka"""

    if file[-4:] == ".pdb":
        structure = parser.PDBparser().parse(file)
        return structure
    
    elif file[-4:] == ".pka":
       pkas = parser.read_pka()
       return pkas

    else:
        raise ValueError("File extention not recognized by prodes")
    
