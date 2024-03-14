import numpy as np
from prodes.calculations import grid_wizard as gw
from prodes.io.parser import PDBparser

structure = PDBparser().parse("tests/data/1GDW_h.pdb")
atoms = structure.heavy_atoms[:40]
