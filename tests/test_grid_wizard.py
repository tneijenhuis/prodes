from prodes.io.parser import PDBparser

structure = PDBparser().parse("tests/data/1GDW_h.pdb.zip")
atoms = structure.heavy_atoms[:40]
