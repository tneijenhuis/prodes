import numpy as np


class Chain:

    def __init__(self, name, structure):
        self.residues = np.empty([0])
        self.atoms = np.empty([0])
        self.name = name
        self.structure = structure
        self._x = None
        self._y = None
        self._z = None

    @property
    def x(self):
        """returns the x coordinate of the com"""

        if self._x:
            return self._x
        else:
            x = 0.0
            for atom in self.heavy_atoms:
                x += atom.x
            self._x = x / len(self.heavy_atoms)
            return self._x

    @property
    def y(self):
        """returns the y coordinate of the com"""

        if self._y:
            return self._y

        else:
            y = 0.0
            for atom in self.heavy_atoms:
                y += atom.y
            self._y = y / len(self.heavy_atoms)
            return self._y

    @property
    def z(self):
        """Returns the z coordinate of the com"""

        if self._z:
            return self._z
        else:
            z = 0.0
            for atom in self.heavy_atoms:
                z += atom.z
            self._z = z / len(self.heavy_atoms)
            return self._z

    @property
    def mw(self):
        """Returns the molecular mass"""

        total_mass = 0.0
        for residue in self.residues:
            total_mass += residue.mass
# Add the C terminal O to the masses 
        total_mass +=  18.01524
                
        return round(total_mass, 2)

    @property
    def heavy_atoms(self):
        """Returns a list containing all heavy atoms within the structure"""

        return np.array([atom for residue in self.residues for atom in residue.heavy_atoms])
