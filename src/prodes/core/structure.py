
import numpy as np


class Structure:

    def __init__(self, name=None):
        self.name = name
        self.chains = np.empty([0])
        self.residues = np.empty([0])
        self.atoms = np.empty([0])
        self._x = None
        self._y = None
        self._z = None
        self.surface_done = False

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

    def isoelectric_point(self, decimals=3):
        """estimates the isoelectic point of the protein"""

        ph_upper = 14
        ph_lower = 0
        for _ in range(100000):

            diff = ph_upper - ph_lower
            if round(diff, decimals) == 0:
                break

            ph = (ph_upper+ph_lower)/2
            if self.charge(ph, formal=False) < 0:
                ph_upper = ph
            else:
                ph_lower = ph
        return round(ph, decimals)

    @property
    def heavy_atoms(self):
        """Returns a list containing all heavy atoms within the structure"""

        return np.array([atom for residue in self.residues for atom in residue.heavy_atoms])

    @property
    def furthest_heavy_atom(self):
        """returns the heavy atom which is furthest from the protein centroid"""
        from prodes.calculations.standard_equations import distance

        longest_distance = 0.0
        furthest = None
        for atom in self.heavy_atoms:
            dist = distance(atom, self)
            if dist > longest_distance:
                longest_distance = dist
                furthest = atom

        return furthest

    def surface_area(self, probe_r=1.4):

        if self.surface_done is False:
            raise NameError("First calculate the SASA before calling the surface area")
        
        total_area = 0
        for atom in self.heavy_atoms:
            total_area += atom.surface_area(probe_r=probe_r)

        return total_area

    def charge(self, ph=7, formal=True):
        """Cacluclates the charge of a protein"""
        charge = 0
        for residue in self.residues:
            charge += residue.charge(ph, formal)
        return charge

    def dipole(self, ph=7, formal=True):
        """Calculates the depole moment using atom partial charges"""

        from prodes.core.point import Point
        from prodes.calculations.standard_equations import distance

        debye = 4.803
        dipole_vector = Point(0, 0, 0)
        centre = Point(0, 0, 0)

        for atom in [atoms for residues in self.residues for atoms in residues.charged_atoms(ph)]:

            dipole_vector.x += (atom.x - self.x) * atom.charge(ph, formal)
            dipole_vector.y += (atom.y - self.y) * atom.charge(ph, formal)
            dipole_vector.z += (atom.z - self.z) * atom.charge(ph, formal)

        dipole_vector.x *= debye
        dipole_vector.y *= debye
        dipole_vector.z *= debye
        dipole = distance(dipole_vector, centre)

        return round(dipole, 3)


    def count_residues_on_surf(self, cutoff=0.20):
        """counts the number of residues on the protein surface using the RSA"""

        from prodes.data import all_residues
        protein_aa = all_residues.keys()

        residue_count = {}
        for residue in protein_aa:
            residue_count[residue] = 0

        for residue in self.residues:
            if residue.rsa >= cutoff:
                residue_count[residue.name] += 1

        if sum(residue_count.values()) == 0:
            raise NameError("No surface has been defined")

        return residue_count

    def residue_surf_fractions(self, r_probe):
        """Calculate the residue surface fractions"""

        from prodes.data import all_residues
        protein_aa = all_residues().keys()
        residue_areas = {}
        for residue in protein_aa:
            residue_areas[residue] = 0

        for residue in self.residues:
            residue_areas[residue.name] += residue.area(r_probe)

        residue_fractions = {}

        for residue in residue_areas:
            residue_fractions[residue] = round(residue_areas[residue]/self.surface_area(r_probe), 3)

        return residue_fractions

    def rotate(self, x_degree=0, z_degree=0):
        """Rotates the structure around the x and z axis"""

        from math import cos, sin, radians

        x_degree = radians(x_degree)
        z_degree = radians(z_degree)
        for atom in self.atoms:
            x = atom.x - self.x
            y = atom.y - self.y
            z = atom.z - self.z
# Rotate over x

            y_rotated = y*cos(x_degree) - z*sin(x_degree)
            z_rotated = y*sin(x_degree) + z*cos(x_degree)

# Rotate over z
            x_rotated = x*cos(z_degree) - y_rotated*sin(z_degree)
            y_rotated = x*sin(z_degree) + y_rotated*cos(z_degree)

            atom.x = x_rotated + self.x
            atom.y = y_rotated + self.y
            atom.z = z_rotated + self.z

    def redo_pkas(self, pka_dict):
        """Takes a dictionary of {residues_numbers : pka}.
        Changes the pka all residues in the dict"""

        for residue in self.residues:
            if residue.number in pka_dict:
                pka = pka_dict[residue.number]
                residue._pka = pka
