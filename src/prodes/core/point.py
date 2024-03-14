from dataclasses import dataclass
import numpy as np


@dataclass
class Point:
    """Represents a point in space, with coordinates x,y,z"""

    x: float
    y: float
    z: float


class Surface_point():
    """Class for points which represent the surface sphere of an atom"""
    def __init__(self, x, y, z, atom):
        self.x = x
        self.y = y
        self.z = z
        self.atom = atom


class Property_point():
    """Class which represents the the points onto which properties are projected"""

    def __init__(self, x, y, z, ep=None, lipo=None):
        self.x = x
        self.y = y
        self.z = z
        self.__ep = ep
        self.__lipo = lipo

    @property
    def ep(self):
        """electrostatic potential"""
        if self.__ep is not None:
            return self.__ep

        else:
            raise NameError("ep is not yet defined, try set_ep before calling ep")

    @property
    def lipo(self):
        """projected lipophilicity"""
        if self.__lipo is not None:
            return self.__lipo

        else:
            raise NameError("lipo is not yet defined, try set_lipo before calling lipo")

    def set_ep(self, atoms, ph=7, formal=True, cutoff=10000):
        """Projects the partial charges of an array of atoms onto the point
        -------
        Args :
            Atoms : array of Atom objects"""

        from prodes.calculations.standard_equations import distance
        from prodes.calculations import distance_functions

        ep = 0
        for atom in np.array([atom for atom in atoms if atom.charge(ph=ph, formal=formal) != 0]):

            dist = distance(self, atom)
            if dist <= cutoff:
                charge = distance_functions.atom_charge_coulomb(atom.charge(ph=ph, formal=formal))

                ep_cont = distance_functions.charge_simple(charge, dist*10**-10, 4)
                ep += ep_cont

        self.__ep = round(ep, 2)

    def set_lipo(self, atoms, cutoff=10, scale="mj_scaled"):
        """Projects the lipophilysity of an array of atoms onto the point
        -------
        Arguments:
            Atoms : array of Atom objects"""

        from math import exp

        from prodes.data import hydrophobic_scale
        from prodes.calculations.standard_equations import distance

        lipophilicity = 0

        scale_dict = hydrophobic_scale(scale)
        for atom in np.array([atom for atom in atoms if atom.element != "H"]):
            dist = distance(self, atom)
            if dist <= cutoff:
                if atom.name != "OXT":
                    lipophilicity += scale_dict[atom.residue_name] * exp(-dist)
                else:
                    lipophilicity += 1 * exp(-dist)

        self.__lipo = lipophilicity


class Cell():
    """Class which represents the the points onto which properties are projected"""

    def __init__(self, size, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.size = size
        self.content = np.empty([0])
        self.empty = True

    def filtered_content(self, *content_filter):
        """returns an array with filtered objects"""
        return np.array([component for component in self.content if type(component).__name__ in content_filter])
        
    def add_content(self, content):
        """Adds something to the cell instance"""
        self.content = np.concatenate([self.content, np.array([content])])
        self.empty = False
