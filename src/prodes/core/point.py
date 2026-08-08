from dataclasses import dataclass

import numpy as np


@dataclass
class Point:
    """Represents a point in space, with coordinates x,y,z"""

    x: float
    y: float
    z: float


class Surface_point:
    """Class for points which represent the surface sphere of an atom"""

    def __init__(self, x, y, z, atom):
        self.x = x
        self.y = y
        self.z = z
        self.atom = atom


class Property_point:
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

    def set_values(self, ep=None, lipo=None):
        """Writes an already calculated electrostatic potential and/or lipophilicity onto the point.

        Used to reassemble results computed in worker processes: a forked worker
        mutates its own copy of the point, so the parent writes the returned
        values onto the original object itself.
        """

        if ep is not None:
            self.__ep = ep

        if lipo is not None:
            self.__lipo = lipo

    def set_ep(self, atoms, ph=7, formal=True, cutoff=10000):
        """Projects the partial charges of an array of atoms onto the point
        -------
        Args :
            Atoms : array of Atom objects"""

        charged = [atom for atom in atoms if atom.charge(ph=ph, formal=formal) != 0]
        if not charged:
            self.__ep = 0
            return

        coords = np.array([[a.x, a.y, a.z] for a in charged])
        charges = np.array([a.charge(ph=ph, formal=formal) for a in charged])
        point_coord = np.array([self.x, self.y, self.z])

        dists = np.sqrt(np.sum((coords - point_coord) ** 2, axis=1))
        mask = dists <= cutoff
        if not np.any(mask):
            self.__ep = 0
            return

        coulomb_charges = charges[mask] * 1.6e-19
        dists_m = dists[mask] * 1e-10
        absolute_permittivity = 8.854e-12
        ep = np.sum(coulomb_charges / (4 * absolute_permittivity * dists_m * 4 * np.pi))
        self.__ep = round(float(ep), 2)

    def set_lipo(self, atoms, cutoff=10, scale="mj_scaled"):
        """Projects the lipophilysity of an array of atoms onto the point
        -------
        Arguments:
            Atoms : array of Atom objects"""

        from prodes.data import hydrophobic_scale

        non_h = [atom for atom in atoms if atom.element != "H"]
        if not non_h:
            self.__lipo = 0
            return

        scale_dict = hydrophobic_scale(scale)
        coords = np.array([[a.x, a.y, a.z] for a in non_h])
        point_coord = np.array([self.x, self.y, self.z])

        dists = np.sqrt(np.sum((coords - point_coord) ** 2, axis=1))
        mask = dists <= cutoff
        if not np.any(mask):
            self.__lipo = 0
            return

        dists_masked = dists[mask]
        atoms_masked = [a for a, m in zip(non_h, mask, strict=True) if m]
        contributions = np.array([1.0 if a.name == "OXT" else scale_dict[a.residue_name] for a in atoms_masked])
        lipophilicity = float(np.sum(contributions * np.exp(-dists_masked)))
        self.__lipo = lipophilicity


class Cell:
    """Class which represents the the points onto which properties are projected"""

    def __init__(self, size, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.size = size
        self._content_list = []
        self._content_array = None
        self.empty = True

    @property
    def content(self):
        if self._content_array is None:
            self._content_array = np.array(self._content_list)
        return self._content_array

    def filtered_content(self, *content_filter):
        """returns an array with filtered objects"""
        return np.array([component for component in self.content if type(component).__name__ in content_filter])

    def add_content(self, content):
        """Adds something to the cell instance"""
        self._content_list.append(content)
        self._content_array = None
        self.empty = False
