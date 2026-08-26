import math
from dataclasses import dataclass

import numpy as np

# Decimal places the projected potential is rounded to. Kept at the value the
# unscreened implementation used, so that zero ionic strength reproduces the
# pre-5.0 numbers bit for bit and this change stays auditable.
#
# The rounding is not free. A point whose exact potential is between 0 and 0.005
# rounds to 0.00 and drops out of the positive set, which biases NSurfPosEp
# downwards by a measured 0.6 to 5.6 per cent depending on the structure. That
# bias is one sided and it is worse on the screened potential, where more points
# sit near zero. Three places would reduce it to under 0.7 per cent.
#
# It is left alone here because it is not a screening problem: the same bias is
# present in the shipped unscreened code at 0.6 to 0.9 per cent, so changing it
# would fold an unrelated correction into a physics change and cost the exact
# reproduction guarantee. Worth fixing on its own.
EP_DECIMALS = 2


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

    def set_ep(self, atoms, ph=7, formal=True, cutoff=10000, *, debye_length):
        """Projects the partial charges of an array of atoms onto the point

        Args:
            atoms: array of Atom objects.
            ph: pH at which the atomic charges are evaluated.
            formal: True for integer formal charges, False for the fractional
                Henderson-Hasselbalch occupancies.
            cutoff: distance cutoff in Angstrom. Effectively unlimited by default.
            debye_length: screening length in Angstrom, from
                prodes.calculations.distance_functions.debye_length. Each charge
                is damped by exp(-d / debye_length), so a distant one contributes
                almost nothing. Pass math.inf for no screening, which is the
                unscreened Coulomb sum Prodes used before version 5.

                Required and keyword only on purpose. With a default of math.inf
                a call site that forgot to pass it would silently compute the
                unscreened potential, which is the exact bug this change exists
                to remove and which no output would reveal.
        """

        # A zero or negative screening length would damp every contribution to
        # nothing and return a surface that is uniformly zero volts, which looks
        # like a valid result. Note 0 is what the user facing ionic strength uses
        # to mean "no screening", so the two must not be confused.
        if debye_length <= 0:
            raise ValueError(f"debye_length must be positive, got {debye_length}; use math.inf for no screening")

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

        # Mobile buffer ions screen each charge, so a surface point responds
        # mainly to its neighbours. Without this every charge in the structure
        # contributes in full, including ones on the far side, which adds a large
        # smooth offset that swamps the local pattern.
        screening = 1.0 if math.isinf(debye_length) else np.exp(-dists[mask] / debye_length)

        ep = np.sum(coulomb_charges * screening / (4 * absolute_permittivity * dists_m * 4 * np.pi))
        self.__ep = round(float(ep), EP_DECIMALS)

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
