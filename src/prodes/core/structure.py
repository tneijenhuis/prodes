import logging

import numpy as np

logger = logging.getLogger(__name__)


class Structure:

    def __init__(self, name=None):
        self.name = name
        self.chains = np.empty([0])
        self.residues = np.empty([0])
        self.atoms = np.empty([0])
        self._x = None
        self._y = None
        self._z = None
        self._heavy_atoms = None
        self.surface_done = False
        # Filled in by prodes.calculations.disulfides.assign_disulfides, which
        # the parser calls on every structure it reads.
        self.disulfides = []

    def _compute_centroid(self):
        coords = np.array([[a.x, a.y, a.z] for a in self.heavy_atoms])
        self._x, self._y, self._z = coords.mean(axis=0)

    @property
    def x(self):
        """returns the x coordinate of the com"""
        if self._x is None:
            self._compute_centroid()
        return self._x

    @property
    def y(self):
        """returns the y coordinate of the com"""
        if self._y is None:
            self._compute_centroid()
        return self._y

    @property
    def z(self):
        """Returns the z coordinate of the com"""
        if self._z is None:
            self._compute_centroid()
        return self._z

    @property
    def mw(self):
        """Returns the molecular mass"""

        total_mass = 0.0
        for residue in self.residues:
            total_mass += residue.mass
        # Add the C terminal O to the masses
        total_mass += 18.01524

        return round(total_mass, 2)

    def isoelectric_point(self, decimals=3):
        """estimates the isoelectic point of the protein"""

        ph_upper = 14
        ph_lower = 0
        for _ in range(100000):

            diff = ph_upper - ph_lower
            if round(diff, decimals) == 0:
                break

            ph = (ph_upper + ph_lower) / 2
            if self.charge(ph, formal=False) < 0:
                ph_upper = ph
            else:
                ph_lower = ph
        return round(ph, decimals)

    @property
    def heavy_atoms(self):
        """Returns a list containing all heavy atoms within the structure"""
        if self._heavy_atoms is None:
            self._heavy_atoms = np.array([atom for residue in self.residues for atom in residue.heavy_atoms])
        return self._heavy_atoms

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

        from prodes.calculations.standard_equations import distance
        from prodes.core.point import Point

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
            residue_fractions[residue] = round(residue_areas[residue] / self.surface_area(r_probe), 3)

        return residue_fractions

    def rotate(self, x_degree=0, z_degree=0):
        """Rotates the structure around the x and z axis"""

        from math import cos, radians, sin

        x_degree = radians(x_degree)
        z_degree = radians(z_degree)
        for atom in self.atoms:
            x = atom.x - self.x
            y = atom.y - self.y
            z = atom.z - self.z
            # Rotate over x

            y_rotated = y * cos(x_degree) - z * sin(x_degree)
            z_rotated = y * sin(x_degree) + z * cos(x_degree)

            # Rotate over z
            x_rotated = x * cos(z_degree) - y_rotated * sin(z_degree)
            y_rotated = x * sin(z_degree) + y_rotated * cos(z_degree)

            atom.x = x_rotated + self.x
            atom.y = y_rotated + self.y
            atom.z = z_rotated + self.z

    def titratable_groups(self):
        """Returns every (residue, group key) pair in the structure that could take a predicted pKa.

        A disulfide-bonded cysteine is not among them: it has no titratable side
        chain, so a pKa predictor is right to say nothing about it.
        """

        return [(residue, key) for residue in self.residues for group in (residue.pkas or []) for key in group]

    def redo_pkas(self, pka_dict):
        """Applies predicted pKa values from a file, one titratable group at a time.

        Takes {residue number: [{group: pka}]}, as read by
        prodes.io.parser.read_pka, where the group is a three letter residue
        name for a side chain or "N+" or "C-" for a terminus.

        Values are merged into the residue's own groups rather than replacing
        its list wholesale. Replacing it meant that a file naming a residue's
        side chain but not its terminus deleted the terminal pKa, after which
        the terminal atom silently took its charge from the side-chain value.
        Groups the file does not mention keep their default, which is what the
        README has always promised.

        Two kinds of entry are dropped rather than applied:

        - a side-chain pKa for a cysteine that is in a disulfide, because the
          group it predicts does not exist. PROPKA agrees, reporting such a
          cysteine as 99.99, its marker for "does not titrate"; a real value
          from another predictor is logged at info level before it is discarded.
        - a pKa for a group the residue does not have at all, for instance a
          serine pKa from H++. Such an entry used to reach charged_atoms and
          raise a TypeError there.

        Groups the file leaves out are counted and reported. They keep the
        textbook value for their residue type, which is usually not what the
        absence was meant to convey.
        """

        from prodes.io.pka_converter import PROPKA_NOT_TITRATABLE

        applied = set()
        for residue in self.residues:
            for group in pka_dict.get(residue.number, []):
                for key, pka in group.items():

                    if residue.set_group_pka(key, pka):
                        applied.add((id(residue), key))

                    elif residue.disulfide_partner is not None and key == residue.name:
                        if pka != PROPKA_NOT_TITRATABLE:
                            logger.info(
                                "the pKa file gives %s %s a pKa of %s, but it is bonded into a disulfide and cannot titrate; ignoring the predicted value",
                                residue.name,
                                residue.number,
                                pka,
                            )

                    else:
                        logger.warning("the pKa file gives %s %s a %s pKa, which that residue does not have; ignoring it", residue.name, residue.number, key)

        missing = [(residue, key) for residue, key in self.titratable_groups() if (id(residue), key) not in applied]
        if missing:
            named = ", ".join(f"{key} {residue.number}" for residue, key in missing[:5])
            logger.warning(
                "%d titratable groups are not in the pKa file and keep the textbook value for their residue type: %s%s",
                len(missing),
                named,
                ", ..." if len(missing) > 5 else "",
            )
