import numpy as np


class Residue:

    def __init__(self, name, structure, number=0, chain=None, _pka=None, terminus=None):
        self.atoms = np.empty([0])
        self.name = name
        self.number = number
        self.structure = structure
        self.chain = chain
        self._pka = _pka
        self.terminus = terminus
        self._heavy_atoms = None
        self.disulfide_partner = None

    @property
    def heavy_atoms(self):
        """Returns a list containing all heavy atoms within the structure"""
        if self._heavy_atoms is None:
            self._heavy_atoms = np.array([atoms for atoms in self.atoms if atoms.element != "H"])
        return self._heavy_atoms

    @property
    def protons(self):
        """Returns a list containing all protons within the residue"""

        return [atoms for atoms in self.atoms if atoms.element == "H"]

    @property
    def pkas(self):
        """The titratable groups of this residue, as [{group: pka}], or None if it has none.

        The side chain comes first and the terminus second, but nothing should
        rely on that: read a value with group_pka or side_chain_pka instead of
        by position, or a residue whose side chain is missing hands out its
        terminal pKa in place of one it does not have.
        """

        if self._pka is None:
            from prodes import data

            pka = data.residue_data(self.name)["pka"]
            pkas = []
            # A cysteine in a disulfide has no thiol proton and does not titrate
            # at all, so it gets no side-chain entry, the same as a residue that
            # was never titratable. Doing it here rather than in charge() means
            # every consumer of pkas agrees about it.
            if pka is not None and self.disulfide_partner is None:
                pkas.append({self.name: pka})
            if self == self.chain.residues[0]:
                pkas.append({"N+": data.n_term_pka()})
            elif self == self.chain.residues[-1]:
                pkas.append({"C-": data.c_term_pka()})
            if len(pkas) > 0:
                self._pka = pkas

        return self._pka

    def group_pka(self, key):
        """Returns the pKa of one titratable group, or None if this residue has no such group.

        The key is a three-letter residue name for a side chain, or "N+" or
        "C-" for a terminus.
        """

        for group in self.pkas or []:
            if key in group:
                return group[key]

        return None

    @property
    def side_chain_pka(self):
        """Returns the pKa of the titratable side chain, or None if there is none.

        None means the group does not exist, which covers both a residue type
        that never titrates and a cysteine that is bonded into a disulfide.
        """

        return self.group_pka(self.name)

    def set_group_pka(self, key, pka):
        """Replaces the pKa of one group, leaving the residue's other groups alone.

        Returns True when the residue has that group and the value was applied,
        and False when it does not, which is how a caller can tell that a
        predicted value has nowhere to go.
        """

        for group in self.pkas or []:
            if key in group:
                group[key] = pka
                return True

        return False

    def mark_disulfide_partner(self, other):
        """Records which cysteine this one's SG is bonded to, or None for none.

        Clearing the cached pKas keeps the residue correct even when something
        has already read them, so detection does not have to be the first thing
        that touches a freshly parsed structure, and re-running it does not
        leave a stale partner behind.
        """

        self.disulfide_partner = other
        self._pka = None

    @property
    def mass(self):
        """returns the mass"""

        from prodes.data import residue_data

        return residue_data(self.name)["mass"]

    def charge(self, ph, formal=True):
        """Calculates the charge of the residue"""

        charge = 0
        # charged_atoms = self.charged_atoms(ph)
        for atom in self.heavy_atoms:
            charge += atom.charge(ph, formal)
        return charge

    def charged_atoms(self, ph):
        """searches for the charged atom of the residue"""

        from prodes.calculations.standard_equations import neg_charge, pos_charge
        from prodes.data import residue_data

        charged_atoms = []
        if self.pkas:

            for pka in self.pkas:
                pka_value = list(pka.values())[0]
                key = list(pka.keys())[0]
                if key in ["GLU", "ASP", "TYR", "CYS", "C-", "-C"]:
                    charge = neg_charge(pka_value, ph)
                    if key == "C-" or key == "-C":
                        for atom in self.atoms:
                            if atom.name == self.terminus:
                                charged_atoms.append(atom)
                                atom._charge = charge
                    else:
                        for atom in self.atoms:
                            charged = residue_data(key)["charged_atoms"]
                            if atom.name in charged:
                                charged_atoms.append(atom)
                                atom._charge = charge / len(charged)
                else:
                    charge = pos_charge(pka_value, ph)
                    if key == "-N" or key == "N+":
                        for atom in self.atoms:
                            if atom.name == self.terminus:
                                charged_atoms.append(atom)
                                atom._charge = charge
                    else:
                        for atom in self.atoms:
                            charged = residue_data(key)["charged_atoms"]
                            if atom.name in charged:
                                charged_atoms.append(atom)
                                atom._charge = charge / len(charged)
        return charged_atoms

    def area(self, probe_r=1.4):
        """Calculates the rsa of a residue"""

        if self.structure.surface_done is False or self.structure is None:
            raise NameError("First calculate the SASA before calling the surface area")

        total_area = 0
        for atom in self.heavy_atoms:
            total_area += atom.surface_area(probe_r=probe_r)

        return total_area

    def rsa(self, probe_r=1.4):
        """Calculate the RSA based on the surface area"""

        if self.structure.surface_done is False or self.structure is None:
            raise NameError("First calculate the SASA before calling the surface area")

        else:
            from prodes.data import residue_data

            gly_x_gly = residue_data["gly_x_gly"]
            area = self.area(probe_r)
            return area / gly_x_gly
