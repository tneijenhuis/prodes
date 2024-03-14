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

    @property
    def heavy_atoms(self):
        """Returns a list containing all heavy atoms within the structure"""

        return np.array([atoms for atoms in self.atoms if atoms.element != "H"])

    @property
    def protons(self):
        """Returns a list containing all protons within the residue"""

        return [atoms for atoms in self.atoms if atoms.element == "H"]

    @property
    def pkas(self):
        if self._pka is None:
            from prodes import data 
            pka = data.residue_data(self.name)["pka"]
            pkas = []
            if pka is not None:
                pkas.append({self.name: pka})
            if self == self.chain.residues[0]:
                pkas.append({"N+": data.n_term_pka()})
            elif self == self.chain.residues[-1]:
                pkas.append({"C-": data.c_term_pka()})
            if len(pkas) > 0:
                self._pka = pkas

        return self._pka

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

        from prodes.data import residue_data
        from prodes.calculations.standard_equations import pos_charge, neg_charge
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
                                atom._charge = charge/len(charged)
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
                                atom._charge = charge/len(charged)
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
