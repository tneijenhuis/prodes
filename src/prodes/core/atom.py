from prodes import data

class Atom:
    """Contains all information of each atom from the PDB file"""

    def __init__(self, identifier, name, residue_name, chain_name, residue_number, x, y, z,
                 segment_id="", element=None, structure=None, chain=None,
                 residue=None, _charge=None):

        self.identifier = identifier
        self.name = name
        self.residue_name = residue_name
        self.chain_name = chain_name
        self.residue_number = residue_number
        self.x = x
        self.y = y
        self.z = z
        self.segment_id = segment_id
        self.element = element
        self.structure = structure
        self.chain = chain
        self.residue = residue
        self.lipo = None
        self.cloud = []
        self.cell = None
        self._charge = _charge

    @property
    def radius(self):
        """Will find the vdw radius of the atom.
        If atom is not in the list of standard elements, none will be returened"""
        
        if self.element is None:
            element = self.name[0]

        else:
            element = self.element

        if element != "H":
            if self.name in data.residue_data(self.residue_name)["aromatic_carbons"]:
                return data.vdw_radius("Cr")
            else:
                return data.vdw_radius(element)
        else:
            return None

    def charge(self, ph=7, formal=True):
        """Returns the charge of a atom"""

        from prodes import data
        from prodes.calculations.standard_equations import pos_charge, neg_charge
        
        charge = 0

        residue_data = data.residue_data(self.residue_name)
        potential_charge = residue_data["potential_charge"]
        if potential_charge is not None:
        
            if self.name in residue_data["charged_atoms"]:
                if formal is True:
                    if potential_charge > 0:
                        if list(self.residue.pkas[0].values())[0] > ph:
                            charge = 1/len(residue_data["charged_atoms"])
                        else:
                            charge = 0
                    
                    elif potential_charge < 0:
                        if list(self.residue.pkas[0].values())[0] < ph:
                            charge = -1/len(residue_data["charged_atoms"])
                        else:
                            charge = 0

                else: 
                    if potential_charge > 0:
                        charge = pos_charge(list(self.residue.pkas[0].values())[0], ph)/len(residue_data["charged_atoms"])
                    
                    else:
                        charge = neg_charge(list(self.residue.pkas[0].values())[0], ph)/len(residue_data["charged_atoms"])
        
        if self.name == self.residue.terminus:

            if formal is True:
                if self.name == "N":
                    pka = list(self.residue.pkas[-1].values())[0]
                    if pka > ph:
                        charge = 1

                    else:
                        charge = 0

                elif self.name == "C":
                    pka = list(self.residue.pkas[-1].values())[0]
                    if pka < ph:
                        charge = -1

                    else:
                        charge = 0

            else:
                if self.name == "N":
                    charge = pos_charge(list(self.residue.pkas[-1].values())[0], ph)
                
                elif self.name == "C":
                    charge = neg_charge(list(self.residue.pkas[-1].values())[0], ph)

        return charge
        # if self._charge is None:

        #     charged_atoms = self.residue.charged_atoms(ph)

        #     if self.name not in charged_atoms:
        #         self._charge = 0
        
        # return self._charge
        # if exact:
        #     return self._charge
        
        # elif self._charge != 0:
        #     from prodes.data import residue_data


        #     if self.name == self.residue.terminus:           
        #         if self._charge > 0.01:
        #             return 1
        #         elif self._charge < -0.01:
        #             return -1
        #         else:
        #             return 0
        #     else:
        #         charged_sidechain_atoms = len(residue_data(self.residue_name)["charged_atoms"])
        #         if self._charge*charged_sidechain_atoms > 0.01:
        #             return 1/charged_sidechain_atoms
        #         elif self._charge*charged_sidechain_atoms < -0.01:
        #             return -1/charged_sidechain_atoms
        #         else:
        #             return 0
        # else:
        #     return 0


    def max_area(self, probe_r=1.4):
        """calculates the maximal surface area of an atom"""

        from math import pi

        r = self.radius + probe_r
        return (4 * pi * (r**2))

    def max_n_points(self, probe_r=1.4, points_per_a= 2):
        area = self.max_area(probe_r=probe_r)
        return int(area*points_per_a)

    def surface_area(self, probe_r=1.4, points_per_a=2):
        """calculates the surface area"""
        
        if self.structure.surface_done is False or self.structure == None:
            raise NameError("First calculate the SASA before calling the surface area")

        else:
            n_surf_points = len(self.cloud)
            max_area = self.max_area(probe_r=probe_r)
            max_n_points = self.max_n_points(probe_r=probe_r, points_per_a=points_per_a) 
            return max_area/max_n_points*n_surf_points
