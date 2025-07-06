import numpy as np
from prodes.core.structure import Structure
from prodes.core.chain import Chain
from prodes.core.residue import Residue
from prodes.core.atom import Atom


class PDBparser:

    def parse(self, file, identifier="ATOM"):
        """Parses pdb files"""

        self.identifier = identifier
        name = file.split("/")[-1].split(".")[0]

        if file[-4:] == ".pdb":
            with open(file) as pdb:
                structure = self._read_pdb(pdb, name)

        else:
            raise ValueError("File extention not recognized by parser")

        return structure

    def _read_pdb(self, pdb, name):
        """Function takes a PDB formatted file and returns a Structure object which contains all atom information of the file"""

        current_structure = self._create_structure(name)

        for line in pdb:
            identifier = line[0:6].strip()
            if identifier == self.identifier:
                information = self._read_line(line)
                current_atom = Atom(identifier, *information)
                self._add_atom(current_atom, current_structure)

        current_structure.residues[-1].terminus = "C"
        return current_structure

    def _create_structure(self, name):
        """initializes the parser"""

        self.current_chain = ""
        self.current_residue = ""
        self.residue_numbers = {}
        self.chain_names = []
        return Structure(name)

    def _read_line(self, line):

        name = line[12:17].strip()
        residue_name = line[17:20].strip()
        chain_name = line[21].strip()
        residue_number = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        segid = line[72:76].strip()
        element = line[76:78].strip()

        return name, residue_name, chain_name, residue_number, x, y, z, segid, element

    def _add_atom(self, atom, structure):
        """Adds a atom object to a structure object"""

        addarray = np.array([atom])
        structure.atoms = np.concatenate((structure.atoms, addarray))
        atom.structure = structure

# Find the chain of the atom
        if atom.chain_name not in self.chain_names:
            if len(self.chain_names) != 0:
                self.current_chain.residues[-1].terminus = "C"
            self.current_chain = self._make_new_chain(structure, atom)

        elif atom.chain_name != self.current_chain.name:
            self.current_chain = structure.chains[self.chain_names.index(atom.chain_name)]

        atom.chain = self.current_chain
        self._add_atom_to_chain(atom, self.current_chain)

# Find the residue of the atom
        if atom.residue_number not in self.residue_numbers[self.current_chain.name]:
            self.current_residue = self._make_new_residue(structure, atom)

        atom.residue = self.current_residue
        self._add_atom_to_residue(atom, self.current_residue)

    def _make_new_chain(self, structure, atom):
        """Makes a new chain object"""

        self.current_chain = Chain(atom.chain_name, atom.structure)
        self._add_chain_to_structure(self.current_chain, structure)
        self.chain_names.append(atom.chain_name)
        self.residue_numbers[self.current_chain.name] = []
        return self.current_chain

    def _add_atom_to_chain(self, atom, chain):
        """Adds an atom object to a chain object"""

        addarray = np.array([atom])
        chain.atoms = np.concatenate([chain.atoms, addarray])

    def _add_chain_to_structure(self, chain, structure):
        """Adds a chain object to a structure object"""

        addarray = np.array([chain])
        structure.chains = np.concatenate((structure.chains, addarray))

    def _make_new_residue(self, structure, atom):
        """Makes a new residue object"""

        self.current_residue = Residue(atom.residue_name, atom.structure, atom.residue_number, atom.chain)
        if len(self.current_chain.residues) == 0:
            self.current_residue.terminus = "N"
        self._add_residue_to_structure(self.current_residue, structure)
        self.residue_numbers[atom.chain_name].append(atom.residue_number)
        self._add_residue_to_chain(self.current_residue, self.current_chain)
        return self.current_residue

    def _add_residue_to_structure(self, residue, structure):
        """"adds residue to a structure"""

        addarray = np.array([residue])
        structure.residues = np.concatenate((structure.residues, addarray))

    def _add_atom_to_residue(self, atom, residue):
        """adds atom to a residue"""

        addarray = np.array([atom])
        residue.atoms = np.concatenate([residue.atoms, addarray])

    def _add_residue_to_chain(self, residue, chain):
        "adds residue to a chain"

        addarray = np.array([residue])
        chain.residues = np.concatenate([chain.residues, addarray])


def read_pka(file):
    """Reads a pka json file and will return a dict containing residue numbers and pkas
    pka json files can be generated using the pka_converter in prodes.io"""

    import json
   
    with open(file) as f:
        pka_dict = {}
        pkas = json.loads(f.read())
        for residue, pka_list in pkas.items():
            residue = int(residue)
            pka_dict[residue] = []
            for pka in pka_list:
                for identifier, pka_value in pka.items():
                    pka_dict[residue].append({identifier: float(pka_value)})
        
    return pka_dict    

def write_pdb(structure, filename, chain = "all"):

    """
    takes a Structure object to generate a PDB formatted file using the Atoms information
    
    arguments:
    structure = Structure object
    filename = PATH + name of the to be generated file
    """

    from prodes import data

    if ".pdb" not in filename:
        print("can only make files witb pdb extention")
        
    else:
        with open(filename, 'w') as f:
            if chain == "all":
                atoms = structure.atoms

            else:
                atoms = np.empty([0])
                for struct_chain in structure.chains:
                    if struct_chain.name in chain:
                        atoms = np.concatenate((atoms, struct_chain.atoms))
            viable_residues = data.all_residues().keys()
            atom_nmbr = 0
            for atom in atoms:
                
                if atom_nmbr > 0:
                    if col4.strip() in viable_residues and atom.residue_name in viable_residues and col5.strip() != atom.chain_name:
                        
                        atom_nmbr += 1
                        
                        col2 = str(atom_nmbr)
                        for i in range((5 - len(col2))):
                            col2 = " " + col2
                        
                        f.write("TER   {}     {}{}{}\n".format(col2, col4, col5, col6))
                    elif col4.strip() in viable_residues and atom.residue_name not in viable_residues:

                        atom_nmbr += 1
                        
                        col2 = str(atom_nmbr)
                        for i in range((5 - len(col2))):
                            col2 = " " + col2
                        
                        f.write("TER   {}     {}{}{}\n".format(col2, col4, col5, col6))

                atom_nmbr +=1
                
                col1 = atom.identifier
                for i in range((6 - len(col1))):
                    col1 = col1 + " "
                
                col2 = str(atom_nmbr)
                for i in range((5 - len(col2))):
                    col2 = " " + col2
                
                col3 = atom.name   
                for i in range((5 - len(col3))):
                    if i < 2:
                        col3 = " " + col3
                    else:
                        col3 = col3 + " "
                
                col4 = atom.residue_name
                for i in range((5 - len(col4))):
                    if i < 1:
                        col4 = " " + col4
                    else:
                        col4 = col4 + " "
                
                col5 = atom.chain_name

                col6 = str(atom.residue_name)
                for i in range((4 - len(col6))):
                    col6 = " " + col6
                
                col7 = str(atom.x)
                if len(col7) > 11:
                    col7 = col7[:7]
                for i in range((11 - len(col7))):
                    col7 = " " + col7

                col8 = str(atom.y)
                if len(col8) > 8:
                    col8 = col8[:7]
                for i in range(8 - len(col8)):
                    col8 = " " + col8
                
                col9 = str(atom.z)
                if len(col9) > 8:
                    col9 = col9[:7]
                for i in range(8 - len(col9)):
                    col9 = " " + col9 
                
                col10 = ""
                for i in range(6-len(col10)):
                    col10 = " " + col10
                
                col11 = ""
                for i in range(6-len(col11)):
                    col11 = " " + col11
                
                col12 = atom.segment_id
                for i in range(4-len(atom.segment_id)):
                    col12 = " " + col12

                col13 = atom.element
                for i in range(2-len(col13)):
                    col13 = " " + col13

                f.write("{}{}{}{}{}{} {}{}{}{}{}      {}{}\n".format(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13))
            f.write("END")


class Builder():

    def build_dummy_atom(self, x, y, z, chain_name = "A", ep = "", size = ""):
        identifier = "ATOM"
        name = "X"
        residue_name = "DUM"
        chain_name = chain_name
        residue_number = 0 
        occupancy = ""
        segid = ""
        temperature_factor = ""
        element = "X"
        

        atom = Atom(identifier, name, residue_name, chain_name, residue_number, x, y, z, occupancy,
                temperature_factor, segid, element)

        return atom
