from rdkit import Chem
import copy
from label import Label

class MarkInChI(object):

    def reorder(self, inchi, inchiplus):

        # This function reorders the substituents
        # Reorder Zz atoms
        #Get number of Zz atoms
        index = inchi.find("Te")+2
        char = inchi[index]
        no = 0
        no_atoms = ""
        while char != "/":
            no_atoms += str(char)
            index += 1
            char = inchi[index]
        if no_atoms != "":  # more than one Zz atom in the formula
            no = int(no_atoms)
        # Get the rank of each Zz atom
        line = inchi.split("/")[2][1:]  # canonical numbering string
        list_of_no = []
        for n in line.split("-"):
            # produce an ordered list canonical numbers
            if "(" in n:
                n_list = n.replace("(", "-").replace(")", "-").split("-")
                list_of_no += n_list
            else:
                list_of_no.append(n)
        list_int = set(map(int, list_of_no))  # convert list of numbers to int
        start_zz = max(list_int)-no+1  # lowest canonical number of a Zz atom
        #Scan the line of numbers to define the order of Zz atoms
        ordered_list = []
        no_list = list(range(start_zz, max(list_int)+1))
        for atom_no in list_of_no:
            if int(atom_no) in no_list and int(atom_no) not in ordered_list:
                ordered_list.append(int(atom_no))
        # Rearrange the inchiplus to include the substituents for each Zz
        # atom in the right order
        inchiplus_new = inchiplus.copy()
        i = 0
        for number in ordered_list:
            inchiplus_new[i] = inchiplus[number-start_zz]
            i += 1
        return inchiplus_new

    def replacement(self, main_mol, substituent):

        # This function performs replacements on atoms "-"
        # substituent = rank-atom@replacement or rankH-replacement
        No_H = True  # True if replacement atom is not H
        new_mol = copy.deepcopy(main_mol)
        id = substituent.split("-")
        if "H" not in id[0]:
            rank, atom, replacement = tuple(id[0]) + tuple(id[1].split("@"))
            replacement = Chem.MolFromSmiles(replacement).GetAtoms()[0]
        else:
            rank, replacement = tuple(id)
            # find the atom from the rank
            atom = self.label.find_atom(rank[:-1], self.inchi.split("/")[1])
            if replacement != "H":
                # convert replacement to an Atom() object
                replacement = Chem.MolFromSmiles(replacement).GetAtoms()[0]
            else:
                No_H = False
        is_aromatic = False  # false if atom being replaced is not aromatic
        # get atomic mass for the atom being replaced
        table = Chem.GetPeriodicTable()
        atomic_mass = int(table.GetMostCommonIsotopeMass(atom))
        if No_H:
            # normal replacement no H
            if "H" not in rank:
                # get isotopic label of atom being replaced
                num = int(self.ranks[rank])
                for mol_atom in new_mol.GetAtoms():
                    cn1 = mol_atom.GetIsotope() == num
                    cn2 = mol_atom.GetSymbol() == str.upper(atom)
                    if cn1 and cn2:
                        # Replace the labelled atom
                        is_aromatic = mol_atom.GetIsAromatic()
                        replacement.SetIsAromatic(is_aromatic)
                        rwmol = Chem.RWMol(new_mol)
                        idx = mol_atom.GetIdx()
                        rwmol.ReplaceAtom(idx, replacement)
                        new_mol = copy.deepcopy(rwmol.GetMol())
            else:
                # get isotopic label of atom being replaced
                num = int(self.ranks[rank[:-1]])
                rwmol = Chem.RWMol(new_mol)
                new_rwmol = copy.deepcopy(rwmol)
                for mol_atom in rwmol.GetAtoms():
                    cn1 = mol_atom.GetIsotope() == num
                    cn2 = mol_atom.GetSymbol() == str.upper(atom)
                    if cn1 and cn2:
                        # Replace the labelled atom to og and add new atom
                        iso_num = self.label.find_isotope(self.inchi, rank[:-1])
                        print(f"iso_num: {iso_num}")
                        idx = mol_atom.GetIdx()
                        is_aromatic = mol_atom.GetIsAromatic()
                        add_index = rwmol.GetNumAtoms()
                        new_rwmol.GetAtoms()[idx].SetIsotope(iso_num)
                        new_rwmol.AddAtom(replacement)
                        single = Chem.rdchem.BondType.SINGLE
                        new_rwmol.AddBond(idx, add_index, order=single)
                        new_mol = copy.deepcopy(new_rwmol.GetMol())
        return copy.deepcopy(new_mol)

    def replace(self, main_mol, substituent):

        # This function replaces undefined atoms
        # Case not hydrogen
        new_mol = copy.deepcopy(main_mol)
        final_mol = None
        if substituent != "" and substituent != "H":
            sub_inchi = ""
            sub_mol = None
            sub_index = 0
            main_index = 0
            # We assume that the markinchi will be either:
            # 1- one atom that is convertable to a smiles and then to
            # an Atom()
            # 2- more than one atom written in the format of a "part-inchi"
            if substituent.find("/") != -1:
                # case 2
                sub_inchi = "InChI=1B/"+substituent  # convert to inchi
                sub_mol = Chem.rdinchi.InchiToMol(sub_inchi)[0]
                # Get number of atoms in the substituent after substracting
                # 1 for the removed Zz atom
                num = sub_mol.GetNumAtoms()-1
                # case 1: substituent doesn't start with "Te"
                # we get index of atom connected to it going backwards
                # and delete "Te" (see get_index function in label.py)
                if sub_mol.GetAtoms()[0].GetSymbol() != "Te":
                    sub_mol, sub_index = self.label.get_index(sub_mol)
                # case 2: substituent starts with "Te"
                # just delete "Te" atom and sub_index is left as zero
                else:
                    new_submol = self.label.delete_zz(sub_mol)
                # same cases as in sub_mol
                if main_mol.GetAtoms()[0].GetSymbol() != "Te":
                    new_mol, main_index = self.label.get_index(main_mol)
                    main_index += num
                else:
                    new_mol = self.label.delete_zz(main_mol)
                    main_index = num
            else:
                # case 1: atom
                # convert it to Mol() rather than Atom() to preserve
                # its properties
                sub_mol = Chem.MolFromSmiles(substituent)
                if main_mol.GetAtoms()[0].GetSymbol() != "Te":
                    new_mol, main_index = self.label.get_index(main_mol)
                    main_index += 1
                else:
                    new_mol = self.label.delete_zz(main_mol)
                    main_index = 1
            # connect main core to substituent through a single bond
            combo = Chem.CombineMols(sub_mol, new_mol)
            edcombo = Chem.EditableMol(combo)
            print(f"edcombo: {Chem.MolToSmiles(edcombo.GetMol())}")
            print(f"sub_index: {sub_index}, main_index: {main_index}")
            single = Chem.rdchem.BondType.SINGLE
            edcombo.AddBond(sub_index, main_index, order=single)
            final_mol = edcombo.GetMol()
        else:
            # substituent is hydrogen: just delete "Te" in the core
            final_mol = self.label.delete_zz(new_mol)
        return final_mol


    def run(self, inchiplus, main_mol):

        # Algorithm
        # Get first list of substituents and substitute in order
        inchiplus_item = inchiplus.pop(0)  # first R group or change
        # Treat attachment cases
        attachments = []
        # molecule will help us understand if we have replacement "-" or
        # variable attachment "," case
        molecule = inchiplus_item.split("!")[0].split("/")[0]
        if "," in molecule:
            # Get a clean list of attachments e.g. ["1H", "2H", "4H"...etc]
            # or ["1", "2", "4"...etc]
            attachments = molecule.split(",")
            attachments[-1] = attachments[-1].split("-")[0]  # remove suffix
            index = inchiplus_item.find("-")
            suffix = inchiplus_item[(index+1):]  # part after "-"
            suffix_list = suffix.split("!")  # R group replacements
            inchiplus_item = ""
            for number in attachments:
                if "H" in number:  # v. attachment of hydrogen case
                    main_atom = ""
                else:
                    # v. attachment of different atoms Example 3.ii
                    main_atom = suffix_list[0]+"@"
                # split substituent to several sub_substiutions and put it
                # back in inchiplus_item
                for suff in suffix_list:
                    sub = number + "-" + main_atom + suff
                    inchiplus_item += (sub+"!")
            inchiplus_item = inchiplus_item[:-1]
        else:
            # only atom replacement case
            if "-" in molecule:
                # rank = canonical inchi number of atom to be replaced
                rank = molecule.split("-")[0]
                index = inchiplus_item.find("-")
                suffix = inchiplus_item[index+1:]  # part after "-"
                suffix_list = suffix.split("!")  # list of atoms
                atom = suffix_list[0]  # atom to be replaced
                inchiplus_item = ""
                for suff in suffix_list:
                    # e.g. 1-C@N (change C with canonical no. 1 to N)
                    sub = rank + "-" + atom + "@" + suff
                    inchiplus_item += (sub+"!")
                inchiplus_item = inchiplus_item[:-1]  # remove "!" at the end
        # split instances of R group separated by "!"
        grouplist=inchiplus_item.split("!")
        print(f"inchiplus: {inchiplus}")
        print(f"grouplist: {grouplist}")
        for substituent in grouplist:
            if substituent == "H":  # make H implicit
                substituent = ""
            new_mol = copy.deepcopy(main_mol)
            new_inchi = ""
            if "-" in substituent.split("/")[0]:
                # replace normal atom "-"
                new_mol = self.replacement(new_mol, substituent)
            else:
                # replace Zz atom
                new_mol = self.replace(new_mol, substituent)
            if len(inchiplus) == 0:  # if finished substitutions
                # produce inchi
                new_mol = copy.deepcopy(self.label.sanitize_labels(self.ranks,
                                                                  self.inchi,
                                                                  new_mol))
                new_inchi = Chem.MolToInchi(self.label.sanitize(new_mol))
                # if it is a new inchi then put it in the list
                if new_inchi not in self.list_of_inchi:
                    self.list_of_inchi.append(new_inchi)
            else:
                # There is still more substitutions then continue
                # with current mol and inchiplus
                new_inchiplus = inchiplus.copy()
                self.run(new_inchiplus, new_mol)
        return

if __name__=="__main__":
    # Only run the code below if this class is run directly by python not gui
    # Correct InChI syntax to work with rdkit
    inchi_obj = MarkInChI()
    inchi = input("Please enter the MarkInChI: ")
    inchi = inchi.replace("MarkInChI", "InChI", 1)
    inchi = inchi.replace("Zz", "Te")
    # Check it is actually a markinchi
    if inchi.find("<>") == -1:
        print("Not a MarkInChI")
    else:
        inchiplus = inchi.split("<>")
        # Get main inchi and substituents
        inchiplus_item = inchiplus.pop(0)  # inchiplus_item = main_inchi
        inchi_obj.inchi = inchiplus_item  # store original inchi
        inchi_obj.list_of_inchi = []  # list of produced single inchis
        # create an instance of the labelling class
        inchi_obj.label = label = Label()
        # isotopically label the inchi where replacements will occur
        inchiplus_item, inchi_obj.ranks = label.label_inchi(inchiplus_item,
                                                            inchiplus)
        if len(inchiplus) > 1 and inchi.find("Te") != -1:  # needs reordering
            inchiplus = inchi_obj.reorder(inchiplus_item, inchiplus)
        # Convert main inchi to mol
        main_mol = Chem.rdinchi.InchiToMol(inchiplus_item)[0]
        # Sanitize (only in rdkit)
        new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(main_mol))
        # run alogrithm and then print the resulted list of single inchis
        inchi_obj.run(inchiplus, new_mol)
        print(inchi_obj.list_of_inchi)
