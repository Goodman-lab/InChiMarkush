import sys, os
from rdkit import Chem
import copy
from label import Label
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'markmol2markinchi'))
from zz_convert import zz_convert

class MarkInChI(object):

    def __init__(self, inchi):

        zz = zz_convert()
        inchi = inchi.replace("MarkInChI", "InChI", 1)
        # Check it is actually a markinchi
        if inchi.find("<M>") == -1:
            print("Not a MarkInChI")
        else:
            inchiplus = inchi.split("<M>")
            # Get main inchi and substituents
            inchiplus_item = inchiplus.pop(0)  # inchiplus_item = main_inchi
            if inchiplus_item.find("Zz") != -1:
                inchiplus_item = zz.zz_to_te(inchiplus_item)
            # TODO: use zz_to_te instead
            for sub in range(0, len(inchiplus)):
                inchiplus[sub] = inchiplus[sub].replace("Zz", "Te")
            self.inchi = inchiplus_item  # store original inchi
            self.list_of_inchi = []  # list of produced single inchis
            # create an instance of the labelling class
            self.label = label = Label()
            # isotopically label the inchi where replacements will occur
            inchiplus_item, self.ranks = label.label_inchi(inchiplus_item,
                                                                inchiplus)
            print(f"after labelling: {inchiplus_item}")
            print(f"ranks: {self.ranks}")
            # Convert main inchi to mol
            print(f"core_inchi: {inchiplus_item}")
            main_mol = Chem.rdinchi.InchiToMol(inchiplus_item)[0]
            # Sanitize (only in rdkit)
            new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(main_mol))
            # run alogrithm and then print the resulted list of single inchis
            self.run_count = 0
            self.run(inchiplus, new_mol)
            print(self.list_of_inchi)
            print(f"Number of inchi produced: {len(self.list_of_inchi)}")

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
                print(f"ranks: {self.ranks}")
                print(f"num: {num}")
                print(f"symbol: {atom}")
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
                # smiles = Chem.MolToSmiles(new_mol)
                # smiles.replace("[C]","C")
                # new_mol = Chem.MolFromSmiles(smiles)
                rwmol = Chem.RWMol(new_mol)
                new_rwmol = copy.deepcopy(rwmol)
                for mol_atom in rwmol.GetAtoms():
                    cn1 = mol_atom.GetIsotope() == num
                    cn2 = mol_atom.GetSymbol() == str.upper(atom)
                    if cn1 and cn2:
                        # Replace the labelled atom to og and add new atom
                        iso_num = self.label.find_isotope(self.inchi, rank[:-1])
                        idx = mol_atom.GetIdx()
                        is_aromatic = mol_atom.GetIsAromatic()
                        add_index = rwmol.GetNumAtoms()
                        new_rwmol.GetAtoms()[idx].SetIsotope(iso_num)
                        new_rwmol.AddAtom(replacement)
                        single = Chem.rdchem.BondType.SINGLE
                        new_rwmol.AddBond(idx, add_index, order=single)
                        new_mol = copy.deepcopy(new_rwmol.GetMol())
        print(f"new_mol: {Chem.MolToSmiles(new_mol)}")
        new_mol = self.label.sanitize_charges(new_mol)
        return copy.deepcopy(new_mol)

    def replace(self, main_mol, substituent):
        # This function replaces undefined atoms
        # Case not hydrogen
        final_mol = None
        if substituent != "" and substituent != "H":
            sub_inchi = ""
            sub_mol = None
            # We assume that the markinchi will be either:
            # 1- one atom that is convertable to a smiles and then to
            # an Atom()
            # 2- more than one atom written in the format of a "part-inchi"
            if substituent.find("/") != -1:
                # case 2
                sub_inchi = "InChI=1B/"+substituent  # convert to inchi
                sub_mol = Chem.rdinchi.InchiToMol(sub_inchi)[0]
                final_mol = self.label.combine(main_mol, sub_mol)
            else:
                # case 1: atom
                # convert it to Mol() rather than Atom() to preserve
                # its properties
                sub_mol = Chem.MolFromSmiles(substituent)
                final_mol = self.label.combine(main_mol, sub_mol)
        else:
            # substituent is hydrogen: just delete "Te" in the core
            final_mol = self.label.delete_zz(main_mol)
        return final_mol

    def run(self, inchiplus, main_mol):

        # Algorithm
        # Get first list of substituents and substitute in order
        self.run_count += 1
        print(f"inchiplus: {inchiplus} - {self.run_count}")
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
        print(f"grouplist: {grouplist}")
        for substituent in grouplist:
            if substituent == "H":  # make H implicit
                substituent = ""
            new_mol = copy.deepcopy(main_mol)
            new_inchi = ""
            print(f"subsitutent: {substituent}")
            print(f"before substitution: {Chem.MolToSmiles(new_mol)}")
            if "-" in substituent.split("/")[0]:
                # replace normal atom "-"
                new_mol = self.replacement(new_mol, substituent)
                print(f"after substitution: {Chem.MolToSmiles(new_mol)}")
            else:
                # replace Zz atom
                new_mol = self.replace(new_mol, substituent)
                print(f"after substitution: {Chem.MolToSmiles(new_mol)}")
            if len(inchiplus) == 0:  # if finished substitutions
                # produce inchi
                new_mol = copy.deepcopy(self.label.sanitize_labels(self.ranks,
                                                                  self.inchi,
                                                                  new_mol))
                new_inchi = Chem.MolToInchi(self.label.sanitize(new_mol))
                print(f"new_inchi: {new_inchi}")
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
    inchi = input("Please enter the MarkInChI: ")
    inchi_obj = MarkInChI(inchi)
