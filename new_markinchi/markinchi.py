from rdkit import Chem

class MarkInChI(object):

    def find_atom(self, rank, formula):

        # This function finds an atom given its canonical
        # position
        atoms = []
        dict = {}
        atom = ""
        no = ""
        index = 0
        for char in formula:
            if char.isupper():
                if atom not in dict.keys() and atom != "":
                    dict[atom] = "1"
                no = ""
                atom = char
            if char.islower():
                atom += char
            if not char.isalpha():
                no += char
                dict[atom] = no
        else:
            if no == "" and atom not in dict.keys():
                dict[atom]="1"
        for key in dict.keys():
            if key != "H":
                i = int(dict[key])
                atoms += [key]*i
        return atoms[int(rank)-1]

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
        if no_atoms != "":
            no = int(no_atoms)
        # Get the rank of each Zz atom
        line = inchi.split("/")[2].replace("c", "")
        list_of_no = []
        for n in line.split("-"):
            if "(" in n:
                n_list = n.replace("(", "-").replace(")", "-").split("-")
                list_of_no += n_list
            else:
                list_of_no.append(n)
        list_int = set(map(int, list_of_no))
        start_zz = max(list_int)-no+1
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

    def replacement(self, smiles_core, substituent):

        # This function performs replacements on atoms "-"
        iso_num = 0
        id = substituent.split("-")
        if "H" not in substituent:
            rank, main_atom, atom = tuple(id[0]) + tuple(id[1].split("@"))
        else:
            rank, atom = tuple(id)
        inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles_core))
        # Label the atom
        if inchi.find("/i") == -1:  # no isotopic layer
            if inchi.find("/f") == -1:  # no fixed layer
                if inchi.find("/r") == -1:  # no organometallic layer
                    inchi += "/i"+rank.split("H")[0]+"+50"  # case 1
                else:
                    # case 3
                    isotope = "/i"+rank.split("H")[0]+"+50"
                    index = inchi.find("/r")
                    inchi = inchi[:index]+isotope+inchi[index:]
            else:  # there is fixed layer
                # case 2
                isotope = "/i"+rank.split("H")[0]+"+50"
                index = inchi.find("/f")
                inchi = inchi[:index]+isotope+inchi[index:]
        else:
            index = inchi.find("/i")
            indexf = inchi[index+2:].find("/")
            sub = ""
            if indexf == -1:
                # isotopic layer is last layer
                sub = ","+inchi[index+2:]
            else:
                # isotopic layer not last layer
                sub = ","+inchi[index+2:index+2+indexf]
            for i in range(0, len(sub)):
                if sub[i] == ",":
                    num = sub[i:].split(",")[1].split("+")[0]
                    if num == rank:
                        iso_num = int(sub[i:].split(",")[1].split("+")[0])
                        index2 = sub[i+1:].find(",")
                        if index2 != -1:
                            isotope = rank.split("H")[0]+"+50"
                            inchi = inchi[:index+1+i]+isotope+inchi[index2+i:]
                            break
                        else:
                            isotope = rank.split("H")[0]+"+50"
                            inchi = inchi[:index+1+i]+isotope+inchi[index+2+len(sub):]
                            break
                    else:
                        if int(num) > int(rank.split("H")[0]):
                            isotope = rank.split("H")[0]+"+50,"
                            inchi = inchi[:index+2+i]+isotope+inchi[:index+2+i]
                            break
            else:  # if all numbers are lower we just add to the end of sub
                isotope = ","+rank.split("H")[0]+"+50"
                sub += isotope
                inchi = inchi[:index]+sub+inchi[index+2+indexf:]
        # At this point, the atom is labelled
        # Now we should add a substituent if there is hydrogen
        lab_smiles = Chem.MolToSmiles(Chem.rdinchi.InchiToMol(inchi)[0])
        atom1 = self.find_atom(rank.split("H")[0], inchi.split("/")[1])
        table = Chem.GetPeriodicTable()
        atomic_mass = int(table.GetMostCommonIsotopeMass(atom1))
        # determine if atom is lower_case or upper_case
        upper = not lab_smiles.find(str(atomic_mass+50)+atom1) == -1
        if rank.find("H") != -1:
            sub_new = ""
            sub_old = ""
            index = 0
            if upper:
                sub_new = str(atomic_mass+iso_num)+atom1
                sub_old = str(atomic_mass+50)+atom1
                index = inchi.find(str(atomic_mass+50)+atom1)
            else:
                sub_new = str(atomic_mass+iso_num)+str.lower(atom1)
                sub_old = str(atomic_mass+50)+str.lower(atom1)
                index = lab_smiles.find(str(atomic_mass+50)+str.lower(atom1))
            lab_smiles = lab_smiles.replace(sub_old, sub_new)
            # remove their isotomic mass from a smiles if they exist
            lab_smiles = Chem.MolToSmiles(Chem.rdinchi.InchiToMol(
                                          Chem.MolToInchi(Chem.MolFromSmiles(lab_smiles)))[0])
            index2 = lab_smiles[index:].find("]")
            lab_smiles = lab_smiles[:index+index2+1]+"("+atom+")"+lab_smiles[index+index2+1:]

        else:
            if not upper:
                atom = str.lower(atom)
                atom1 = str.lower(atom1)

            index = lab_smiles.find("["+str(atomic_mass+50)+atom1)
            index2 = lab_smiles[index:].find("]")
            lab_smiles = lab_smiles[:index]+atom+lab_smiles[index+index2+1:]

        """
        line = inchi.split("/")[2][1:]
        list_of_no = []
        for n in line.split("-"):
            if "(" in n:
                n_list = n.replace("(", "-").replace(")", "-").split("-")
                list_of_no += n_list
            else:
                list_of_no.append(n)
        i = 1  # position in smiles
        for number in list_of_no:
            if rank[0] == number:
                break
            i += 1
        r = 0
        count = 0
        if "H" in rank:
            atom = self.find_atom(rank[0], inchi.split("/")[1])
            for char in smiles_core:
                if char.casefold() == atom.casefold():
                    r += 1
                    if r == i:
                        sc = smiles_core
                        smiles_core = sc[:count]+"("+main_atom+")"+sc[count:]
                count += 1
        else:
            for char in smiles_core:
                if char.casefold() == main_atom.casefold():
                    r += 1
                    if r == i:
                        if char.islower():
                            sc = smiles_core
                            smiles_core = sc[:count]+atom.lower()+sc[count+1:]
                        else:
                            sc = smiles_core
                            smiles_core = sc[:count]+atom.upper()+sc[count+1:]
                count += 1
        """
        return lab_smiles


    def replace(self, smiles_core, substituent):

        # This function replaces undefined atoms
        # Case not hydrogen
        if substituent != "":
            sub_inchi = "InChI=1B/"+substituent
            sub_mol = Chem.rdinchi.InchiToMol(sub_inchi)
            sub_smiles = Chem.MolToSmiles(sub_mol[0])
            sub_smiles = sub_smiles.replace("[Te]", "", 1)
            sub_smiles = sub_smiles.replace("[",
                                            "").replace("]", "")
            smiles_core = smiles_core.replace("[Te]", sub_smiles, 1)
        else:
            smiles_core = smiles_core.replace("[Te]", "")
        smiles_core = smiles_core.replace("=", "").replace("()", "")
        return smiles_core


    def run(self, inchiplus, smiles_core):

        self.no_run += 1
        # Algorithm
        # Get first list of substituents and substitute in order
        inchiplus_item = inchiplus.pop(0)
        # Treat attachment cases
        attachments = []
        molecule = inchiplus_item.split("!")[0].split("/")[0]
        if "," in molecule:
            attachments = molecule.split(",")
            attachments[-1] = attachments[-1].split("-")[0]
            index = inchiplus_item.find("-")
            suffix = inchiplus_item[(index+1):]
            suffix_list = suffix.split("!")
            inchiplus_item = ""
            for number in attachments:
                if "H" in number:
                    main_atom = ""
                else:
                    main_atom = suffix_list[0]+"@"
                for suff in suffix_list:
                    sub = number + "-" + main_atom + suff
                    inchiplus_item += (sub+"!")
            inchiplus_item = inchiplus_item[:-1]
        else:
            if "-" in molecule:
                rank = molecule.split("-")[0]
                index = inchiplus_item.find("-")
                suffix = inchiplus_item[index+1:]
                suffix_list = suffix.split("!")
                atom = suffix_list[0]
                inchiplus_item = ""
                for suff in suffix_list:
                    sub = rank + "-" + atom + "@" + suff
                    inchiplus_item += (sub+"!")
                inchiplus_item = inchiplus_item[:-1]
        grouplist=inchiplus_item.split("!")
        for substituent in grouplist:
            if substituent == "H":
                substituent = ""
            new_smiles=""
            if "-" in substituent.split("/")[0]:
                new_smiles = self.replacement(smiles_core, substituent)
                print(f"Number of runs: {self.no_run}")
                print(f"new_smiles: {new_smiles}")
            else:
                new_smiles=self.replace(smiles_core, substituent)
                print(f"Number of runs: {self.no_run}")
                print(f"new_smiles: {new_smiles}")
            if len(inchiplus) == 0:
                new_inchi=Chem.MolToInchi(Chem.MolFromSmiles(new_smiles))
                if new_inchi not in self.list_of_inchi:
                    self.list_of_inchi.append(new_inchi)
                    print(new_inchi)
            else:
                new_inchiplus = inchiplus.copy()
                self.run(new_inchiplus, new_smiles)
        return


#if True:# __name__=="__main__":
if __name__=="__main__":
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
        inchiplus_item = inchiplus.pop(0)
        if len(inchiplus) > 1 and inchi.find("Te") != -1:  # need reordering
            inchiplus = inchi_obj.reorder(inchiplus_item, inchiplus)
        # Convert main inchi to smiles and run algorithm
        print(inchiplus_item)
        molfrominchi = Chem.rdinchi.InchiToMol(inchiplus_item)
        smiles_core=Chem.MolToSmiles(molfrominchi[0])
        inchi_obj.list_of_inchi = []
        inchi_obj.no_run = 0
        inchi_obj.run(inchiplus, smiles_core)
        print(inchi_obj.list_of_inchi)
