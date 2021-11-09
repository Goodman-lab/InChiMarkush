from rdkit import Chem
import copy

# This module will be used to label atoms as isotopes.
# This will also be able to find the position of an atom given its label
# This will also contain other helpful functions

class Label(object):
    def label(self, inchi, rank, lab):
        # Label the atom
        if inchi.find("/i") == -1:  # no isotopic layer
            if inchi.find("/f") == -1:  # no fixed layer
                if inchi.find("/r") == -1:  # no organometallic layer
                    # just add label to the end of inchi
                    inchi += "/i"+rank.split("H")[0]+"+"+lab
                else:
                    # no /i, no /f but contains /r
                    isotope = "/i"+rank.split("H")[0]+"+"+lab
                    index = inchi.find("/r")
                    inchi = inchi[:index]+isotope+inchi[index:]
            else:
                # no /i but there is /f or /r or both
                isotope = "/i"+rank.split("H")[0]+"+"+lab
                index = inchi.find("/f")
                inchi = inchi[:index]+isotope+inchi[index:]
        else:
            # there is isotopic layer in inchi
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
                # scan the isotopic layer to add the label at the right atom
                if sub[i] == ",":
                    num = sub[i:].split(",")[1].split("+")[0]
                    if num == rank:  # if atom is labelled
                        iso_num = int(sub[i:].split(",")[1].split("+")[0])
                        index2 = sub[i+1:].find(",")
                        if index2 != -1:
                            isotope = rank.split("H")[0]+"+"+lab
                            inchi = inchi[:index+2+i]+isotope+inchi[index2+i:]
                            break
                        else:
                            isotope = rank.split("H")[0]+"+"+lab
                            inchi = inchi[:index+2+i]+isotope+inchi[index+1+len(sub):]
                            break
                    else:  # if atom is not labelled
                        if int(num) > int(rank.split("H")[0]):
                            isotope = rank.split("H")[0]+"+"+lab+","
                            inchi = inchi[:index+2+i]+isotope+inchi[index+2+i:]
                            break
            else:
                # if all numbers are lower we just add to the end of
                # the isotopic layer
                isotope = ","+rank.split("H")[0]+"+"+lab
                sub += isotope
                if indexf == -1:
                    # case: no other layer after isotopic layer
                    inchi = inchi[:index+2]+sub[1:]
                else:
                    # case: there is a layer after isotopic layer
                    inchi = inchi[:index+2]+sub[1:]+inchi[index+2+indexf:]
        return inchi

    def label_inchi(self, inchi, inchiplus):

        """ This function is used to obtain canonical labelling for the
            atoms that will be replaced as the inchi canonical labelling
            is going to change after every substitution. Each atom will
            be labelled by a label where:
            label = rank + 10 + atomic mass of the atom
            10 is added to not get an isotope that exists for low ranks and
            using adding rank will make it canonical.
            This function returns the labelled inchi and a dictionary ranks
            so labells can be obtained easily. """

        ranks = {}  # dict {rank:label}
        new_inchi = inchi
        for sub in inchiplus:
            molecule = sub.split("!")[0].split("/")[0]
            if molecule.find("-") != -1:
                # substituent contain replacements
                if molecule.find(",") != -1:
                    # substituent contain variable attachments
                    for mini_sub in molecule.split(","):
                        # get ranks of atoms to be labelled
                        num = mini_sub.split("-")[0].split("H")[0]
                        table = Chem.GetPeriodicTable()
                        atom = self.find_atom(num, inchi.split("/")[1])
                        atomic_mass = int(table.GetMostCommonIsotopeMass(atom))
                        if num not in ranks.keys():
                            # add rank to dictionary ranks if doesn't exist
                            ranks[num] = int(num)+10+atomic_mass
                            new_inchi = self.label(new_inchi, num, str(int(num)+10))
                            print(new_inchi)
                        else:
                            raise RuntimeError("Same atom replaced more than once")
                else:
                    # only replacements with no variable attachments
                    num = molecule.split("-")[0].split("H")[0]
                    table = Chem.GetPeriodicTable()
                    atom = self.find_atom(num, inchi.split("/")[1])
                    atomic_mass = int(table.GetMostCommonIsotopeMass(atom))
                    if num not in ranks.keys():
                        ranks[num] = int(num)+10+atomic_mass
                        new_inchi = self.label(inchi, num, str(int(num)+10))
                    else:
                        raise RuntimeError("Same atom replaced more than once")
        return new_inchi, ranks

    def delete_zz(self, mol):

        # This function delete the first instance of Zz and return the new mol
        index = 0
        new_mol = copy.deepcopy(self.sanitize(mol))
        for atom in Chem.RWMol(new_mol).GetAtoms():
            if atom.GetSymbol == "Te":
                index = atom.GetIdx()
                break
        rwmol = Chem.RWMol(mol)
        rwmol.RemoveAtom(index)
        new_mol = rwmol.GetMol()
        return copy.deepcopy(new_mol)

    def get_index(self, mol):

        # This function gets the index of the atom connected to first Te
        # and deletes Te
        index = 0
        index1 = 1
        new_mol = copy.deepcopy(self.sanitize(mol))
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "Te":
                index1 = atom.GetIdx()
                for i in reversed(range(0, index1)):
                    if mol.GetBondBetweenAtoms(i, index1) != None:
                        index = i
        rwmol = Chem.RWMol(mol)
        rwmol.RemoveAtom(index1)
        new_mol = rwmol.GetMol()
        return copy.deepcopy(new_mol), index
    def find_atom(self, rank, formula):

        # This function finds an atom given its canonical
        # position and the chemical formula from the inchi
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
    def sanitize(self, mol):
        # there is some issues with RDKIT mol produced from inchi, while the
        # mol produced from SMILES seem to be fine.
        return copy.deepcopy(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))

    def find_isotope(self, inchi, rank):

        # this function finds the isotopic label of an atom in an inchi
        result = 0
        index = inchi.find("/i")
        if index != -1:
            indexf = inchi[index+2:].find("/")
            sub = ""
            if indexf == -1:
                sub = inchi[index+2:]
            else:
                sub = inchi[index+2:indexf]
            for part in sub.split(","):
                rank2, iso = tuple(part.split("+"))
                if rank2 == rank:
                    result = iso
        return result
    def sanitize_labels(self, ranks, inchi, mol):

        """ This function resets the fake isotopic labels of the atoms
            of an inchi to what they were before labelling"""

        rwmol = Chem.RWMol(mol)
        new_mol = copy.deepcopy(rwmol)
        for rank in ranks.keys():
            num = ranks[rank]
            atom = self.find_atom(rank, inchi.split("/")[1])
            iso_num = self.find_isotope(inchi, rank)
            for mol_atom in new_mol.GetAtoms():
                cn1 = mol_atom.GetIsotope() == num
                cn2 = mol_atom.GetSymbol() == str.upper(atom)
                idx = mol_atom.GetIdx()
                if cn1 and cn2:
                    new_mol.GetAtoms()[idx].SetIsotope(int(iso_num))
        return copy.deepcopy(new_mol.GetMol())
