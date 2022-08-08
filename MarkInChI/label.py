import sys

from rdkit import Chem
import copy
import numpy
from helper import helper

# This module will be used to label atoms as isotopes.
# This will also be able to find the position of an atom given its label
# This will also contain other helpful functions

class Label(object):

    def __init__(self):
        self.helper = helper()
    def label(self, inchi, rank, lab):
        # rank, lab are strings
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
            so labels can be obtained easily. """

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
                        else:
                            raise RuntimeError("Same atom replaced more than once")
                else:
                    # only replacements with no variable attachments
                    num = molecule.split("-")[0].split("H")[0]
                    print(f"num: {num}")
                    table = Chem.GetPeriodicTable()
                    atom = self.find_atom(num, inchi.split("/")[1])
                    atomic_mass = int(table.GetMostCommonIsotopeMass(atom))
                    if num not in ranks.keys():
                        ranks[num] = int(num)+10+atomic_mass
                        new_inchi = self.label(new_inchi, num, str(int(num)+10))
                        print(f"new_inchi: {new_inchi}")
                    else:
                        raise RuntimeError("Same atom replaced more than once")
        # Now we label Te atoms:
        if new_inchi.find("Te") != -1:
            no, start_zz = self.helper.zz_no(new_inchi)
            for i in range(0, no):
                new_inchi = self.label(new_inchi, str(i+start_zz),
                                       str(128+10+i+1))
        return new_inchi, ranks

    def delete_zz(self, mol):

        # This function delete the instance of Zz with the lowest rank
        # and returns the new mol
        index = 0
        min_rank = 10000
        new_mol = copy.deepcopy(self.sanitize(mol))
        rwmol = Chem.RWMol(new_mol)
        for atom in rwmol.GetAtoms():
            cn1 = atom.GetIsotope() < min_rank
            cn2 = atom.GetSymbol() == "Te"
            if  cn1 and cn2:
                min_rank = atom.GetIsotope()
                index = atom.GetIdx()
        rwmol.RemoveAtom(index)
        new_mol = rwmol.GetMol()
        return copy.deepcopy(self.sanitize(new_mol))

    def get_index(self, mol, label):

        # This function gets the index of the atom connected to Te with the
        # lowest rank and deletes Te
        index = 0
        index1 = 1
        new_mol = copy.deepcopy(self.sanitize(mol))
        min_rank = 100000000
        rwmol = Chem.RWMol(mol)
        for atom in rwmol.GetAtoms():
            if atom.GetSymbol() == "Te" and atom.GetIsotope() < min_rank:
                min_rank = atom.GetIsotope()
                index1 = atom.GetIdx()
                num = rwmol.GetNumAtoms()-1
                for i in reversed(range(0, num)):
                    if mol.GetBondBetweenAtoms(i, index1) != None:
                        index = i
        pre_label = rwmol.GetAtoms()[index].GetIsotope()
        symbol = rwmol.GetAtoms()[index].GetSymbol()
        # label connected atom by prelabel+label
        post_label = 0
        if pre_label == 0:
            table = Chem.GetPeriodicTable()
            atomic_mass = int(table.GetMostCommonIsotopeMass(symbol))
            post_label = atomic_mass+label
        else:
            post_label = pre_label+label
        rwmol.GetAtoms()[index].SetIsotope(post_label)
        rwmol.RemoveAtom(index1)
        new_mol = rwmol.GetMol()
        new_index = 0
        for atom in new_mol.GetAtoms():
            if atom.GetSymbol == symbol and atom.GetIsotope() == post_label:
                new_index = atom.GetIdx()
        return copy.deepcopy(new_mol), post_label

    def combine(self, main_mol, sub_mol, num = None):

        print("BAF0")
        # label 30 for sub_mol, and label 35 for main_mol
        sub_label = 0
        if sub_mol.GetNumAtoms() > 1:
            print("BAF1")
            sub_mol, sub_label = self.get_index(sub_mol, 30)
        else:
            print("BAF2")
            pre_label = sub_mol.GetAtoms()[0].GetIsotope()
            post_label = 0
            if pre_label == 0:
                print("BAF3")
                table = Chem.GetPeriodicTable()
                atom = sub_mol.GetAtoms()[0].GetSymbol()
                atomic_mass = int(table.GetMostCommonIsotopeMass(atom))
                post_label = atomic_mass+30
            else:
                post_label = pre_label+30
            print("BAF4")
            sub_mol.GetAtoms()[0].SetIsotope(post_label)
            sub_label = post_label

        if num != None:
            main_label = num
            new_mol = main_mol
        else:
            new_mol, main_label = self.get_index(main_mol, 35)
        combo = Chem.CombineMols(sub_mol, new_mol)
        edcombo = Chem.EditableMol(combo)
        single = Chem.rdchem.BondType.SINGLE
        sub_index = 0
        main_index = 0
        back_mol = edcombo.GetMol()
        for atom in back_mol.GetAtoms():
            if atom.GetIsotope() == sub_label:
                sub_index = atom.GetIdx()
                atom.SetIsotope(sub_label-30)
            if atom.GetIsotope() == main_label:
                main_index = atom.GetIdx()
                if num == None:
                    atom.SetIsotope(main_label-35)
                else:
                    atom.SetIsotope(0)
        edcombo = Chem.EditableMol(back_mol)
        edcombo.AddBond(sub_index, main_index, order=single)
        final_mol = self.sanitize(edcombo.GetMol())
        return copy.deepcopy(final_mol)


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
        # This also sanitize isotopic labels eg. [ch2] -> c
        new_mol = copy.deepcopy(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
        for atom in new_mol.GetAtoms():
            table = Chem.GetPeriodicTable()
            symbol = atom.GetSymbol()
            atomic_mass = int(table.GetMostCommonIsotopeMass(symbol))
            if atom.GetIsotope() == atomic_mass:
                atom.SetIsotope(0)
        return new_mol

    def sanitize_charges(self, mol):


        """
        list_of_indices = []
        item = "["
        for index, elem in enumerate(list(smiles)):
            if elem == item:
                list_of_indices.append(item)
        i = 0
        """

        mol.UpdatePropertyCache(strict=False)
        table = Chem.GetPeriodicTable()
        for atom in mol.GetAtoms():
            atomic_number = atom.GetAtomicNum()
            valence = atom.GetTotalValence()
            default_valence = table.GetDefaultValence(atomic_number)
            electrons = table.GetNOuterElecs(atomic_number)
            if valence != default_valence:
                print(f"Warning: atom {atom.GetSymbol} not in default valence")
            num = (electrons-valence)
            charge = int((num%2)*numpy.sign(num))
            if atom.GetSymbol() != "C":
                atom.SetFormalCharge(charge)
        return copy.deepcopy(mol)

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
