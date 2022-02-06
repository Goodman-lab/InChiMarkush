import copy
from rdkit import Chem
from zz_convert import zz_convert

class markmol(object):

    def convert(self, content):

        # This function converts a normal markush SDF file to RDKIT compatible
        # SDF file. Input: list. Output: list.
        new_content = self.replaceR(content)
        new_content = self.delete(new_content)
        new_content = self.add(new_content)
        new_content = self.label(new_content)
        return copy.deepcopy(new_content)

    def replaceR(self, content):

        # This function replaces R groups 'R#' in the file by Te atoms 'Te'.
        #
        ctab = 0
        new_content = []
        replace = True
        for line in content:
            if replace:
                new_line = line.replace("R#", "Te")
            else:
                # delete after $RGP
                new_line = " "
                replace = True
                self.ctabs.append(ctab)
            if line.find("$END CTAB") != -1:
                ctab += 1
            if line.find("$RGP") != -1:
                replace = False
            if line.find("M  RGP") != -1:
                # get SN of R gourps
                parts = line.split("RGP")[1].split()
                for i in range(1, len(parts), 2):
                    self.Rpositions.append(parts[i])
            if line.find("M  APO") != -1:
                # find points of connections for each substituent
                connection = line.split()[3]
                self.connections.append(connection)

            new_content.append(new_line)
        for line in new_content: print(line)
        return copy.deepcopy(new_content)

    def delete(self, content):

        new_content = []
        for line in content:
            # conditions a & b
            a = (not line.isspace()) and line[0] != "$"
            b = not (line[0] == "M" and (line.find("END") == -1))
            if a and b:
                new_content.append(line)
        for line in new_content: print(line)
        return copy.deepcopy(new_content)

    def add(self, content):

        content[0] = "\n\n\n"
        new_content = []
        add_line = "$$$$\n\n\n\n"
        i = 0
        count = 1
        for line in content:
            if i > 0 and count == int(self.connections[i-1]):
                # must be 2D structure
                if line.find("0.0000") != -1:
                    index = line.find("0.0000")
                    new_line = line[:index+11]+"8"+line[index+12:]
                    new_content.append(new_line)
                else:
                    new_content.append(line)
            else:
                new_content.append(line)
            if line[0] == "M" and line.find("END") != -1:
                count = 1
                i += 1
                new_content.append(add_line)
            if line.find(".") != -1:
                count += 1

        for line in new_content: print(line)
        return copy.deepcopy(new_content)

    def label(self, content):

        n = len(self.Rpositions)
        p = str(n)
        iso_line = f"M  ISO"
        iso_line += (3-len(p))*" "+p
        i = 1
        print(f"Rpositions: {self.Rpositions}")
        for p in self.Rpositions:
            iso_line += (4-len(p))*" "+p+" "+str(130+i)
            i += 1
        iso_line += "\n"
        new_content = []
        replace = True
        for line in content:
            if line.find("M  END") != -1 and replace:
                new_content.append(iso_line)
                replace = False
            new_content.append(line)
        for line in new_content: print(line)
        return copy.deepcopy(new_content)

    def relabel_core(self, inchi):

        # for now we just delete /i layer
        index = inchi.find("/i")
        new_inchi = ""
        if index != -1:
            suf = ""
            if inchi[index+2:].find("/") != -1:
                suf = "/"+"/".join(inchi[index+2:].split("/")[1:])
            new_inchi = inchi[:index]+suf
        else:
            new_inchi = inchi
        return new_inchi

    def relabel_sub(self, mol):

        # returns inchi for > 1 main atom and smiles otherwise
        sub_inchi = ""
        for atom in mol.GetAtoms():
            if atom.GetIsotope() > 0:
                idx = atom.GetIdx()
                rwmol = Chem.RWMol(mol)
                single = Chem.rdchem.BondType.SINGLE
                rwmol.GetAtoms()[idx].SetIsotope(0)
                if mol.GetNumAtoms() > 1:
                    rwmol.AddAtom(Chem.Atom("Te"))
                    rwmol.AddBond(idx, mol.GetNumAtoms(), order = single)
                    sub_inchi = Chem.MolToInchi(rwmol)
                    # check sub_inchi is valid
                    if Chem.rdinchi.InchiToMol(sub_inchi) == None:
                        raise RuntimeError("Sub_inchi invalid")
                    sub_inchi = "/".join(sub_inchi.split("/")[1:])
                else:
                    sub_inchi = Chem.MolToSmiles(rwmol)
        return sub_inchi
    def produce_markinchi(self):

        zz = zz_convert()
        core_mol = self.core_mol
        core_inchi = Chem.MolToInchi(core_mol)
        order = []
        index = core_inchi.find("/i")
        label_part = ""
        if index != -1:
            indexf = core_inchi[index+2:].find("/")
            if indexf != -1:
                label_part = core_inchi[index+2:index+2+indexf]
            else:
                label_part = core_inchi[index+2:]
        print(label_part)
        for part in label_part.split(","):
            order.append(str(int(part.split("+")[1])-2))
        print(order)
        print(core_inchi)
        mark_inchi = self.relabel_core(zz.te_to_zz(core_inchi))
        Rsubstituents = self.Rsubstituents
        for num in order:
            subs = Rsubstituents[int(num)-1]
            print(f"len(subs): {len(subs)}")
            mark_inchi += "<>"
            for sub in subs:
                print(2)
                sub_inchi = self.relabel_sub(sub)
                if sub_inchi.find("Te") != -1:
                    sub_inchi = zz.te_to_zz("InChI=1B/"+sub_inchi)
                    sub_inchi = "/".join(sub_inchi.split("/")[1:])
                print(f"sub_inchi: {sub_inchi}")
                mark_inchi += sub_inchi+"!"
            mark_inchi = mark_inchi[:-1]
        mark_inchi = mark_inchi.replace("InChI=1S/", "MarkInChI=1B/")
        mark_inchi = mark_inchi.replace("[HH]", "H")
        return mark_inchi



if __name__=="__main__":
    name = input("Enter the file name with the extension:")
    file = open(name, "r")
    content = file.readlines()
    mark_obj = markmol()
    zz = zz_convert()
    mark_obj.Rpositions = []  # number of each Te atom
    mark_obj.Rsubstituents = []
    mark_obj.ctabs = []  # no. of each substituent after $RGP
    mark_obj.connections = []
    content = mark_obj.convert(content)
    new_name = name.split(".")[0]+"_RDKIT.sdf"
    new_file = open(new_name, "w")
    new_file.writelines(content)
    new_file.close()
    supply = Chem.SDMolSupplier(new_name)
    substituents = []
    for mol in supply:
        if mol is not None:
            print(Chem.MolToSmiles(mol))
            substituents.append(mol)
    mark_obj.core_mol = copy.deepcopy(supply[0])
    i = 1
    ctabs = mark_obj.ctabs
    while i < len(ctabs):
        mark_obj.Rsubstituents.append(substituents[int(ctabs[i-1]):int(ctabs[i])])
        i += 1
    mark_obj.Rsubstituents.append(substituents[int(ctabs[i-1]):])
    print(mark_obj.Rsubstituents)
    print(mark_obj.produce_markinchi())
    file.close()
