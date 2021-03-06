import copy
from rdkit import Chem
from zz_convert import zz_convert
from label import Label

class markmol(object):

    """ This class is used for the conversion of a markush mol to a markush
        inchi."""

    def convert(self, content):

        # This function converts a normal markush SDF file to RDKIT compatible
        # SDF file. Input: list. Output: list.
        new_content = self.replaceR(content)
        new_content = self.var_attach(new_content)
        new_content = self.delete(new_content)
        new_content = self.add(new_content)
        new_content = self.label(new_content)
        return copy.deepcopy(new_content)

    def replaceR(self, content):

        # This function 1. replaces R groups 'R#' in the file by Te atoms 'Te'.
        # 2. gets list of atoms indices.
        # 3. get indices of each R group (e.g. R1 is atom 9, etc.)
        # 4. get the point of connections of R group substituents
        # 5. delete the number after "$RGP" line
        # 6. get the number of substituents per R group (using "ctab")
        ctab = 0
        new_content = []
        replace = True
        list_of_atoms = {}
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
            if line.find("M  ALS") != -1:
                index = line.split()[2]
                atom_list = line.split()[5:]
                list_of_atoms[index] = atom_list
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
        self.list_of_atoms = list_of_atoms
        return copy.deepcopy(new_content)

    def delete(self, content):

        # This function deletes all lines that are:
        # 1. empty space lines (a).
        # 2. start with '$' signs (a).
        # 3. all that start "M" that are not "M END" (b).
        # 4. 'F' atom list line that comes after the bonds' block (c).

        new_content = []
        for line in content:
            # conditions a & b & c (see above)
            a = (not line.isspace()) and line[0] != "$"
            b = not (line[0] == "M" and (line.find("END") == -1))
            c = not ((line.find("F") != -1) and (((len(line)-15)) % 4 == 0))
            if a and b and c:
                new_content.append(line)
        return copy.deepcopy(new_content)

    def add(self, content):

        # This function adds 3 empty lines before each substituent and 4 $ signs
        # (to make its own mol).
        # labels connection points on substituents.
        content[0] = "\n\n\n"
        new_content = []
        add_line = "$$$$\n\n\n\n"
        i = 0
        count = 1
        atom_line = False
        for line in content:
            # label the 'connection' atom on the substituent with isotopic
            # label 8.
            if i > 0 and count == int(self.connections[i-1]):
                # must be mol V2000
                if atom_line:
                    if len(line) > 68:
                        new_line = line[:35]+"8"+line[36:]
                        new_content.append(new_line)
                    else:
                        atom_line = False
                        new_content.append(line)
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
            if line.find("999 V2000") != -1:
                atom_line = True
        return copy.deepcopy(new_content)

    #### variable attachments
    def var_attach(self, content):

        # find the indices of the 'empty' atoms
        # find the attachments for each bond, store that information,
        # and link it to the right atom.
        # find the indices of the atoms bound to the 'empty' atoms
        atom_inds = []  # indices of empty atoms
        attachments = []
        bonds = []  # store the block of bonds
        real_bonds = {}
        atom_line = False
        i = 0
        bond_line = False
        for line in content:
            if atom_line:
                i += 1
                # get atom symbols in core
                if len(line) > 68:
                    if len(line[31:34].split()) > 0:
                        self.atom_symbols[i] = line[31:34].split()[0]
                    else:
                        self.atom_symbols[i] = "empty"
                    if line[31] == " ":
                        atom_inds.append(str(i))  # store empty atom index
                else:
                    # atoms' block is finished and bonds' block begins
                    atom_line = False
                    bond_line = True
            if bond_line:
                # store bonds
                if line[0] != "M":
                    # bond
                    bonds.append(line)
                else:
                    bond_line = False
            if line.find("M  SAL") != -1:
                # get variable attachment points
                attachments.append(line.split()[4:])
            if line[0] == "M" and line.find("END") != -1:
                # end for loop after core is finished
                break
            if line.find("999 V2000") != -1:
                atom_line = True

        # store the bond between the empty atom and the real attachment atom.
        for b in bonds:
            for ind in atom_inds:
                if b.split()[0] == ind:
                    real_bonds[ind] = b.split()[1]
                else:
                    if b.split()[1] == ind:
                        real_bonds[ind] = b.split()[0]

        # modify
        # substract the number of 'empty atoms'*2 from the number of atoms
        # AND the number of bonds on the top, delete 'empty atoms' and atoms
        # that are bound to them (variable-bond atoms), and corresponding bonds
        atom_ids = {}  # dict of number:atom.
        no = len(atom_inds)
        new_content = []
        allow_change = True
        atom_line = False
        to_be_deleted = list(real_bonds.values())+list(real_bonds.keys())
        count = 0
        file_count = 0
        no_atoms = 0
        for line in content:
            if atom_line:
                count += 1
                if len(line) > 68:
                    delete_line = False
                    list_count = 0
                    for d in to_be_deleted:
                        if str(count) == d:
                            delete_line = True
                            # get all real atoms
                            if 2*list_count < len(to_be_deleted):
                                atom_ids[count] = line[31:34].split()[0]
                        list_count += 1

                    if not delete_line:
                        new_line = line
                        for at in self.list_of_atoms.keys():
                            # change the L atom line to the first instance atom
                            # in the atoms' block.
                            if count == int(at):
                                sym = self.list_of_atoms[at][0]
                                nn = len(sym)-1
                                new_line = line.replace(("L"+nn*" "), sym)
                                self.atom_symbols[count] = sym
                        new_content.append(new_line)
                else:
                    atom_line = False
                    new_content.append(line)
            else:
                if line.find("999 V2000") != -1 and allow_change:
                    # adjust top core mol statement (e.g. no of atoms, etc.)
                    atom_line = True
                    no_atoms = int(line.split()[0])-2*no
                    no_bonds = int(line.split()[1])-no
                    atom_part = (3-len(str(no_atoms)))*" "+str(no_atoms)
                    bond_part = (3-len(str(no_bonds)))*" "+str(no_bonds)
                    print(f"beginning: {3-len(str(no_atoms))}")
                    new_line = atom_part+bond_part+"  "+line[8:]
                    allow_change = False
                    new_content.append(new_line)
                else:
                    if not atom_line and not allow_change:
                        if len(line) == 22:
                            delete_line = False
                            for dd in to_be_deleted:
                                cn1 = line.split()[0] == dd
                                cn2 = line.split()[1] == dd
                                if cn1 or cn2:
                                    delete_line = True
                            if not delete_line:
                                new_content.append(line)
                        else:
                            new_content.append(line)
                    else:
                        new_content.append(line)
            if line.find ("M  END") != -1:
                # when the core mol finishes append the rest of the mol file
                # and stop exit the for loop.
                new_content += content[file_count+1:]
                break
            file_count += 1
        print(f"atom_inds: {atom_inds}")
        print(f"attachments: {attachments}")
        print(f"real_bonds: {real_bonds}")
        print(f"atom_ids: {atom_ids}")
        print(f"atom_symbols: {self.atom_symbols}")
        self.attachments = attachments
        self.attach_ids = atom_ids
        self.no_atoms = no_atoms
        #new_content = self.label_attachments(copy.deepcopy(new_content))
        return copy.deepcopy(new_content)

    def label_attachments(self, content):

        # label all attachments (5+rank in mol file)
        all_attachments = []
        for att in self.attachments: all_attachments += att
        new_content = []
        file_count = 0
        atom_count = 0
        atom_line = False
        for line in content:
            if atom_line:
                if len(line) > 68:
                    atom_count += 1
                    new_line = line
                    for attach in all_attachments:
                        if str(atom_count) == attach:
                            sub_label = str(int(attach)+5)
                            label = sub_label+(3-len(sub_label))*" "
                            new_line = line[:35]+label+line[38:]
                    new_content.append(new_line)
                else:
                    # atom_line = False
                    new_content += content[file_count:]
                    break
            else:
                new_content.append(line)
            if line.find("999 V2000") != -1:
                atom_line = True
            file_count += 1
        return copy.deepcopy(new_content)

    def label(self, content):

        attach_list = []
        for l in self.attachments:
            attach_list += l
        R = self.Rpositions
        print(f"Rpositions: {R}")
        L = list(self.list_of_atoms.keys())
        atoms = R + L + list(self.attach_ids.keys()) + attach_list
        # all atoms that need labelling now form a set "all_atoms"
        all_atoms = []
        for atom in atoms:
            if atom not in all_atoms:
                all_atoms.append(int(atom))
        # all_atoms is int list
        all_atoms.sort()
        print(f"all_atoms: {all_atoms}")
        i = 0
        new_all_atoms = copy.deepcopy(all_atoms)
        # new_all_atoms contains all atoms other than the empty and real atoms
        print(f"no_atoms: {self.no_atoms}")
        for at in all_atoms:
            if int(at) > self.no_atoms:
                new_all_atoms= all_atoms[:i]
                break
            i += 1
        all_atoms = new_all_atoms
        print(f"all_atoms: {all_atoms}")
        n = len(all_atoms)
        p = str(n)
        iso_line = f"M  ISO"
        iso_line += (3-len(p))*" "+p
        for p in all_atoms:
            # label
            a = self.atom_symbols[p]
            table = Chem.GetPeriodicTable()
            atomic_mass = int(table.GetMostCommonIsotopeMass(a))
            label = str(atomic_mass + 5 + p)
            iso_line += (4-len(str(p)))*" "+str(p)+(4-len(label))*" "+label
        iso_line += "\n"
        new_content = []
        replace = True
        for line in content:
            if line.find("M  END") != -1 and replace:
                new_content.append(iso_line)
                replace = False
            new_content.append(line)
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

        # returns inchi for subs with more than 1 main atom and smiles otherwise
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

        # R groups
        zz = zz_convert()
        core_mol = self.core_mol
        core_inchi = Chem.MolToInchi(core_mol)
        print(f"core_inchi: {core_inchi}")
        order = []
        replace_order = {}
        index = core_inchi.find("/i")
        label_part = ""
        if index != -1:
            indexf = core_inchi[index+2:].find("/")
            if indexf != -1:
                label_part = core_inchi[index+2:index+2+indexf]
            else:
                label_part = core_inchi[index+2:]
        print(f"label_part: {label_part}")
        canonical_dict = {} # mol_label:inchi_label
        for part in label_part.split(","):
            canonical_dict[part.split("+")[1]] = part.split("+")[0]
            rank = part.split("+")[0]
            formula = core_inchi.split("/")[1]
            symbol = self.help_label.find_atom(rank, formula)
            if symbol == "Te":
                order.append(str(int(part.split("+")[1])-6))
            else:
                if str(int(part.split("+")[1])-5) in self.list_of_atoms.keys():
                    replace_order[(str(int(part.split("+")[1])-5))] = str(int(part.split("+")[0]))
        mark_inchi = ""
        if core_inchi.find("Te") != -1:
            mark_inchi = self.relabel_core(zz.te_to_zz(core_inchi))
            print(f"core_inchi_2: {core_inchi}")
        else:
            mark_inchi = self.relabel_core(core_inchi)
            print(f"core_inchi_2: {core_inchi}")
        Rsubstituents = self.Rsubstituents
        for num in order:
            ind = self.Rpositions.index(num)
            subs = Rsubstituents[ind]
            sub_inchis = []
            mark_inchi += "<M>"
            for sub in subs:
                sub_inchi = self.relabel_sub(sub)
                if sub_inchi.find("Te") != -1:
                    sub_inchi = zz.te_to_zz("InChI=1B/"+sub_inchi)
                    sub_inchi = "/".join(sub_inchi.split("/")[1:])
                sub_inchis.append(sub_inchi)
            sub_inchis.sort()
            for sub_inchi1 in sub_inchis:
                mark_inchi += sub_inchi1+"!"
            mark_inchi = mark_inchi[:-1]
        mark_inchi = mark_inchi.replace("InChI=1S/", "MarkInChI=1B/")
        mark_inchi = mark_inchi.replace("[HH]", "H")
        #### list of atoms
        part = ""
        print(f"replace_order: {replace_order}")
        for ind in replace_order.keys():
            print(self.list_of_atoms)
            atom_list = self.list_of_atoms[ind]
            atom_list.sort()
            part += "<M>"
            lab = replace_order[ind]
            part += lab+"-"
            for atom in atom_list:
                part += atom +"!"
            part = part[:-1]
        print(f"list of atoms part: {part}")
        mark_inchi += part
        #### variable attachments
        print(f"mark_inchi: {mark_inchi}")
        var_part = ""
        atom_ids = self.attach_ids
        for i in range(0, len(list(atom_ids.keys()))):
            mol_rank = list(atom_ids.keys())[i]
            symbol = atom_ids[mol_rank]
            subs = []
            sub_inchis = []
            var_part += "<M>"
            for mi in self.attachments[i]:
                print(f"self.attachments: {self.attachments}")
                other_symbol = self.atom_symbols[int(mi)]
                print(f"other_symbol: {other_symbol}")
                no = 0
                if other_symbol == "Te":
                    no = 6
                else:
                    no = 5
                new_attach = canonical_dict[str(int(mi)+no)]
                var_part += new_attach+"H"+","
            var_part = var_part[:-1]
            var_part += "-"
            if symbol == "Te":
                ind = self.Rpositions.index(str(mol_rank))
                subs = self.Rsubstituents[ind]
                print(f"subs: {subs}")
                sub_inchi = ""
                for sub in subs:
                    sub_inchi = self.relabel_sub(sub)
                    print(f"sub_inchi: {sub_inchi}")
                    if sub_inchi.find("Te") != -1:
                        sub_inchi = zz.te_to_zz("InChI=1B/"+sub_inchi)
                        sub_inchi = "/".join(sub_inchi.split("/")[1:])
                    sub_inchis.append(sub_inchi)
                sub_inchis.sort()
                for sub_inchi1 in sub_inchis:
                    var_part += sub_inchi1+"!"
                var_part = var_part[:-1]
            else:
                # check if it is a list of atoms
                if str(mol_rank) in self.list_of_atoms.keys():
                    atoms = self.list_of_atoms[str(mol_rank)]
                    atoms.sort()
                    for atom in atoms:
                        var_part += atom+"!"
                    var_part = var_part[:-1]
                else:
                    var_part += symbol
                # check if it is a normal atom substituent
        mark_inchi += var_part
        return mark_inchi

if __name__=="__main__":
    name = input("Enter the file name with the extension:")
    file = open(name, "r")
    content = file.readlines()
    mark_obj = markmol()
    zz = zz_convert()
    mark_obj.help_label = Label()
    mark_obj.atom_symbols = {}
    mark_obj.Rpositions = []  # number of each Te atom
    mark_obj.Rsubstituents = []
    mark_obj.ctabs = []  # no. of each substituent after $RGP
    mark_obj.connections = []
    mark_obj.list_of_atoms = {}
    mark_obj.attachments = []
    mark_obj.attach_ids = {}
    mark_obj.no_atoms = 0
    content = mark_obj.convert(content)
    new_name = name.split(".")[0]+"_RDKIT.sdf"
    new_file = open(new_name, "w")
    new_file.writelines(content)
    new_file.close()
    supply = Chem.SDMolSupplier(new_name)
    substituents = []
    for mol in supply:
        if mol is not None:
            substituents.append(mol)
    mark_obj.core_mol = copy.deepcopy(supply[0])
    i = 1
    ctabs = mark_obj.ctabs
    while i < len(ctabs):
        mark_obj.Rsubstituents.append(substituents[int(ctabs[i-1]):int(ctabs[i])])
        i += 1
    mark_obj.Rsubstituents.append(substituents[int(ctabs[i-1]):])
    print(mark_obj.produce_markinchi())
    file.close()
