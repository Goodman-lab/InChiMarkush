from rdkit import Chem

class helper(object):

  def zz_no(self, inchi):

    #Get number of Zz atoms
    index = inchi.find("Te")+2
    char = inchi[index]
    no = 1
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

    return no, start_zz
