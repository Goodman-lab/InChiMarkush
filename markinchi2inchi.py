#################################################################
#
# source activate my-rdkit-env
#
################################################################

from rdkit import Chem


def run(inchiplus, grouplist):
    inchiplus=sys.argv[1].split("<>")
    grouplist=inchiplus[1].split("!")

    molfrominchi = Chem.rdinchi.InchiToMol(inchiplus[0])
    smiles_core=Chem.MolToSmiles(molfrominchi[0])

    for substituent in grouplist:
        if substituent == "H":
            substituent = ""
        new_smiles=smiles_core.replace("[Te]", substituent)

        new_inchi=Chem.MolToInchi(Chem.MolFromSmiles(new_smiles))
        print(new_inchi)



if __name__=="__main__":
    from sys import argv

    if argv[1].find("Te") == -1:
        print("No tellurium")
    elif argv[1].find("<>") == -1:
        print("No substituents")
    else:
        inchiplus = argv[1].split("<>")
        grouplist = inchiplus[1].split("!")
        run(inchiplus, grouplist)
