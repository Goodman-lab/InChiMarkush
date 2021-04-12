#################################################################
#
# source activate my-rdkit-env
#
################################################################

from rdkit import Chem


def run(inchiplus, smiles_core):
  #print("=============")
  #print(inchiplus)
  #print(smiles_core)
  inchiplus_item = inchiplus.pop(0)
  grouplist=inchiplus_item.split("!")
  for substituent in grouplist:
    if substituent == "H":
      substituent = ""
    new_smiles=smiles_core.replace("[TeH]", substituent, 1)
    if new_smiles.find("()") >= 0:
      new_smiles=smiles_core.replace("[TeH]", substituent, 1).replace("()","")
    if len(inchiplus) == 0:
      new_inchi=Chem.MolToInchi(Chem.MolFromSmiles(new_smiles))
      #print(new_smiles)
      print(new_inchi)
    else:
      #print(inchiplus)
      #print(new_smiles)
      new_inchiplus = inchiplus.copy()
      run(new_inchiplus, new_smiles)
      #print("back from recursion")
      #print(substituent,new_smiles,inchiplus)
  return


#if True:# __name__=="__main__":
if __name__=="__main__":
    from sys import argv

    if argv[1].find("Te") == -1:
        print("No tellurium")
    elif argv[1].find("<>") == -1:
        print("No substituents")
    else:
        inchiplus = argv[1].split("<>")
        #grouplist = inchiplus[1].split("!")
        #run(inchiplus, grouplist)
        inchiplus_item = inchiplus.pop(0)
        molfrominchi = Chem.rdinchi.InchiToMol(inchiplus_item)
        smiles_core=Chem.MolToSmiles(molfrominchi[0])
        #print(inchiplus_item)
        #print(inchiplus)
        #print(smiles_core)
        run(inchiplus, smiles_core)

