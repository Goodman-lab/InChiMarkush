#################################################################
#
# source activate my-rdkit-env
#
import sys
from rdkit import Chem
#################################################################
#print(sys.argv)


#m = Chem.MolFromMolFile(sys.argv[1])
#mu = Chem.MolFromMolBlock(sys.argv[1], sanitize = False)
#print("read molfile: "+sys.argv[1])
#inchi = Chem.inchi.MolToInchi(m)
#print(inchi)
#print(Chem.MolToSmiles(m))

if sys.argv[1].find("Te") == -1:
  print("No tellurium")
  exit()
if sys.argv[1].find("<>") == -1:
  print("No substituents")
  exit()


inchiplus=sys.argv[1].split("<>")
grouplist=inchiplus[1].split("!")
#print(molfrominchi)

#print(inchiplus)
#print(grouplist)



molfrominchi = Chem.rdinchi.InchiToMol(inchiplus[0])
smiles_core=Chem.MolToSmiles(molfrominchi[0])
#print(smiles_core)

for substituent in grouplist:
  #print(substituent)
  if substituent == "H":
    substituent = ""
  new_smiles=smiles_core.replace("[Te]", substituent)
  #print(new_smiles)
  new_inchi=Chem.MolToInchi(Chem.MolFromSmiles(new_smiles))
  print(new_inchi)

