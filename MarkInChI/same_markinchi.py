# switched off function for comparing numbers of non-H atoms
# useful only if some atoms are replaced and MarkInChI is generated from different structures
# should not happen because MarkInChI should be canonical

import sys, os
from markinchi import MarkInChI

class CompareMarkInChI(object):

    def __init__(self, inchi1, inchi2):

        # First check if they are MarkInChIs
        if ("<M>" not in inchi1 or "MarkInChI" not in inchi1) or ("<M>" not in inchi2 or "MarkInChI" not in inchi2):
            print("One or both inputs not MarkInChI")
            sys.exit()

        # Compare lengths of the MarkInChIs
        if len(inchi1) == len(inchi2):

            # Compare the MarkInChI string
            if inchi1 == inchi2:
                print("Equivalent (exactly the same)")
                sys.exit()

            self.get_parts(inchi1, inchi2)
            subgroup = 0
            # Compare if substituents and grouplists are the same
            # (if same core InChI, but more substituents, MarkInChIs are not equivalent)
            if self.inchi_main1 == self.inchi_main2:
                if len(self.grouplist1) == len(self.grouplist2):
                    while subgroup < len(self.grouplist1_sorted):
                        if len(self.grouplist1_sorted[subgroup]) != len(self.grouplist2_sorted[subgroup]):
                            print("Not equivalent (same core, different number of substituents)")
                            sys.exit()

                        subgroup = subgroup + 1

                    if self.grouplist1_sorted == self.grouplist2_sorted:
                        print("Equivalent (same core, same substituents)")
                    else:
                        print("Not equivalent (same core, different substituents)")
                sys.exit()
        else:
            # Check if the molecular formulas differ if they don't contain Zz
            molformula1 = inchi1.split("/")[1]
            molformula2 = inchi2.split("/")[1]
            molformulas = molformula1 + molformula2

            if molformula1 != molformula2:
                if "Zz" not in molformulas:
                    print("Not equivalent (different molecular formulas)")
                    sys.exit()

            self.get_parts(inchi1, inchi2)
            # Think this through if it really cannot be truth
            if self.inchi_main1 == self.inchi_main2:
                print("Not equivalent (different length, same core)")
                sys.exit()

        # Suppress printing
        sys.stdout = open(os.devnull, 'w')

        markinchi1 = MarkInChI(inchi1).list_of_inchi
        markinchi2 = MarkInChI(inchi2).list_of_inchi

        #Enable printing
        sys.stdout = sys.__stdout__

        markinchisorted1 = sorted(markinchi1)
        markinchisorted2 = sorted(markinchi2)

        # Comparing MarkInChIs using list of InChIs if none of the above criteria was able to determine it
        if markinchisorted1 == markinchisorted2:
            print("Equivalent via list of InChIs")
        else:
            print("Not equivalent via list of InChIs")


    def get_parts(self, inchi1, inchi2):
        # Obtaining core InChIs
        inchi_parts1 = inchi1.split("<M>")
        inchi_parts2 = inchi2.split("<M>")

        self.inchi_main1 = inchi_parts1.pop(0)
        self.inchi_main2 = inchi_parts2.pop(0)

        # Compare substituents (make grouplist consisting of subgrouplists)

        self.grouplist1 = []
        self.grouplist2 = []
        sub1 = 0
        sub2 = 0
        while sub1 < len(inchi_parts1):
            subgrouplist1 = sorted(inchi_parts1[sub1].split("!"))
            self.grouplist1.append(subgrouplist1)
            sub1 = sub1 + 1

        while sub2 < len(inchi_parts2):
            subgrouplist2 = sorted(inchi_parts2[sub2].split("!"))
            self.grouplist2.append(subgrouplist2)
            sub2 = sub2 + 1

        self.grouplist1_sorted = sorted(self.grouplist1)
        self.grouplist2_sorted = sorted(self.grouplist2)

        return

if __name__ == "__main__":
    inchi1 = input("Please enter the first MarkInChI: ")
    inchi2 = input("Please enter the second MarkInChI: ")
    Comp = CompareMarkInChI(inchi1, inchi2)