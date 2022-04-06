from markinchi import MarkInChI
class compare_markinchi(object):
    def __init__(self, inchi1, inchi2):
        markinchi1 = MarkInChI(inchi1).list_of_inchi
        markinchi2 = MarkInChI(inchi2).list_of_inchi
        set_markinchi1 = []
        set_markinchi2 = []
        for sub_inchi in markinchi1:
            if sub_inchi not in set_markinchi1:
                set_markinchi1.append(sub_inchi)
        for sub_inchi in markinchi2:
            if sub_inchi not in set_markinchi2:
                set_markinchi2.append(sub_inchi)
        intersection = []
        for inchi_1 in set_markinchi1:
            for inchi_2 in set_markinchi2:
                if inchi_1 == inchi_2:
                    intersection.append(inchi_1)
                    break
        print(f"intersection: {intersection}")

if __name__ == "__main__":
    inchi1 = input("Please enter the first markinchi: ")
    inchi2 = input("Please enter the second markinchi: ")
    Comp = compare_markinchi(inchi1, inchi2)
