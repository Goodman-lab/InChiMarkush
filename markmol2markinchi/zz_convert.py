import sys, os
from rdkit import Chem

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'new_markinchi'))
from helper import helper

class zz_convert(object):

    def __init__(self):
        self.helper = helper()

    def zz_to_te(self, inchi):

        # TODO: replace Zz with Te
        new_inchi = inchi.replace("Zz", "Te")

        # TODO: determine ranking of zz atoms
        no, start_zz = self.helper.zz_no(new_inchi)

        # TODO: add one hydrogen to each Te using the hydrogen layer
        index = new_inchi.find("/h")
        h_layer = ""
        if index != -1:
            # find "H," or "H/" and see the number before
            h_layer = new_inchi.split("/h")[1].split("/")[0]
            if h_layer.find("H,") != -1:
                in_index = h_layer.find("H,")
                _index = h_layer[:in_index].rfind("-")
                com_index = h_layer[:in_index].rfind(",")
                number = ""
                if _index != -1:
                    if com_index > _index:
                        number = h_layer[com_index+1:in_index]
                    else:
                        number = h_layer[_index+1:in_index]
                else:
                    if com_index != -1:
                        number = h_layer[com_index+1:in_index]
                    else:
                        number = h_layer[:in_index]

                h_layer = self.h_operations(number, start_zz, no, in_index,
                                            _index, h_layer)
            else:
                if h_layer[-1] == "H":
                    in_index = h_layer.find("H")
                    _index = h_layer[:in_index].rfind("-")
                    com_index = h_layer[:in_index].rfind(",")
                    number = ""
                    if _index != -1:
                        if com_index > _index:
                            number = h_layer[com_index+1:in_index]
                        else:
                            number = h_layer[_index+1:in_index]
                    else:
                        if com_index != -1:
                            number = h_layer[com_index+1:in_index]
                        else:
                            number = h_layer[:in_index]

                    h_layer = self.h_operations(number, start_zz, no, in_index,
                                                _index, h_layer)
                else:
                    if no > 1:
                        h_layer = str(start_zz)+"-"+str(start_zz+no-1)+","+h_layer
                    else:
                        h_layer = str(start_zz)+","+h_layer

            if new_inchi.split("/h")[1].find("/") != -1:
                suf = "/"+"/".join(new_inchi.split("/h")[1].split("/")[1:])
                new_inchi = new_inchi[:index+2]+h_layer+suf
            else:
                new_inchi = new_inchi[:index+2]+h_layer
        else:
            if no > 1:
                h_layer = "/h"+str(start_zz)+"-"+str(start_zz+no-1)+"H"
            else:
                h_layer = "/h"+str(start_zz)+"H"
            if new_inchi.split("/c")[1].find("/") != -1:
                index1 = new_inchi.find("/c")
                indexf = new_inchi[index1+1:].find("/")
                indexf += index1+1
                new_inchi = new_inchi[:indexf]+h_layer+new_inchi[indexf:]
            else:
                new_inchi += h_layer
        formula = new_inchi.split("/")[1]
        index0 = formula.find("H")
        num = no
        if index0 != -1:
            indexf = 1
            for char in formula[index0+1:]:
                if char.isalpha():
                    break
                indexf += 1
            if indexf == 1:
                num += 1
            else:
                num += int(formula[index0+1:index0+indexf])
            formula = formula[:index0+1]+str(num)+formula[index0+indexf:]
            new_inchi = new_inchi.split("/")[0]+"/"+formula+"/"+"/".join(new_inchi.split("/")[2:])
        else:
            indexc = formula.find("C")
            indexf = indexc+1
            for char in formula[indxc+1:]:
                if char.is_alpha():
                    break
                indexf += 1
            sub = "H"+str(no)
            formula = formula[:indexf]+sub+formula[indexf:]
        return new_inchi


    def h_operations(self, number, start_zz, no, in_index, _index, h_layer):

        if int(number) == start_zz - 1:
            # case: highest number that has one hydrogen is right before
            # the lowest Zz atom number
            com_index = h_layer.rfind(",")
            if _index != -1:
                if _index > com_index:
                    h_layer = h_layer[:_index]+str(start_zz+no-1)+h_layer[in_index:]
                else:
                    suf = number+"-"+str(start_zz+no-1)+h_layer[in_index:]
                    h_layer = h_layer[:com_index+1]+suf
            else:
                if com_index != -1:
                    suf = number+"-"+str(start_zz+no-1)+h_layer[in_index:]
                    h_layer = h_layer[:com_index+1]+suf
                else:
                    h_layer = number+"-"+str(start_zz+no-1)+h_layer[in_index:]
        else:
            if no > 1:
                # More than 1 Zz (use "-")
                sub = ","+str(start_zz)+"-"+str(start_zz+no-1)
                h_layer = h_layer[:in_index]+sub+h_layer[in_index:]
            else:
                # only one Zz so no need for ("-")
                h_layer = h_layer[:in_index]+","+str(start_zz)+h_layer[in_index:]
        return h_layer

    def te_to_zz(self, inchi):

        # Bug: need to remove hydrogens from the formula
        # TODO: replace Te with Zz
        new_inchi = inchi.replace("Te", "Zz")

        # TODO: determine ranking of zz atoms
        no, start_zz = self.helper.zz_no(inchi)

        # TODO: remove one hydrogen from each Te using the hydrogen layer
        index = new_inchi.find("/h")
        h_layer = ""
        if index != -1:
            h_layer = new_inchi.split("/h")[1].split("/")[0]
            if h_layer.find("H,") != -1:
                in_index = h_layer.find("H,")
                _index = h_layer[:in_index].rfind("-")
                com_index = h_layer[:in_index].rfind(",")
                number = ""
                number0 = ""
                if _index != -1:
                    if com_index > _index:
                        number = h_layer[com_index+1:in_index]
                        if int(number) == start_zz:
                            h_layer = h_layer[:com_index]+h_layer[in_index:]
                    else:
                        number = h_layer[_index+1:in_index]
                        if com_index != -1:
                            number0 = h_layer[com_index+1:_index]
                        else:
                            number0 = h_layer[:_index]
                        if int(number) > start_zz:
                            # TODO: find the number before "-" to do comparisons
                            if int(number0) < start_zz-1:
                                h_layer = h_layer[:_index+1]+str(start_zz-1)+h_layer[in_index:]
                            elif (int(number0) == start_zz-1):
                                h_layer = h_layer[:_index]+h_layer[in_index:]
                        elif (int(number)== start_zz):
                            if int(number0) < start_zz-1:
                                h_layer = h_layer[:_index+1]+str(start_zz-1)+h_layer[in_index:]
                            elif (int(number0) == start_zz):
                                if com_index != -1:
                                    h_layer = h_layer[:com_index]+h_layer[in_index:]
                                else:
                                    h_layer = h_layer[in_index:]
                            elif (int(number0) == start_zz-1):
                                h_layer = h_layer[:_index]+h_layer[in_index:]
                else:
                    if com_index != -1:
                        number = h_layer[com_index+1:in_index]
                        if int(number) == start_zz:
                            h_layer = h_layer[:com_index]+h_layer[in_index:]
                    else:
                        number = h_layer[:in_index]
                        if int(number) == start_zz:
                            h_layer = h_layer[in_index+2:]
            else:
                if h_layer[-1] == "H":
                    in_index = h_layer.find("H")
                    _index = h_layer[:in_index].rfind("-")
                    com_index = h_layer[:in_index].rfind(",")
                    number = ""
                    number0 = ""
                    if _index != -1:
                        if com_index > _index:
                            number = h_layer[com_index+1:in_index]
                            if int(number) == start_zz:
                                h_layer = h_layer[:com_index]+h_layer[in_index:]
                        else:
                            number = h_layer[_index+1:in_index]
                            if com_index != -1:
                                number0 = h_layer[com_index+1:_index]
                            else:
                                number0 = h_layer[:_index]
                            if int(number) > start_zz:
                                # TODO: find the number before "-" to do comparisons
                                if int(number0) < start_zz-1:
                                    h_layer = h_layer[:_index+1]+str(start_zz-1)+h_layer[in_index:]
                                elif (int(number0) == start_zz):
                                    if com_index != -1:
                                        h_layer = h_layer[:com_index]+h_layer[in_index:]
                                    else:
                                        h_layer = h_layer[in_index:]
                                elif (int(number0) == start_zz-1):
                                    h_layer = h_layer[:_index]+h_layer[in_index:]
                            elif (int(number)== start_zz):
                                if int(number0) < start_zz-1:
                                    h_layer = h_layer[:_index+1]+str(start_zz-1)+h_layer[in_index:]
                                elif (int(number0) == start_zz-1):
                                    h_layer = h_layer[:_index]+h_layer[in_index:]
                    else:
                        if com_index != -1:
                            number = h_layer[com_index+1:in_index]
                            if int(number) == start_zz:
                                h_layer = h_layer[:com_index]+h_layer[in_index:]
                        else:
                            number = h_layer[:in_index]
                            if int(number) == start_zz:
                                h_layer = ""
            if new_inchi.split("/h")[1].find("/") != -1:
                suf = "/"+"/".join(new_inchi.split("/h")[1].split("/")[1:])
                new_inchi = new_inchi[:index+2]+h_layer+suf
            else:
                new_inchi = new_inchi[:index+2]+h_layer
        formula = new_inchi.split("/")[1]
        index0 = formula.find("H")
        num = no*-1
        if index0 != -1:
            indexf = 1
            for char in formula[index0+1:]:
                if char.isalpha():
                    break
                indexf += 1
            if indexf == 1:
                num += 1
            else:
                num += int(formula[index0+1:index0+indexf])
            formula = formula[:index0+1]+str(num)+formula[index0+indexf:]
            new_inchi = new_inchi.split("/")[0]+"/"+formula+"/"+"/".join(new_inchi.split("/")[2:])
        else:
            indexc = formula.find("C")
            indexf = indexc+1
            for char in formula[indxc+1:]:
                if char.is_alpha():
                    break
                indexf += 1
            sub = "H"+str(no)
            formula = formula[:indexf]+sub+formula[indexf:]
        return new_inchi




if __name__ == "__main__":
    inchi = input("Please enter the inchi with Zz:")
    zz = zz_convert()
    new_inchi = zz.zz_to_te(inchi)
    print(new_inchi)
