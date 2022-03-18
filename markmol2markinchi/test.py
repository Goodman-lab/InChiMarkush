from markmol import markmol
mark_obj = markmol()
file = open("ex5-2000.sdf", 'r')
content = file.readlines()
mark_obj.var_attach(content)
