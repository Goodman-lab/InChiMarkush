# MarkInChI GUI

This folder contains a python
script, markinchi.py, to convert a MarkInChI
to a list of inchis and plot them. It also has another python script, markmol.py, that can produce markinchi strings from markush 2000V sdf files. Moreover, there is a python gui script, markinchi_gui.py, and a .exe application, markinchi_gui.exe, that can be built using the batch file, PI.bat, in a virtual environment that contains the rdkit package. This version
covers all 5 examples described in the
MarkInChI and VInChI Working Group meeting
on 17th of May 2021 (see in resources) and which you can also find in the 'Examples' folder. The whole
code is documented and you can make use of some
presentations and .pdf documents provided in
the resources section below.

## Current Limitations

This version was only tested on the cases covered by the 5 examples in the 17th May meeting document (see resource 1). Thus, there is a number of limitations including but not restricted to:

1- Stereochemistry: this version of code doesn't treat the stereochemistry of newly formed stereocentres and leaves them undefined.

2- Terminal R groups: this version only covers terminal R groups and doesn't cover R groups with more than one connection. This means that ring variations, nested R groups, and any representation that requires more than one bond to the R group can't be represented.

3- Number of atoms in MarkInChI Core: as isotopic labelling is used to introduce a new canonical labelling system to the MarkInChI core, there is a limit to the inchi canonical rank of the atom being replaced that it doesn't exceed 116.

4- Markmol doesn't support isotopically labelled molecules: this is the part of the application that converts a markush 2000v .sdf file to a markinchi.

## Changes to the documentation described in the 17th May 2021 working group meeting file

Some changes have been applied to make MarkInChI more compatible and efficient:

1- Use of <M> instead of <>: this was applied to make markinchi compatible with RInChI.

2- Not including charges in lists of atoms: e.g. in example 5 '8-C!N+', now the programme can identify charges automatically and there is no need to include them i.e. '8-C!N+' becomes '8-C!N'. This is a more useful representation as charges cannot be included in lists of atoms in MarvinSketch.

## Final Notes on the EXE Application

The application can be built using the PyInstaller package and the software makes extensive use of the Python RDKIT package. MarvinSketch, a product of ChemAxon, is used to produce the markush sdf files in the examples. A copy of the rdkit license can be found in the main folder.

## Resources

1- [17th May 2021 Working Group Meeting](https://drive.google.com/file/d/14VPgQNHCs5_X2_AXaWxGXCOD12z79k_z/view?usp=sharing)

2- [Presentation 1 - Abdullah Kattineh](https://docs.google.com/presentation/d/1F631duOvL39CWStcOzdq5riRriocgaHC/edit?usp=sharing&ouid=105229675019634902210&rtpof=true&sd=true)

3- [Presentation 2 - Abdullah Kattineh](https://docs.google.com/presentation/d/1Ml4-kUDWRmLlXuE9YN66a9nvKDNoP7H2/edit?usp=sharing&ouid=105229675019634902210&rtpof=true&sd=true)

4- [Presentation 3 - Abdullah Kattineh](https://docs.google.com/presentation/d/1Q45tEMzSBWpYQxvWDJojaQH_Q-P0I0BW/edit?usp=sharing&ouid=104156901695306507330&rtpof=true&sd=true)

5- [Presentation 4 - Abdullah Kattineh](https://docs.google.com/presentation/d/17Cg3LWVf9LIutk8ik0HVcpAGzmqwZ3xk/edit?usp=sharing&ouid=104156901695306507330&rtpof=true&sd=true)
