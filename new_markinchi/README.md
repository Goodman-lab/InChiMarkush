# MarkInChI to List of InChI

This folder contains both a python
script and a python gui to convert MarkInChI
to a list of inchis and view them. This version
covers all 5 examples described in the
MarkInChI and VInChI Working Group meeting
on 17th of May 2021 (see in resources). The whole
code is documented and you can make use of some
presentations and .pdf documents provided in
the resources section below.

## Current Limitations

This version was only tested on the cases covered by the 5 examples in the 17th May meeting document (see resource 1). Thus, there is a number of limitations including but not restricted to:

1- Stereochemistry: this version of code doesn't treat the stereochemistry of newly formed centres and leaves them undefined.

2- Nested variable groups: this version doesn't cover a variable group inside another variable group.

3- Number of atoms in MarkInChI Core: as isotopic labelling is used to introduce a new canonical labelling system to the MarkInChI core, there is a limit to the rank of the atom being replaced that it doesn't exceed 116.

## Resources

1- [17th May 2021 Working Group Meeting](https://drive.google.com/file/d/14VPgQNHCs5_X2_AXaWxGXCOD12z79k_z/view?usp=sharing)

2- [Presentation 1 - Abdullah Kattineh](https://docs.google.com/presentation/d/1F631duOvL39CWStcOzdq5riRriocgaHC/edit?usp=sharing&ouid=105229675019634902210&rtpof=true&sd=true)

3- [Presentation 2 - Abdullah Kattineh](https://docs.google.com/presentation/d/1Ml4-kUDWRmLlXuE9YN66a9nvKDNoP7H2/edit?usp=sharing&ouid=105229675019634902210&rtpof=true&sd=true)
