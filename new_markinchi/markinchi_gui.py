import os, sys
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import ANCHOR, PhotoImage, NW, Scrollbar, END
from PIL import ImageTk,Image
import pygubu
from tkinter.filedialog import askopenfile  # convert mol to markinchi
import copy
from rdkit import Chem
from rdkit.Chem import Draw
from markinchi import MarkInChI
from label import Label
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'markmol2markinchi'))
from zz_convert import zz_convert
from markmol import markmol

PROJECT_PATH = os.path.abspath(os.path.dirname(__file__))
PROJECT_UI = os.path.join(PROJECT_PATH, "markinchi_gui.ui")

class MarkinchiGuiApp(object):
    def __init__(self, master=None):

        # Get the pygubu builder that will build tkinter widgets from
        # XML file
        self.builder = builder = pygubu.Builder()
        builder.add_resource_path(PROJECT_PATH)
        builder.add_from_file(PROJECT_UI)
        self.mainwindow = builder.get_object('frame1', master)
        self.master = master
        master.resizable(False, False)  # Make the program window not resizable
        master.title("MarkInChI to List of InChIs")  # Set the title of the program
        # build menu
        my_menu = builder.get_object('menu1')
        master.config(menu=my_menu)
        # build entry
        self.entry = builder.get_object('entry1')
        # Build the canvases from the XML file.
        self.my_canvas = builder.get_object('canvas1')
        self.my_canvas.configure(bg='white')
        self.mol_canvas = builder.get_object('canvas2')
        self.mol_canvas.grid_forget()


        # Make window scrollable
        # Build scrollabrs from XML
        my_scrollbar = builder.get_object('scrollbar1')
        h_scrollbar = builder.get_object('scrollbar3')

        """ Put the scrollbar on the y axis of the
            canvas and activate it in the canvax"""

        # Create ANOTHER Frame INSIDE the Canvas
        self.second_frame = builder.get_object('frame2')
        # Configure and set the scrollbars on the canvas
        my_scrollbar.configure(command=self.my_canvas.yview)
        h_scrollbar.configure(command=self.my_canvas.xview)
        self.my_canvas.configure(xscrollcommand=h_scrollbar.set,
                                 yscrollcommand=my_scrollbar.set)
        # Bind the canvas to the second frame
        sx = self.my_canvas
        self.second_frame.bind('<Configure>', lambda e:
                               sx.configure(scrollregion=sx.bbox("all")))

        # Add that New frame To a Window In The Canvas and set its position
        self.my_canvas.create_window((100, 0), window=self.second_frame,
                                     anchor="nw")
        # Bind to the mousewheel when cursor is on an object
        self.second_frame.bind('<Enter>', self._bound_to_mousewheel)
        self.second_frame.bind('<Leave>', self._unbound_to_mousewheel)
        builder.connect_callbacks(self)  # bind layout objects to functions

    def _unbound_to_mousewheel(self, event):

        # Unbinds object when cursor is not on it
        self.my_canvas.unbind_all("<MouseWheel>")

    def _bound_to_mousewheel(self, event):

        # Binds object when cursor is on it
        self.my_canvas.bind_all("<MouseWheel>", self._on_mousewheel)

    def _on_mousewheel(self, event):

        # Binds mousewheel to the vertical scrollbar
        self.my_canvas.yview_scroll(int(-1*(event.delta/120)), "units")

    def convert(self):

        """ Command of the convert button. It converts markinchi to a list
            of single inchis and then displays them in the listbox"""

        entry = self.entry
        self.listbox = self.builder.get_object('listbox1')
        self.listbox.delete(0, END)
        markinchi = entry.get()
        list_of_inchi = self.get_inchis(markinchi)
        i = 1
        print(f"Number of inchi produced: {len(list_of_inchi)}")
        for inchi in list_of_inchi:
            self.listbox.insert(i, inchi)
            i += 1


    def get_inchis(self, inchi):

        """ This function converts a markinchi to a list of inchis"""

        inchi_obj = MarkInChI()
        zz = zz_convert()
        inchi = inchi.replace("MarkInChI", "InChI", 1)
        #inchi = inchi.replace("Zz", "Te")
        # Check it is actually a markinchi
        if inchi.find("<>") == -1:
            # TODO: flag the user
            print("Not a MarkInChI")
        else:
            inchiplus = inchi.split("<>")
            # Get main inchi and substituents
            inchiplus_item = inchiplus.pop(0)
            if inchiplus_item.find("Zz") != -1:
                inchiplus_item = zz.zz_to_te(inchiplus_item)

            for sub in range(0, len(inchiplus)):
                inchiplus[sub] = inchiplus[sub].replace("Zz", "Te")
            inchi_obj.inchi = inchiplus_item
            inchi_obj.list_of_inchi = []
            inchi_obj.no_run = 0
            inchi_obj.label = label = Label()
            inchiplus_item, inchi_obj.ranks = label.label_inchi(inchiplus_item,
                                                                 inchiplus)
            # Convert main inchi to mol and run algorithm
            main_mol = Chem.rdinchi.InchiToMol(inchiplus_item)[0]
            # Sanitize (only in rdkit)
            new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(main_mol))
            inchi_obj.run(inchiplus, new_mol)
        return copy.deepcopy(inchi_obj.list_of_inchi)

    def plot(self):

        """ This function is a command of the plot button. It plots
            the single inchi selected in the listbox in a canvas
            below the button """

        inchi = self.listbox.get(ANCHOR)
        mol = Chem.rdinchi.InchiToMol(inchi)[0]
        Draw.MolToFile(mol, "mol.png")
        self.img = ImageTk.PhotoImage(Image.open("mol.png"))
        self.mol_canvas.create_image(20,20, anchor=NW, image=self.img)
        self.mol_canvas.image = self.img
        self.mol_canvas.grid(row=5, column=0, columnspan=2)

    def open(self):

        self.entry.delete(1, 1000)
        path = os.getcwd()
        file = askopenfile(mode='r', initialdir=path, defaultextension='.sdf',
                           filetypes=[("SDF file", ".sdf")])
        if file is not None:
            name = os.path.basename(file.name)
            content = file.readlines()
            mark_obj = markmol()
            zz = zz_convert()
            mark_obj.Rpositions = []
            mark_obj.Rsubstituents = []
            mark_obj.ctabs = []
            mark_obj.connections = []
            content = mark_obj.convert(content)
            new_name = name.split(".")[0]+"_RDKIT.sdf"
            new_file = open(new_name, "w")
            new_file.writelines(content)
            new_file.close()
            supply = Chem.SDMolSupplier(new_name)
            substituents = []
            for mol in supply:
                if mol is not None:
                    substituents.append(mol)
            mark_obj.core_mol = copy.deepcopy(supply[0])
            i = 1
            ctabs = mark_obj.ctabs
            while i < len(mark_obj.ctabs):
                mark_obj.Rsubstituents.append(substituents[int(ctabs[i-1]):int(ctabs[i])])
                i += 1
            mark_obj.Rsubstituents.append(substituents[int(ctabs[i-1]):])
            print(mark_obj.produce_markinchi())
            self.entry.insert(0, mark_obj.produce_markinchi())
            file.close()

    def quit(self):

        self.master.quit()

    def run(self):
        self.mainwindow.mainloop()  # show window


if __name__ == '__main__':
    # start the app
    root = tk.Tk()
    app = MarkinchiGuiApp(root)
    app.run()
