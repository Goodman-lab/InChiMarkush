import os, sys, traceback
import tkinter as tk
import tkinter.ttk as ttk
# Import pygubu objects necessary for pyinstaller to work
import pygubu.builder.tkstdwidgets
import pygubu.builder.ttkstdwidgets
import pygubu.builder.widgets.scrollbarhelper
import pygubu.builder.widgets.tkscrollbarhelper
from tkinter import messagebox, ANCHOR, PhotoImage, NW, Scrollbar, END
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

try:
    CDIR = os.path.abspath(os.path.dirname(__file__))
except NameError:
    CDIR = os.path.dirname(os.path.abspath(sys.argv[0]))
PROJECT_UI = os.path.join(CDIR, "markinchi_gui.ui")

class MarkinchiGuiApp(object):
    def __init__(self, master=None):

        # Get the pygubu builder that will build tkinter widgets from
        # XML file
        self.builder = builder = pygubu.Builder()
        builder.add_resource_path(CDIR)
        builder.add_from_file(PROJECT_UI)
        self.mainwindow = builder.get_object('frame1', master)
        self.master = master
        master.resizable(False, False)  # Make the program window not resizable
        master.title("MarkInChI to List of InChIs")  # Set the title of the program
        # build menu
        my_menu = builder.get_object('menu1')
        master.config(menu=my_menu)
        # build entry and listbox
        self.entry = builder.get_object('entry1')
        self.listbox = self.builder.get_object('listbox1')
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
        self.listbox.delete(0, END)
        markinchi = entry.get()
        c_a = markinchi.find("<>") != -1
        c_b = markinchi.find("MarkInChI") == -1
        if markinchi.find("<M>") == -1 or c_a or c_b:
            # Error Message
            msg = "Please enter a valid MarkInChI with correct separators <M>"
            messagebox.showerror("Error", msg)
        else:
            try:
                list_of_inchi = self.get_inchis(markinchi)
                i = 1
                print(f"Number of inchi produced: {len(list_of_inchi)}")
                for inchi in list_of_inchi:
                    self.listbox.insert(i, inchi)
                    i += 1
            except AttributeError as msg:
                ex_type, ex_value, ex_traceback = sys.exc_info()
                trace_back = traceback.extract_tb(ex_traceback)
                msg1 = "please see error message below:\n"
                messagebox.showerror("Error",
                                     msg1+repr(msg)+"\n"+str(trace_back))


    def get_inchis(self, inchi):

        """ This function converts a markinchi to a list of inchis"""

        inchi_obj = MarkInChI(inchi)
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

        self.entry.delete(0, END)
        self.listbox.delete(0, END)
        path = os.getcwd()
        file = askopenfile(mode='r', initialdir=path, defaultextension='.sdf',
                           filetypes=[("SDF file", ".sdf")])
        if file is not None:
            name = os.path.basename(file.name)
            content = file.readlines()
            mark_obj = markmol()
            zz = zz_convert()
            mark_obj.help_label = Label()
            mark_obj.atom_symbols = {}
            mark_obj.Rpositions = []  # number of each Te atom
            mark_obj.Rsubstituents = []
            mark_obj.ctabs = []  # no. of each substituent after $RGP
            mark_obj.connections = []
            mark_obj.list_of_atoms = {}
            mark_obj.attachments = []
            mark_obj.attach_ids = {}
            mark_obj.no_atoms = 0
            mark_obj.attach_reordered = False
            mark_obj.main_dict_renumber = {}
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
            print(f"ctabs: {ctabs}")
            while i < len(ctabs) - 1:
                mark_obj.Rsubstituents.append(substituents[int(ctabs[i-1]):int(ctabs[i])])
                i += 1
            else:
                if len(ctabs) > 0:
                    mark_obj.Rsubstituents.append(substituents[int(ctabs[i-1]):])
            divide_subst = len(mark_obj.subst_order)
            R_subst = mark_obj.Rsubstituents[:divide_subst]
            other_subst = mark_obj.Rsubstituents[divide_subst:]
            subst_dict = dict(zip(mark_obj.subst_order, R_subst))
            subst_dict = dict(sorted(subst_dict.items(), key=lambda item: item[0]))
            mark_obj.Rsubstituents = list(subst_dict.values()) + other_subst

            try:
                markinchi_final = mark_obj.produce_markinchi()
                print(markinchi_final)
                self.entry.insert(0, markinchi_final)
            except IndexError:
                msg = "Please enter a valid V2000 Markush sdf file"
                messagebox.showerror("Error", msg)
            file.close()
        else:
            msg = "Please choose a file"
            messagebox.showerror("Error", msg)

    def quit(self):

        self.master.quit()

    def run(self):
        self.mainwindow.mainloop()  # show window


if __name__ == '__main__':
    # start the app
    root = tk.Tk()
    app = MarkinchiGuiApp(root)
    app.run()
