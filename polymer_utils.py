# Finished class definitions and functions for polymer loading utilities

from email.policy import default
import sys
from rdkit import Chem
from pathlib import Path
from random import randint
from copy import deepcopy
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.utils import OpenEyeToolkitWrapper
# from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem
from collections import defaultdict
import py3Dmol
from IPython.utils import io
import ipywidgets as widgets
from IPython.display import Javascript, display, clear_output

class ParsePolymer:
    instances = {}
    def __init__(self, file):
        import time
        # store an instances id to differentiate instances later
        self.instance_id = int(time.time() * 1000.0)
        self.__class__.instances[self.instance_id] = self
        print(self.instance_id)
        # store and load file 
        self.file = file
        self.error_status = False
        # make sure file exists
        path = Path(file)
        if not path.exists():
            print(f"file path {str(path)} does not exist, returning from class init")
            self.error_status = True
            return

        # initialize other class vars 
        self.assigned_atoms = set()
        self.unassigned_atoms = set()
        self.assigned_bonds = set()
        self.unassigned_bonds = set()
        self.monomer_ids = None
        self.sampled_molecule = None
        self.sample_seed = None
        self.monomer_context_ids = None
        # load the molecule 
        rdmol = Chem.rdmolfiles.MolFromPDBFile(self.file, removeHs=False)
        for atom in rdmol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
            atom.SetProp("atomLabel", atom.GetSymbol() + f"{atom.GetAtomMapNum()}")

        self.full_molecule = rdmol
        self.n_atoms = rdmol.GetNumAtoms()
        self.atoms_str = None
        self.atoms_serial_str = None
        self.bonds_str = None
        self.double_bond_list = []

    def sample_molecule(self, seed=-1, size=80):
        # if seed == -1, when the function will chose a random seed location in the molecule
        # if seed is outside the bounds of the molecule, an error message will play

        # select a seed location if not specified
        if seed == -1:
            seed = randint(0, self.n_atoms - 1)
        elif seed >= self.n_atoms:
            print(f"seed ({seed}) is too large for given molecule of size ({self.n_atoms}) -> (index starts from 0)")
        elif seed < -1:
            print(f"seed ({seed}) cannot be negative")

        # use the seed location and select a subset of the molecule 
        # until the size is met or an entire subsection of the molecule is selected
        def get_substructure_ids(rdmol, seed, size):
            found = [] # atom ids found that we need not go over again
            active = [seed] # atom ids where the algorithm is currently centered 
            loop_condition = True
            while len(active) != 0 and loop_condition:
                for active_idx in active:
                    active_atom = rdmol.GetAtomWithIdx(active_idx)
                    for neighbor in active_atom.GetNeighbors():
                        idx = neighbor.GetIdx()
                        if (idx in found) or (idx in active):
                            continue
                        else:
                            active.append(idx)
                    active.remove(active_idx)
                    found.append(active_idx)
                    if len(found) >= size:
                        loop_condition = False
                        break
            return found

        rd_ids = get_substructure_ids(self.full_molecule, seed, size)
        if len(rd_ids) < size:
            print(f"molecule substructure is smaller than the specificed size of {size}.")

        # smarts solution using rdkits smarts function
        sub_smarts = Chem.MolFragmentToSmarts(self.full_molecule, atomsToUse = rd_ids)
        sub_rdmol = Chem.rdmolfiles.MolFromSmarts(sub_smarts)

        # output the seed location and the selected subset rdmolecule 
        self.sample_seed = seed
        self.sampled_molecule = sub_rdmol

    def generate_ids(self, seed, exclusive_wall_ids=[], inclusive_wall_ids=[]):
        wall_ids = [*exclusive_wall_ids, *inclusive_wall_ids]
        found = [] # atom ids found that we need not go over again
        active = [seed] # atom ids where the algorithm is currently centered 
        n = 0
        while len(active) != 0 and n < 2*self.n_atoms:
            for active_idx in active:
                active_atom = self.full_molecule.GetAtomWithIdx(active_idx)
                for neighbor in active_atom.GetNeighbors():
                    idx = neighbor.GetIdx()
                    if (idx in found) or (idx in active):
                        continue
                    elif (idx in wall_ids):
                        if idx in inclusive_wall_ids:
                            found.append(idx)
                        continue
                    else:
                        active.append(idx)
                active.remove(active_idx)
                found.append(active_idx)
        return found

    def test_polymer_load(self, visualize=False, **kwargs):
        # executes from_pdb using the given biopolymer substructure library and records how many
        # atoms and bonds have chemical information assigned.
        # returns: original molecule ids that are assigned and have assigned bonds 

        mol = Molecule.from_pdb(self.file)
        for atom in mol.atoms:
            if atom.metadata['already_matched']:
                self.assigned_atoms.add(atom.molecule_atom_index)
            else:
                self.unassigned_atoms.add(atom.molecule_atom_index)
        for bond in mol.bonds:
            # check for assigned bonds 
            if bond.bond_order in [1,2,3]:
                self.assigned_bonds.add((bond.atom1_index, bond.atom2_index))
            else:
                self.unassigned_bonds.add((bond.atom1_index, bond.atom2_index))
                
        # print info on the polymer loading
        print(f"number of atoms assigned: {len(self.assigned_atoms)}")
        print(f"number of bonds assigned: {len(self.assigned_bonds)}")
        print(f"number of atoms not assigned: {len(self.unassigned_atoms)}")
        print(f"number of bonds not assigned: {len(self.unassigned_bonds)}")
        if visualize:
            view = None
            return view
        else: 
            return mol
    
    def visualize(self, mode="rdkit", width=500, height=300, show_all_hydrogens=True, highlight_ids=[]):
        # change rdkit map nums to zero
        if not isinstance(self.sampled_molecule, Chem.rdchem.Mol):
            rdmol_new = deepcopy(self.full_molecule)
        else:
            rdmol_new = deepcopy(self.sampled_molecule)
        
        if mode == "rdkit":
            from IPython.display import SVG
            from rdkit.Chem.Draw import (  # type: ignore[import]
                rdDepictor,
                rdMolDraw2D,
            )
            from rdkit.Chem.rdmolops import RemoveHs  # type: ignore[import]
            [atom.SetAtomMapNum(0) for atom in rdmol_new.GetAtoms()]

            width = 500 if width is None else width
            height = 300 if height is None else height

            if not show_all_hydrogens:
                # updateExplicitCount: Keep a record of the hydrogens we remove.
                # This is used in visualization to distinguish eg radicals from normal species
                rdmol_new = RemoveHs(rdmol_new, updateExplicitCount=True)

            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(rdmol_new)
            rdmol_new = rdMolDraw2D.PrepareMolForDrawing(rdmol_new)

            drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
            # drawer.drawOptions().addAtomIndices = True
            # drawer.drawOptions().addLabels = True
            drawer.drawOptions().includeAtomTags = True
            drawer.drawOptions().maxFontSize = 10
            drawer.drawOptions().minFontSize = 7
            # drawer.drawOptions().circleAtoms = True
            if len(highlight_ids) > 0:
                drawer.DrawMolecule(rdmol_new, highlightAtoms=highlight_ids)
            else:
                drawer.DrawMolecule(rdmol_new)
            drawer.FinishDrawing()

            return SVG(drawer.GetDrawingText())
        elif mode == "nglview":
            # use nglview for this
            import nglview as nv

            # first, assign all atoms to be highlighted to their own ligand group
            grey_highlighted_info = Chem.AtomPDBResidueInfo()
            grey_highlighted_info.SetResidueName("asddsf") #completely random name that must be > 5 letters for some reason? 

            red_highlighted_info = Chem.AtomPDBResidueInfo()
            red_highlighted_info.SetResidueName("zxcvbn") #completely random name that must be > 5 letters for some reason? 

            unhighlighted_info = Chem.AtomPDBResidueInfo()
            unhighlighted_info.SetResidueName("sfdfsg") #completely random name that must be > 5 letters for some reason? 
            for atom in rdmol_new.GetAtoms():
                if atom.GetIdx() in highlight_ids:
                    atom.SetMonomerInfo(red_highlighted_info)
                else:
                    atom.SetMonomerInfo(unhighlighted_info)

            # next, construct the viewer         
            view = nv.show_rdkit(rdmol_new)

            # grey highlights
            view.add_spacefill(".ZXC", opacity=0.4, color="red")
            view.add_spacefill(".ZXC and _H", opacity=0.4, color="white")

            return view 
        elif mode == "py3Dmol":
            # for rdmol_new, assign new residue names and atoms names to map between
            # rdkit and the py3Dmol later
            # atom names go from H1 to H99, or C1 to C99, etc
            # residue names go from AAA to ZZZ
            # this should allow for something like 1382400 hydrogens,carbons,etc
            element_counts = defaultdict(int)
            # dictionary for maping map to the rdmol later 
            map_dict = dict()
            for atom in rdmol_new.GetAtoms():
                element_count = element_counts[atom.GetAtomicNum()]
                res_num = element_count/99
                atom_num = (element_count%99) + 1
                res_string = f"{chr(int(res_num/26/26%26)+65)}{chr(int(res_num/26%26)+65)}{chr(int(res_num%26)+65)}"
                atom.GetPDBResidueInfo().SetName(f"{atom.GetSymbol()}{atom_num:02d}")
                atom.GetPDBResidueInfo().SetResidueName(res_string)
                atom.GetPDBResidueInfo().SetResidueNumber(int(res_num) + 1)
                element_counts[atom.GetAtomicNum()] += 1
                map_dict[f"{int(res_num) + 1}{atom.GetSymbol()}{atom_num:02d}"] = atom.GetAtomMapNum()

            # create pdb block and view object  
            pdb_block = Chem.rdmolfiles.MolToPDBBlock(rdmol_new)
            view = py3Dmol.view(width=width, height=height)
            view.addModel(pdb_block,'pdb', {"keepH": True})
            view.setStyle({"model": -1}, {"stick": {}})
            original_clicker_code = '''function(atom,viewer,event,container) {
                   if(!atom.label) {
                        atom.label = viewer.addLabel(atom.atom+":"+atom.serial,{position: atom, backgroundColor: 'mintcream', fontColor:'black', backgroundOpacity: "0.3"});
                        console.log(Object.getOwnPropertyNames(atom))
                        console.log(atom.serial)
                        viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFF0000}})
                        
                        atom_data = document.getElementById("atoms_of_class_%d");
                        var a = atom_data.dataset.atoms;
                        var s = atom_data.dataset.atoms_serial;

                        var arr = a.split(',');
                        arr.push(atom.resi.toString() + atom.atom);
                        atom_data.dataset.atoms = arr;

                        var arr_s = s.split(',');
                        arr_s.push(atom.serial);
                        atom_data.dataset.atoms_serial = arr_s;
                   }
                   else {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                        viewer.setStyle({"serial": atom.serial}, {"stick": {"colorscheme": "default"}})
                        
                        atom_data = document.getElementById("atoms_of_class_%d");
                        var a = atom_data.dataset.atoms;
                        var s = atom_data.dataset.atoms_serial;
                        var arr = a.split(',');
                        var arr_s = s.split(',');
                        
                        var atom_str = atom.resi.toString() + atom.atom;
                        arr = arr.filter(function(item) {
                            return item !== atom_str
                        })
                        atom_data.dataset.atoms = arr;

                        arr_s = arr_s.filter(function(item) {
                            return Number(item) !== atom.serial
                        })
                        atom_data.dataset.atoms_serial = arr_s;
                   }
                   viewer.render();}''' % (self.instance_id, self.instance_id)
            view.setClickable({},True,original_clicker_code)
            # create an html div element to store atom indices
            code = """
                    var instance_id = "atoms_of_class_%d";
                    var element = document.getElementById(instance_id);
                    if(element){
                        element.dataset.atoms = []
                        element.dataset.atoms_serial = []
                    } else {
                        var atom = document.createElement("div");
                        atom.innerHTML = "";
                        atom.id = instance_id
                        atom.dataset.atoms = [];
                        atom.dataset.atoms_serial = [];
                        document.body.appendChild(atom);
                    }
                    """ % (self.instance_id)
            display(Javascript(code))
            # button to export atoms
            clear_button = widgets.Button(description="Clear", button_style="danger")
            button1 = widgets.Button(description="[#1] Finalize Selection")
            button2 = widgets.Button(description="[#2] Print Selection")
            double_bonds_button = widgets.Button(description="Choose Double bonds")
            output = widgets.Output()
            
            def clear_selection(b):
                code = """
                        var atom_data = document.getElementById("atoms_of_class_%d");
                        atom_data.dataset.atoms = []
                        """ % (self.instance_id)
                display(Javascript(code))
                clear_output(wait=False)
                view.setStyle({"model": -1}, {"stick": {"colorscheme": "default"}})
                view.removeAllLabels()
                display(view, widgets.HBox((clear_button, button1, button2,double_bonds_button)), output, display_id=True)
            def finalize_selection(b):
                global atoms_str
                atoms_str = None
                code = """
                        var atom_data = document.getElementById("atoms_of_class_%d");
                        var a = JSON.stringify(atom_data.dataset.atoms);
                        var s = JSON.stringify(atom_data.dataset.atoms_serial);
                        Jupyter.notebook.kernel.execute("atoms_str = " + a);
                        Jupyter.notebook.kernel.execute("ParsePolymer.instances[%d].atoms_str = " + a);
                        Jupyter.notebook.kernel.execute("ParsePolymer.instances[%d].atoms_serial_str = " + s);
                        console.log("exporting atom serial ids");
                        """ % (self.instance_id, self.instance_id, self.instance_id)
                display(Javascript(code))
                return
            def print_selection(b):
                string = self.atoms_str
                if string == '':
                    print("must select some atoms before printing")
                    return
                atom_map_ids = [map_dict[i] for i in string.split(",") if i != '']
                print(atom_map_ids)
                self.monomer_ids = atom_map_ids
                new_smarts, smarts_block = self.generate_smarts_entry()
                print(smarts_block)
            def double_bonds_click(b):
                if b.description == "Choose Double bonds": # initiate "choosing double bonds mode"
                    labels_string = self.atoms_str
                    serial_string = self.atoms_serial_str
                    if serial_string == '' or serial_string == None:
                        print("must select some atoms before attempting to set double bonds")
                        return
                    b.description = "Set Double bonds"
                    b.button_style="warning"
                    # get old red atoms and labels
                    atom_serial_ids = [int(i) for i in serial_string.split(",") if i != '']
                    atom_labels = [i for i in labels_string.split(",") if i != '']
                    # create new html element that will store which double bonds are clicked
                    code = """
                            var instance_id = "atoms_of_class_%d";
                            var element = document.getElementById(instance_id);
                            if(element){
                                element.dataset.double_bonds = [];
                            } else {
                                var atom = document.createElement("div");
                                atom.innerHTML = "";
                                atom.id = instance_id
                                atom.dataset.atoms = [];
                                atom.dataset.double_bonds=[];
                                document.body.appendChild(atom);
                            }
                            """ % (self.instance_id)
                    display(Javascript(code))
                    # set new clickable for the view object and reload
                    code = '''function(atom,viewer,event,container) {
                        var atom_data = document.getElementById("atoms_of_class_%d");
                        var s = atom_data.dataset.atoms_serial;
                        var arr_s = s.split(',');
                        arr_s = arr_s.filter(function(item) {
                            return Number(item) == atom.serial
                        })
                        if (!!arr_s.length) 
                        {
                            if(atom.color == 0xFFFF00) {
                                    viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFF0000}});
                                    atom.color = 0xFF0000;

                                    var a = atom_data.dataset.double_bonds;
                                    var arr = a.split(',');
                                    
                                    var atom_str = atom.resi.toString() + atom.atom + ":" + atom.serial + ":" + atom.x + ":" + atom.y + ":" + atom.z;
                                    arr = arr.filter(function(item) {
                                        return item !== atom_str;
                                    })
                                    atom_data.dataset.double_bonds = arr;
                            }
                            else {
                                    viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFFFF00}});
                                    atom.color = 0xFFFF00;

                                    var a = atom_data.dataset.double_bonds;
                                    var arr = a.split(',');
                                    arr.push(atom.resi.toString() + atom.atom + ":" + atom.serial + ":" + atom.x + ":" + atom.y + ":" + atom.z);
                                    atom_data.dataset.double_bonds = arr;      
                            }
                        }
                        viewer.render();}''' % (self.instance_id)
                    view.setClickable({},True,code)
                    # set all of the old labels and red colors and reload
                    for label,serial in zip(atom_labels,atom_serial_ids):
                        #TODO: improve this:
                        symbol = label
                        for i in range(0, len(label)):
                            character = label[i]
                            if not character.isdigit():
                                break
                            symbol = symbol[1:]
                        view.setStyle({"serial": serial}, {"stick": {"color": 0xFF0000}})
                        view.addLabel(f"{symbol}:{serial}",{'backgroundColor': 'mintcream', 'fontColor':'black', 'backgroundOpacity': 0.3},{"serial":serial})
                        
                    clear_output(wait=False)
                    display(view, widgets.HBox((clear_button, button1, button2,double_bonds_button)), output, display_id=True)

                elif b.description == "Set Double bonds": # user has chosen double bonds which now need to be set
                    b.description = "Update Molecule"
                    b.button_style="warning"
                    # get the atoms that have been highlighted
                    code = """
                            var atom_data = document.getElementById("atoms_of_class_%d");
                            var b = JSON.stringify(atom_data.dataset.double_bonds);
                            Jupyter.notebook.kernel.execute("ParsePolymer.instances[%d].bonds_str = " + b);
                            console.log("exporting atom serial ids");
                            """ % (self.instance_id, self.instance_id)
                    display(Javascript(code))
                elif b.description == "Update Molecule":
                    b.description = "Choose Double bonds"
                    b.button_style=""
                    # resets all colors and labels and represents double bond in the figure
                    labels_string = self.atoms_str
                    serial_string = self.atoms_serial_str
                    # get old red atoms and labels
                    atom_serial_ids = [int(i) for i in serial_string.split(",") if i != '']
                    atom_labels = [i for i in labels_string.split(",") if i != '']
                    # reset the double bonds info in the html
                    code = """
                            var instance_id = "atoms_of_class_%d";
                            var element = document.getElementById(instance_id);
                            if(element){
                                element.dataset.double_bonds = [];
                            } else {
                                var atom = document.createElement("div");
                                atom.innerHTML = "";
                                atom.id = instance_id
                                atom.dataset.atoms = [];
                                atom.dataset.double_bonds=[];
                                document.body.appendChild(atom);
                            }
                            """ % (self.instance_id)
                    display(Javascript(code))
                    # original code for the clicker
                    view.setClickable({},True,original_clicker_code)
                    # set all of the old labels and red colors and reload as was done before
                    for label,serial in zip(atom_labels,atom_serial_ids):
                        #TODO: improve this:
                        symbol = label
                        for i in range(0, len(label)):
                            character = label[i]
                            if not character.isdigit():
                                break
                            symbol = symbol[1:]
                        view.setStyle({"serial": serial}, {"stick": {"color": 0xFF0000}})
                        view.addLabel(f"{symbol}:{serial}",{'backgroundColor': 'mintcream', 'fontColor':'black', 'backgroundOpacity': 0.3},{"serial":serial})
                    # get bond information
                    bonds_list = self.bonds_str.split(",")
                    bonds_list = [b for b in bonds_list if b != ""]
                    if len(bonds_list) > 2:
                        print("can only select two atoms at a time when specifying double bonds")
                        return
                    bond1 = bonds_list[0].split(":")
                    bond2 = bonds_list[1].split(":")
                    label1 = bond1[0]
                    label2 = bond2[0]
                    x1 = float(bond1[2])
                    y1 = float(bond1[3])
                    z1 = float(bond1[4])
                    x2 = float(bond2[2])
                    y2 = float(bond2[3])
                    z2 = float(bond2[4])

                    view.addCylinder({"start":{"x":x1,"y":y1,"z":z1},
                                      "end":{"x":x2,"y":y2,"z":z2},
                                      "radius": 0.35,
                                      "fromCap":0,
                                      "toCap":0,
                                      "dashed": True,
                                      "opacity": "0.6"})

                    clear_output(wait=False)
                    print(x1, y1, z1, x2, y2, z2)
                    display(view, widgets.HBox((clear_button, button1, button2,double_bonds_button)), output, display_id=True)
                    # store double bond info in self.double_bond_list
                    self.double_bond_list.append(tuple(sorted([map_dict[label1],map_dict[label2]])))

            clear_button.on_click(clear_selection)
            button1.on_click(finalize_selection)
            button2.on_click(print_selection)
            double_bonds_button.on_click(double_bonds_click)

            return display(view, widgets.HBox((clear_button, button1, button2,double_bonds_button)), output, display_id=True)
 
    def generate_monomer_smarts(self, assign_map_ids_to=[]):
        # if assign_map_ids_to is empty, assume map ids to monomer 
        if len(self.monomer_ids) == 0:
            print("must specify monomer ids to generate library entry, monomer smarts, or library charges")
            return
        # assume either Molecule or RDKit molecule
        if len(assign_map_ids_to) == 0:
            assign_map_ids_to = self.monomer_ids
        if self.monomer_context_ids == None:
            smarts_ids = self.monomer_ids
        else:
            smarts_ids = self.monomer_context_ids
        if self.sampled_molecule != None:
            rdmol_copy = deepcopy(self.sampled_molecule)
        else:
            rdmol_copy = deepcopy(self.full_molecule)
        # map from atom map numbers to molecule indices
        smarts_mapped_ids = []
        # important that this is done in a separate loop
        for a in rdmol_copy.GetAtoms():
            if a.GetAtomMapNum() in smarts_ids: # smarts_ids are given as map numbers
                    smarts_mapped_ids.append(a.GetIdx()) # and must be convered to indices
                    # ^^this has no effect for unsampled molecules^^
        if len(self.double_bond_list) > 0:
            for bond in rdmol_copy.GetBonds():
                bond_ids = tuple(sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]))
                if bond_ids in self.double_bond_list:
                    bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
        n=1
        for a in rdmol_copy.GetAtoms():
            if a.GetAtomMapNum() in smarts_ids:
                if a.GetAtomMapNum() in assign_map_ids_to:
                    a.SetAtomMapNum(n)
                    n += 1
                else:
                    a.SetAtomMapNum(0)
            
        print(smarts_mapped_ids)
        smarts = Chem.MolFragmentToSmarts(rdmol_copy, atomsToUse = smarts_mapped_ids)

        return smarts, rdmol_copy
        
    def _smarts_string_cleaner(self, smarts):
        atomic_nums = {"#1": "H",
                    "#6": "C",
                    "#7": "N",
                    "#8": "O",
                    "#15": "P",
                    "#16": "S"
                    }
        OneLetterSymbols = ["H", "C", "c", "N", "O", "P", "S", "F"]
        TwoLetterSymbols = ["Cl"]
        # give atomic numbers there corresponding symbols if needed
        n = 1
        i = 0
        searching = True
        element_counts = defaultdict(int)
        elements_list = []
        while searching:
            select = smarts[i:i+3]
            if select[0] == "#":
                # replace atomic number identifier with map indexed symbol
                if select[2] in "1234567890": # a two digit atomic number
                    atomic_num_string = select[:3] # select all two digits
                    symbol = atomic_nums[atomic_num_string]
                    new_smarts = smarts[:i] + symbol + smarts[i+3:]
                    smarts = new_smarts
                else:
                    atomic_num_string = select[:2] # select only the single digit
                    symbol = atomic_nums[atomic_num_string]
                    new_smarts = smarts[:i] + symbol + smarts[i+2:]
                    smarts = new_smarts
                # attempt to add a symbol now indexed based on fellow atoms with the same symbol
                element_counts[symbol] += 1
                elements_list.append(symbol + str(element_counts[symbol]))
            elif select[:2] in TwoLetterSymbols:
                symbol = select
                element_counts[symbol] += 1
                elements_list.append(symbol + str(element_counts[symbol]))
            elif select[0] in OneLetterSymbols:
                symbol = select[0]
                element_counts[symbol] += 1
                elements_list.append(symbol + str(element_counts[symbol]))
            if (i + 1) == len(smarts):
                break
            i += 1
        # assign map numbers starting from the first atom
        block = "\"" + smarts + "\": ["
        for element in elements_list[:-1]:
            block += "\n\t\""
            block += element
            block += "\","
        block += "\n\t\""
        block += elements_list[-1]
        block += "\"\n]"
        return smarts, block

    def generate_smarts_entry(self):
        smarts,_ = self.generate_monomer_smarts()
        new_smarts, smarts_block = self._smarts_string_cleaner(smarts)
        return new_smarts, smarts_block

    def generate_library_charges_am1(self, toolkit_registry=OpenEyeToolkitWrapper()):
        # creates charges for a monomer within a context
        # the monomer must be a subset molecule of the context 
        # custom_charges is in the form of a dict, where the keys are the map
        # indexes of the monomer and the values are the charges for atoms with
        # that index
        # smarts, sub_rdmol = self.generate_monomer_smarts()

        mol = self.test_polymer_load()
        rdmol = mol.to_rdkit()
        n=1
        for a in rdmol.GetAtoms():
            if a.GetIdx() in self.monomer_ids:
                a.SetAtomMapNum(n)
                n += 1
            else:
                a.SetAtomMapNum(0)
        if self.monomer_context_ids == None:
            print("monomer context must be specified")
            return
        sub_smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse = self.monomer_context_ids)
        sub_rdmol = Chem.rdmolfiles.MolFromSmarts(sub_smarts)

        # add some random chains to the sub_rdmol to both fix any
        # neighbor spots where the molecule was cut but also to 
        # mitigate any major charge shifts due to neighbors. However, 
        # if neighbors significantly affect charges, then the appropriate 
        # molecule_context_ids should be set. 

        mw = Chem.RWMol(sub_rdmol)
        # manually modify the rdmol depending on the valencies 
        #!!!!!!!!!THIS ONLY WORKS NOW FOR CARBON ATOM LINKS !!!!!!!!!!!!!!!!
        def incomplete_atom(atom):
            is_incomplete = False
            if atom.GetAtomicNum() == 6:
                if (len([n for n in atom.GetNeighbors()]) + len([n for n in atom.GetBonds() if n.GetBondType() == Chem.rdchem.BondType.DOUBLE])) < 4:
                    is_incomplete = True
            return is_incomplete

        for atom in sub_rdmol.GetAtoms():
            if incomplete_atom(atom):
                # add a carbon end group
                C_idx = mw.AddAtom(Chem.Atom(6))
                H1_idx = mw.AddAtom(Chem.Atom(1))
                H2_idx = mw.AddAtom(Chem.Atom(1))
                H3_idx = mw.AddAtom(Chem.Atom(1))
                mw.AddBond(atom.GetIdx(), C_idx, Chem.rdchem.BondType.SINGLE)
                mw.AddBond(C_idx, H1_idx, Chem.rdchem.BondType.SINGLE)
                mw.AddBond(C_idx, H2_idx, Chem.rdchem.BondType.SINGLE)
                mw.AddBond(C_idx, H3_idx, Chem.rdchem.BondType.SINGLE)
        sub_rdmol = Chem.rdchem.Mol(mw)
        for atom in sub_rdmol.GetAtoms():
            if incomplete_atom(atom):
                # add hydrogen
                H1_idx = mw.AddAtom(Chem.Atom(1))
                mw.AddBond(atom.GetIdx(), H1_idx, Chem.rdchem.BondType.SINGLE)

        charged_mol = Molecule.from_rdkit(Chem.rdchem.Mol(mw), allow_undefined_stereo=True, hydrogens_are_explicit=True)

        charged_mol.assign_partial_charges(
            partial_charge_method="am1bcc", toolkit_registry=toolkit_registry
        )

        for off_atom, rdatom in zip(charged_mol.atoms, mw.GetAtoms()):
            if rdatom.GetAtomMapNum() != 0:
                off_atom.metadata.update({"map_num": rdatom.GetAtomMapNum()})
            else:
                off_atom.metadata.update({"map_num": 0})

        return sub_smarts, charged_mol

    def generate_library_charge_entry(self, custom_charges=dict()):
        if len(custom_charges) == 0:
            smarts, charged_mol = self.generate_library_charges_am1()
            # unpack the charges 
            charges = charged_mol.partial_charges
            custom_charges = dict()
            for atom, charge in zip(charged_mol.atoms, charges):
                map_num = atom.metadata["map_num"]
                if map_num != 0:
                    custom_charges[map_num] = charge
        # create a library charge entry 
        entry_string = "<LibraryCharge smirks=\"" + smarts + "\""
        for map_num, charge in custom_charges.items():
            charge_string = f"{charge:.8g}"
            charge_string = charge_string.replace("e", "* elementary_charge\"")
            entry_string = entry_string + f" charge{map_num}=\"" + charge_string
        entry_string = entry_string + "></LibraryCharge>"
        return entry_string

if __name__ == "__main__":
    file = "openff_polymer_testing/polymer_examples/rdkit_simple_polymers/naturalrubber.pdb"
    p = ParsePolymer(file)
    # p.monomer_ids = p.generate_ids(102, [79], [])
    # p.monomer_context_ids = p.generate_ids(102, [79], [])
    # entry = p.generate_library_charge_entry()
    # print(entry)

    p.visualize(mode="py3Dmol")
        