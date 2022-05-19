import sys
from pyrsistent import T
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
import time
import os
import tempfile

class ChemistryEngine:
    instances = {}
    def __init__(self, file):
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
        self.selected_monomers = None # dictinoary of smarts of previously found monomers
        self.sampled_molecule = None # used when sampling a subset
        self.sample_seed = None      # of the molecule
        # load the molecule 
        suffix = path.suffix.lower()
        if 'pdb' in suffix:
            rdmol = Chem.rdmolfiles.MolFromPDBFile(self.file, removeHs=False)
        elif 'sdf' in suffix:
            suppl = Chem.rdmolfiles.SDMolSupplier(self.file, removeHs=False)
            rdmol = suppl[0]
        for atom in rdmol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
            atom.SetProp("atomLabel", atom.GetSymbol() + f"{atom.GetAtomMapNum()}")

        self.full_molecule = rdmol
        self.n_atoms = rdmol.GetNumAtoms()
        self.atoms_str = None
        self.atoms_serial_str = None
        self.bonds_str = None
        self.double_bond_list = []

        self.map_molecule = None

    def generate_monomer_smarts(self, monomer_ids, monomer_context_ids=[], assign_map_ids_to=[]):
        # if assign_map_ids_to is empty, assume map ids to monomer 
        if len(monomer_ids) == 0:
            print("must specify monomer ids to generate library entry, monomer smarts, or library charges")
            return
        # assume either Molecule or RDKit molecule
        if len(assign_map_ids_to) == 0:
            assign_map_ids_to = monomer_ids
        if len(monomer_context_ids) == 0:
            smarts_ids = monomer_ids
        else:
            smarts_ids = monomer_context_ids
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
                else:
                    atomic_num_string = select[:2] # select only the single digit
                    symbol = atomic_nums[atomic_num_string]
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
        return block

    def generate_smarts_entry(self, monomer_ids):
        smarts,_ = self.generate_monomer_smarts(monomer_ids)
        smarts_block = self._smarts_string_cleaner(smarts)
        return smarts, smarts_block

    def test_polymer_load(self):
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
        return self.assigned_atoms, self.assigned_bonds, self.unassigned_atoms, self.unassigned_bonds

    def create_map_dict(self):
        map_molecule = deepcopy(self.full_molecule)
        # this should allow for something like 1382400 hydrogens,carbons,etc
        element_counts = defaultdict(int)
        # dictionary for maping map to the rdmol later 
        map_dict = dict()
        for atom in map_molecule.GetAtoms():
            element_count = element_counts[atom.GetAtomicNum()]
            res_num = element_count/99
            atom_num = (element_count%99) + 1
            res_string = f"{chr(int(res_num/26/26%26)+65)}{chr(int(res_num/26%26)+65)}{chr(int(res_num%26)+65)}"
            atom.GetPDBResidueInfo().SetName(f"{atom.GetSymbol()}{atom_num:02d}")
            atom.GetPDBResidueInfo().SetResidueName(res_string)
            atom.GetPDBResidueInfo().SetResidueNumber(int(res_num) + 1)
            atom.GetPDBResidueInfo().SetSerialNumber(int(res_num)) # zero indexed
            element_counts[atom.GetAtomicNum()] += 1
            map_dict[f"{int(res_num) + 1}{atom.GetSymbol()}{atom_num:02d}"] = atom.GetAtomMapNum()
        pdb_block = Chem.rdmolfiles.MolToPDBBlock(map_molecule)
        self.map_molecule = map_molecule # used to convert between mapping and atomMapNum
        return pdb_block, map_dict

    def _ids_to_map_ids(self, ids):
        # takes as input the ids of a molecule and outputs 
        # a list of the map numbers of those atoms
        # **will only become useful when I start sampling molecules
        if isinstance(ids, list):
            map_nums = []
            for atom in self.full_molecule.GetAtoms():
                if atom.GetIdx() in ids:
                    map_nums.append(atom.GetAtomMapNum())
            return map_nums
        elif isinstance(ids, int):
            for atom in self.full_molecule.GetAtoms():
                if atom.GetIdx() == ids:
                    return atom.GetAtomMapNum()
            return None
        return None
    
    def symbol_from_map_id(self, map_id):
        for atom in self.full_molecule.GetAtoms():
            if atom.GetAtomMapNum() == map_id:
                return atom.GetSymbol()
        return None

    def get_file_block(self, format, additional_specs={}):
        # format: 'sdf' or 'pdb' at the moment
        # additional_specs: dictionary of chemistry specs indexed by atom map number
        # example to get an sdf with added double bond:
        # block = self.get_file_block('sdf', {'double': (21,23)})
        mol_copy = deepcopy(self.full_molecule)
        # modify the molecule:
        doubles = additional_specs.get('triple', default=[])
        sorted_doubles = [tuple(sorted(t)) for t in doubles]

    T   triples = additional_specs.get('triple', default=[])
        sorted_triples = [tuple(sorted(t)) for t in triples]

        for spec_name, map_ids in additional_specs.items():
            if spec_name in ['double', 'triple']:
                map_ids = tuple(sorted(map_ids))
                for bond in mol_copy.GetBonds():
                    bond_ids = tuple(sorted(tuple(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())))
                    if map_ids == bond_ids:
                        if spec_name == 'double':
                            bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
            
        if format == 'pdb':
            block, map_dict = self.create_map_dict()
            return block
        elif format == 'sdf':
            with tempfile.TemporaryDirectory() as tmpdir:
                prev_dir = os.getcwd()
                os.chdir(os.path.abspath(tmpdir))
                writer = Chem.rdmolfiles.SDWriter('molecule.sdf')
                writer.write(mol_copy)

                with open("molecule.sdf", "r") as file:
                    sdf_block = file.read()
                os.chdir(prev_dir)
            return sdf_block
        else:
            print("no valid format specified")


class PolymerVisualizer3D:
    instances = {}
    def __init__(self, chemistry_engine):
        self.instance_id = int(time.time() * 1000.0)
        self.__class__.instances[self.instance_id] = self
        self.chemistry_engine = chemistry_engine
        pdb_block, map_dict = chemistry_engine.create_map_dict()
        self.map_dict = map_dict
        self.width = 800
        self.height = 500
        self.view = self._view_from_block(pdb_block, 'pdb')
        self.view.setStyle({"model": -1}, {"stick": {}})
        # buttons
        clear_button = widgets.Button(description="Clear", button_style="danger")
        finalize_button = widgets.Button(description="Finalize Selection")
        print_button = widgets.Button(description="Print Selection", disabled=True)
        double_bonds_button = widgets.Button(description="Make Double bonds", disabled=True)
        triple_bonds_button = widgets.Button(description="Make Triple bonds", disabled=True)
        formal_charge_button = widgets.Button(description="assign formal charges", disabled=True)
        formal_charge_menu = widgets.Dropdown(
                                                        options=['-3','-2','-1','0','1','2','3'],
                                                        value='0',
                                                        description='',
                                                        disabled=True,
                                                    )
        test_load_button = widgets.Button(description="Test Load", disabled=True)
        order_substructure_button = widgets.Button(description="Reorder Entries", disabled=True)
        delete_substructure_button = widgets.Button(description="Delete Entry", disabled=True)

        clear_button.on_click(self._button_clear_selection)
        finalize_button.on_click(self._button_finalize_selection)
        print_button.on_click(self._button_print_selection)
        double_bonds_button.on_click(self._button_assign_double_bonds)
        triple_bonds_button.on_click(self._button_assign_triple_bonds)
        formal_charge_button.on_click(self._button_assign_formal_charge)
        test_load_button.on_click(self._button_test_load)
        order_substructure_button.on_click(self._button_order_substructure_smarts)
        delete_substructure_button.on_click(self._button_delete_substructure_smarts)

        # any references to any widgets should only be made to this dict 
        self.widgets = {'clear': clear_button,
                        'finalize': finalize_button,
                        'print': print_button,
                        'double': double_bonds_button,
                        'triple': triple_bonds_button,
                        'charge_button': formal_charge_button,
                        'charge_menu': formal_charge_menu,
                        'test_load': test_load_button,
                        'order_smarts': order_substructure_button,
                        'delete_smarts': delete_substructure_button}

        self.highlights = {} # dict of atoms and their highlighted colors
        self.spheres = {} # dict of which atoms have spheres 
        self.selected_atoms = defaultdict(list)
        self.click_mode = "select_monomer"
        self.valid_click_modes = ['select_monomer', 'do_nothing', 'select_double', 'select_triple', 'select_charged']

    # the only function the user should ever have to call themselves 
    def show(self):
        self._reset_view()
        return

    def __repr__(self):
        print("in __repr__()")
        self.show()
        return f"file: {self.chemistry_engine.file}"

    def _view_from_block(self, block, format):
        view = py3Dmol.view(width=self.width, height=self.height)
        view.addModel(block, format, {"keepH": True})
        return view

    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #-------------------------------Clickables-------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _set_clickable(self, new_mode):
        if new_mode not in self.valid_click_modes:
            print("internal: new_mode is not valid")
            return
        self.click_mode = new_mode
        if new_mode == "select_monomer":
            selection = {"model": -1}
            code = '''function(atom,viewer,event,container) {
                    if(!atom.style['clicksphere']) {
                            viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFF0000}, "clicksphere": {"radius": 0.25}})
                            
                            var atom_identifier = atom.resi.toString() + atom.atom;
                            var serial = atom.serial;
                            Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._store_selected_atom(" + serial + ")");
                    }
                    else {
                            viewer.setStyle({"serial": atom.serial}, {"stick": {"colorscheme": "default"}})
                        
                            var atom_identifier = atom.resi.toString() + atom.atom;
                            var serial = atom.serial;
                            Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._remove_selected_atom(" + serial + ")");
                    }
                    viewer.render();}''' % (self.instance_id, self.instance_id)
        elif new_mode == 'do_nothing':
            selection = {"model": -1}
            code = '''function(atom,viewer,event,container) {
                    void(0)
                    }'''
        elif new_mode in ["select_double", "select_triple"]:
            selection = {"serial": self.selected_atoms['select_monomer']}
            code = '''function(atom,viewer,event,container) {
                        if(atom.color == 0xFFFF00) {
                                viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFF0000}});
                                atom.color = 0xFF0000;

                                var serial = atom.serial;
                                Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._remove_selected_atom(" + serial + ")");
                        }
                        else {
                                viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFFFF00}});
                                atom.color = 0xFFFF00;

                                var serial = atom.serial;
                                Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._store_selected_atom(" + serial + ")");     
                        }
                        viewer.render();}''' % (self.instance_id, self.instance_id)
        elif new_mode == "select_charged":
            selection = {"serial": self.selected_atoms['select_monomer']}
            return #not yet implemented
        else:
            print("internal error in _set_clickable")
            return

        self.view.setClickable(selection,True,code)
    
    def _store_selected_atom(self, serial):
        print(f"in the store with serial number: {serial}")
        data = self.selected_atoms[self.click_mode]
        if serial not in data:
            self.selected_atoms[self.click_mode].append(serial)
        else:
            print("error: serial number already in self.selected_atoms")
    def _remove_selected_atom(self, serial):
        data = self.selected_atoms[self.click_mode]
        if serial in data:
            self.selected_atoms[self.click_mode].remove(serial)
        else:
            print("error: serial not in self.selected_atoms")
    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #-------------------------------Hoverables-------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _set_hoverable(self, mode):
        if mode not in self.valid_click_modes:
            print("internal: new_mode is not valid")
            return
        func1 = '''
            function(atom,viewer,event,container) {
                    var id = atom.serial - 1;
                    if(!atom.label) {
                        atom.label = viewer.addLabel(atom.elem+":"+id,{position: atom, backgroundColor: 'mintcream', fontColor:'black', backgroundOpacity: "0.3"});
                    }
                }
            '''
        func2 = '''
                function(atom,viewer,event,container) {
                        if(atom.label) {
                            viewer.removeLabel(atom.label);
                            delete atom.label;
                        }
                    }
            '''
        if mode != 'select_monomer':
            selection = {'serial': self.selected_atoms['select_monomer'], 'invert': True}
        else:
            selection = {'model': -1}
        self.view.setHoverable(selection, True, func1, func2)
        self.view.setHoverDuration(100)

    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #------------------------------View Actions------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _reload_view(self, center_to_monomer=False):
        # called anytime the view is changed
        clear_output(wait=False)
        display(self.view, self.buttons)
    
    def _reset_view(self):
        # removes shapes and resets all clickables/views
        self.view.removeAllShapes()
        self.view.setStyle({"model": -1}, {"stick": {"colorscheme": "default"}})
        self.view.removeAllLabels()
        self.selected_atoms = defaultdict(list)
        self._set_clickable('select_monomer')
        self._set_hoverable('select_monomer')
        self._reset_buttons()
        self._load_buttons()
        self._reload_view()
    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #------------------------------Button Utils------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _load_buttons(self):
        # must be called any time buttons are changed
        c1 = widgets.VBox((self.widgets['clear'], self.widgets['finalize']))
        c2 = widgets.VBox((self.widgets['print'], self.widgets['order_smarts']))
        c3 = widgets.VBox((self.widgets['test_load'], self.widgets['delete_smarts']))
        c4 = widgets.VBox((self.widgets['double'], self.widgets['triple']))
        c5 = widgets.VBox((self.widgets['charge_button'], self.widgets['charge_menu']))

        self.buttons = widgets.HBox((c1,c2,c3,c4,c5))

    def _enable_all_buttons(self):
        for id, widget in self.widgets.items():
            widget.disabled = False

    def _enable_button_with_id(self, id):
        self.widgets[id].disabled = False

    def _disable_button_with_id(self, id):
        self.widgets[id].disabled = True

    def _disable_all_buttons_except(self, exceptions):
        for id, widget in self.widgets.items():
            if id not in [*exceptions, 'clear']:
                widget.disabled = True

    def _reset_buttons(self):
        # first, set all styles back to default
        for id, widget in self.widgets.items():
            if id not in ['clear', 'charge_menu']:
                widget.button_style = ''
        # then, disable buttons that are disabled on startup
        self._disable_all_buttons_except(['clear', 'finalize'])

    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #-------------------------Button Click Actions---------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _button_clear_selection(self, b):
        self._reset_buttons()

    def _button_finalize_selection(self, b):
        # cycles between "Finalize Selection" and "Edit Selection"
        if b.description == "Finalize Selection":
            b.description = "Edit Selection"
            b.button_style = "warning"
            self._enable_all_buttons()
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
            # make sure all info on highlights and labels are retained
            self.view.center({"serial": self.selected_atoms['select_monomer']})
            self.view.setStyle({"serial": self.selected_atoms['select_monomer']}, {"stick": {"colorscheme": "default"}, "clicksphere": {"radius": 0.25}})
            self.view.setStyle({"serial": self.selected_atoms['select_monomer'], "invert": True}, {"line": {}})
            self._reload_view()
        elif b.description == "Edit Selection":
            b.description = "Finalize Selection"
            b.button_style = ""
            self._reset_buttons()
            self._set_clickable('select_monomer')
            self._set_hoverable('select_monomer')
            # only keep info on what monomer atoms are selected
            # but discard info on double, triple, formal charges
            monomer_atoms = self.selected_atoms[self.click_mode]
            self.selected_atoms = defaultdict(list)
            self.selected_atoms[self.click_mode] = monomer_atoms
            self.view.center({"serial": self.selected_atoms['select_monomer']})
            self.view.setStyle({"serial": self.selected_atoms['select_monomer']}, {"stick": {"color": 0xFF0000}, "clicksphere": {"radius": 0.25}})
            self.view.setStyle({"serial": self.selected_atoms['select_monomer'], "invert": True}, {"stick": {"colorscheme": "default"}})
            self._reload_view()
        else:
            print("internal: no action for _button_finalize_selection!")

    def _button_print_selection(self, b):
        ids = [(i - 1) for i in self.selected_atoms['select_monomer']]
        map_ids = self.chemistry_engine._ids_to_map_ids(ids)
        smarts, library_entry = self.chemistry_engine.generate_smarts_entry(map_ids)
        print(library_entry)
        
    def _button_assign_double_bonds(self, b):
        # cycles between
        print("double bonds!")
    def _button_assign_triple_bonds(self, b):
        print("triple bonds!")
    def _button_assign_formal_charge(self, b):
        print("charges!")
    def _button_test_load(self, b):
        print("test load!")
    def _button_order_substructure_smarts(self, b):
        print("reorder!")
    def _button_delete_substructure_smarts(self, b):
        print("delete!")

if __name__ == "__main__":
    file = "openff_polymer_testing/polymer_examples/rdkit_simple_polymers/naturalrubber.pdb"
    p = ChemistryEngine(file)
    # monomer_ids = [102, 103, 105, 96, 92, 93, 95, 104, 94, 97, 99, 101, 98, 100]
    # smarts, block = p.generate_smarts_entry(monomer_ids)
    block = p.get_file_block('pdb')
    print(block)

    # p.visualize(mode="py3Dmol")
# what we want:
# user opens a polymer sdf or pdb and loads it into an rdkit molecule
# user selects a monomer and "saves" it to the chemistry engine
# with a required string input that identifies the monomer
# once saved, saved monomers can be "reordered" with some function
# that does that through a list of monomer strings or a list of integers
# another reordering method would be to have a bunch of drop-down menus 

# allow a test load of the molecule that reloads the view with all unselected
# atoms identified with a red circle (which the user can get rid of if they want
# with a second click of the test load button)

# steps:
# (1) user calls ChemistryEngine with a file which creates an rdkit molecule 
# and all class info which will store monomer information, 
# -> any call to ChemistryEngine is made using the AtomMapNum of the molecule
#   -> atom maping between 3dmol.js and rdkit is handled by PolymerVisualizer3D
# (2) user calls PolymerVisualizer3D which creates the view object and displays the 
# widget 
# -> PolymerVisualizer acts as the visual intermediate between the user and rdkit molecule 

# how to map between 3Dmol.js and PolymerVisualizer -> probably use residue info like now
# residue# + ; + atom_name. ChemistryEngine will return a dictionary that maps between this "key"
# and the atommapnumber, which can later be used to call ChemistryEngine 

