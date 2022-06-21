from multiprocessing import connection
from ntpath import join
import sys
from openff.toolkit.utils.toolkit_registry import ToolkitRegistry
from pyparsing import ParseSyntaxException
# from pyrsistent import T
from rdkit import Chem
from pathlib import Path
from random import randint
from copy import deepcopy
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.utils import OpenEyeToolkitWrapper
# from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem
from collections import defaultdict, OrderedDict
from collections.abc import Iterable
import py3Dmol
from IPython.utils import io
import ipywidgets as widgets
from IPython.display import Javascript, display, clear_output
import time
import os
import tempfile
import json
import itertools

import networkx as nx
from networkx.algorithms import isomorphism
from openmm import unit as openmm_unit
from openmm.app import PDBFile
from rdkit import Chem

class ChemistryEngine:
    instances = {}
    def __init__(self, file, name=""):
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
        if name == "":
            self.name = path.stem
        else:
            self.name = name
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

    def generate_monomer_smarts(self, monomer_ids, monomer_context_ids=[], assign_map_ids_to=[], additional_specs={}):
        # if assign_map_ids_to is empty, assume map ids to monomer 
        if len(monomer_ids) == 0:
            print("must specify monomer ids to generate library entry, monomer smarts, or library charges")
            return
        # assume either Molecule or RDKit molecule
        if 'captured' in additional_specs.keys():
            assign_map_ids_to = additional_specs['captured']

        # if len(assign_map_ids_to) == 0:
        #     assign_map_ids_to = monomer_ids
        if len(monomer_context_ids) == 0:
            smarts_ids = monomer_ids
        else:
            smarts_ids = monomer_context_ids
        if self.sampled_molecule != None:
            rdmol_copy = deepcopy(self.sampled_molecule)
        else:
            rdmol_copy = deepcopy(self.full_molecule)
        # modify the molecule:
        doubles = additional_specs.get('double', [])
        sorted_doubles = [tuple(sorted(t)) for t in doubles]

        triples = additional_specs.get('triple', [])
        sorted_triples = [tuple(sorted(t)) for t in triples]

        for bond in rdmol_copy.GetBonds():
            bond_ids = tuple(sorted([bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()]))
            if bond_ids in sorted_doubles:
                bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
            elif bond_ids in sorted_triples:
                bond.SetBondType(Chem.rdchem.BondType.TRIPLE)

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
                    "#16": "S",
                    "#17": "Cl"
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
        return block, elements_list

    def generate_smarts_entry(self, monomer_ids, additional_specs={}):
        smarts,_ = self.generate_monomer_smarts(monomer_ids, additional_specs=additional_specs)
        smarts_block, elements_list = self._smarts_string_cleaner(smarts)
        return smarts, smarts_block, elements_list

    def test_polymer_load(self, substructure_lib = "", verbose=True):
        # executes from_pdb using the given biopolymer substructure library and records how many
        # atoms and bonds have chemical information assigned.
        # returns: original molecule ids that are assigned and have assigned bonds 
        self.assigned_atoms = set()
        self.unassigned_atoms = set()
        self.assigned_bonds = set()
        self.unassigned_bonds = set()

        mol = Molecule.from_pdb(self.file, substructure_lib)
        chemical_info = {"double": [], "triple": []}
        for atom in mol.atoms:
            if atom.metadata['already_matched']:
                self.assigned_atoms.add(atom.molecule_atom_index)
            else:
                self.unassigned_atoms.add(atom.molecule_atom_index)
        for bond in mol.bonds:
            # check for assigned bonds 
            if bond.bond_order in [1,2,3]:
                self.assigned_bonds.add((bond.atom1_index, bond.atom2_index))
                if bond.bond_order == 2:
                    chemical_info["double"].append((bond.atom1_index, bond.atom2_index))
                if bond.bond_order == 3:
                    chemical_info["triple"].append((bond.atom1_index, bond.atom2_index))
            else:
                self.unassigned_bonds.add((bond.atom1_index, bond.atom2_index))
                
        # print info on the polymer loading
        if verbose:
            print(f"number of atoms assigned: {len(self.assigned_atoms)}")
            print(f"number of bonds assigned: {len(self.assigned_bonds)}")
            print(f"number of atoms not assigned: {len(self.unassigned_atoms)}")
            print(f"number of bonds not assigned: {len(self.unassigned_bonds)}")
        return self.assigned_atoms, self.assigned_bonds, self.unassigned_atoms, self.unassigned_bonds, chemical_info

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
        doubles = additional_specs.get('double', [])
        sorted_doubles = [tuple(sorted(t)) for t in doubles]

        triples = additional_specs.get('triple', [])
        sorted_triples = [tuple(sorted(t)) for t in triples]

        for bond in mol_copy.GetBonds():
            bond_ids = tuple(sorted([bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()]))
            if bond_ids in sorted_doubles:
                bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
            elif bond_ids in sorted_triples:
                bond.SetBondType(Chem.rdchem.BondType.TRIPLE)
            
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

class SubstructureEngine:
    # a standalone class that deals with inputing and outputing user defined monomers 
    # either in sdf or smarts format
    def __init__(self):
        # a monomer is composed of a name, an rdmol, and a dictionary of metadata
        self.substructures = OrderedDict()
        self.bond_substructures = OrderedDict()

    #input
    def add_substructures_as_sdf(self, name, sdf_file):
        names_listed = isinstance(name, list)
        if names_listed:
            check = [isinstance(a, str) for a in name]
            if not all(check):
                print("all names must be of type <str>")
        name_string = isinstance(name, str)
        if not (names_listed or name_string):
            print("no valid input for naming")
            return
        # read the file
        rd_reader = Chem.rdmolfiles.SDMolSupplier(sdf_file, removeHs=False)
        i = 0
        substructure_list = []
        for rdmol in rd_reader:
            # first, get name
            if names_listed:
                try:
                    n = name[i]
                except IndexError:
                    print("not enough names given")
                    return
            elif name_string:
                if i == 0:
                    n = name
                else:
                    n = name + str(i)
            if n in self.substructures.keys():
                print("duplicate name, continuing")
                continue

            # for each, read the molecule and get custom metadata of where to cut molecule
            has_exclude = rdmol.HasProp("exclude")
            has_include = rdmol.HasProp("include")
            ids_to_include = []
            if has_exclude and has_include:
                print("error in molecule reading: can specify atoms in exclude or include, but not both")
                return
            elif has_exclude or has_include:
                if has_exclude:
                    cuts = rdmol.GetProp("exclude")
                elif has_include:
                    cuts = rdmol.GetProp("include")
                cut_objs = [eval(a) for a in cuts.split(",") if a != ""]
                cut_ints = []
                for obj in cut_objs:
                    if isinstance(obj, Iterable):
                        ints = [int(a) for a in obj]
                        cut_ints = cut_ints + ints
                    else:
                        cut_ints.append(int(obj))
                if has_exclude:
                    # must get all ids that are not in cut_ints:
                    for atom in rdmol.GetAtoms():
                        idx = atom.GetIdx()
                        if idx not in cut_ints:
                            ids_to_include.append(idx)
                elif has_include:
                    ids_to_include = cut_ints
            else:
                for atom in rdmol.GetAtoms():
                    ids_to_include.append(atom.GetIdx())
            # now store info
            substructure = Substructure(n)
            substructure.atoms = ids_to_include
            substructure.rdmol = rdmol
            substructure.smarts = Chem.rdmolfiles.MolToSmarts(rdmol)
            substructure_list.append(substructure)
            i += 1
        for substructure in substructure_list:
            self.substructures[substructure.name] = substructure
        return
    def add_substructures_as_smarts(self, name, smarts, include=[], exclude=[], bond=False):
        has_exclude = bool(len(exclude))
        has_include = bool(len(include))
        ids_to_include = []
        if has_exclude and has_include:
            print("error in molecule reading: can specify atoms in exclude or include, but not both")
            return
        rdmol = Chem.rdmolfiles.MolFromSmarts(smarts)
        if has_exclude:
            # must get all ids that are not in cut_ints:
            for atom in rdmol.GetAtoms():
                idx = atom.GetIdx()
                if idx not in exclude:
                    ids_to_include.append(idx)
        elif has_include:
            ids_to_include = include
        else:
            for atom in rdmol.GetAtoms():
                ids_to_include.append(atom.GetIdx())
        substructure = Substructure(name)
        substructure.atoms = ids_to_include
        substructure.rdmol = rdmol
        substructure.smarts = smarts
        if bond:
            self.bond_substructures[name] = substructure
        else:
            self.substructures[name] = substructure
        return
    def substructure_name_exists(self, name):
        for substructure in self.substructures.values():
            if name == substructure.name:
                return True
        return False
    #output
    def get_substructure_smarts_fragment(self, name):
        substructure = self.substructures[name]
        smarts = Chem.MolFragmentToSmarts(substructure.rdmol, atomsToUse=substructure.atoms)
        return smarts
    def output_substructures_json(self, file_name):
        json_dict = dict()
        for name, substructure in [*self.substructures.items(), *self.bond_substructures.items()]:
            # smarts, elements_list
            # smarts = self.get_substructure_smarts_fragment(name)
            smarts = substructure.smarts
            elements_list = []
            symbol_nums = defaultdict(int)
            for atom in substructure.rdmol.GetAtoms():
                if atom.GetAtomMapNum() > 0:
                    symbol = atom.GetSymbol()
                    symbol_nums[symbol] += 1
                    symbol_num = symbol_nums[symbol]
                    symbol = symbol + str(symbol_num)
                    elements_list.append(symbol)

            if name in json_dict.keys():
                entry = json_dict[substructure.name]
                if smarts in entry.keys():
                    print("duplicate entries")
                else:
                    entry[smarts] = elements_list
                json_dict[substructure.name] = entry
            else:
                json_dict[substructure.name] = {smarts: elements_list}

        with open(file_name, "w") as file:
            json.dump(json_dict, file, indent=4)

class NetworkxMonomerEngine(SubstructureEngine):
    def __init__(self):
        import networkx as nx
        from networkx.algorithms import isomorphism
        from openmm import unit as openmm_unit
        from openmm.app import PDBFile
        from rdkit import Chem

        self.substructures = defaultdict()
        self.bond_substructures = OrderedDict()
        self.monomers = defaultdict()

    def get_substructures_from_pdb(self, name, smarts, pdb_file, caps):
        # TODO: functioanlity for multiple pdb files

        # finds all substructures that match with a given monomer includes 
        # specific noncaptured atoms for each substructure
        rdmol = Chem.rdmolfiles.MolFromSmarts(smarts)
        # remove the caps:
        atoms_to_use = []
        for atom in rdmol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in caps:
                atoms_to_use.append(idx)
        monomer_smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse = atoms_to_use)
        monomer_unit = Chem.rdmolfiles.MolFromSmarts(monomer_smarts)
        polymer = Chem.rdmolfiles.MolFromPDBFile(pdb_file, removeHs=False)
        substructures = defaultdict()
        
        is_isomorphic, isomorphisms = self._get_isomorphisms(monomer_unit, pdb_file)
        captured_isomorphisms = list(isomorphisms)
        if not is_isomorphic:
            print("monomer not found in pdb, no substructures generated")
            return
        monomer = Monomer(name)
        monomer.atoms = atoms_to_use
        monomer.caps = caps
        monomer.rdmol = rdmol # original monomer with caps attached, not monomer_unit
        monomer.isomorphisms = captured_isomorphisms

        self.monomers[name] = monomer
        mapped_atoms_dict, mapped_atoms_set = self._get_existing_isomorphic_atoms()
        unique_substructure_count = 0
        unique_bond_count = 0
        for captured_mapping in captured_isomorphisms:
            open_atoms = dict() #atoms with an open connection to the rest of the polymer
            substructure_ids = set(captured_mapping.keys()) # start with the initial isomorphism
            captured_ids = set(captured_mapping.keys())
            original_ids = set(captured_mapping.keys())
            for joint_idx in captured_mapping.keys():
                atom = polymer.GetAtomWithIdx(joint_idx)
                for neighbor in atom.GetNeighbors():
                    idx = neighbor.GetIdx()
                    if idx not in captured_mapping.keys():
                        open_atoms[idx] = joint_idx
            for open_atom, joint_atom in open_atoms.items():
                # extend from the end of the isomorphism until the chain ends or another isomorphism is reached. 
                # use a BFS implementation with (for now) 3 conditions for when to exit:
                # 1) the end of the graph is reached (usually just a hydrogen) in a reasonable number of nodes (7)
                # 2) an isomorphism is adjacent, in which case find the entirety of the neighboring monomer
                queue = [open_atom]
                visited_atoms = set()
                searching = True
                isomorphism_adjacent = False
                end_reached = False
                while searching:
                    current_id = queue[0]
                    current_atom = polymer.GetAtomWithIdx(current_id)
                    for neighbor in current_atom.GetNeighbors():
                        idx = neighbor.GetIdx()
                        if idx not in visited_atoms and idx not in substructure_ids:
                            queue.append(current_id)
                    queue.pop(0)
                    visited_atoms.add(current_id)
                    # now test for if the search can stop
                    # condition 1: the end of the graph is more than 7 atoms away:
                    if len(visited_atoms) > 7:
                        searching = False
                    # condition 2: the end of the graph is reached
                    elif len(queue) == 0:
                        searching = False
                        end_reached = True
                    # condition 3: an isomorphism is adjacenet 
                    elif current_id in mapped_atoms_set:
                        searching = False
                        isomorphism_adjacent = True
                if isomorphism_adjacent:
                    # find which monomer unit is adjacent and use that monomer unit as the non-captured atoms
                    # additionally, a unique "bond_substructure" must be created that represents bonds between monomers 
                    most_unique_adjacency = set()
                    union_length = 999999
                    for monomer in self.monomers.values():
                        isomorphisms = monomer.isomorphisms
                        for mapping in isomorphisms:
                            ids = set(mapping.keys())
                            if open_atom in ids:
                                if len(ids.union(original_ids)) < union_length:
                                    most_unique_adjacency = ids
                    if most_unique_adjacency == None:
                        print("error 3")
                    # for now, assume the adjacent monomer is directly adjacent
                    substructure_ids = substructure_ids | most_unique_adjacency

                    # finally, create a bond substructure between the most unique adjacency and original_ids
                    # TODO: how to not have to assume that these are single bonds!
                    polymer_copy = deepcopy(polymer)
                    # unspecified bonds except for the inter-monomer bond
                    for atom in polymer_copy.GetAtoms():
                        atom.SetAtomMapNum(0)
                    bond = polymer_copy.GetBondBetweenAtoms(open_atom, joint_atom)
                    bond.GetBeginAtom().SetAtomMapNum(1)
                    bond.GetEndAtom().SetAtomMapNum(1)

                    bond_smarts = Chem.MolFragmentToSmarts(polymer_copy, atomsToUse = (most_unique_adjacency|original_ids))
                    bond_smarts = bond_smarts.replace(".", "~")
                    query_mol = Chem.rdmolfiles.MolFromSmarts(bond_smarts)
                    is_new, isomorphic_bond = self._substructure_is_new(query_mol, bond=True)
                    # check for uniqueness. replace smaller bonds with larger bonds if isomorphic
                    if is_new:
                        unique_bond_count += 1
                        # is_new, isomorphic_bond = self._substructure_is_new(query_mol)
                        self.add_substructures_as_smarts(f"{name}_bond_{unique_bond_count}", bond_smarts, bond=True)
                    else:
                        old_bond = isomorphic_bond.rdmol
                        if query_mol.GetNumAtoms() > old_bond.GetNumAtoms():
                            self.bond_substructures.pop(isomorphic_substructure.name)
                            self.add_substructures_as_smarts(f"{name}_bond_{unique_bond_count}", bond_smarts, bond=True)

                elif end_reached:
                    substructure_ids = substructure_ids | set(visited_atoms)
                    captured_ids = captured_ids | set(visited_atoms)
                else: # either the end of the molecule or no isomorphisms in site (for now, this is bad)
                    substructure_ids = substructure_ids | set(visited_atoms) # ids of captured + non-captured atoms
                    # captured_ids defined above 
            # now substructure_ids and captured_ids should be updated to include larger context

            polymer_copy = deepcopy(polymer)
            n = 1
            # assign chemical information and atom map nums
            for atom in polymer_copy.GetAtoms():
                if atom.GetIdx() in captured_ids:
                    atom.SetAtomMapNum(n)
                    n += 1
                    if atom.GetIdx() in captured_mapping.keys():
                        substructure_atom = monomer_unit.GetAtomWithIdx(captured_mapping[atom.GetIdx()])
                        atom.SetFormalCharge(substructure_atom.GetFormalCharge())
            for bond in polymer_copy.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if a1 in captured_mapping.keys() and a2 in captured_mapping.keys():
                    bond.SetBondType(monomer_unit.GetBondBetweenAtoms(captured_mapping[a1], captured_mapping[a2]).GetBondType())
                else:
                    bond.SetBondType(Chem.rdchem.BondType.UNSPECIFIED)
            substructure_smarts = Chem.MolFragmentToSmarts(polymer_copy, atomsToUse = substructure_ids)
            query_mol = Chem.rdmolfiles.MolFromSmarts(substructure_smarts)
            is_new, isomorphic_substructure = self._substructure_is_new(query_mol)
            # check for uniqueness. replace smaller substructures with larger substructures if isomorphic
            if is_new:
                unique_substructure_count += 1
                # is_new, isomorphic_substructure = self._substructure_is_new(query_mol)
                self.add_substructures_as_smarts(f"{name}_{unique_substructure_count}", substructure_smarts)
            else:
                old_substructure = isomorphic_substructure.rdmol
                if query_mol.GetNumAtoms() > old_substructure.GetNumAtoms():
                    self.substructures.pop(isomorphic_substructure.name)
                    unique_substructure_count += 1
                    self.add_substructures_as_smarts(f"{name}_{unique_substructure_count}", substructure_smarts)
        return
    
    def get_substructures_from_connection_dict(self, name, smarts, connection_dict, caps, custom_end_groups=defaultdict(set), inter_monomer_bond_type=Chem.BondType.SINGLE):
        # test case #1: {1:2}, rubber
        # finds substructs and places caps into sets identified by what map number the cap is attached to
        rdmol = Chem.rdmolfiles.MolFromSmarts(smarts)
        terminal_groups = custom_end_groups
        atoms_to_use = []
        for atom in rdmol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in caps:
                atoms_to_use.append(idx)
        # now find capping groups
        
        for atom in rdmol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            idx = atom.GetIdx()
            if map_num > 0:
                # start searching for the end of the monomer in the direction of the capping atoms
                terminal_group_ids = set()
                first_cap = None
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in caps:
                        first_cap = neighbor.GetIdx()
                        break
                searching_list = [first_cap]
                terminal_group_ids = set([first_cap])
                while len(searching_list) != 0:
                    current_id = searching_list.pop(0)
                    if current_id in caps:
                        terminal_group_ids.add(current_id)
                    current_atom = rdmol.GetAtomWithIdx(current_id)
                    for neighbor in current_atom.GetNeighbors():
                        neighbor_id = neighbor.GetIdx()
                        if neighbor_id not in terminal_group_ids and neighbor_id in caps:
                            searching_list.append(neighbor_id)
                terminal_groups[map_num] = terminal_group_ids
        # now to assemble all combinations of monomers and end_groups
        # this combines information from the connection_dict and terminal_groups 
        # dictionaries. 
        side_groups = defaultdict(list)
        for head, tails in connection_dict.items():
            # tails can be a list for multiple tail groups
            if isinstance(tails, list):
                for tail in tails:
                    side_groups[head].append(tail)
                    if head in terminal_groups.keys():
                        side_groups[head].append(terminal_groups[head])
                    side_groups[tail].append(head)
                    if tail in terminal_groups.keys():
                        side_groups[tail].append(terminal_groups[tail])
            else:
                tail = tails
                side_groups[head].append(tail)
                if head in terminal_groups.keys():
                    side_groups[head].append(terminal_groups[head])
                side_groups[tail].append(head)
                if tail in terminal_groups.keys():
                    side_groups[tail].append(terminal_groups[tail])
        # make dict into a list
        groups_list = []
        for head, tails in side_groups.items():
            sub_list = []
            for tail in tails:
                sub_list.append(tuple([head, tail]))
            groups_list.append(sub_list)
        # now iterate through all multiples of the side_groups to create substructures
        unique_substructure_count = 1
        for iters in itertools.product(*groups_list):
            # for each side group, build the substructure
            monomer_smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse = atoms_to_use)
            sub_rdmol = Chem.rdmolfiles.MolFromSmarts(monomer_smarts)
            # identify the middle atoms for capturing atoms later:
            for atom in sub_rdmol.GetAtoms():
                atom.SetBoolProp("captured", True)
            for side_group in iters:
                head, tail = side_group
                if isinstance(tail, int):
                    middle_monomer_unit = True
                else:
                    middle_monomer_unit = False
                # add the tail to the mutable rdmol
                if middle_monomer_unit:
                    middle_rdmol = Chem.rdmolfiles.MolFromSmarts(monomer_smarts) # the neighboring rdmol
                    for atom in middle_rdmol.GetAtoms():
                        atom.SetAtomMapNum(-1 * atom.GetAtomMapNum()) # neg to differentiate head from tail
                        atom.SetBoolProp("captured", False)
                    # connect middle_rdmol to mutable_rdmol at the connection point map nums
                    # find start and end point for bond:
                    sub_rdmol = Chem.CombineMols(sub_rdmol, middle_rdmol)
                    # add bond based on connection dict
                    for atom in sub_rdmol.GetAtoms():
                        map_num = atom.GetAtomMapNum()
                        if map_num == head:
                            start_id = atom.GetIdx()
                        if map_num == -tail:
                            end_id = atom.GetIdx()
                    editable_mol = Chem.EditableMol(sub_rdmol)
                    editable_mol.AddBond(start_id, end_id ,order=Chem.rdchem.BondType.SINGLE)
                    sub_rdmol = editable_mol.GetMol()
                    for atom in sub_rdmol.GetAtoms():
                        if atom.GetAtomMapNum() < 0:
                            atom.SetAtomMapNum(0)
                else:
                    # add the end group
                    # first, find which cap is "open", ie. bonded to the atom with >0 AtomMapNum
                    open_cap = 0
                    found = False
                    for cap_id in tail:
                        atom = rdmol.GetAtomWithIdx(cap_id)
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomMapNum() > 0:
                                open_cap = cap_id
                                found = True
                                break
                        if found: # break early
                            break
                    cap_rdmol = deepcopy(rdmol)
                    for atom in cap_rdmol.GetAtoms():
                        if atom.GetIdx() == open_cap:
                            atom.SetAtomMapNum(1)
                        else:
                            atom.SetAtomMapNum(0)
                    cap_smarts = Chem.rdmolfiles.MolFragmentToSmarts(cap_rdmol, atomsToUse = list(tail))
                    cap_rdmol = Chem.rdmolfiles.MolFromSmarts(cap_smarts)
                    for atom in cap_rdmol.GetAtoms():
                        atom.SetAtomMapNum(-1 * atom.GetAtomMapNum())
                        atom.SetBoolProp("captured", True)
                    sub_rdmol = Chem.CombineMols(sub_rdmol, cap_rdmol)
                    for atom in sub_rdmol.GetAtoms():
                        map_num = atom.GetAtomMapNum()
                        if map_num == head:
                            start_id = atom.GetIdx()
                        if map_num == -1:
                            end_id = atom.GetIdx()
                    editable_mol = Chem.EditableMol(sub_rdmol)
                    editable_mol.AddBond(start_id, end_id ,order=Chem.rdchem.BondType.SINGLE)
                    sub_rdmol = editable_mol.GetMol()
                    for atom in sub_rdmol.GetAtoms():
                        if atom.GetAtomMapNum() < 0:
                            atom.SetAtomMapNum(0)
            # create substructure entry
            n = 1
            for atom in sub_rdmol.GetAtoms():
                if atom.GetBoolProp("captured"):
                    atom.SetAtomMapNum(n)
                    n += 1
            substructure_smarts = Chem.rdmolfiles.MolToSmarts(sub_rdmol)
            print(substructure_smarts)
            self.add_substructures_as_smarts(f"{name}_{unique_substructure_count}", substructure_smarts)  
            unique_substructure_count += 1
        # iterate through once more to create bond smarts
        bond_entries = []
        for start, end in connection_dict.items():
            bond = {start, end}
            if bond not in bond_entries:
                bond_entries.append(bond)
        unique_bond_count = 1
        for bond in bond_entries:
            start, end = bond
            monomer_smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse = atoms_to_use)
            sub_rdmol1 = Chem.rdmolfiles.MolFromSmarts(monomer_smarts)
            sub_rdmol2 = Chem.rdmolfiles.MolFromSmarts(monomer_smarts)
            for atom in sub_rdmol2.GetAtoms():
                atom.SetAtomMapNum(-1 * atom.GetAtomMapNum())
            sub_rdmol = Chem.CombineMols(sub_rdmol1, sub_rdmol2)
            # add bond based on connection dict
            for atom in sub_rdmol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num == start:
                    start_id = atom.GetIdx()
                if map_num == -end:
                    end_id = atom.GetIdx()
            editable_mol = Chem.EditableMol(sub_rdmol)
            editable_mol.AddBond(start_id, end_id, order=inter_monomer_bond_type)
            sub_rdmol = editable_mol.GetMol()
            for atom in sub_rdmol.GetAtoms():
                if atom.GetIdx() in [start_id, end_id]:
                    atom.SetAtomMapNum(1)
                else:
                    atom.SetAtomMapNum(0)
            substructure_smarts = Chem.rdmolfiles.MolToSmarts(sub_rdmol)
            print(substructure_smarts)
            self.add_substructures_as_smarts(f"{name}_bond_{unique_bond_count}", substructure_smarts, bond=True)
            unique_bond_count += 1

    def build_polymer_from_connection_dict(self, name, smarts, connection_dict, caps, custom_end_groups=defaultdict(set), length=10):
        return

    def _get_existing_isomorphic_atoms(self):
        # returns a list of atom ids that have already been mapped
        mapped_atoms_dict = defaultdict()
        mapped_atoms_list = []
        
        for monomer in self.monomers.values():
            isomorphic_atoms_list = []
            for mapping in monomer.isomorphisms:
                isomorphic_atoms_list = isomorphic_atoms_list + list(mapping.keys())
            isomorphic_atoms = set(isomorphic_atoms_list)
            mapped_atoms_dict[monomer.name] = isomorphic_atoms
            mapped_atoms_list = mapped_atoms_list + list(isomorphic_atoms)
        mapped_atoms_set = set(mapped_atoms_list)

        return mapped_atoms_dict, mapped_atoms_set

    def _get_isomorphisms(self, query, structure):
        # returns an isomorphism map using networkx from query to structure where 
        # both query and structure are rdkit molecules 
        def _rdmol_to_networkx(rdmol):
            _bondtypes = {
                # 0: Chem.BondType.AROMATIC,
                Chem.BondType.SINGLE: 1,
                Chem.BondType.AROMATIC: 1.5,
                Chem.BondType.DOUBLE: 2,
                Chem.BondType.TRIPLE: 3,
                Chem.BondType.QUADRUPLE: 4,
                Chem.BondType.QUINTUPLE: 5,
                Chem.BondType.HEXTUPLE: 6,
            }
            rdmol_G = nx.Graph()
            n_hydrogens = [0] * rdmol.GetNumAtoms()
            for atom in rdmol.GetAtoms():
                atomic_number = atom.GetAtomicNum()
                # Assign sequential negative numbers as atomic numbers for hydrogens attached to the same heavy atom.
                # We do the same to hydrogens in the protein graph. This makes it so we
                # don't have to deal with redundant self-symmetric matches.
                if atomic_number == 1:
                    heavy_atom_idx = atom.GetNeighbors()[0].GetIdx()
                    n_hydrogens[heavy_atom_idx] += 1
                    atomic_number = -1 * n_hydrogens[heavy_atom_idx]
                
                rdmol_G.add_node(
                    atom.GetIdx(),
                    atomic_number=atomic_number,
                    formal_charge=atom.GetFormalCharge(),
                    map_num=atom.GetAtomMapNum()
                )
                # These substructures (and only these substructures) should be able to overlap previous matches.
                # They handle bonds between substructures.
            for bond in rdmol.GetBonds():
                bond_type = bond.GetBondType()

                # All bonds in the graph should have been explicitly assigned by this point.
                if bond_type == Chem.rdchem.BondType.UNSPECIFIED:
                    raise Exception
                    # bond_type = Chem.rdchem.BondType.SINGLE
                    # bond_type = Chem.rdchem.BondType.AROMATIC
                    # bond_type = Chem.rdchem.BondType.ONEANDAHALF
                rdmol_G.add_edge(
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    bond_order=_bondtypes[bond_type],
                )
            return rdmol_G

        def _openmm_topology_to_networkx(openmm_topology):
            """
            Construct an OpenFF Topology object from an OpenMM Topology object.

            Parameters
            ----------
            substructure_library : dict{str:list[str, list[str]]}
                A dictionary of substructures. substructure_library[aa_name] = list[tagged SMARTS, list[atom_names]]
            openmm_topology : openmm.app.Topology
                An OpenMM Topology object

            Returns
            -------
            omm_topology_G : networkx graph
                A networkX graph representation of the openmm topology with chemical information added from the
                substructure dictionary. Atoms are nodes and bonds are edges.
                Nodes (atoms) have attributes for `atomic_number` (int) and `formal_charge` (int).
                Edges (bonds) have attributes for `bond_order` (Chem.rdchem.BondType).
                Any edges that are not assgined a bond order will have the value Chem.rdchem.BondType.UNSPECIFIED
                and should be considered an error.
            """
            import networkx as nx

            omm_topology_G = nx.Graph()
            for atom in openmm_topology.atoms():
                omm_topology_G.add_node(
                    atom.index,
                    atomic_number=atom.element.atomic_number,
                    formal_charge=0.0,
                    atom_name=atom.name,
                    residue_name=atom.residue.name,
                    residue_number=atom.residue.index,
                )

            n_hydrogens = [0] * openmm_topology.getNumAtoms()
            for bond in openmm_topology.bonds():
                omm_topology_G.add_edge(
                    bond.atom1.index,
                    bond.atom2.index,
                    bond_order=Chem.rdchem.BondType.UNSPECIFIED,  # bond.order
                )
                # omm_topology_G.add_edge(
                #     bond.atom1.index,
                #     bond.atom2.index,
                #     bond_order=1,  # quick fix for now
                # )
                # Assign sequential negative numbers as atomic numbers for hydrogens attached to the same heavy atom.
                # We do the same to the substructure templates that are used for matching. This saves runtime because
                # it removes redundant self-symmetric matches.
                if bond.atom1.element.atomic_number == 1:
                    h_index = bond.atom1.index
                    heavy_atom_index = bond.atom2.index
                    n_hydrogens[heavy_atom_index] += 1
                    omm_topology_G.nodes[h_index]["atomic_number"] = (
                        -1 * n_hydrogens[heavy_atom_index]
                    )
                if bond.atom2.element.atomic_number == 1:
                    h_index = bond.atom2.index
                    heavy_atom_index = bond.atom1.index
                    n_hydrogens[heavy_atom_index] += 1
                    omm_topology_G.nodes[h_index]["atomic_number"] = (
                        -1 * n_hydrogens[heavy_atom_index]
                    )

            # Try matching this substructure to the whole molecule graph
            # node_match = isomorphism.categorical_node_match(
            #     ["atomic_number", "already_matched"], [-100, False]
            # )
            return omm_topology_G

        def node_match(data1, data2):
            if data1.get("atomic_number", -100) == data2.get("atomic_number", -100):
                if data1.get("atom_map", 0) == data2.get("atom_map", 0):
                    # return true if atomic numbers match on non_captured atoms
                    return True
                else:
                    # else, captured atoms must have maching values for "already_matched"
                    if data1.get("already_matched", False) == data2.get("already_matched", False):
                        return True
                    else:
                        return False
            else:
                return False

        if isinstance(structure, str):
            if ".pdb" in structure or ".PDB" in structure:
                from openmm.app import PDBFile
                pdb = PDBFile(structure)

                rdmol_G = _rdmol_to_networkx(query)
                omm_topology_G = _openmm_topology_to_networkx(pdb.topology)
                GM = isomorphism.GraphMatcher(
                    omm_topology_G, rdmol_G, node_match=node_match
                )
                return GM.subgraph_is_isomorphic(), GM.subgraph_isomorphisms_iter()
            else:
                return -1, -1
        elif isinstance(structure, Chem.rdchem.Mol):
            rdmol_G = _rdmol_to_networkx(query)
            structure_G = _rdmol_to_networkx(structure)
            GM = isomorphism.GraphMatcher(
                structure_G, rdmol_G, node_match=node_match
            )
            return GM.subgraph_is_isomorphic(), GM.subgraph_isomorphisms_iter()
        else:
            return -1, -1

    def _substructure_is_new(self, query_mol, bond=False):
        if bond:
            substructures = self.bond_substructures
        else:
            substructures = self.substructures

        is_unique = True
        isomorphic_substructure = None
        for sub in substructures.values():
            sub_mol = sub.rdmol
            is_isomorphic_forwards, isomorphisms = self._get_isomorphisms(sub_mol, query_mol)
            is_isomorphic_backwards, isomorphisms = self._get_isomorphisms(query_mol, sub_mol)
            if is_isomorphic_forwards or is_isomorphic_backwards:
                is_unique = False
                isomorphic_substructure = sub
                break
        return is_unique, isomorphic_substructure

class PolymerNetworkRepresentation():
    # stores polymer connectivity information as a Networkx graph
    def __init__(self):
        self.graph = nx.Graph()
        self.length = 0
        self.open_nodes = []
    
    def add_monomer_randomly(self, bond_info = {}):
        # places a new monomer onto a random dangling edge, if possible
        found_connection = False
        if self.length == 0:
            found_connection = True
            self.graph.add_node(
                self.length,
                open_bonds = list(bond_info.keys()),
                terminal_group = False
            )
            self.open_nodes.append(self.length)
            self.length += 1
            
        else:
            # must now find an open edge that is able to accept a bond from the incoming monomer 
            for open_node in self.open_nodes:
                open_bonds = self.graph.nodes[open_node]['open_bonds']
                new_monomer_attachment = -1
                existing_attachment = -1
                for bond_key, attachment_values in bond_info.items():
                    if isinstance(attachment_values, set):
                        # a set means that this bond is an end group
                        continue
                    if isinstance(attachment_values, int):
                        attachment_values = [attachment_values]
                    matches = set(attachment_values) & set(open_bonds)
                    if len(matches) > 0:
                        # choose the first match
                        match = matches.pop()
                        # so, open_bonds contains match, and bond_info looks like: {bond_key: [..., match, ...]}
                        new_monomer_attachment = bond_key
                        existing_attachment = match
                    else:
                        continue
                if new_monomer_attachment == -1:
                    continue
                # first, remove existing_attachment from open_bonds attribute of open_node
                open_bonds.remove(match)
                if len(open_bonds) == 0:
                    self.open_nodes.remove(open_node)
                self.graph.nodes[open_node]['open_bonds'] = open_bonds
                # next, create a new node 
                new_node_open_bonds = list(bond_info.keys()).remove(new_monomer_attachment)
                self.graph.add_node(
                    self.length,
                    open_bonds = new_node_open_bonds,
                    terminal_group = False
                )
                if len(new_node_open_bonds) > 0:
                    self.open_nodes.append(self.length)
                # finally, add an edge between the nodes
                self.graph.add_edge(
                    open_node,
                    self.length,
                    bond_info = {open_node: existing_attachment, self.length: new_monomer_attachment}
                    )
                added_monomer_id = self.length
                self.length += 1

                # if there exist any edge groups, add them here
                for bond_key, attachment_values in bond_info.items():
                    if isinstance(attachment_values, set):
                        self.graph.add_node(
                            self.length,
                            open_bonds = [],
                            terminal_group = True,
                            terminal_group_ids = attachment_values
                        )
                        self.graph.add_edge(
                            added_monomer_id,
                            self.length,
                            bond_info = {added_monomer_id: bond_key, self.length: attachment_values}
                        )
                        self.length += 1
                found_connection = True
                
        return found_connection

    def complete_monomer(self, terminal_group_id_dict):
        # finds dangling edges and turns them into end groups 
        for open_node in self.open_nodes:
            open_bonds = self.graph.nodes[open_node]['open_bonds']
            for open_bond in open_bonds:
                # add a new node and edge
                if open_bond in terminal_group_id_dict.keys():
                    open_bond_terminal_group = terminal_group_id_dict[open_bond]
                    self.graph.add_node(
                        self.length,
                        open_bonds = [],
                        terminal_group = True,
                        terminal_group_ids = open_bond_terminal_group
                    )
                    self.graph.add_edge(
                        open_node,
                        self.length,
                        bond_info = {open_node: open_bond, self.length: open_bond_terminal_group}
                    )
                    self.length += 1
                else:
                    print(f"missing terminal group information for bond number: {open_bond}")
                    continue
        return

    def to_adjacency_matrix(self):
        return
    def from_adjacency_matrix(self, matrix):
        return

class BruteForceParameterEstimator(NetworkxMonomerEngine):
    # calculates the charges and parameters from an entire pdb file loaded with substructures
    # and estimates parameter info for a set of the given substructures
    def __init__(self, pdb_file):
        self.substructures = defaultdict()
        self.bond_substructures = defaultdict()
        self.substructure_properties = defaultdict()
        self.monomers = defaultdict()
        self.pdb_file = pdb_file
        self.monomer_charge_maps = defaultdict()

    def get_elf10_charge_estimates(self):
        # computes elf10 charges for the entire pdb and averages charges to give library charges
        # for each of the substructures
        if len(self.substructures) == 0:
            print("need existing substructures to estimate charges")
            return None, None
        self.output_substructures_json("substructures.json")
        mol = Molecule.from_pdb(self.pdb_file, "substructures.json")
        mol = [mol, *mol.enumerate_stereoisomers(max_isomers=2, toolkit_registry=OpenEyeToolkitWrapper())][-1]
        mol.assign_partial_charges(partial_charge_method='am1bccelf10', toolkit_registry=OpenEyeToolkitWrapper())
        struct_rdmol = mol.to_rdkit()
        # now, group partial charges by isomorphisms
        monomer_charge_maps = dict() # a dict of dicts
        # structure: {<monomer name>: {0:<charge>, 1:<charge>, 2:<charge>,...}, <monomer_name>: {...}, ...}
        for substructure in self.substructures.values():
            sub_rdmol = substructure.rdmol
            is_isomorphic, isomorphisms = self._get_isomorphisms(sub_rdmol, struct_rdmol)
            captured_isomorphisms = list(isomorphisms)
            charge_dict = defaultdict(list)
            for captured_mapping in captured_isomorphisms:
                for struct_idx, sub_idx in captured_mapping.items():
                    charge_dict[sub_idx].append(mol.partial_charges[struct_idx])
            # average charges that map to the same substructure atoms
            averaged_charge_dict = defaultdict(int)
            for id, charges in charge_dict.items():
                avg = sum(charges)/len(charges)
                averaged_charge_dict[id] = avg
            monomer_charge_maps[substructure.name] = averaged_charge_dict
        self.monomer_charge_maps = monomer_charge_maps
        return monomer_charge_maps

    def generate_library_charge_entry(self, name):
        substruct = self.substructures[name]
        charges_dict = dict()
        if len(self.monomer_charge_maps) == 0:
            self.monomer_charge_maps = self.get_elf10_charge_estimates() # charges map to ids, not atom map nums
        custom_charges = self.monomer_charge_maps[name]
        for atom_id in custom_charges.keys():
            atom = substruct.rdmol.GetAtomWithIdx(atom_id)
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                charges_dict[map_num] = custom_charges[atom_id]

        # create a library charge entry 
        entry_string = "<LibraryCharge smirks=\"" + substruct.smarts + "\""
        for map_num, charge in charges_dict.items():
            charge_string = f"{charge:.8g}"
            charge_string = charge_string.replace(" ", " * ")
            entry_string = entry_string + f" charge{map_num}=\"" + charge_string + "\""
        entry_string = entry_string + "></LibraryCharge>"
        return entry_string

class Substructure:
    def __init__(self, name):
        self.name = name
        self.smarts = ""
        self.atoms = []
        self.rdmol = None

class Monomer:
    # utility class to store monomer info
    def __init__(self, name):
        self.name = name
        self.atoms = []
        self.chemical_info = defaultdict(list)
        self.rdmol = None
        self.caps = []
        self.isomorphisms = None

class PolymerVisualizer3D:
    instances = {}
    def __init__(self, chemistry_engine):
        self.name = chemistry_engine.name
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
        name_input_box = widgets.Text(
                                        value='',
                                        placeholder='Monomer Name',
                                        description='',
                                        disabled=True,
                                        layout = widgets.Layout(width='150px')
                                    )
        save_monomer_button = widgets.Button(description="Save Monomer", disabled=True)
        print_button = widgets.Button(description="Print Selection", disabled=True)
        double_bonds_button = widgets.Button(description="Make Double bonds", disabled=True)
        triple_bonds_button = widgets.Button(description="Make Triple bonds", disabled=True)
        formal_charge_button = widgets.Button(description="assign formal charges", disabled=True)
        formal_charge_menu = widgets.Dropdown(
                                                        options=['-3','-2','-1','0','1','2','3'],
                                                        value='0',
                                                        description='',
                                                        disabled=True,
                                                        layout = widgets.Layout(width='50px')
                                                    )
        test_load_button = widgets.Button(description="Test Load", disabled=True)
        order_substructure_button = widgets.Button(description="Reorder (Push to Top)", disabled=True)
        order_tags = widgets.ToggleButtons(
                                        options=[''],
                                        description='',
                                        disabled=True,
                                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                        tooltips=[''],
                                        layout = widgets.Layout(width='100px')
                                    )
        delete_substructure_button = widgets.Button(description="Delete Entry", disabled=True)
        delete_menu = widgets.Dropdown(
                                            options=['Select Monomer Name'],
                                            value='Select Monomer Name',
                                            description='',
                                            disabled=True,
                                        )
        capture_atoms_button = widgets.Button(description="Assign Captured Atoms", disabled=True)

        clear_button.on_click(self._button_clear_selection)
        finalize_button.on_click(self._button_finalize_selection)
        save_monomer_button.on_click(self._button_save_monomer)
        print_button.on_click(self._button_print_selection)
        double_bonds_button.on_click(self._button_assign_double_bonds)
        triple_bonds_button.on_click(self._button_assign_triple_bonds)
        formal_charge_button.on_click(self._button_assign_formal_charge)
        test_load_button.on_click(self._button_test_load)
        order_substructure_button.on_click(self._button_order_substructure_smarts)
        delete_substructure_button.on_click(self._button_delete_substructure_smarts)
        capture_atoms_button.on_click(self._button_capture_atoms)


        # any references to any widgets should only be made to this dict 
        self.widgets = {'clear': clear_button,
                        'finalize': finalize_button,
                        'name_input': name_input_box,
                        'save_monomer': save_monomer_button,
                        'print': print_button,
                        'double': double_bonds_button,
                        'triple': triple_bonds_button,
                        'capture': capture_atoms_button,
                        'charge_button': formal_charge_button,
                        'charge_menu': formal_charge_menu,
                        'test_load': test_load_button,
                        'order_smarts': order_substructure_button,
                        'order_tags': order_tags,
                        'delete_smarts': delete_substructure_button,
                        'delete_menu': delete_menu}

        self.highlights = {} # dict of atoms and their highlighted colors
        self.spheres = {} # dict of which atoms have spheres 
        self.selected_atoms = defaultdict(list)
        self.chemical_info = defaultdict(list)
        self.monomers = OrderedDict()
        self.click_mode = "select_monomer"
        self.valid_click_modes = ['select_monomer', 'do_nothing', 'select_double', 'select_triple', 'select_charged', 'select_captured']

    # the only function the user should ever have to call themselves 
    def show(self):
        self._reset_view()
        return

    def __repr__(self):
        print("in __repr__()")
        self.show()
        return f"file: {self.chemistry_engine.file}"

    def _view_from_block(self, block, format, modify=True):
        view = py3Dmol.view(width=self.width, height=self.height)
        view.addModel(block, format, {"keepH": True})
        if modify:
            self.index_mode = format
        return view

    def get_selected_atoms(self, click_mode):
        data = self.selected_atoms[click_mode]
        if self.index_mode == 'pdb':
            return [(i+1) for i in data]
        else:
            return data

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
                            
                            var serial = atom.serial;
                            Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._store_selected_atom(" + serial + ")");
                    }
                    else {
                            viewer.setStyle({"serial": atom.serial}, {"stick": {"colorscheme": "default"}})
                        
                            var serial = atom.serial;
                            Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._remove_selected_atom(" + serial + ")");
                    }
                    viewer.render();}''' % (self.instance_id, self.instance_id)
        elif new_mode == 'do_nothing':
            selection = {"model": -1}
            code = '''function(atom,viewer,event,container) {
                    void(0)
                    }'''
        elif new_mode in ["select_double", "select_triple", "select_captured"]:
            selection = {"serial": self.get_selected_atoms('select_monomer')}
            # selection = {"model": -1}
            code = '''function(atom,viewer,event,container) {
                        if(!atom.style['clicksphere']) {
                                viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFFFF00}, "clicksphere": {"radius": 0.25}})
            
                                var serial = atom.serial;
                                Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._store_selected_atom(" + serial + ")");
                        }
                        else {
                                viewer.setStyle({"serial": atom.serial}, {"stick": {"colorscheme": "default"}})
     
                                var serial = atom.serial;
                                Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._remove_selected_atom(" + serial + ")");
                        }
                        viewer.render();}''' % (self.instance_id, self.instance_id)
        elif new_mode == "select_charged":
            selection = {"serial": self.get_selected_atoms('select_monomer')}
            return #not yet implemented
        else:
            print("internal error in _set_clickable")
            return

        self.view.setClickable(selection,True,code)
    
    def _store_selected_atom(self, serial):
        if self.index_mode == 'pdb':
            serial = serial - 1
            print(serial)
        elif self.index_mode == 'sdf':
            pass
        print(f"in the store with serial number: {serial}")
        data = self.selected_atoms[self.click_mode]
        if serial not in data:
            self.selected_atoms[self.click_mode].append(serial)
        else:
            print("error: serial number already in self.selected_atoms")
    def _remove_selected_atom(self, serial):
        if self.index_mode == 'pdb':
            serial = serial - 1
        elif self.index_mode == 'sdf':
            pass
        data = self.selected_atoms[self.click_mode]
        if serial in data:
            self.selected_atoms[self.click_mode].remove(serial)
        else:
            print("error: serial not in self.selected_atoms")

    def _atom_id_to_serial(self, id):
        if self.index_mode == 'pdb':
            return id + 1
        else:
            return id
    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #-------------------------------Hoverables-------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _set_hoverable(self, mode):
        if mode not in self.valid_click_modes:
            print("internal: new_mode is not valid")
            return
        func1 = '''
            function(atom,viewer,event,container) {
                    if(!atom.label) {
                        atom.label = viewer.addLabel(atom.elem+":"+ atom.serial,{position: atom, backgroundColor: 'mintcream', fontColor:'black', backgroundOpacity: "0.3"});
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
        # if mode != 'select_monomer':
        #     selection = {'serial': self.get_selected_atoms('select_monomer'), 'invert': True}
        # else:
        #     selection = {'model': -1}
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
        c2 = widgets.VBox((self.widgets['print'], self.widgets['test_load']))
        c3 = widgets.VBox((self.widgets['name_input'], self.widgets['save_monomer']))
        c4 = widgets.VBox((self.widgets['double'], self.widgets['triple'], self.widgets['capture']))
        c5 = widgets.VBox((self.widgets['charge_button'], self.widgets['charge_menu']))

        monomer_selection_buttons = widgets.HBox((c1,c2,c3,c4,c5))

        r1 = widgets.HBox((self.widgets['delete_smarts'], self.widgets['delete_menu']))
        r2 = widgets.HBox((self.widgets['order_smarts'], self.widgets['order_tags']))
        
        self.buttons = widgets.VBox((monomer_selection_buttons, r1, r2))

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
            else:
                widget.disabled = False

    def _reset_buttons(self):
        # first, set all styles back to default
        for id, widget in self.widgets.items():
            if id not in ['clear', 'charge_menu']:
                widget.button_style = ''
        # then, disable buttons that are disabled on startup
        enabled_buttons = ['clear', 'finalize']
        # if there are 1 or more monomers, enable print and delete buttons
        if len(self.monomers) > 0:
            enabled_buttons = enabled_buttons + ['print', 'delete_smarts', 'delete_menu', 'test_load']
        # if there are 2 or more, then reordering functionality becomes useful
            if len(self.monomers) > 1:
                enabled_buttons = enabled_buttons + ['order_smarts', 'order_tags']
        
        self._disable_all_buttons_except(enabled_buttons)
        # save changes to appropriate widgets menus
        names = self.monomers.keys()
        self.widgets['order_tags'].options = list(names)
        self.widgets['delete_menu'].options = list(names) + ['Select Monomer Name']
        self.widgets['delete_menu'].value = 'Select Monomer Name'

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
            for button in ['order_tags', 'order_smarts', 'delete_smarts', 'delete_menu', 'test_load', ]:
                self._disable_button_with_id(button)
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
            # make sure all info on highlights and labels are retained
            self.view.center({"serial": self.get_selected_atoms('select_monomer')})
            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer')}, {"stick": {"colorscheme": "default"}})
            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer'), "invert": True}, {"line": {}})
            self._reload_view()
        elif b.description == "Edit Selection":
            b.description = "Finalize Selection"
            b.button_style = ""
            # reload pdb in case an sdf was used
            pdb_block = self.chemistry_engine.get_file_block('pdb', {})
            self.view = self._view_from_block(pdb_block, 'pdb')
            self._reset_buttons()
            self._set_clickable('select_monomer')
            self._set_hoverable('select_monomer')
            # only keep info on what monomer atoms are selected
            # but discard info on double, triple, formal charges
            monomer_atoms = self.selected_atoms['select_monomer']
            self.selected_atoms = defaultdict(list)
            self.selected_atoms[self.click_mode] = monomer_atoms
            self.chemical_info = defaultdict(list)
            self.view.center({"serial": self.get_selected_atoms('select_monomer')})
            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer')}, {"stick": {"color": 0xFF0000}, "clicksphere": {"radius": 0.25}})
            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer'), "invert": True}, {"stick": {"colorscheme": "default"}})
            self._reload_view()
        else:
            print("internal: no action for _button_finalize_selection!")

    def _button_save_monomer(self, b):
        pass 
        name = self.widgets['name_input'].value
        if name == "":
            print("must specify monomer name")
            return
        
        monomer = Monomer(name)
        atoms = deepcopy(self.selected_atoms['select_monomer'])
        info = deepcopy(self.chemical_info)
        monomer.atoms = atoms
        monomer.chemical_info = info
        self.monomers[name] = monomer
        
        self._button_finalize_selection(self.widgets['finalize']) # same actions as clicking "Edit Selection"

    def _button_print_selection(self, b):
        # if we are in the chemistry assignment mode
        if self.widgets['finalize'].description == 'Edit Selection':
            ids = self.selected_atoms['select_monomer']
            map_ids = self.chemistry_engine._ids_to_map_ids(ids)
            smarts, library_entry, _ = self.chemistry_engine.generate_smarts_entry(map_ids, self.chemical_info)
            print(library_entry)
        else:
            # print all monomers
            json_dict = dict()
            for name, monomer in self.monomers.items():
                map_ids = self.chemistry_engine._ids_to_map_ids(monomer.atoms)
                print(map_ids)
                smarts,_, elements_list = self.chemistry_engine.generate_smarts_entry(map_ids, monomer.chemical_info)
                if monomer.name in json_dict.keys():
                    entry = json_dict[monomer.name]
                    if smarts in entry.keys():
                        print("duplicate entries")
                    else:
                        entry[smarts] = elements_list
                    json_dict[monomer.name] = entry
                else:
                    json_dict[monomer.name] = {smarts: elements_list}
            print(json.dumps(json_dict, indent=4))
        
    def _button_assign_double_bonds(self, b):
        # cycles between "Make Double bonds" and "Commit Double Bond"
        if b.description == "Make Double bonds":
            b.description = "Commit Double Bond"
            b.button_style = "warning"
            self._set_clickable('select_double')
            self._set_hoverable('select_double')
            self.view.center({"serial": self.get_selected_atoms('select_monomer')})
            self._reload_view()
        elif b.description == "Commit Double Bond":
            b.description = "Make Double bonds"
            b.button_style = ""
            
            # load the new double bond into the view with the chemistry engine + sdf
            # check that only two atoms were selected
            bond = self.selected_atoms['select_double']
            self.selected_atoms['select_double'] = []
            if len(bond) == 2:
                new_bond = tuple(bond) # sdf format indexed from 0
                bonds = self.chemical_info['double'] # does default dict handle this? TODO
                bonds.append(new_bond)
                self.chemical_info['double'] = bonds
                sdf_block = self.chemistry_engine.get_file_block('sdf', self.chemical_info)
                self.view = self._view_from_block(sdf_block, 'sdf')
            else:
                print("must select only 2 atoms")
                self.selected_atoms[self.click_mode] = []

            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer')}, {"stick": {"colorscheme": "default"}})
            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer'), "invert": True}, {"line": {}})

            self.view.center({"serial": self.get_selected_atoms('select_monomer')})
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
            self._reload_view()
        else:
            print("internal: no action for _button_assign_double_bonds!")
    def _button_assign_triple_bonds(self, b):
        # cycles between "Make Triple bonds" and "Commit Triple Bond"
        if b.description == "Make Triple bonds":
            b.description = "Commit Triple Bond"
            b.button_style = "warning"
            self._set_clickable('select_triple')
            self._set_hoverable('select_triple')
            self.view.center({"serial": self.get_selected_atoms('select_monomer')})
            self._reload_view()
        elif b.description == "Commit Triple Bond":
            b.description = "Make Triple bonds"
            b.button_style = ""
            
            # load the new triple bond into the view with the chemistry engine + sdf
            # check that only two atoms were selected
            bond = self.selected_atoms['select_triple']
            self.selected_atoms['select_triple'] = []
            if len(bond) == 2:
                new_bond = tuple(bond) # sdf format indexed from 0
                bonds = self.chemical_info['triple'] # does default dict handle this? TODO
                bonds.append(new_bond)
                self.chemical_info['triple'] = bonds
                sdf_block = self.chemistry_engine.get_file_block('sdf', self.chemical_info)
                self.view = self._view_from_block(sdf_block, 'sdf')
            else:
                print("must select only 2 atoms")
                self.selected_atoms[self.click_mode] = []

            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer')}, {"stick": {"colorscheme": "default"}})
            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer'), "invert": True}, {"line": {}})

            self.view.center({"serial": self.get_selected_atoms('select_monomer')})
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
            self._reload_view()
        else:
            print("internal: no action for _button_assign_triple_bonds!")
    def _button_capture_atoms(self, b):
        # cycles between "Make Double bonds" and "Commit Double Bond"
        if b.description == "Assign Captured Atoms":
            b.description = "Commit Captured Atoms"
            b.button_style = "warning"
            self._set_clickable('select_captured')
            self._set_hoverable('select_captured')
            self.view.center({"serial": self.get_selected_atoms('select_monomer')})
            self._reload_view()
        elif b.description == "Commit Captured Atoms":
            b.description = "Assign Captured Atoms"
            b.button_style = ""
            
            # load the new double bond into the view with the chemistry engine + sdf
            # check that only two atoms were selected
            self.chemical_info['captured'] = self.selected_atoms['select_captured']
            self.selected_atoms['select_captured'] = []

            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer')}, {"stick": {"colorscheme": "default"}})
            self.view.setStyle({"serial": self.get_selected_atoms('select_monomer'), "invert": True}, {"line": {}})

            self.view.center({"serial": self.get_selected_atoms('select_monomer')})
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
            self._reload_view()
        else:
            print("internal: no action for _button_capture_atoms!")

    def _button_assign_formal_charge(self, b):
        print("charges!")
    def _button_test_load(self, b):
        json_dict = dict()
        for name, monomer in self.monomers.items():
            map_ids = self.chemistry_engine._ids_to_map_ids(monomer.atoms)
            smarts,_, elements_list = self.chemistry_engine.generate_smarts_entry(map_ids, monomer.chemical_info)
            if monomer.name in json_dict.keys():
                entry = json_dict[monomer.name]
                if smarts in entry.keys():
                    print("duplicate entries")
                else:
                    entry[smarts] = elements_list
                json_dict[monomer.name] = entry
            else:
                json_dict[monomer.name] = {smarts: elements_list}

        with open(f"substructures_{self.name}.json", "w") as file:
            json.dump(json_dict, file, indent=4)
        assigned_atoms, assigned_bonds, unassigned_atoms, unassigned_bonds, chemical_info = self.chemistry_engine.test_polymer_load(f"substructures_{self.name}.json")
        # get new view with new chemical information
        sdf_block = self.chemistry_engine.get_file_block('sdf', chemical_info)
        self.view = self._view_from_block(sdf_block, 'sdf')

        # add spheres to see where there are unassigned atoms
        for atom in unassigned_atoms:
            self.view.addSphere({"center":{"serial": self._atom_id_to_serial(atom)},
                                "radius": 0.75,
                                "color": 'red',
                                "opacity": "0.3"})
        self.view.setStyle({"model": -1}, {"stick": {}})
        self._reload_view()
        print("test load")
        print(chemical_info)
    def _button_order_substructure_smarts(self, b):
        name = self.widgets['order_tags'].value
        options = list(self.widgets['order_tags'].options)
        if name in options and options == list(self.monomers.keys()):
            options.remove(name)
        else:
            print("error in order_substructure_smarts")
        options = [name] + options
        self.widgets['order_tags'].options = options
        self.widgets['delete_menu'].options = options + ['Select Monomer Name']
        self.widgets['order_tags'].value = None
        # reorder Monomers OrderedDict:
        self.monomers.move_to_end(name, last=False)

    def _button_delete_substructure_smarts(self, b):
        name = self.widgets['delete_menu'].value
        if name == 'Select Monomer Name':
            print('first select Monomer name from menu, then delete')
        elif name in self.monomers.keys():
            self.monomers.pop(name)
        else:
            print("internal error in delete_substructure_smarts")
        self._reset_buttons()

    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #---------------------------Visualization Tools--------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def test_load(self, substructures_json, verbose=False):
        # functionally identical to _button_test_load() but instead returns
        # a separate view (not self.view) and reads from a user_defined json
        # file
        assigned_atoms, assigned_bonds, unassigned_atoms, unassigned_bonds, chemical_info = self.chemistry_engine.test_polymer_load(substructures_json)
        # get new view with new chemical information
        sdf_block = self.chemistry_engine.get_file_block('sdf', chemical_info)
        view = self._view_from_block(sdf_block, 'sdf', modify=False)
        original_mode = self.index_mode
        self.index_mode = 'sdf'
        # add spheres to see where there are unassigned atoms
        for atom in unassigned_atoms:
            view.addSphere({"center":{"serial": self._atom_id_to_serial(atom)},
                                "radius": 0.75,
                                "color": 'red',
                                "opacity": "0.3"})
        view.setStyle({"model": -1}, {"stick": {}})
        self.index_mode = original_mode
        if verbose == False:
            clear_output(wait=False)
        display(view)


if __name__ == "__main__":
    # file = "openff_polymer_testing/polymer_examples/rdkit_simple_polymers/naturalrubber.pdb"
    # p = ChemistryEngine(file)
    # monomer_ids = [5, 4, 13, 6, 8, 7, 3, 12, 9, 11, 10, 0]
    # smarts, block = p.generate_smarts_entry(monomer_ids, {'double': [(3, 4)]})
    # # block = p.get_file_block('sdf', {"double": [(0,1)]})
    # print(block)

    # os.chdir("openff_polymer_testing")
    # file = "polymer_examples/rdkit_simple_polymers/naturalrubber.pdb"
    # engine = ChemistryEngine(file, name="order_dependent")
    # viz = PolymerVisualizer3D(engine)
    # viz.test_load("substructures_order_dependent.json", verbose=True)

    # p.visualize(mode="py3Dmol")

    # engine1 = MonomerInputEngine()
    # # test using PVC:
    # engine1.add_monomer_as_smarts('end1', "[#6]-[#6](-[#17])(-[H])-[#6:1](-[#6:2](-[#17:3])(-[H:4])-[H:7])(-[H:5])-[H:6]")
    # engine1.add_monomer_as_smarts('end2', "[#6:1](-[#6](-[H])(-[#6])-[H])(-[H:2])(-[#17:3])-[#6:4](-[H:5])(-[H:6])-[H:7]")
    # engine1.add_monomer_as_smarts('middle', "[#6]-[#6](-[#6:1](-[H:2])(-[H:3])-[#6:4](-[#17:5])(-[H:6])-[#6](-[#6])(-[H])-[H])(-[H])-[#17]")
    # engine1.output_substructures_json("smarts_test.json")
    # # write some sdf files and use those to input monomers
    # # writing using existing rdmols for convenience
    # # with Chem.SDWriter("test_monomers.sdf") as w:
    # #     for name, monomer in engine1.monomers.items():
    # #         rdmol = monomer.rdmol
    # #         w.write(rdmol)
    # # now try reading in sdfile
    # engine2 = MonomerInputEngine()
    # engine2.add_monomers_as_sdf(['end1', 'end2', 'middle'], "test_monomers.sdf")
    # engine2.output_substructures_json("sdf_test.json")

    smarts = "[C:1](-[H])(-[H])(-[H])-[O]-[C:2](-[H])(-[H])(-[H])"
    pdb_file = "openff_polymer_testing/polymer_examples/rdkit_simple_polymers/PEO.pdb"
    substruct_file = "automatic_PEO_substructures.json"
    charges_file_name = "PEO_charges.txt"

    engine = NetworkxMonomerEngine()
    engine.get_substructures_from_connection_dict("PEO", smarts, {1:2}, [1,7])
    engine.output_substructures_json(substruct_file)

# takeaways from 6/9 meeting:
# need to go from "monomer info" to a fully parameterized pdb file with substructures, and the 
# substructures are the main focus here! So, how does a user input information? 
# say the user gave an sdf of the intermediate monomer with the maximum number 
# of connections->I could find where that structure exists which would include recognizing 
# if another monomer exists on either end or if it terminates. 

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

