from rdkit import Chem
import json
from collections import defaultdict
from copy import deepcopy
import sys, getopt
import itertools

class SubstructureGenerator:
    def __init__(self):
        self.monomers = defaultdict()

    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #                      READERS
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

    def add_monomer_as_smarts_fragment(self, smarts, name, caps=[], add_caps_from_discarded_ids=True):
        rdmol = Chem.rdmolfiles.MolFromSmarts(smarts)
        if caps == []: # attempt to find caps from MapNums if not specified
            for atom in rdmol.GetAtoms():
                if atom.GetAtomMapNum() == 0:
                    caps.append(atom.GetIdx())
        if caps == []:
            return -1
        print(caps)
        # remove caps and store monomer + caps in Monomer instance
        atoms_to_use = [] # opposte of caps list
        atom_map_index = 1
        for atom in rdmol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in caps:
                atoms_to_use.append(idx)
                atom.SetAtomMapNum(atom_map_index)
                atom_map_index += 1
            else: 
                atom.SetAtomMapNum(0)

        # caps are found by finding all connected fragments with ids in "caps"
        cap_fragments = []
        if add_caps_from_discarded_ids:
            for atom in rdmol.GetAtoms():
                id = atom.GetIdx()
                if id not in caps:
                    continue
                # if neighbor in existing fragment, add id to that fragment
                connected_fragments = []
                for neighbor in atom.GetNeighbors():
                    n_id = neighbor.GetIdx()
                    for cap_frag, fragment_id in zip(cap_fragments, range(0, len(cap_fragments))):
                        if n_id in cap_frag:
                            connected_fragments.append(fragment_id)
                
                if len(connected_fragments) > 0:
                    new_fragment = {id}
                    for frag_id in sorted(connected_fragments, reverse=True):
                        frag = cap_fragments.pop(frag_id)
                        new_fragment = new_fragment | frag
                    cap_fragments.append(new_fragment)
                else:
                    cap_fragments.append(set([id]))

        # find where each monomer connects to the cap, and replace with wildtype atoms 
        # but with the same bond order 
        mid_rdmol = deepcopy(rdmol)
        for atom in mid_rdmol.GetAtoms():
            if atom.GetIdx() in caps:
                continue
            # search atom neighbors for connections to capped atoms
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in caps:
                    neighbor.SetAtomicNum(0)
                    neighbor.SetQuery(Chem.AtomFromSmarts("[*]")) # must reset query for wildtype atom to be printed
                    atoms_to_use.append(neighbor.GetIdx())
                
        monomer_smarts = Chem.MolFragmentToSmarts(mid_rdmol, atomsToUse = atoms_to_use)
        self.add_monomer(name, monomer_smarts)
        # finally, for each connected cap, find where the cap connects to the monomer
        for cap_fragment in cap_fragments:
            cap_ids_to_use = list(cap_fragment)
            cap_rdmol = deepcopy(rdmol)
            for atom in cap_rdmol.GetAtoms():
                if atom.GetIdx() not in cap_fragment:
                    continue
                # search atom neighbors for connections to monomer_atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomMapNum() > 0: # if part of main body of monomer
                        neighbor.SetAtomicNum(0)
                        neighbor.SetQuery(Chem.AtomFromSmarts("[*]")) # must reset query for wildtype atom to be printed
                        # neighbor should retain its original atomMapNum
                        cap_ids_to_use.append(neighbor.GetIdx())
            cap_smarts = Chem.MolFragmentToSmarts(cap_rdmol, atomsToUse = cap_ids_to_use)
            self.add_monomer_cap(name, cap_smarts)
            
    def add_monomer_as_sdf_fragment(self, sdf_file, name, caps=[], add_caps_from_discarded_ids=True):
        # identical in function to `self.add_monomer_as_smarts_fragment`, but reads monomer 
        # from an sdf file instead of a smarts. Can handle sdf with multiple molecules, but 
        # must specify a list of names or a single name string that will be enumerated from 1,2,3,...
        suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
        if isinstance(name, list): # a list of names
            name_index = 0
            for rdmol in suppl:
                try:
                    name = name[name_index]
                except IndexError:
                    print("not enough names for SDFFile, exiting")
                    return
                name_index += 1
                smarts = Chem.rdmolfiles.MolToSmarts(rdmol)
                self.add_monomer_as_smarts_fragment(smarts, name, add_caps_from_discarded_ids=add_caps_from_discarded_ids,  caps=[])
        elif isinstance(name, str): # a string that can be enumerated <str> + <name_index>
            name_index = 1
            for rdmol in suppl:
                new_name = name + str(name_index)
                name_index += 1 
                smarts = Chem.rdmolfiles.MolToSmarts(rdmol)
                print(smarts)
                self.add_monomer_as_smarts_fragment(smarts, new_name, caps=[], add_caps_from_discarded_ids=add_caps_from_discarded_ids)
        else:
            print(f"cannot handle name(s) of type: {type(name)}")

    def read_monomer_info_dict(self, monomer_dict):
        # reads monomer_dict of the same format as that returns from `self.get_monomer_info_dict`
        if "monomers" not in monomer_dict.keys():
            print(f"no \"monomers\" field present in monomer_dict out of the given keys ({monomer_dict.keys()})")
        monomers = monomer_dict['monomers']
        for name, monomer_smarts in monomers.items():
            self.add_monomer(name, monomer_smarts)
        if "caps" in monomer_dict.keys():
            for name, caps_list in monomer_dict["caps"].items():
                for cap in caps_list:
                    self.add_monomer_cap(name, cap)
        # check for miscellaneous keys
        for key in monomer_dict.keys():
            if key not in ["monomers", "caps"]:
                print(f"malformated input: key \"{key}\" is not supported, skipped")
    
    def input_monomer_info_json(self, file_name):
        with open(file_name, "r") as file:
            data = json.load(file)
            self.read_monomer_info_dict(data)

    def add_monomer(self, name, smarts):
        # adds monomer information
        # name: any string to specify the samrts (must be unique)
        # smarts: smarts of the monomer unit as it would appear in the middle of 
        #         a polymer. Wildtype atoms ("[*]") represent the beginning of a 
        #         neighboring monomer. MUST CONTAIN EXPLICIT HYDROGENS
        #         ex) [*]-[C:1](-[H])(-[H])-[C:2](-[H])(-[H])-[*]
        #               -> use AtomMapNums ([C:1]) to assign caps later 
        if name in self.monomers.keys():
            print("monomer name already exists")
            return
        monomer = Monomer(name)
        monomer.caps = []
        monomer.smarts = smarts

        self.monomers[name] = monomer
        pass

    def add_monomer_cap(self, name, smarts):
        # adds cap to an already existing monomer 
        # name: name of existing monomer
        # smarts: smarts of the cap to add to the monomer. Wildtype atoms ("[*]")
        #         represent where the cap connects to the monomer, which must contain
        #         the atomMapNum of the monomer atom to connect to. MUST CONTAIN EXPLICIT
        #         HYDROGENS
        #           ex) for [*]-[C:1](-[H])(-[H])-[C:2](-[H])(-[H])-[*], a valid cap would be
        #           represented as [*:1]-[H] (note how [*:1] corresponds to [C:1])
        #           or maybe [*:2]-[C](-[H])(-[H])(-[H])(note how [*:2] corresponds to [C:2])
        if name not in self.monomers.keys():
            print("monomer name does not exist")
            return
        monomer = self.monomers[name]
        # TODO: add some validation here for edge cases
        monomer.caps.append(smarts)
        pass

    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #                        WRITERS
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

    def _enumerate_substructures_with_caps(self, monomer_smarts, monomer_caps, remove_complete_substructures=True):
        # for a named substructure, returns all possible combinations
        # of caps and inter-monomer bonds. 
        # monomer_smarts: smarts of the monomer, with wildtype atoms (required) repressenting 
        #                 intermonomer_bonds and any possible connection to a cap
        # monomer_caps: list of smarts with atom_mapped wildtype atoms representing a connection
        #               to a atom. The AtomMapNum of the wildtype atom must match the AtomMapNum
        #               of the atom represented by the wildtype atom
        if monomer_caps == []:
            return [monomer_smarts]

        rdmol = Chem.MolFromSmarts(monomer_smarts)
        cap_groups = defaultdict(list)
        # the default "cap" is the intermonomer bond represented with wildtype atoms
        for atom in rdmol.GetAtoms():
            if atom.GetAtomicNum() == 0: # found a wildtype atom
                atom_id = atom.GetIdx()
                rdmol_copy = deepcopy(rdmol)
                atom_copy = rdmol_copy.GetAtomWithIdx(atom_id)
                # should only have one neighbor
                neighbor = atom_copy.GetNeighbors()[0] #should only have one neighbor
                atom_map_num = neighbor.GetAtomMapNum()
                if atom_map_num == 0: # if connected to intermonomer bond, else, error
                    print("ill-formated monomer: inter-monomer connection without atom map number")
                    return []
                      
                neighbor.SetAtomicNum(0)
                neighbor.SetQuery(Chem.AtomFromSmarts("[*]"))
                neighbor.SetAtomMapNum(atom_map_num)
                ids_to_include = [atom.GetIdx(), neighbor.GetIdx()]
                cap_smarts = Chem.MolFragmentToSmarts(rdmol_copy, atomsToUse = ids_to_include)
                cap_groups[atom_map_num].append(cap_smarts)
        # also add end-of-polymer caps that are specified
        for cap_smarts in monomer_caps:
            cap_rdmol = Chem.MolFromSmarts(cap_smarts)
            # find with atom on the monomer the cap bonds to
            for atom in cap_rdmol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    cap_groups[atom.GetAtomMapNum()].append(cap_smarts)
                    break
        # enumerate all combinations of the cap_groups
        cap_groups = list(cap_groups.values())
        enumerated_substructures = []
        for iters in itertools.product(*cap_groups): # enumerate all combinations of cap_groups
            # first, remove all wildtype atoms
            nonwildtype_atoms = []
            for atom in rdmol.GetAtoms():
                if atom.GetAtomicNum() > 0:
                    nonwildtype_atoms.append(atom.GetIdx())
            substructure_smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse = nonwildtype_atoms)
            substructure = Chem.MolFromSmarts(substructure_smarts)
            # attach caps to the substructure (except for intermonomder bonds, which are already present)
            for cap_smarts in iters:
                #TODO: better way of doing this
                # attach the cap smarts to the substructure
                cap_rdmol = Chem.MolFromSmarts(cap_smarts)
                # remove attachment point
                attachment_atom_map_num = -1
                cap_start_id = -1
                cap_ids_to_include = []
                inter_monomer_bond_order = None
                for atom in cap_rdmol.GetAtoms():
                    if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0:
                        attachment_atom_map_num = atom.GetAtomMapNum()
                        # should have only one neighbor
                        neighbor = atom.GetNeighbors()[0]
                        cap_start_id = neighbor.GetIdx()
                        bond = cap_rdmol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        inter_monomer_bond_order = bond.GetBondType()
                    else:
                        cap_ids_to_include.append(atom.GetIdx())
                # use atom MapNums to identify attachment points
                for atom in cap_rdmol.GetAtoms():
                    if atom.GetIdx() == cap_start_id:
                        atom.SetAtomMapNum(1)
                    else:
                        atom.SetAtomMapNum(0)
                cap_fragment = Chem.MolFragmentToSmarts(cap_rdmol, atomsToUse=cap_ids_to_include)
                cap_rdmol = Chem.MolFromSmarts(cap_fragment)
                for atom in cap_rdmol.GetAtoms():
                    atom.SetAtomMapNum(-atom.GetAtomMapNum())
                # find where to attach the cap_fragment using atom map num
                new_substructure = Chem.CombineMols(substructure, cap_rdmol)
                # find start and end ids for the bond
                start = -1
                end = -1
                for atom in new_substructure.GetAtoms():
                    if atom.GetAtomMapNum() < 0:
                        start = atom.GetIdx()
                    elif atom.GetAtomMapNum() == attachment_atom_map_num:
                        end = atom.GetIdx()
                editable_mol = Chem.EditableMol(new_substructure)
                editable_mol.AddBond(start, end, order = inter_monomer_bond_order)
                substructure = editable_mol.GetMol()
            if remove_complete_substructures:
                incomplete_substructure = False
                for atom in substructure.GetAtoms():
                    if atom.GetAtomicNum() == 0:
                        incomplete_substructure = True
                if not incomplete_substructure:
                    continue
            # AtomMapNums no longer serve any purpose:
            for atom in substructure.GetAtoms():
                atom.SetAtomMapNum(0)
            enumerated_substructures.append(Chem.MolToSmarts(substructure))
        return enumerated_substructures


    def get_monomer_info_dict(self):
        monomer_dict = {"monomers": defaultdict(str), "caps": defaultdict(list)}
        for name, monomer in self.monomers.items():
            monomer_dict["monomers"][name] = monomer.smarts
            monomer_dict["caps"][name] = monomer.caps

        return monomer_dict 

    def output_monomer_info_json(self, file_name):
        monomer_dict = self.get_monomer_info_dict()
        with open(file_name, "w") as file:
            json.dump(monomer_dict, file, indent=4)


    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #                   ERROR HANDLING
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def check_monomer(rdmol):
        # rdmol: Smarts or Chem.Mol object
        if isinstance(rdmol, str):
            rdmol = Chem.MolToSmarts(rdmol)

        return
    def check_cap(rdmol):
        if isinstance(rdmol, str):
            rdmol = Chem.MolToSmarts(rdmol)
        return
    
class Monomer:
    # utility class to store monomer info
    def __init__(self, name):
        self.name = name
        self.smarts = ""
        self.caps = []

class InvalidRDMol(Exception):
    # raised during validation and error checking to catch invalid inputted/generated monomers
    def __init__(self, message, smarts=""):
        self.message = message
        self.smarts = smarts
        super().__init__(self.message)
    
    def __str__(self):
        if self.smarts == "":
            return f"rdmol found to be invalid: {self.message}"
        else:
            return f"rdmol with smarts \"{self.smarts}\" found to be invalid: {self.message}"

class InvalidMonomer(InvalidRDMol):
    pass
class InvalidCap(InvalidRDMol):
    pass 

if __name__ == "__main__":
    args = sys.argv
    name = None
    file_location = None
    json_file_name = "substrucures.json"
    print_dict = False
    if len(args) > 1:
        opts, args = getopt.getopt(args[1:], "f:n:o:p", ['file=', 'name='])
        for opt, arg in opts:
            if opt in ["-f", "--file"]:
                file_location = arg
            elif opt in ["-n", "--name"]:
                name = arg
            elif opt == "-p":
                print_dict = True
            elif opt == "-o":
                json_file_name = arg
    else:
        print("please specify sdf file location")
        sys.exit()
    print(name)
    print(file_location)
    print(json_file_name)
    engine = SubstructureGenerator()
    engine.add_monomer_as_sdf_fragment(file_location, name)
    engine.output_monomer_info_json(json_file_name)
    if print_dict:
        print(engine.get_monomer_info_dict())
    