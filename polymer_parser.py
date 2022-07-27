from rdkit import Chem
import json
from collections import defaultdict
from copy import deepcopy

class MonomerEngine:
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
        suppl = Chem.SDMolSupplier(sdf_file)
        if isinstance(name, list): # a list of names
            name_index = 0
            for rdmol in suppl:
                try:
                    name = name[name_index]
                except IndexError:
                    print("not enough names for SDFFile, exiting")
                    return
                smarts = Chem.rdmolfiles.MolToSmarts(rdmol)
                self.add_monomer_as_smarts_fragment(smarts, name, caps, add_caps_from_discarded_ids)
        elif isinstance(name, str): # a string that can be enumerated <str> + <name_index>
            name_index = 1
            for rdmol in suppl:
                name = name + str(name_index) 
                smarts = Chem.rdmolfiles.MolToSmarts(rdmol)
                self.add_monomer_as_smarts_fragment(smarts, name, caps, add_caps_from_discarded_ids)
        else:
            print(f"cannot handle name(s) of type: {type(name)}")

    def read_monomer_info_dict(self, monomer_dict):
        # reads monomer_dict of the same format as that returns from `self.get_monomer_info_dict`
        if "monomers" not in monomer_dict.keys():
            print(f"no \"monomers\" field present in monomer_dict out of the given keys ({monomer_dict.keys()})")
        monomers = monomer_dict['monomers']
        for name, monomer_smarts in monomers:
            self.add_monomer(name, monomer_smarts)
        if "caps" in monomer_dict.keys():
            for name, caps_list in monomer_dict["caps"]:
                for cap in caps_list:
                    self.add_monomer_cap(name, cap)
        # check for miscellaneous keys
        for key in monomer_dict.keys():
            if key not in ["monomers", "caps"]:
                print(f"malformated input: key \"{key}\" is not supported, skipped")
    
    def input_monomer_info_json(self, file_name):
        with open(file_name, "r") as file:
            data = json.load(file_name)
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

class Monomer:
    # utility class to store monomer info
    def __init__(self, name):
        self.name = name
        self.smarts = ""
        self.caps = []

if __name__ == "__main__":
    from openff.toolkit.topology.molecule import Molecule

    # pdb_file = "openff_polymer_testing/polymer_examples/rdkit_simple_polymers/PEG_PLGA_heteropolymer.pdb"
    # json_file = "PEG_PLGA_substructures.json"
    # monomer_info = {"PEG": ("[C](-[H])(-[H])(-[H])-[O]-[C](-[H])(-[H])-[C](-[H])(-[H])-[O](-[H])", [0,1,2,3,11,12]),
    #                 "PLGA1": ("[H]-[O]-[C](=[O])-[C](-[H])(-[C](-[H])(-[H])(-[H]))-[H]",[0,10]),
    #                 "PLGA2": ("[H]-[O]-[C](=[O])-[C](-[H])(-[H])-[H]", [0,7])}
    # engine = MonomerEngine()
    # for name, substructure_and_caps in monomer_info.items():
    #     smarts, caps = substructure_and_caps
    #     engine.add_monomer_as_smarts_fragment(smarts, name, caps)
    # engine.output_monomer_info_json(json_file)

    # mol = Molecule().from_pdb_and_monomer_info(pdb_file, json_file)

    pdb_file = "openff_polymer_testing/polymer_examples/amino_acids/T4-protein.pdb"
    json_file = "T4_protein_substructures.json"
    monomer_info = {
        "ALA": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#1:8])-[#6:9](=[#8:10])-[O]-[H]", []),
        "ARG": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#7:14](-[#1:15])-[#6:16](-[#7:17](-[#1:18])-[#1:19])=[#7&+:20](-[#1:21])-[#1:22])-[#6:23](=[#8:24])-[O]-[H]", []),
        "ASH": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#8:9])-[#8:10]-[#1:11])-[#6:12](-[#8:13])-[O]-[H]", []),
        "ASN": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](=[#8:9])-[#7:10](-[#1:11])-[#1:12])-[#6:13](=[#8:14])-[O]-[H]", []),
        "ASP": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](=[#8:9])-[#8:10])-[#6:11](=[#8:12])-[O]-[H]", []),
        "CYS": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#16:8]-[#1:9])-[#6:10](=[#8:11])-[O]-[H]", []),
        "GLH": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](-[#8:12])-[#8:13]-[#1:14])-[#6:15](-[#8:16])-[O]-[H]", []),
        "GLN": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](=[#8:12])-[#7:13](-[#1:14])-[#1:15])-[#6:16](=[#8:17])-[O]-[H]", []),
        "GLU": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](=[#8:12])-[#8:13])-[#6:14](=[#8:15])-[O]-[H]", []),
        "GLY": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#1:5])-[#6:6](=[#8:7])-[O]-[H]", []),
        "HID": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1-[#7:9](-[#1:10])-[#6:11](-[#1:12])-[#7:13]-[#6:14]-1-[#1:15])-[#6:16](-[#8:17])-[O]-[H]", []),
        "HIE": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1-[#7:9]-[#6:10](-[#1:11])-[#7:12](-[#1:13])-[#6:14]-1-[#1:15])-[#6:16](-[#8:17])-[O]-[H]", []),
        "HIP": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1-[#7:9](-[#1:10])-[#6:11](-[#1:12])-[#7:13](-[#1:14])-[#6:15]-1-[#1:16])-[#6:17](-[#8:18])-[O]-[H]", []),
        "ILE": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6@@:5](-[#1:6])(-[#6:7](-[#1:8])(-[#1:9])-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#6:14](-[#1:15])(-[#1:16])-[#1:17])-[#6:18](=[#8:19])-[O]-[H]", []),
        "LEU": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#6:10](-[#1:11])(-[#1:12])-[#1:13])-[#6:14](-[#1:15])(-[#1:16])-[#1:17])-[#6:18](=[#8:19])-[O]-[H]", []),
        "LYN": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#6:14](-[#1:15])(-[#1:16])-[#7:17](-[#1:18])-[#1:19])-[#6:20](-[#8:21])-[O]-[H]", []),
        "LYS": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#6:14](-[#1:15])(-[#1:16])-[#7&+:17](-[#1:18])(-[#1:19])-[#1:20])-[#6:21](=[#8:22])-[O]-[H]", []),
        "MET": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#16:11]-[#6:12](-[#1:13])(-[#1:14])-[#1:15])-[#6:16](=[#8:17])-[O]-[H]", []),
        "PHE": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1:[#6:9](-[#1:10]):[#6:11](-[#1:12]):[#6:13](-[#1:14]):[#6:15](-[#1:16]):[#6:17]:1-[#1:18])-[#6:19](=[#8:20])-[O]-[H]", []),
        "PRO": ("[H]-[#7:1]1-[#6:2](-[#1:3])(-[#1:4])-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6@@:11]-1(-[#1:12])-[#6:13](=[#8:14])-[O]-[H]", []),
        "SER": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#8:8]-[#1:9])-[#6:10](=[#8:11])-[O]-[H]", []),
        "THR": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6@@:5](-[#1:6])(-[#6:7](-[#1:8])(-[#1:9])-[#1:10])-[#8:11]-[#1:12])-[#6:13](=[#8:14])-[O]-[H]", []),
        "TRP": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1=[#6:9](-[#1:10])-[#7:11](-[#1:12])-[#6:13]2:[#6:14](-[#1:15]):[#6:16](-[#1:17]):[#6:18](-[#1:19]):[#6:20](-[#1:21]):[#6:22]:2-1)-[#6:23](=[#8:24])-[O]-[H]", []),
        "TYR": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1:[#6:9](-[#1:10]):[#6:11](-[#1:12]):[#6:13](-[#8:14]-[#1:15]):[#6:16](-[#1:17]):[#6:18]:1-[#1:19])-[#6:20](=[#8:21])-[O]-[H]", []),
        "VAL": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#6:7](-[#1:8])(-[#1:9])-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#1:14])-[#6:15](=[#8:16])-[O]-[H]", []),
    }
    engine = MonomerEngine()
    for name, substructure_and_caps in monomer_info.items():
        smarts, caps = substructure_and_caps
        engine.add_monomer_as_smarts_fragment(smarts, name, caps)
    engine.output_monomer_info_json(json_file)

    mol = Molecule().from_pdb_and_monomer_info(pdb_file, json_file)

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

    # smarts = "[C:1](-[H])(-[H])(-[H])-[O]-[C:2](-[H])(-[H])(-[H])"

    # engine = NetworkxMonomerEngine()
    # rdmol, graph = engine.build_random_polymer_from_connection_dict("PEO", smarts, {1:2}, [1,7], length=4)

    # smarts = "[C:1](-[H])(-[H])(-[H])-[O]-[C:2](-[H])(-[H])(-[H])"
    # pdb_file = "openff_polymer_testing/polymer_examples/rdkit_simple_polymers/PEO.pdb"
    # substruct_file = "automatic_PEO_substructures.json"
    # charges_file_name = "PEO_charges.txt"

    # engine = NetworkxMonomerEngine()

    # engine.get_substructures_from_connections("PEO", smarts, {1:2}, [1,7])
    # print(polymer)
    # [print(a.GetSymbol()) for a in polymer.GetAtoms()]
    # [print(b.GetBondType()) for b in polymer.GetBonds()]

    # engine = PolymerNetworkRepresentation()
    # for i in range(0, 40):
    #     engine.add_monomer_randomly("PAMAM", bond_info = [(1,2), (2,1), (2,1)])
    # engine.complete_monomer({1: {2}, 2: {9}})
    # graph = engine.graph

    # labels_dict = {}
    # for node in graph.nodes:
    #     if graph.nodes[node]["terminal_group"] == True:
    #         labels_dict[node] = graph.nodes[node]["terminal_group_ids"]
    #     else:
    #         labels_dict[node] = f"{node}, {graph.nodes[node]['open_bonds']}"
    # nx.draw(graph, pos = nx.kamada_kawai_layout(graph), with_labels=True, labels=labels_dict)
    # plt.savefig("bla.png")

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

