from rdkit import Chem
import json
from collections import defaultdict
from copy import deepcopy
import sys, getopt
from substructure_generator import *
from openff.toolkit.topology import Topology, Molecule

if __name__ == "__main__":

    pdb_file = "polymer_examples/compatible_pdbs/simple_polymers/PEO_PLGA.pdb"
    json_file = "PEO_PLGA.json"
    monomer_info = {"PEG": ("[H:1]-[C:2](-[C:3](-[O:4]-[H:5])(-[H:6])-[H:7])(-[H:8])-[H:9]", [0,4]),
                    "PLGA1": ("[H]-[O]-[C](=[O])-[C](-[H])(-[C](-[H])(-[H])(-[H]))-[H]",[0,10]),
                    "PLGA2": ("[H]-[O]-[C](=[O])-[C](-[H])(-[H])-[H]", [0,7])}
    engine = SubstructureGenerator()
    for name, substructure_and_caps in monomer_info.items():
        smarts, caps = substructure_and_caps
        engine.add_monomer_as_smarts_fragment(smarts, name, caps)
    engine.output_monomer_info_json(json_file)

    mol = Topology.from_pdb_and_monomer_info(pdb_file, json_file)

    # pdb_file = "openff_polymer_testing/polymer_examples/amino_acids/T4-protein.pdb"
    # json_file = "T4_protein_substructures.json"
    # monomer_info = {
    #     "ALA": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#1:8])-[#6:9](=[#8:10])-[O]-[H]", []),
    #     "ARG": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#7:14](-[#1:15])-[#6:16](-[#7:17](-[#1:18])-[#1:19])=[#7&+:20](-[#1:21])-[#1:22])-[#6:23](=[#8:24])-[O]-[H]", []),
    #     "ASH": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#8:9])-[#8:10]-[#1:11])-[#6:12](-[#8:13])-[O]-[H]", []),
    #     "ASN": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](=[#8:9])-[#7:10](-[#1:11])-[#1:12])-[#6:13](=[#8:14])-[O]-[H]", []),
    #     "ASP": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](=[#8:9])-[#8:10])-[#6:11](=[#8:12])-[O]-[H]", []),
    #     "CYS": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#16:8]-[#1:9])-[#6:10](=[#8:11])-[O]-[H]", []),
    #     "GLH": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](-[#8:12])-[#8:13]-[#1:14])-[#6:15](-[#8:16])-[O]-[H]", []),
    #     "GLN": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](=[#8:12])-[#7:13](-[#1:14])-[#1:15])-[#6:16](=[#8:17])-[O]-[H]", []),
    #     "GLU": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](=[#8:12])-[#8:13])-[#6:14](=[#8:15])-[O]-[H]", []),
    #     "GLY": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#1:5])-[#6:6](=[#8:7])-[O]-[H]", []),
    #     "HID": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1-[#7:9](-[#1:10])-[#6:11](-[#1:12])-[#7:13]-[#6:14]-1-[#1:15])-[#6:16](-[#8:17])-[O]-[H]", []),
    #     "HIE": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1-[#7:9]-[#6:10](-[#1:11])-[#7:12](-[#1:13])-[#6:14]-1-[#1:15])-[#6:16](-[#8:17])-[O]-[H]", []),
    #     "HIP": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1-[#7:9](-[#1:10])-[#6:11](-[#1:12])-[#7:13](-[#1:14])-[#6:15]-1-[#1:16])-[#6:17](-[#8:18])-[O]-[H]", []),
    #     "ILE": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6@@:5](-[#1:6])(-[#6:7](-[#1:8])(-[#1:9])-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#6:14](-[#1:15])(-[#1:16])-[#1:17])-[#6:18](=[#8:19])-[O]-[H]", []),
    #     "LEU": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#6:10](-[#1:11])(-[#1:12])-[#1:13])-[#6:14](-[#1:15])(-[#1:16])-[#1:17])-[#6:18](=[#8:19])-[O]-[H]", []),
    #     "LYN": ("[H]-[#7:1](-[#1:2])-[#6:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#6:14](-[#1:15])(-[#1:16])-[#7:17](-[#1:18])-[#1:19])-[#6:20](-[#8:21])-[O]-[H]", []),
    #     "LYS": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#6:14](-[#1:15])(-[#1:16])-[#7&+:17](-[#1:18])(-[#1:19])-[#1:20])-[#6:21](=[#8:22])-[O]-[H]", []),
    #     "MET": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#16:11]-[#6:12](-[#1:13])(-[#1:14])-[#1:15])-[#6:16](=[#8:17])-[O]-[H]", []),
    #     "PHE": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1:[#6:9](-[#1:10]):[#6:11](-[#1:12]):[#6:13](-[#1:14]):[#6:15](-[#1:16]):[#6:17]:1-[#1:18])-[#6:19](=[#8:20])-[O]-[H]", []),
    #     "PRO": ("[H]-[#7:1]1-[#6:2](-[#1:3])(-[#1:4])-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#6@@:11]-1(-[#1:12])-[#6:13](=[#8:14])-[O]-[H]", []),
    #     "SER": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#8:8]-[#1:9])-[#6:10](=[#8:11])-[O]-[H]", []),
    #     "THR": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6@@:5](-[#1:6])(-[#6:7](-[#1:8])(-[#1:9])-[#1:10])-[#8:11]-[#1:12])-[#6:13](=[#8:14])-[O]-[H]", []),
    #     "TRP": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1=[#6:9](-[#1:10])-[#7:11](-[#1:12])-[#6:13]2:[#6:14](-[#1:15]):[#6:16](-[#1:17]):[#6:18](-[#1:19]):[#6:20](-[#1:21]):[#6:22]:2-1)-[#6:23](=[#8:24])-[O]-[H]", []),
    #     "TYR": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8]1:[#6:9](-[#1:10]):[#6:11](-[#1:12]):[#6:13](-[#8:14]-[#1:15]):[#6:16](-[#1:17]):[#6:18]:1-[#1:19])-[#6:20](=[#8:21])-[O]-[H]", []),
    #     "VAL": ("[H]-[#7:1](-[#1:2])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#6:7](-[#1:8])(-[#1:9])-[#1:10])-[#6:11](-[#1:12])(-[#1:13])-[#1:14])-[#6:15](=[#8:16])-[O]-[H]", []),
    # }
    # engine = SubstructureGenerator()
    # for name, substructure_and_caps in monomer_info.items():
    #     smarts, caps = substructure_and_caps
    #     engine.add_monomer_as_smarts_fragment(smarts, name, caps)
    # engine.add_monomer("MET_special", "[H]-[#7:1](-[#1:2])(-[#1:0])-[#6@@:3](-[#1:4])(-[#6:5](-[#1:6])(-[#1:7])-[#6:8](-[#1:9])(-[#1:10])-[#16:11]-[#6:12](-[#1:13])(-[#1:14])-[#1:15])-[#6:16](=[#8:17])-*")
    # engine.output_monomer_info_json(json_file)

    # assigned_atoms = set()
    # unassigned_atoms = set()
    # assigned_bonds = set()
    # unassigned_bonds = set()

    # mol,  = Molecule.from_pdb_and_monomer_info(pdb_file, json_file)
    # chemical_info = {"double": [], "triple": []}
    # isomorphisms = defaultdict(list)
    # for atom in mol.atoms:
    #     if atom.metadata['already_matched']:
    #         assigned_atoms.add(atom.molecule_atom_index)
    #         isomorphisms[atom.metadata['residue_name']].append(atom.molecule_atom_index)
    #     else:
    #         unassigned_atoms.add(atom.molecule_atom_index)
    # for bond in mol.bonds:
    #     # check for assigned bonds 
    #     if bond.bond_order in [1,1.5,2,3]:
    #         assigned_bonds.add((bond.atom1_index, bond.atom2_index))
    #         if bond.bond_order == 2:
    #             chemical_info["double"].append((bond.atom1_index, bond.atom2_index))
    #         if bond.bond_order == 3:
    #             chemical_info["triple"].append((bond.atom1_index, bond.atom2_index))
    #     else:
    #         unassigned_bonds.add((bond.atom1_index, bond.atom2_index))
        
    # print(f"number of atoms assigned: {len(assigned_atoms)}")
    # print(f"number of bonds assigned: {len(assigned_bonds)}")
    # print(f"number of atoms not assigned: {len(unassigned_atoms)}")
    # print(f"number of bonds not assigned: {len(unassigned_bonds)}")

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

