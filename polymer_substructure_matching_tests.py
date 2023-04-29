from rdkit import Chem
from openmm.app import PDBFile
from networkx.algorithms import isomorphism
import networkx as nx

def _fuzzy_query(query):
    """return a copy of Query which is less specific:
    - ignore aromaticity and hybridization of atoms (i.e. [#6] not C)
    - ignore bond orders
    - ignore formal charges
    """
    from rdkit import Chem

    # it's tricky from the Python API to properly edit queries,
    # but you can do SetQuery on Atoms/Bonds to edit them quite powerfully
    generic = Chem.MolFromSmarts("**")
    generic_bond = generic.GetBondWithIdx(0)
    # N.B. This isn't likely to be an active
    generic_mol = (
        Chem.MolFromSmarts(  # TODO: optimisation, create this once somewhere
            "".join("[#{}]".format(i + 1) for i in range(112))
        )
    )

    fuzzy = Chem.Mol(query)
    for a in fuzzy.GetAtoms():
        a.SetFormalCharge(0)
        a.SetQuery(
            generic_mol.GetAtomWithIdx(a.GetAtomicNum() - 1)
        )  # i.e. H looks up atom 0 in our generic mol
        a.SetNoImplicit(True)
    for b in fuzzy.GetBonds():
        b.SetIsAromatic(False)
        b.SetBondType(Chem.rdchem.BondType.SINGLE)
        b.SetQuery(generic_bond)

    return fuzzy

def _get_connectivity_from_openmm_top(omm_top):

    # convert openmm topology to rdkit Molecule
    # all bonds initially SINGLE, all charge initially neutral
    rwmol = Chem.RWMol()
    for atom in omm_top.atoms():
        idx = rwmol.AddAtom(Chem.Atom(atom.element.atomic_number))
        res = Chem.AtomPDBResidueInfo()
        res.SetResidueName(atom.residue.name)
        res.SetResidueNumber(int(atom.residue.id))
        res.SetChainId(atom.residue.chain.id)
        rwatom = rwmol.GetAtomWithIdx(idx)
        rwatom.SetPDBResidueInfo(res)
    # we're fully explicit
    for atom in rwmol.GetAtoms():
        atom.SetNoImplicit(True)
    for bond in omm_top.bonds():
        rwmol.AddBond(bond[0].index, bond[1].index, Chem.BondType.SINGLE)

    return rwmol

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

def _get_isomorphisms(query, structure):
    # returns an isomorphism map using networkx from query to structure where 
    # both query and structure are rdkit molecules 

    def node_match(data1, data2):
        if data1.get("atomic_number", -100) == 0 or data2.get("atomic_number", -100) == 0:
            return True
        elif data1.get("atomic_number", -100) == data2.get("atomic_number", -100):
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

    # if isinstance(structure, str):
    #     if ".pdb" in structure or ".PDB" in structure:
    #         rdmol_G = _rdmol_to_networkx(query)
    #         omm_topology_G = _pdb_to_networkx(structure)
    #         GM = isomorphism.GraphMatcher(
    #             omm_topology_G, rdmol_G, node_match=node_match
    #         )
    #         return GM.subgraph_is_isomorphic(), GM.subgraph_isomorphisms_iter()
    #     else:
    #         return -1, -1
    if isinstance(structure, Chem.rdchem.Mol):
        rdmol_G = _rdmol_to_networkx(query)
        structure_G = _rdmol_to_networkx(structure)
        GM = isomorphism.GraphMatcher(
            structure_G, rdmol_G, node_match=node_match
        )
        return GM.subgraph_is_isomorphic(), GM.subgraph_isomorphisms_iter()
    elif isinstance(structure, nx.classes.graph.Graph):
        rdmol_G = _rdmol_to_networkx(query)
        GM = isomorphism.GraphMatcher(
            structure, rdmol_G, node_match=node_match
        )
        return GM.subgraph_is_isomorphic(), GM.subgraph_isomorphisms_iter()
    else:
        return -1, -1

substructs = ["[N:1]([C@:2]([C:3](=[O:4])[O:9][H:22])([C@:5]([C:6]([C:8]([H:19])([H:20])[H:21])([H:14])[H:15])([C:7]([H:16])([H:17])[H:18])[H:13])[H:12])([H:10])[H:11]", 
              "[N:1]([C@:2]([C:3](=[O:4])[O:9][H:21])([C@:5]([C:6]([C:8]([H:18])([H:19])[H:20])([H:13])[H:14])([C:7]([H:15])([H:16])[H:17])[H:12])[H:11])[H:10]",
              "[N:1]([C@:2]([C:3](=[O:4])[O-:9])([C@:5]([C:6]([C:8]([H:18])([H:19])[H:20])([H:13])[H:14])([C:7]([H:15])([H:16])[H:17])[H:12])[H:11])[H:10]",
              "[N:1]([C@:2]([C:3]=[O:4])([C@:5]([C:6]([C:8]([H:17])([H:18])[H:19])([H:12])[H:13])([C:7]([H:14])([H:15])[H:16])[H:11])[H:10])[H:9]",
              "[N:1]([C@:2]([C:3](=[O:4])[O-:9])([C@:5]([C:6]([C:8]([H:19])([H:20])[H:21])([H:14])[H:15])([C:7]([H:16])([H:17])[H:18])[H:13])[H:12])([H:10])[H:11]",
              "[N:1]([C@:2]([C:3]=[O:4])([C@:5]([C:6]([C:8]([H:18])([H:19])[H:20])([H:13])[H:14])([C:7]([H:15])([H:16])[H:17])[H:12])[H:11])([H:9])[H:10]",
              "[N:1]([C@:2]([C:3](=[O:4])[O:9][H:20])([C@:5]([C:6]([C:8]([H:17])([H:18])[H:19])([H:12])[H:13])([C:7]([H:14])([H:15])[H:16])[H:11])[H:10])([H:21])[H:22]",
              "[N:1]([C@:2]([C:3](=[O:4])[O-:9])([C@:5]([C:6]([C:8]([H:17])([H:18])[H:19])([H:12])[H:13])([C:7]([H:14])([H:15])[H:16])[H:11])[H:10])([H:20])[H:21]",
              "[N:1]([C@:2]([C:3]=[O:4])([C@:5]([C:6]([C:8]([H:16])([H:17])[H:18])([H:11])[H:12])([C:7]([H:13])([H:14])[H:15])[H:10])[H:9])([H:19])[H:20]",
              "[N+:1]([C@:2]([C:3](=[O:4])[O-:9])([C@:5]([C:6]([C:8]([H:17])([H:18])[H:19])([H:12])[H:13])([C:7]([H:14])([H:15])[H:16])[H:11])[H:10])([H:20])([H:21])[H:22]",
              "[N+:1]([C@:2]([C:3]=[O:4])([C@:5]([C:6]([C:8]([H:16])([H:17])[H:18])([H:11])[H:12])([C:7]([H:13])([H:14])[H:15])[H:10])[H:9])([H:19])([H:20])[H:21]"]

top = PDBFile("openff_polymer_testing/6cww.pdb").topology
rdmol = _get_connectivity_from_openmm_top(top)

for atom in rdmol.GetAtoms():
    atom.SetNoImplicit(True)
for bond in rdmol.GetBonds():
    bond.SetBondType(Chem.BondType.SINGLE)

for smarts in substructs:
    fuzzy = _fuzzy_query(Chem.MolFromSmarts(smarts))
    print(Chem.MolToSmiles(fuzzy))

    matches = sorted([set(match) for match in rdmol.GetSubstructMatches(fuzzy, maxMatches=1000)])
    for match in matches:
        ids = list(match)
        found_smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse=ids, isomericSmarts=True)
        print(found_smarts)
        break
    break

# for smarts in substructs:
#     query = Chem.MolFromSmarts(smarts)
#     is_isomorphic, isomorphisms = _get_isomorphisms(query, rdmol)
#     matches = sorted([set(match.keys()) for match in isomorphisms])
#     for match in matches:
#         ids = list(match)
#         found_smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse=ids, isomericSmarts=True)
#         print(found_smarts)
#         break




