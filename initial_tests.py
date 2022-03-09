from openff.toolkit.topology.molecule import FrozenMolecule, Molecule
# from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem

def get_smarts_substructure_from_rdkit_idx(mol, ids):
    rdmol = mol.to_rdkit()
    smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse = ids)
    return smarts, rdmol
def get_atoms_between(mol, wall_ids, seed):
    # a fun little function that gets all atom ids until it hits "walls" specified by the user
    if len(wall_ids) < 2:
        print("must provide more two or more atom ids")
        return -1
    found = [] # atom ids found that we need not go over again
    active = [seed] # atom ids where the algorithm is currently centered 
    n = 0
    rdmol = mol.to_rdkit() #makes sure ids are the rdkit ids
    while len(active) != 0 and n < 2*mol.n_atoms:
        for active_idx in active:
            active_atom = rdmol.GetAtomWithIdx(active_idx)
            for neighbor in active_atom.GetNeighbors():
                idx = neighbor.GetIdx()
                if (idx in found) or (idx in active):
                    continue
                elif (idx in wall_ids):
                    found.append(idx)
                else:
                    active.append(idx)
            active.remove(active_idx)
            found.append(active_idx)
    return found

mol = Molecule.from_pdb('openff_polymer_testing/polymer_examples/rdkit_simple_polymers/naturalrubber.pdb')
# mol.generate_conformers()
# confs = mol._conformers
# mol._conformers = [mol._conformers[0]]

# [print(f"{a.metadata['dict']}") for a in mol.atoms[:20]]
# [print(f"{a.molecule_atom_index} {a.atomic_number}") for a in mol.atoms]

# since not all of the residues are assigned and some negative hydrogen atomic number exist 
assigned_atoms = set()
for atom in mol.atoms:
    if atom.atomic_number < 0:
        atom._atomic_number = -atom.atomic_number
    dictionary = eval(atom.metadata['dict']) # this is very not good
    # basically this is a very 
    try:
        test = dictionary['already_matched']
        assigned_atoms.add(atom.molecule_atom_index)
    except Exception:
        pass
print(assigned_atoms)

# now, molecules are ready to be visualized 


# rdmol = mol.to_rdkit()
# [print(f"{a.GetIdx()} {a.GetAtomicNum()}") for a in rdmol.GetAtoms()]





