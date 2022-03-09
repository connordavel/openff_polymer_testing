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

# -----------------------------------------------------------------
# trying to get a molecule to recognize one of its substructures
# -----------------------------------------------------------------

# mol = Molecule.from_smiles('CC(=O)N[C@H](C)C(=O)N[C@H](Cc1ccc(OP([O-])([O-])(=O))cc1)C(=O)NC')
# mol.generate_conformers()
# confs = mol._conformers
# mol._conformers = [mol._conformers[0]]


# selected_atoms = get_atoms_between(mol, [19, 22, 8], 19)
# print(sorted(selected_atoms))
# smarts, rdmol = get_smarts_substructure_from_rdkit_idx(mol, ids=selected_atoms)
# print(smarts)


# off_matches = mol.chemical_environment_matches(smarts)
# print(off_matches)

# qmol = Chem.MolFromSmarts(smarts)
# basic_matches = rdmol.GetSubstructMatches(qmol, useChirality=True)
# print(basic_matches)


# --------------------------------------------------------------
# trying to get a substructure to recognize itself v.2
# --------------------------------------------------------------
# smarts = "[N:1]([C@:2]([C:3](=[O:4])[O:9][H:17])([C:5]([C:6](=[O:7])[N:8]([H:15])[H:16])([H:13])[H:14])[H:12])([H:10])[H:11]"
# print(smarts)
# rdmol = Chem.MolFromSmarts(smarts)
# mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
# matches = mol.chemical_environment_matches(smarts)
# print(matches)

# -------------------------------------------------------------------------
# trying to recognize a simple polymer from the substructure dictionary
# -------------------------------------------------------------------------

poly = Molecule.from_smiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
print("matches from chemical_envirnoment_matches:")
matches = poly.chemical_environment_matches('[#6:1]-[#6:2]')
print(matches)
print("results of perceive_residues")
print(poly.perceive_residues())
# see if it worked
[print(a.metadata) for a in poly.atoms if a.atomic_number == 6]


