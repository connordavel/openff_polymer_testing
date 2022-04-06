from openff.toolkit.topology.molecule import FrozenMolecule, Molecule
# from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem
from create_smarts_substructure import *
import os 

os.chdir("openff_polymer_testing")
file = 'polymer_examples/rdkit_simple_polymers/polyamide_equil.pdb'
mol = Molecule.from_pdb(file)
# rdmol = Chem.rdmolfiles.MolFromPDBFile(file)
assigned_atoms = set()
for atom in mol.atoms:
    dictionary = eval(atom.metadata['dict']) # this is very not good
    try:
        test = dictionary['already_matched']
        assigned_atoms.add(atom.molecule_atom_index)
    except Exception:
        pass
print(len(assigned_atoms))
print(mol.n_atoms)
# print(rdmol.GetNumAtoms())