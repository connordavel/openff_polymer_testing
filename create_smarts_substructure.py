# function that takes a random subset of a pdb molecule and returns an rdkit molecule
# useful when you have a really large molecule and only want to look at a few chains at a time

from pathlib import Path
import sys
from pygments import highlight
from rdkit import Chem
from random import randint
from copy import deepcopy
from openff.toolkit.topology.molecule import FrozenMolecule, Molecule
# from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem

def get_smarts_substructure_from_rdkit_idx(mol, ids):
    # assume either Molecule or RDKit molecule
    if type(mol) == Molecule:
        rdmol = mol.to_rdkit()
    else:
        rdmol = mol
    atom_map_nums = [a.GetAtomMapNum() for a in rdmol.GetAtoms()]
    if atom_map_nums == ([0] * rdmol.GetNumAtoms()): # if empty atom map numbers
        n = 1
        for a in rdmol.GetAtoms():
            if a.GetIdx() in ids:
                a.SetAtomMapNum(n)
                n += 1
        # [a.SetAtomMapNum(a.GetIdx() + 1) for a in rdmol.GetAtoms()]
    smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse = ids)

    return smarts, rdmol

def get_atoms_between_rdmol(rdmol, exclusive_wall_ids, seed, inclusive_wall_ids=[]):
    # a fun little function that gets all atom ids until it hits "walls" specified by the user
    wall_ids = [*exclusive_wall_ids, *inclusive_wall_ids]
    found = [] # atom ids found that we need not go over again
    active = [seed] # atom ids where the algorithm is currently centered 
    n = 0
    while len(active) != 0 and n < 2*rdmol.GetNumAtoms():
        for active_idx in active:
            active_atom = rdmol.GetAtomWithIdx(active_idx)
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

def sample_molecule(pdbfile, subset_size, seed = -1):
    # if seed == -1, when the function will chose a random seed location in the molecule
    # if seed is outside the bounds of the molecule, an error message will play
    
    # read in the rdmol:
    rdmol = Chem.rdmolfiles.MolFromPDBFile(pdbfile, removeHs=False)

    # select a seed location if not specified
    n_atoms = rdmol.GetNumAtoms()
    if seed == -1:
        seed = randint(0, n_atoms - 1)
    elif seed >= n_atoms:
        n_atoms = rdmol.GetNumAtoms()
        print(f"seed ({seed}) is too large for given molecule of size ({n_atoms}) -> (index starts from 0)")
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

    rd_ids = get_substructure_ids(rdmol, seed, subset_size)
    if len(rd_ids) < subset_size:
        print(f"molecule substructure is smaller than the specificed size of {subset_size}.")
    # there is no clean way to get a chemically correct substructure with these ids, so
    # I will use rdkit's smarts functionallity to select atoms and hopefully sanitize the molecule
    # of most large chemical errors (info stored in pdb_block and sub_rdmol). I will also return
    # a separate rdmol called "manual_substructure" just to see if manually deleting atoms can 
    # actually work 

    # manual solution -> done this way since atom ids tend to change when deleting atoms 
    manual_substructure = deepcopy(rdmol)
    mw = Chem.RWMol(manual_substructure)
    for atom in mw.GetAtoms():
        if atom.GetIdx() in rd_ids:
            continue
        atom.SetAtomicNum(0)
    manual_substructure = None
    manual_substructure = Chem.DeleteSubstructs(mw, Chem.MolFromSmarts('[#0]'))

    # smarts solution using rdkits smarts function
    [atom.SetAtomMapNum(atom.GetIdx()) for atom in rdmol.GetAtoms()]
    sub_smarts = Chem.MolFragmentToSmarts(rdmol, atomsToUse = rd_ids)
    sub_rdmol = Chem.rdmolfiles.MolFromSmarts(sub_smarts)
    pdb_block = Chem.rdmolfiles.MolToPDBBlock(sub_rdmol)

    # output the seed location and the selected subset rdmolecule 

    return sub_rdmol, seed, pdb_block, manual_substructure

def smarts_string_cleaner(smarts):
    from collections import defaultdict
    # for some reason, rdkit likes to return smarts strings with atomic numbers instead of symbols
    # so this is a very messy way to deal with it
    # this also returns a formated string block for insertion into the sybstructure library 
    # honestly, if you trace an error back to this function, best just do it manually 
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

def create_smarts_substructure(pdbfile, 
                               exclusive_wall_ids=[], 
                               seed=-1, 
                               inclusive_wall_ids=[], 
                               sample_size=-1,
                               sample_seed=-1):
    # the culmination of the functions from this file
    # With only a pdbfile location specified, this function will only
    #       return an offmol and rdmol and nothing else
    # exlusive_wall_ids, seed, and inclusive_wall_ids all determine how the function
    #       selects a smarts substructure from the molecule
    # For small pdb files, sample_size and sample_seed need not be used at all
    # If sample size is specified (not -1), "sample_molecule()" will be called
    #       to get a subset of the original molecule. This helps with large molecule files
    # If sample seed is specified, the sample_molecule algoithm will start at a certain
    #       rdkit index 

    
    # sample if needed
    if sample_size != -1:
        rdmol, seed_loc, pdb_block, manual_rdmol = sample_molecule(pdbfile, sample_size, sample_seed)
    else:
        rdmol = Chem.rdmolfiles.MolFromPDBFile(pdbfile, removeHs=False)

    # prevents implicit hydrogens from displaying later on
    for atom in rdmol.GetAtoms():
        atom.SetProp("atomLabel", atom.GetSymbol())

    rdmol_with_original_map_ids = deepcopy(rdmol)
    new_smarts = None
    smarts_block = None
    selected_atoms = []

    if seed != -1:
        selected_atoms = get_atoms_between_rdmol(rdmol, exclusive_wall_ids, seed, inclusive_wall_ids)
        smarts, rdmol_sub = get_smarts_substructure_from_rdkit_idx(rdmol, ids=selected_atoms)
        if "#" in smarts:  # if smarts uses atomic numbers
            new_smarts, smarts_block = smarts_string_cleaner(smarts)
    
    #run test to return atoms that are currently able to be assigned
    assigned_atom_map_ids, _, _ = get_assignments(pdbfile)
    assigned_atoms_ids = set()
    unassigned_atoms_ids = set()
    for atom in rdmol_with_original_map_ids.GetAtoms():
        if atom.GetAtomMapNum() in assigned_atom_map_ids:
            assigned_atoms_ids.add(atom.GetIdx())
        else:
            unassigned_atoms_ids.add(atom.GetIdx())

    return rdmol_with_original_map_ids, new_smarts, smarts_block, list(selected_atoms), assigned_atoms_ids, unassigned_atoms_ids

def get_assignments(pdbfile):
    # returns which atoms where assigned in "from_pdb" successfully
    mol = Molecule.from_pdb(pdbfile)
    assigned_atoms = set()
    unassigned_atoms = set()
    for atom in mol.atoms:
        if atom.metadata['already_matched']:
            assigned_atoms.add(atom.molecule_atom_index)
        else:
            unassigned_atoms.add(atom.molecule_atom_index)
    return assigned_atoms, unassigned_atoms, mol

def rdkit_visualize(
        rdmol,
        width=None,
        height=None,
        show_all_hydrogens=True,
        grey_highlight_ids=[], 
        yellow_highlight_ids=[],
        show_3D = False,
        
    ):
    # a copy of the visualize function from the openff Molecule class
    # avoids having to go through the openforcefield toolkit, which tends
    # to not handle incomplete/nonsensical molecules well 

    from IPython.display import SVG
    from rdkit.Chem.Draw import (  # type: ignore[import]
        rdDepictor,
        rdMolDraw2D,
    )
    from rdkit.Chem.rdmolops import RemoveHs  # type: ignore[import]

    # change rdkit map nums to zero
    rdmol_new = deepcopy(rdmol)

    if show_3D == False:
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
        drawer.drawOptions().addAtomIndices = True
        if len(grey_highlight_ids) > 0:
            drawer.DrawMolecule(rdmol_new, highlightAtoms=grey_highlight_ids)
        else:
            drawer.DrawMolecule(rdmol_new)
        drawer.FinishDrawing()

        return SVG(drawer.GetDrawingText())
    else:
        # use nglview for this
        import nglview as nv

        # first, assign all atoms to be highlighted to their own ligand group
        grey_highlighted_info = Chem.AtomPDBResidueInfo()
        grey_highlighted_info.SetResidueName("asddsf") #completely random name that must be > 5 letters for some reason? 

        yellow_highlighted_info = Chem.AtomPDBResidueInfo()
        yellow_highlighted_info.SetResidueName("zxcvbn") #completely random name that must be > 5 letters for some reason? 

        unhighlighted_info = Chem.AtomPDBResidueInfo()
        unhighlighted_info.SetResidueName("sfdfsg") #completely random name that must be > 5 letters for some reason? 
        for atom in rdmol_new.GetAtoms():
            if atom.GetIdx() in yellow_highlight_ids:
                atom.SetMonomerInfo(yellow_highlighted_info)
            elif atom.GetIdx() in grey_highlight_ids:
                atom.SetMonomerInfo(grey_highlighted_info)
            else:
                atom.SetMonomerInfo(unhighlighted_info)

        # next, construct the viewer         
        view = nv.show_rdkit(rdmol_new)

        # grey highlights
        view.add_spacefill(".ASD", opacity=0.4, color="white")
        view.add_spacefill(".ASD and _H", opacity=0.4, color="red")
        view.add_spacefill(".ASD and _H", opacity=0.4, color="red")

        #yellow highlights
        view.add_spacefill(".ZXC", opacity=0.4, color="yellow")

        return view


if __name__ == "__main__":
    path_str = "polymer_examples/rdkit_simple_polymers/polyethylene.pdb"
    path_loc = Path(path_str)
    if not path_loc.exists():
        path_loc = Path("openff_polymer_testing/" + path_str)
    if not path_loc.exists():
        print("could not find path given")
        sys.exit()

    rdmol, new_smarts, smarts_block, selected_atoms, assigned_atoms, unassigned_atoms = create_smarts_substructure(str(path_loc), 
                                                                                        exclusive_wall_ids=[1], 
                                                                                        seed=20, 
                                                                                        inclusive_wall_ids=[], 
                                                                                        sample_size=70,
                                                                                        sample_seed=-1)

    # view = rdkit_visualize(rdmol, 900, 900, highlight_ids=selected_atoms, show_3D=True)
    # print(view)
    print(len(assigned_atoms))
    print(unassigned_atoms)