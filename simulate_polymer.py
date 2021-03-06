from pathlib import Path
import sys
from pygments import highlight
from rdkit import Chem
from random import randint
from copy import deepcopy
from openff.toolkit.topology.molecule import FrozenMolecule, Molecule
# from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem
import openmm
from openmm.app import PDBFile, Element
from simtk import unit
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
import numpy
import time
import parmed

def simulate_polymer(pdbfile, substructure_file, offxml_file, output):
    # mol should already have one conformer...

    mol_conf = Molecule.from_pdb(pdbfile, substructure_file)
    pdbfile = PDBFile(pdbfile)
    omm_topology = pdbfile.topology

    off_topology = mol_conf.to_topology()
    forcefield = ForceField(offxml_file)
    start = time.time()
    system = forcefield.create_openmm_system(off_topology, allow_nonintegral_charges=True) #WARNING: I have no idea that this means 
    end = time.time()
    difference = end - start
    time_step = 2*unit.femtoseconds  # simulation timestep
    temperature = 1000*unit.kelvin  # simulation temperature
    friction = 1/unit.picosecond  # collision rate
    integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
    simulation = openmm.app.Simulation(omm_topology, system, integrator)
    positions = pdbfile.getPositions() 
    simulation.context.setPositions(positions)

    pdb_reporter = openmm.app.PDBReporter(f'{output}.pdb', 10)
    dcd_reporter = openmm.app.DCDReporter(f'{output}.dcd', 10)
    simulation.reporters.append(pdb_reporter)
    simulation.reporters.append(dcd_reporter)
    
    simulation.minimizeEnergy(maxIterations=1000)
    simulation.step(1000)
    st = simulation.context.getState(getPositions=True, getEnergy=True)
    print(st.getPotentialEnergy())
    # print(st.getPositions())
    unitless_positions = []
    for vec in st.getPositions():
        x = vec.x * 10     # please don't let me forget this
        y = vec.y * 10     # how to do this in a... better... way 
        z = vec.z * 10
        unitless_positions.append([x, y, z])
    unitless_positions = numpy.array(unitless_positions)
    return st, difference

if __name__ == "__main__":
    path_str = "openff_polymer_testing/polymer_examples/rdkit_simple_polymers/PEO.pdb"
    substructure_file = "automatic_PEO_substructures.json"
    path_loc = Path(path_str)
    if not path_loc.exists():
        path_loc = Path("openff_polymer_testing/" + path_str)
    if not path_loc.exists():
        print("could not find path given")
        sys.exit()

    # offxml_file = 'openff_unconstrained_no_library_charges-2.0.0.offxml'
    # st, diff = simulate_polymer(str(path_loc), substructure_file, offxml_file, "PEO_traj_no_library_charges")
    # print(f"time to create openmm system: {diff}")

    offxml_file = 'openff_unconstrained_with_library_charges-2.0.0.offxml'
    st, diff = simulate_polymer(str(path_loc), substructure_file, offxml_file, "PEO_traj_with_library_charges")
    print(f"time to create openmm system: {diff}")

