from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sys


def main(input_file, time_length=20, temperature=300, step_size=0.002, output_file=None):
    num_steps = float(time_length) * 1000 / step_size # convert ns to ps
    
    if output_file is None:
        output_file = input_file.replace('.pdb', '_output.pdb')
    pdb = PDBFile(input_file)
    forcefield = ForceField("amber/protein.ff14SB.xml", 'amber/tip3p.xml')
    # forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=10*angstrom, constraints=HBonds)
    integrator = LangevinIntegrator(temperature*kelvin, 1/picosecond, step_size*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    
    simulation.reporters.append(PDBReporter('output.pdb', 1000))
    simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True,
            potentialEnergy=True, temperature=True))
    simulation.step(num_steps)
    
if __name__ == '__main__':
    input_file = sys.argv[1]
    if len(sys.argv) > 2:
        time_length = sys.argv[2] # nanoseconds
    else:
        time_length = 20
    if len(sys.argv) > 3:
        temperature = sys.argv[3] # Kelvin
    else:
        temperature = 300
    if len(sys.argv) > 4:
        step_size = sys.argv[4]
    else:
        step_size = 0.002
    if len(sys.argv) > 5:
        output_file = sys.argv[5]
    else:
        output_file = None
    
    main(input_file, time_length, temperature, step_size, output_file) 