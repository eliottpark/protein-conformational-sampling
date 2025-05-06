from openmm.app import *
from openmm import *
from simtk.unit import *
import sys
import os
import shutil
import datetime
from multiprocessing import Process


def main(input_file, time_length=20, temperature=300, step_size=0.002, recording_step_interval=10000, output_file=None):
    # Create new subfolder
    file_name = os.path.splitext(os.path.basename(input_file))[0]
    if not os.path.exists('runs'):
        os.makedirs('runs')
    folder_name = os.path.join(os.getcwd(), 'runs', file_name + datetime.datetime.now().strftime("_%Y%m%d_%H%M%S"))
    if not os.path.exists(folder_name):
        os.makedirs(folder_name) 
    # Copy input file to new folder
    new_input_file = os.path.join(folder_name, os.path.basename(input_file))
    shutil.copyfile(input_file, new_input_file)
    input_file = new_input_file
    os.chdir(folder_name)
    # Create new log file
    if os.path.exists('state.log'):
        os.remove('state.log')
    log = open('state.log', 'w')
    # Create new output file
    if output_file is None:
        output_file = new_input_file.replace('.pdb', '_output.pdb')
    if os.path.exists(output_file):
        os.remove(output_file)
    output_trajectory = output_file.replace('.pdb', '_trajectory.dcd')
    if os.path.exists(output_trajectory):
        os.remove(output_trajectory)
    
    num_steps = float(time_length) * 1000 / step_size # convert ns to ps
    print(f"Running simulation for {num_steps} steps")
    print(f"Simulation time length: {time_length} ns")
    print(f"Temperature: {temperature} K")
    print(f"Step size: {step_size} ps")
    print(f"Output file: {output_file}")
    print(f"Output trajectory file: {output_trajectory}")
    print(f"Recording step interval: {recording_step_interval}")
    
    
    pdb = PDBFile(input_file)
    forcefield = ForceField("amber/protein.ff14SB.xml", "tip3p.xml") 
    # forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
            #nonbondedCutoff=10*angstrom, constraints=HBonds)

    integrator = LangevinIntegrator(float(temperature)*kelvin, 1/picosecond, float(step_size)*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    
    simulation.reporters.append(PDBReporter(output_file, recording_step_interval))
    simulation.reporters.append(StateDataReporter(log, recording_step_interval, step=True,
            potentialEnergy=True, temperature=True))
    simulation.reporters.append(DCDReporter(output_trajectory, recording_step_interval))
    # simulation.reporters.append(NetCDFReporter(output_trajectory, recording_step_interval,
    #                                       enforcePeriodicBox=True,
    #                                       coordinates=True,
    #                                       velocities=True))
    simulation.step(num_steps)
    log.close()

if __name__ == "__main__":
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
        recording_step_interval = sys.argv[5]
    else:
        recording_step_interval = 10000
    if len(sys.argv) > 6:
        output_file = sys.argv[6]
    else:
        output_file = None
    # p = Process(
    #     target=main,
    #     args=(input_file, time_length, temperature, step_size, output_file),
    #     daemon=False
    # )
    # p.start()
    # # Now p.pid is the child’s OS process ID
    # print(f"Launched simulation in subprocess PID {p.pid}")
    # exit()
    main(input_file, time_length, temperature, step_size, recording_step_interval, output_file)
