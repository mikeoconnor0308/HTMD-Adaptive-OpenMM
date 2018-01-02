"""
Loads a pdb file and coor file from HTMD, converts it to something sensible to run in OpenMM, and runs it
Outputs a dcd file.
"""
#/usr/bin/env python3
from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import argparse
import time
import sys
import os
import mdtraj
import mdtraj.reporters
from htmd.ui import *


parser = argparse.ArgumentParser(description="Run an OpenMM Simulation")
parser.add_argument('topology', type=str, help='PDB file to run simulation from')
parser.add_argument('coor', type=str, help='AceMD Coordinate file to load coordinates from.')
parser.add_argument('-T', '--temperature', type=float, help='Temperature for simulation', default=300)
parser.add_argument('-t', '--simulation_time', type=float,
                    help='Time (in ns) to run simulation', default=50)
parser.add_argument('-l', '--log_interval', type=float, help='Log interval (in ns)', default=0.1)
parser.add_argument('-ts', '--time_step', type=float, help='Time step (in fs)', default=2.0)
parser.add_argument('-fr', '--friction', type=float, help='Friction coefficent (1/ps)', default=1.0)
parser.add_argument('-o', '--output_path', type=str, help='Path to output trajectory files to. Defaults to path of input file', default="")
args = parser.parse_args()

temp = args.temperature
t_sim = args.simulation_time
fric = args.friction
pdb_str = args.topology

# load the coor file into htmd, if it exists and store as a pdb file.
htmd_mol = Molecule(pdb_str)
htmd_mol.read(args.coor)
htmd_mol.write('input_coor.dcd')

# load the input pdb file in mdtraj, for some reason i can't load the dcd at the same time, as periodic bounds are lost.

traj = mdtraj.load('input.pdb')
# add traj argument to set unit cell information.
topology = traj.topology.to_openmm(traj)

forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')

print("Creating model for pdb {} from frame {}".format(pdb_str, 0))

system = forcefield.createSystem(topology, nonbondedMethod=app.CutoffPeriodic,
                                 nonbondedCutoff=2.0 * unit.nanometers, constraints=app.HBonds)

timestep = args.time_step * unit.femtoseconds
integrator = mm.LangevinIntegrator(temp * unit.kelvin, fric / unit.picoseconds,
                                   timestep)
integrator.setConstraintTolerance(0.00001)

platform = mm.Platform.getPlatformByName('CUDA')
simulation = app.Simulation(topology, system, integrator, platform)

#load positions saved from htmd.
traj = mdtraj.load('input_coor.dcd', top='input.pdb')
simulation.context.setPositions(traj.openmm_positions(0))

if args.output_path != "":
    output_path = args.output_path + "/"
else:
    output_path = os.path.splitext(pdb_str)[0]

simulation.context.setVelocitiesToTemperature(temp * unit.kelvin)

print("Running MD ...")
nsteps = int(t_sim * unit.nanoseconds / timestep)
log_interval = int(args.log_interval * unit.nanoseconds / timestep)
simulation.reporters.append(app.DCDReporter('traj.dcd', log_interval))
simulation.reporters.append(app.StateDataReporter(sys.stdout, log_interval, step=True, time=True, progress=True,
                                                  potentialEnergy=True, temperature=True, remainingTime=True,
                                                  speed=True, totalSteps=nsteps, separator=','))
simulation.step(nsteps)

