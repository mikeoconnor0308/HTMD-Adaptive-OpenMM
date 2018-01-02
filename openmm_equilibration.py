"""
Loads a pdb file (+ dcd file if required) minimizes and equilibrates using Amber10 forcefield & implicit solvent.
It also outputs an AceMD .coor file for compatibility with HTMD's adaptive sampling.
For follow up use in adaptive sampling in HTMD.
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

parser = argparse.ArgumentParser(description="Run an OpenMM Simulation")
parser.add_argument('topology', type=str, help='PDB file to run simulation from')
parser.add_argument('-c', '--coords', type=str, help='DCD file to load coordinates from')
parser.add_argument('-f', '--frame', type=int, help='DCD frame to load coordinates from', default=0)
parser.add_argument('-T', '--temperature', type=float, help='Temperature for simulation', default=300)
parser.add_argument('-teq', '--equilibration_time', type=float,
                    help='Time (in ns) to equilibrate system', default=2)
parser.add_argument('-ts', '--time_step', type=float, help='Time step (in fs)', default=1.0)
parser.add_argument('-fr', '--friction', type=float, help='Friction coefficent (1/ps)', default=1.0)
parser.add_argument('-o', '--output_path', type=str, help='Path to output trajectory files to. Defaults to path of input file', default="")
args = parser.parse_args()

temp = args.temperature
t_equil = args.equilibration_time
fric = args.friction
pdb_str = args.topology

# load pdb from file using MDTraj
if args.coords:
    traj = mdtraj.load_dcd(args.coords, top=pdb_str)
else:
    traj = mdtraj.load(pdb_str)
# add traj argument to set unit cell information.
topology = traj.topology.to_openmm(traj)

forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')

print("Creating model for pdb {} from frame {}".format(pdb_str, args.frame))
modeller = app.Modeller(topology, traj.openmm_positions(args.frame))

system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.CutoffPeriodic,
                                 nonbondedCutoff=2.0 * unit.nanometers, constraints=app.HBonds)

timestep = args.time_step * unit.femtoseconds
integrator = mm.LangevinIntegrator(temp * unit.kelvin, fric / unit.picoseconds,
                                   timestep)
integrator.setConstraintTolerance(0.00001)

platform = mm.Platform.getPlatformByName('CUDA')
simulation = app.Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

print("minimizing...")
simulation.minimizeEnergy()

if args.output_path != "":
    output_path = args.output_path + "/"
else:
    output_path = os.path.splitext(pdb_str)[0]

simulation.context.setVelocitiesToTemperature(temp * unit.kelvin)

# equilibrate
print("Equilibrating ...")
nsteps = int(t_equil * unit.nanoseconds / timestep)
simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, time=True, progress=True,
                                                  potentialEnergy=True, temperature=True, remainingTime=True,
                                                  speed=True, totalSteps=nsteps, separator=','))
simulation.step(nsteps)

positions = simulation.context.getState(getPositions=True).getPositions()
pdb_out = '{}/input.pdb'.format(output_path)
app.PDBFile.writeFile(simulation.topology, positions, open(pdb_out, 'w'))

# HTMD adaptive sampling doesn't work very well with pdb files, so for now we'll generate a .coor file to use.
# load the coor file into htmd, if it exists and store as a pdb file.
from htmd.ui import *
htmd_mol = Molecule(pdb_out)
htmd_mol.write('{}/input.coor'.format(output_path))