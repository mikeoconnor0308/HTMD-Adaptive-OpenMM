#!/bin/bash
python openmm_adaptive_sim.py input.pdb input.coor >log.txt 2>&1
# get rid of the input file created by htmd for reading with mdtraj
rm input_coor.dcd
# need to do this for now, as htmd  appears to only copy xtc files over to data, and i don't know how to fix that.
mdconvert -f traj.dcd -o traj.xtc
# remove the dcd file, as the xtc has everything we need.
rm traj.dcd
