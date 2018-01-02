# Project Title

A set of example scripts for running [High Throughput Molecular Dynamics](https://www.acellera.com/products/high-throughput-molecular-dynamics/) (HTMD) with Adaptive Sampling using OpenMM. 

## Getting Started

### Prerequisites

HTMD is of course required, and is most easily installed through conda
(described [here](https://software.acellera.com//academic-download.html)). The scripts were tested with version 1.11.5, and may require adjustments
for different versions. 

Installing HTMD through conda will automatically install the other dependencies, namely OpenMM and MDTraj. 

### Code Structure

This repository has the following structure:

```
htmd-adaptive-openmm
│   README.md
│   generateOpenMMEquilibrated.sh - Generates equilibrated systems from initial_conditions and places them in generators.
│   openmm_equilibration.py - Runs the equilibration with OpenMM, using implicit solvent.
|   runAdaptiveMD.py - Runs HTMD Adaptive Sampling for the system prepared in generators.
└───initial_conditions
    │   2efv_second_attempt.pdb - Example initial conditions.  
└───generator_stuff
    │   openmm_adaptive_sim.py - Example OpenMM production script using implicit solvent
    │   run.sh - Job script called by HTMD that runs the production MD.
└───generators
    └───0
        |   input.coor - Initial input coordinates using AceMD coordinate file. 
        |   input.pdb - Initial input coordinates in PDB form.
        |   openmm_adaptive_sim.py - Copied from generator_stuff
        |   run.sh - Copied from generator_stuff. 
    .
    .
    .
```

### Generating Initial Conditions

PDB files your system with initial conditions to start sampling from should be placed in the folder `initial_conditions`. Included with the repo
are initial conditions for studying slipknot formation in MJ0366. 

The script `generateOpenMMEquilibrated.sh` can be used to take the starting structures in `initial_conditions`, run equilibration and copy all the 
required files for running adaptive sampling from the equilbrated system into folders within the folder `generators`: 
```
./generateOpenMMEquilibrated.sh
```
The script called the python script `openmm_equilibration.py`, which reads in the starting structures and equilibrates them using the Amber10 forcefield 
and Amber10 implicit solvent. This script should be modified as required for particular systems.  

### Running Adaptive Sampling

To run adaptive sampling after generating initial conditions, the script `runAdaptiveMD.py` can be used as a starting point. 
```
python runAdaptiveMD.py
```

To modify the details of the production simulations run with adaptive sampling, edit the script `generator_stuff/openmm_adaptive_sim.py`. 

## Authors

* **Mike O'Connor** - *Initial work* - [Mike O'Connor](https://github.com/mikeoconnor0308)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

