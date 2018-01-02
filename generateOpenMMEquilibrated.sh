#!/usr/bin/env bash
# quick script for generating initial conditions for adaptive sampling with OpenMM in HTMD.
# loops over each file in initial_conditions, and runs openmm_equilibration.py on it, placing the
# results in generator subdirectories.
mkdir ./generators
i=0
for file in initial_conditions/*; do
    echo "Generating initial conditions for file: ${file}";
    output_dir="./generators/${i}"
    mkdir ${output_dir}
    python3 openmm_equilibration.py ${file} -o ${output_dir}
    cp generator_stuff/* ${output_dir}
    ((i++))
done