#!/bin/bash
#
#PBS -N hmmTraining
#PBS -l walltime=00:02:00
#PBS -l pvmem=1gb
#PBS -m bea
#PBS -M ag568@leicester.ac.uk

# Make python  available
module load python/gcc/36

# Execute the MPI job code
python3 hidden_markov_modeling/analysis.py --config hidden_markov_modeling/config.json
