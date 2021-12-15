#!/bin/bash
# Job array
# MUST MAKE -o directory before submitting.
#SBATCH --job-name=pmssm
#SBATCH -p free
#SBATCH --array=0-246
#SBATCH --cpus-per-task=8
#SBATCH -o /path/for/std/out

omp_threads=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$omp_threads

cd /path/to/src/src

# Run the object file
path/to/src/src/dm_file.o $SLURM_ARRAY_TASK_ID
