#! /bin/bash
#
# Script to run different simmiad simulations with different parameters as a 
# SLURM job array on the VBC high-performance computing cluster.

# SLURM
#SBATCH --job-name=neutral_simulations
#SBATCH --output=simulations/neutral/slurm/%J.out
#SBATCH --error=simulations/neutral/slurm/%J.err
#SBATCH --mem-per-cpu=10GB
#SBATCH --qos=medium
#SBATCH --time=1-12:00:00
#SBATCH --ntasks=1
#SBATCH --array=16,18,49,98

#-255

# ENVIRONMENT #
ml build-env/f2022
ml r/4.1.2-foss-2021b

# Generate R scripts for each parameter combination.
Rscript simulations/neutral/001.sim_generator.R
# Submit each simulation.
FILES=(simulations/neutral/simmiad_scripts/*.R)
srun Rscript ${FILES[$SLURM_ARRAY_TASK_ID]}
