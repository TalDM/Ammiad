#! /bin/bash
#
# Script to run different simmiad simulations with different parameters as a 
# SLURM job array on the VBC high-performance computing cluster.

# SLURM
#SBATCH --job-name=structured_sims
#SBATCH --output=simulations/slurm/%x-%a.out
#SBATCH --error=simulations/slurm/%x-%a.err
#SBATCH --mem-per-cpu=10GB
#SBATCH --qos=long
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-255

# ENVIRONMENT #
ml build-env/f2022
ml r/4.1.2-foss-2021b

# Generate R scripts for each parameter combination.
# If not run already
# Rscript simulations/001.sim_generator.R
# Submit each simulation.
FILES=(simulations/simmiad_scripts/*.R)
srun Rscript ${FILES[$SLURM_ARRAY_TASK_ID]}
