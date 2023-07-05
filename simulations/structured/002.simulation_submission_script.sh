#! /bin/bash
#
# Script to run different simmiad simulations with different parameters as a 
# SLURM job array on the VBC high-performance computing cluster.

# SLURM
#SBATCH --job-name=structured_sims
#SBATCH --output=simulations/structured/slurm/structured_sims%J.out
#SBATCH --error=simulations/structured/slurm/structured_sims%J.err
#SBATCH --mem-per-cpu=10GB
#SBATCH --qos=long
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-255

# ENVIRONMENT #
ml build-env/f2022
ml r/4.1.2-foss-2021b

# Generate R scripts for each parameter combination.
Rscript simulations/structured/001.sim_generator.R
# Submit each simulation.
FILES=(simulations/structured/simmiad_scripts/*.R)
srun Rscript ${FILES[$SLURM_ARRAY_TASK_ID]}
