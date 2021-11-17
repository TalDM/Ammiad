#! /bin/bash
#
# Script to run different simmiad simulations with different parameters as a 
# SLURM job array on a high-performance cluster.

# SLURM
#SBATCH --job-name=simmiad_sims
#SBATCH --output=data_analysis/simulations/003.simmiad_sims.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=thomas.ellis@gmi.oeaw.ac.at
#SBATCH --mem-per-cpu=10GB
#SBATCH --qos=long
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-258

# ENVIRONMENT #
ml build-env/f2020
ml r/3.6.0-foss-2019a

# Generate R scripts for each parameter combination.
Rscript data_analysis/simulations/002.sim_generator.R
# Submit each simulation.
FILES=(data_analysis/simulations/simmiad_scripts/*.R)
srun Rscript ${FILES[$SLURM_ARRAY_TASK_ID]}
