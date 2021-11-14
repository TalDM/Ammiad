#! /bin/bash
#
# Script to run different simmiad simulations with different parameters

# SLURM
#SBATCH --job-name=simmiad_sims
#SBATCH --output=05_results/12_simmiad_from_unique_genotype/003.simmiad_sims.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=thomas.ellis@gmi.oeaw.ac.at
#SBATCH --mem-per-cpu=10GB
#SBATCH --qos=long
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-63

# ENVIRONMENT #
ml build-env/f2020
ml r/3.6.0-foss-2019a

# Generate R scripts for each parameter combination.
Rscript 05_results/12_simmiad_from_unique_genotype/002.sim_generator.R
# Submit each simulation.
FILES=(05_results/12_simmiad_from_unique_genotype/simmiad_scripts/*a02*.R)
srun Rscript ${FILES[$SLURM_ARRAY_TASK_ID]}
