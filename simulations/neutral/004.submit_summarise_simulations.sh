#!/bin/bash

# Submit the BRMS model script to run on the cluster.

#SBATCH --job-name=summarise_sims
#SBATCH --nodes=1
#SBATCH --qos=medium
#SBATCH --time=8:00:00
#SBATCH --mem=20gb
#SBATCH --output=simulations/neutral/slurm/summarise_sims%J.out
#SBATCH --error=simulations/neutral/slurm/summarise_sims%J.err

# ENVIRONMENT #
module load build-env/f2022
module load r/4.0.2-foss-2018b

Rscript simulations/neutral/003.summarise_simulations.R
