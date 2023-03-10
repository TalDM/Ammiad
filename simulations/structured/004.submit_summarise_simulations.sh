#!/bin/bash

# Submit the BRMS model script to run on the cluster.

#SBATCH --job-name=summarise_sims
#SBATCH --nodes=1
#SBATCH --qos=medium
#SBATCH --time=8:00:00
#SBATCH --mem=20gb
#SBATCH --output=simulations/structured/slurm/summarise_sims%J.out
#SBATCH --error=simulations/structured/slurm/summarise_sims%J.err

# ENVIRONMENT #
ml build-env/f2022
ml r/4.1.2-foss-2021b

Rscript simulations/structured/003.summarise_simulations.R
