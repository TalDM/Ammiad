#!/bin/bash

# Submit a script to summarise the simmiad output into a more useable format.
# This is heavy enough lifting to warrant a SLURM job that saves the output.

#SBATCH --job-name=summarise_sims
#SBATCH --nodes=1
#SBATCH --qos=medium
#SBATCH --time=8:00:00
#SBATCH --mem=20gb
#SBATCH --output=simulations/neutral/slurm/summarise_sims%J.out
#SBATCH --error=simulations/neutral/slurm/summarise_sims%J.err

# ENVIRONMENT #
ml build-env/f2022
ml r/4.1.2-foss-2021b

Rscript simulations/neutral/003.summarise_simulations.R
