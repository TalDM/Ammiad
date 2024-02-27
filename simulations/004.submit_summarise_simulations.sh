#!/bin/bash

# Submit the BRMS model script to run on the cluster.

#SBATCH --job-name=summarise_sims
#SBATCH --nodes=1
#SBATCH --qos=rapid
#SBATCH --time=1:00:00
#SBATCH --mem=20gb
#SBATCH --output=simulations/slurm/%a.out
#SBATCH --error=simulations/slurm/%x.err

# ENVIRONMENT #
ml build-env/f2022
ml r/4.1.2-foss-2021b

Rscript simulations/003.summarise_simulations.R
