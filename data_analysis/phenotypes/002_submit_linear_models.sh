#!/bin/bash

# Submit the BRMS model script to run on the cluster.

#SBATCH --job-name=nethouse_models
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --qos=medium
#SBATCH --time=4:00:00
#SBATCH --mem=20gb
#SBATCH --output=nethouse_models_%J.out
#SBATCH --error=nethouse_models_%J.err

# ENVIRONMENT #
module load build-env/f2022
module load r/4.2.0-foss-2021b

Rscript data_analysis/phenotypes/002_linear_models.R
