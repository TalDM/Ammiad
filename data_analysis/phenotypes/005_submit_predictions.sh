#!/bin/bash

# Submit the BRMS model script to run on the cluster.

#SBATCH --job-name=predict_phenotypes
#SBATCH --nodes=1
#SBATCH --qos=medium
#SBATCH --time=4:00:00
#SBATCH --mem=20gb
#SBATCH --output=data_analysis/phenotypes/slurm/predict_phenotypes_%J.out
#SBATCH --error=data_analysis/phenotypes/slurm/predict_phenotypes_%J.err

# ENVIRONMENT #
module load build-env/f2022
module load r/4.2.0-foss-2021b

Rscript data_analysis/phenotypes/004_predict_phenotypes.R
