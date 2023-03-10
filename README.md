<h1>Genetic analysis and simulations on Ammiad wild wheat population</h1>

<h2> Summary </h2>

This repository contains code used for the analysis performed on Genotyping-By-Sequencing data of Ammiad wild wheat collection.

<h2> Data processing </h2>

<h3> Demultiplexing GBS Illumina reads </h3>

A python/2.7 [code](data_processing/GBS_demultiplex.py) was used to demultiplex the raw reads obtained from the Illumina .fastq files.

<h3> Variant calling pipeline </h3>

The UNIX commands used for running variant calling are stored in [variant_calling_pipeline.txt](data_processing/variant_calling_pipeline.txt). The pipeline is using BWA, Samtools, GATK, VCFtools and TASSEL (versions enclosed in the file).

<h3> Calculating pairwise identity and assignment to DGGs </h3>

- Distinct genotype groups (DGGs) were assigned by pairwise genetic identity. 
- To calculate identity and assign to DGGs, use [identity.r](data_processing/identity/identity.r). 
- If you wish to create a summary .csv table with the coordinates, euclidean distances, microhabitats and SNP summary, load the [identity_function.r](data_processing/identity/identity_function.r), using the [identity_function_run_example.r](data_processing/identity/identity_function_run_example.r). 

<h3> Other data-processing scripts </h3>

- `import_field_occupancy_data.R` Prepares data on which DGGs were found in 
which sampling points, when. Returns a tibble called `obs_geno` with columns for
sample ID, position in the transect, DGG, year, habitat, distance (in metres) 
along the transect, and transect ID.
- `import_phenotype_data.R` prepares phenotype data for analysis. Returns a 
tibble giving, sample ID, DGG, experimental block, then columns for phenotypes.
- `ammiad_spatial_clustering.R`: Observed spatial clustering at Ammiad. 
Returns a tibble `obs_di` ('di' is short for 'distance identity') giving 
the spatial bin (in metres) from `cut`, the probability that two sampling 
points are occupied by the same DGG (`mean`), number of pairs of sampling
that could be compared, an integer value for the bin, and transect.
- `ammiad_temporal_stability`: Calculate the probability of plants of the same 
DGG in the same sampling point at harvests separated by different number of 
years for each transect.
Returns a tibble `obs_stability` with columns for transect, number of years
between harvests being compared, the probability that a sampling point is 
occupied by the same DGG in those years (stability), and the number of pairs
of sampling point that could be compared.

<h2> Data analysis </h2>

<h3> Rarefaction analysis </h3>

Rarefaction analysis and plots for Ammiad collection, per year, was generated using [rarefaction.r](data_analysis/rarefaction.r).


<h3> Principal component analysis </h3>

Principal component analysis was done using the SNPRelate library directly from the VCF file. Code for running the analysis and generating the plots can be found in [pca.r](data_analysis/pca.r).

 <h3> Habitat sorting: permutations </h3>

Code to estimate whether sorting by microhabitat can be explained by chance, and
to recreate figure S9, are shown in the directory `habitat_permutations`.
The RMarkdown file imports the functions to do this and creates the plot.

<h3> Individual-based simulations </h3>

This folder contains scripts to run individual-based simulations of populations
evolving under purley neutral demographic processes (i.e. a low level of 
outcrossing, seeds disperse within metres of the mother, but there is no 
selection through survival, fecundity or with respect to microhabitats).
These simulations use library code available as the custom R package `simmiad`
available from [Github](https://github.com/ellisztamas/simmiad)
and
[Zenodo](10.5281/zenodo.4762083).

The folder contains four files to run and analyse the simulations:

- `001.simulation_submission_script.sh`: Master bash script to generate and run
simulations. SLURM message about errors and output for each job are saved to 
a folder `slurm`.
 - `002.sim_generator.R`: R script to generate a folder of scripts to run
 `simmiad` with a unique combination of parameters.
 - `003.summarise_simulations.R`: R script to summarise simulations by 
 calculating means and 95% confidence intervals over simulation replicates. 
Saves an `.RDS` file which can be quickly imported.
 - `004.simmiad_from_one_genotype.Rmd`: RMarkdown file to plot the results.

 This recreates figures 3, S11, S12, S13 and S14.

<h3> Estimating plant density </h3>

Plant density was measured in 2020, in transect A and summarised in [2020_plant_densities.csv](data/2020_plant_densities.csv). Spikes were counted around each sampling point peg (1 meter radius) and between each sampling points (in 2 meter width course). Mean density calculation and plot was done using [plant_density_2020.r](data_analysis/plant_density_2020.r).

 <h3> Phenotypes </h3>
 
Scripts to analyse the phenotype data from the common-garden nethouse experiment
and predict phenotypes in the field:

* `001_linear_models.R`       : Script to fit linear models to each phenotype with BRMS
* `002_submit_linear_models.sh` : Submit the linear models as a SLURM job
* `003_heritability.R`          : Estimate and plot broad-sense heritability
* `004_predict_phenotypes.R`    : Predict phenotypes for 801 samples at Ammiad and calculate variance explained by habitat and year
* `005_submit_predictions.sh`   : Submit the predictions as a SLURM job
* `006_plot_predictions.R`      : Plot variance explained by habitat.
 