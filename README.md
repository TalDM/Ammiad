<h1>Genetic analysis and simulations on Ammiad wild wheat population</h1>

<h2> Summary </h2>

This repository contains code used for the analysis performed on Genotyping-By-Sequencing data of Ammiad wild wheat collection.

<h3> Demultiplexing GBS Illumina reads </h3>

A python/2.7 [code](data_processing/GBS_demultiplex.py) was used to demultiplex the raw reads obtained from the Illumina .fastw files.


<h3> Variant calling pipeline </h3>

The UNIX commands used for running variant calling are stored in [variant_calling_pipeline.txt](variant_calling_pipeline.txt). The pipeline is using BWA, Samtools, GATK, VCFtools and TASSEL (versions enclosed in the file).


<h3> Calculating pairwise identity and assignment to DGGs </h3>

- Distinct genotype groups (DGGs) were assigned by pairwise genetic identity. 
- To calculate identity and assign to DGGs, use [identity.r](identity.r). 
- If you wish to create a summary .csv table with the coordinates, euclidean distances, microhabitats and SNP summary, load the [identity_function.r](identity_function.r), using the [identity_function_run_example.r](identity_function_run_example.r). 


<h3> Rarefaction analysis </h3>

Rarefaction analysis and plots for Ammiad collection, per year, was generated using [rarefaction.r](rarefaction.r).


<h3> Principal component analysis </h3>

Principal component analysis was done using the SNPRelate library directly from the VCF file. Code for running the analysis and generating the plots can be found in [pca.r](pca.r).

<h3> Individual-based simulations </h3>

This folder contains scripts to run individual-based simulations of populations
evolving under purley neutral demographic processes (i.e. a low level of 
outcrossing, seeds disperse within metres of the mother, but there is no 
selection through survival, fecundity or with respect to microhabitats).
These simulations use library code available as the custom R package `simmiad`
available from [Github](https://github.com/ellisztamas/simmiad)
and
[Zenodo](10.5281/zenodo.4762083).

The folder contains six files to run and analyse the simulations:

- `001.simulation_submission_script.sh`: Master bash script to generate and run
simulations
 - `002.sim_generator.R`: R script to generate a folder of scripts to run
 `simmiad` with a unique combination of parameters.
 - `003.simmiad_sims.log`: SLURM log files from each simulation 
 - `004.session_info.txt`: The output of `sessionInfo()` from R showing system
 and package versions.
 - `005.summarise_simulations.R`: R script to summarise simulations by 
 calculating means and 95% confidence intervals over simulation replicates. 
Saves an `.RDS` file which can be quickly imported.
 - `006.simmiad_from_one_genotype.Rmd`: RMarkdown file to plot the results.

 This recreates figures 3, S11, S12, S13 and S14.

 <h3> Habitat sorting: permutations </h3>

Code to estimate whether sorting by microhabitat can be explained by chance, and
to recreate figure S9, are shown in the directory `habitat_permutations`.
The RMarkdown file imports the functions to do this and creates the plot.
