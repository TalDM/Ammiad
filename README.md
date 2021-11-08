<h1>Genetic analysis and simulations on Ammiad wild wheat population</h1>

<h2> Summary </h2>

This repository contains code used for the analysis performed on Genotyping-By-Sequencing data of Ammiad wild wheat collection.


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
