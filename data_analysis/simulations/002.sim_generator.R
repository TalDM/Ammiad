# Tom Ellis, May 2021
# A script for generating a series of new R scripts to run simmiad using a
# range of parameter values so that I don't have to create and maintain each
# manually and make copy-paste-errors.
#
# It will output one R script in the output directory for every
# combination of parameter values.

set.seed(124)

# Directory to save scripts to.
out_dir <- "data_analysis/simulations/simmiad_scripts"
dir.create(out_dir, showWarnings = FALSE)
# Directory to tell the scripts to save the output to.
sim_dir <- "data_analysis/simulations/output/"
dir.create(sim_dir, showWarnings = F)

# Initialise parameters
dispersal <- c(0.5, 0.75, 1, 2)
outcrossing_rate <- c(0.0025, 0.005, 0.01, 0.02)
n_generations <- 2000
n_starting_genotypes <- c(1)
density <- c(1, 2, 3,5)
dormancy <- c(0, 0.1, 0.3, 0.5)
n_sample_points <- 30 
sample_spacing <- 5
nsims <- 100
range_limit= 1.5
stability_years <- "(n_generations-36) : n_generations"

for(d in dispersal){
  for(o in outcrossing_rate){
    for(b in n_starting_genotypes){
      for (a in density){
        for(s in dormancy){
          # File root name for this simulation.
          fname <- paste(
            'd',  sprintf('%03d', d*100),
            '_o', sprintf('%04d', o*10000),
            '_b', sprintf('%02d', b),
            '_a', sprintf('%02d', a),
            '_s', sprintf('%02d', s*10),
            sep='')
          
          # Stuff at the top of the script
          cat("# Script to run simmiad generated by 05_results/11_simmiad_seed_bank/sim_generator.R\n", 
              "# Script will be saved to 05_results/11_simmiad_seed_bank/simmiad_scripts/", fname, "\n\n",
              
              "library('simmiad')\n",
              # "library('dplyr')\n",
              "set.seed(",sample(1:100,1),")\n\n",
              
              "# Simulation parameters",
              "\ndispersal <-", d,
              "\noutcrossing <-", o,
              "\nn_starting_genotypes <-", b,
              "\ndensity <-", a,
              "\ndormancy <-", s,
              "\nn_generations <-", n_generations,
              "\nn_sample_points <-", n_sample_points,
              "\nsample_spacing <-", sample_spacing,
              "\nrange_limit <-", range_limit,
              "\nnsims <-", nsims,
              "\nstability_years <-", stability_years,
              "\nfilename <- '", fname, "'",
              "\n\n",
              
              "# Code to run the simulation\n",
              "output <- simmiad(\n",
              '\tmean_dispersal_distance = dispersal,\n',
              '\toutcrossing_rate = outcrossing,\n',
              '\tn_generations = n_generations,\n',
              '\tn_starting_genotypes = n_starting_genotypes,\n',
              '\tdensity = density,\n',
              '\tdormancy = dormancy,\n',
              '\tn_sample_points = n_sample_points,\n',
              '\tsample_spacing = sample_spacing,\n',
              '\trange_limit = range_limit,\n',
              '\tnsims = nsims,\n',
              '\tprogress = FALSE,\n',
              '\tstability_years = stability_years\n',
              ')\n\n',
              
              '# Save to disk\n',
              'write_simmiad(\n',
              '\toutput,\n',
              '\tdirectory = paste("', sim_dir, '", filename, sep="")\n',
              ')',
              sep="",
              file = paste(out_dir, "/",fname, ".R", sep="")
          )
        }
      }
    }
  }
}

