# Tom Ellis, March 2023
# A script for generating a series of new R scripts to run simmiad using a
# range of parameter values so that I don't have to create and maintain each
# manually and make copy-paste-errors.
# It will output one R script in the output directory for every
# combination of parameter values.

set.seed(102)

# Directory to save scripts to.
out_dir <- "simulations/simmiad_scripts"
dir.create(out_dir, showWarnings = FALSE)
# Directory to tell the scripts to save the output to.
sim_dir <- "simulations/output/"
dir.create(sim_dir, showWarnings = FALSE)
# Directory to save SLURM messages to.
slurm_dir <- "simulations/slurm"
dir.create(slurm_dir, showWarnings = FALSE)

# Initialise parameters
years <- 1984
dispersal <- c(0.25, 0.5, 1, 2)
outcrossing_rate <- c(0.5, 2, 4, 8) / 100
n_generations <- 40
density <- c(1, 2, 3, 5)
dormancy <- c(0, 0.1, 0.3, 0.5)
n_sample_points <- 30 
sample_spacing <- 5
nsims <- 100
range_limit= 2
years_to_sample <- '1:n_generations'

for(d in dispersal){
  for(o in outcrossing_rate){
    for(y in years){
      for (a in density){
        for(s in dormancy){
          # File root name for this simulation.
          fname <- paste(
            'd',  sprintf('%03d', d*100),
            '_o', sprintf('%04d', o*10000),
            '_',y,
            '_a', sprintf('%02d', a),
            '_s', sprintf('%02d', s*10),
            sep='')
          
          # Stuff at the top of the script
          cat(
            "library('simmiad')\n",
            "set.seed(",sample(1:100,1),")\n\n",
            
            "# Import observed genotypes\n",
            'source("data_processing/import_field_occupancy_data.R") \n\n',
            
            "# Simulation parameters",
            "\ndispersal <-", d,
            "\noutcrossing <-", o,
            "\nstarting_genotypes <- obs_geno$DGG[obs_geno$Year ==", y, "]",
            "\ndensity <-", a,
            "\ndormancy <-", s,
            "\nn_generations <-", n_generations,
            "\nn_sample_points <-", n_sample_points,
            "\nsample_spacing <-", sample_spacing,
            "\nrange_limit <-", range_limit,
            "\nnsims <-", nsims,
            "\nyears_to_sample <-", years_to_sample,
            "\nfilename <- '", fname, "'",
            "\nhabitat_labels <- obs_geno$Habitat[obs_geno$Year ==", y, "]",
            "\n\n",
            
            "# Code to run the simulation\n",
            "output <- simmiad(\n",
            '\tmean_dispersal_distance = dispersal,\n',
            '\toutcrossing_rate = outcrossing,\n',
            '\tn_generations = n_generations,\n',
            '\tn_starting_genotypes = starting_genotypes,\n',
            '\tdensity = density,\n',
            '\tdormancy = dormancy,\n',
            '\tn_sample_points = n_sample_points,\n',
            '\tsample_spacing = sample_spacing,\n',
            '\trange_limit = range_limit,\n',
            '\tnsims = nsims,\n',
            '\tprogress = FALSE,\n',
            '\tyears_to_sample = years_to_sample,\n',
            '\tpop_structure = "hardcoded",\n',
            '\thabitat_labels = habitat_labels\n',
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

