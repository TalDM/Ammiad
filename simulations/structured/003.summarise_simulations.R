#' Tom Ellis, 15th April 2021
#' 
#' Import simulation results, and average over replicates.
#' Save the output as and RDS file for quick re-use
library("tidyverse")
base_folder <- "simulations/structured/output/"

files <- c(
  #"clustering.csv",
  "count_NAs.csv",
  "distance_identity.csv",
  #"matching_pairs.csv",
  "n_genotypes.csv",
  "stability.csv"
)

sims <- vector('list', length(files))
names(sims) <- files
for(p in files){
  sims[[p]] <- simmiad::summarise_simmiad(base_folder, p)
}
names(sims) <- tools::file_path_sans_ext(names(sims))

# Give distances for distance_identity
sims$distance_identity$through_time$distance <-  seq(5,145, 5)

sims$stability$through_time <- sims$stability$through_time %>% 
  mutate(
    # generation = rep(1:500, 180),
    outcrossing_rate = as.factor(outcrossing_rate),
    Density = as.factor(density),
    transect = "sim"
  )

saveRDS(sims, "simulations/structured/output/simmiad_summary.Rds")
