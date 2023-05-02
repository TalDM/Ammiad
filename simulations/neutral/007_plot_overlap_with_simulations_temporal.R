#' Overlap in spatial clustering between observed data and populations simulated
#' under neutral demography. Plots show the proportion of replicate simulations
#' where pairs sampling points at different distances were more likely to be
#' occupied by plants of the same DGG than in any observed transect. Panels
#' indicate increasing mean seed dispersal distances (left to right) and
#' increasing outcrossing rates (top to bottom). 'Dormancy' indicates the
#' proportion of seeds drawn from the seed bank in each generation, and
#' 'density' the number of plants per square metre.


library('tidyverse')
library('simmiad')
library('ggpubr')

source("data_processing/import_field_occupancy_data.R")
source('data_processing/ammiad_spatial_clustering.R')
source('data_processing/ammiad_temporal_stability.R')

base_folder <- "simulations/neutral/output/"
folders <- list.dirs(base_folder)

# Return a data.frame of parameter values for each simulation.
parameters <- vector("list", length = length(folders)-1)
for(f in 2:length(folders)){ # starts from 2 because element 1 is the parent directory
  p <-  read.csv(
    paste(folders[f], "/parameters.csv", sep=""),
    header = FALSE,
    col.names = c("parameter", "value")
  ) %>%
    filter(
      parameter %in% c(
        'mean_dispersal_distance', 'outcrossing_rate', 'n_generations', 'n_starting_genotypes', 'density', 'dormancy'
      )
    )
  vals <- p$value
  names(vals) <- p$parameter
  parameters[[f]] <- vals
}
parameters <- parameters %>% 
  do.call(what = 'rbind') %>%
  as.data.frame()

# Set up ------------------------------------------------------
  
  # Get the minimum stability observed in any transect for any pair of years.
  min_obs_stability <- obs_stability %>%
  pivot_wider(names_from = transect, values_from = stability) %>% 
  select(A:D) %>% 
  pmap(.f = min) %>% 
  unlist()
# Loop through each simulation output and count how many simulation replicates
# gave a higher probability of observing identical genotypes for each pair of 
# years.
p_stability <- vector("list", length = length(folders)-1)
pb <- txtProgressBar(min = 2, max = length(folders), style = 3)   
for(f in 2:length(folders)){
  # Open the parameter file
  this_sim <- read.csv(
    paste(folders[f], "/stability.csv", sep=""),
    header = FALSE
  )
  this_sim <- this_sim[,seq(2,36,2)] %>% t() # take the first 36 years and transpose
  # How often are simulations greater than observed?
  this_comparison <- min_obs_stability < this_sim
  p_stability[[f-1]] <- rowMeans(this_comparison, na.rm = T)
  setTxtProgressBar(pb, f)
}
close(pb) # Close the connection
# Make it a tidy table
p_stability <- p_stability %>% 
  do.call(what = "rbind") %>% 
  as.data.frame() %>% 
  `colnames<-`(seq(2,36,2)) %>% 
  cbind(parameters) %>% 
  pivot_longer(`2`:`36`, names_to ="n_years", values_to = "prob") %>% 
  mutate(
    n_years = as.numeric(n_years),
    dormancy = as.factor(dormancy),
    outcrossing_rate = as.factor(outcrossing_rate),
    density = as.factor(density)
  )


# New facet label names for dispersal distances
disp.labs <- paste(c(0.5, 0.75, 1, 2),"m",sep="")
names(disp.labs) <- c(0.5, 0.75, 1, 2)
# New facet label names for outcrossing
out.labs <- paste(c(0.005, 0.02, 0.04, 0.08)*100, "%", sep="")
names(out.labs) <- c(0.005, 0.02, 0.04, 0.08)

plot_stability <- p_stability %>% 
  ggplot(aes( x = n_years, y = prob, colour=density, linetype=dormancy  )  ) + 
  geom_line() +
  theme_bw()+ 
  labs(
    x = "Years between harvests",
    y = "Proportion of replicates",
    linetype="Dormancy",
    color ="Density"
  ) +
  scale_color_brewer(palette="Dark2") +
  facet_grid(
    outcrossing_rate ~ mean_dispersal_distance,
    labeller = labeller(outcrossing_rate = out.labs, mean_dispersal_distance = disp.labs)
  ) 

# ggsave(
#   filename = "simulations/neutral/figures/overlap_with_simulations_temporal.png",
#   device = 'png', width = 16.9, units = "cm"
#   )
