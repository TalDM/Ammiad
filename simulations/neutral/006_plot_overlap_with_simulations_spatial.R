#' Tom Ellis, May 2021
#' Calculate how often simulations overlap with observed data for the four transects.
#' This is done for each pair of distance between sampling points and for each pair of numbers of years.
#' Simulations are compared with the *minimum* clustering/stability among transects for each pair of data,
#' so this reflects overlap with *any* observed transect.

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

# For each distance class (bins of 5m) up to 60m between sampling points,
# calculate the minimum prob of observing the same DGG.
min_obs_clustering <- obs_di %>%
  filter(distance <= 60) %>% 
  select(transect, distance, mean) %>% 
  pivot_wider(names_from = transect, values_from = mean) %>% 
  select(A:D) %>% 
  pmap(.f = min) %>% 
  unlist()
# Loop through each simulation output and count how many simulation replicates
# gave a higher probability of observing identical genotypes for each pair of 
# distances.
p_spatial <- vector("list", length = length(folders)-1 )
pb <- txtProgressBar(min = 2, max = length(folders), style = 3)   
for(f in 2:length(folders)){
  # Open the parameter file
  this_sim <- read.csv(
    paste(folders[f], "/distance_identity.csv", sep=""),
    header = FALSE
    )
  this_sim <- this_sim[,1:12] %>% t() # take only distances <=60 metres and transpose
  # How often are simulations greater than observed?
  this_comparison <- min_obs_clustering < this_sim
  p_spatial[[f]] <- rowMeans(this_comparison, na.rm = T)
  
  setTxtProgressBar(pb, f)
}
close(pb) # Close the connection
p_spatial <- p_spatial %>% 
  do.call(what = "rbind") %>% 
  as.data.frame() %>% 
  `colnames<-`(seq(5,60,5)) %>% 
  cbind(parameters) %>% 
  pivot_longer(`5`:`60`, names_to ="distance", values_to = "prob") %>% 
  mutate(
    distance = as.numeric(distance),
    dormancy = as.factor(dormancy),
    outcrossing_rate = as.factor(outcrossing_rate),
    density = as.factor(density)
  )


# Plot --------------------------------------------------------------------


# New facet label names for dispersal distances
disp.labs <- paste(c(0.5, 0.75, 1, 2),"m",sep="")
names(disp.labs) <- c(0.5, 0.75, 1, 2)
# New facet label names for outcrossing
out.labs <- paste(c(0.005, 0.02, 0.04, 0.08)*100, "%", sep="")
names(out.labs) <- c(0.005, 0.02, 0.04, 0.08)

plot_spatial <- p_spatial %>% 
  ggplot(aes(x = distance, y = prob, colour=density, linetype=dormancy  )  ) +
  geom_line() +
  geom_point() + 
  theme_bw()+ 
  labs(
    x = "Distance between sampling points",
    y = "Proportion replicates",
    linetype="Dormancy",
    color ="Density"
  ) +
  scale_color_brewer(palette="Dark2") +
  facet_grid(
    outcrossing_rate ~ mean_dispersal_distance,
    labeller = labeller(outcrossing_rate = out.labs, mean_dispersal_distance = disp.labs)
  )

# ggsave(
#   filename = "simulations/neutral/figures/overlap_with_simulations_spatial.png",
#   device = 'png', width = 16.9, units = "cm"
# )
