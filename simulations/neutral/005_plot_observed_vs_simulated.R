library('tidyverse')
library('simmiad')
library('ggpubr')

source("data_processing/import_field_occupancy_data.R")
source('data_processing/ammiad_spatial_clustering.R')
source('data_processing/ammiad_temporal_stability.R')
sims <- readRDS("simulations/neutral/output/simmiad_summary.Rds")

sim_curve <- list(
  distance_identity =  sims$distance_identity$through_time %>% 
    filter(
      mean_dispersal_distance == 1,
      outcrossing_rate == 0.04, 
      density == 3,
      dormancy == 0.3
    ),
  stability = sims$stability$through_time %>% 
    filter(
      mean_dispersal_distance == 1,
      outcrossing_rate == 0.04,
      density == 3,
      dormancy == 0.3,
      generation <=36
    ) %>% 
    mutate(stability = mean)
)

p1 <- ggplot(obs_di, aes(x=distance, y = mean, colour=transect)) + 
  geom_line() +
  geom_point() +
  labs(
    x= "Distance (m)",
    y = "Prob. same DGG"
  ) +
  geom_ribbon(
    data=sim_curve$distance_identity, 
    aes(x=distance, ymin=lower, ymax = upper),
    colour = 'grey', fill="grey", alpha =0.5
  ) +
  scale_colour_manual(values = c("#4285f4", "#ea4335", "#fbbc05", "#34a853")) + 
  theme_bw()

p2 <- ggplot(obs_stability, aes(x=n_years, y = stability, colour=transect)) + 
  geom_line() +
  geom_point() +
  labs(
    x= "Number of years",
    y = "Prob. same DGG"
  ) +
  geom_ribbon(
    data=sim_curve$stability, 
    aes(x=generation, ymin=lower, ymax = upper),
    colour = 'grey', fill="grey", alpha = 0.5
  ) +
  scale_colour_manual(values = c("#4285f4", "#ea4335", "#fbbc05", "#34a853")) +
  theme_bw()
