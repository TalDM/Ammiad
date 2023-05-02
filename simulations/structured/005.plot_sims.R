library('tidyverse')
library('simmiad')

source("data_processing/import_field_occupancy_data.R")
source('data_processing/ammiad_spatial_clustering.R')
source('data_processing/ammiad_temporal_stability.R')

# Import simulation data
sim_di_by_year <- readRDS("simulations/structured/output/simmiad_summary.Rds")

# Calculate the distance identity by year, and take the probabilities for 
# sampling points within 5m of one another
# Empty list with an element for each transect.
obs_di_by_year <- vector('list', 4)
names(obs_di_by_year) <- LETTERS[1:4]
for(tr in LETTERS[1:4]){
  d <- obs_geno %>% 
    filter(transect == tr) %>% 
    split(.$Year)
  
  di_all_years <- lapply(d, distance_identities, d$`1984`$dist)
  obs_di_by_year[[tr]] <- data.frame(
    data = 'obs',
    transect = tr,
    harvest = as.numeric(names(d)),
    mean = sapply(di_all_years, function(x) {
      j <- x[x$distances <= 5, ]
      sum(j$matches * j$n, na.rm = T) / sum(j$n)
    })
  )
}
obs_di_by_year <- do.call('rbind', obs_di_by_year)

plot_di_by_year <- obs_di_by_year %>% 
  ggplot(aes( x = harvest, y = mean, colour = transect, group = transect )) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    data = sim_di_by_year %>% 
      filter(
        start_year == '1984',
        dispersal == 1, outcrossing == 0.04, density == 3, dormancy == 0.3), 
    aes( ymin = lower, ymax = upper ),
    colour = 'grey', fill="grey", alpha = 0.5
  ) + 
  lims(
    x = c(1984, 2024)
  ) +
  labs(
    x = "Harvest year",
    y = "Prob. same DGG"
  ) +
  scale_colour_manual(values = c("#4285f4", "#ea4335", "#fbbc05", "#34a853")) +
  theme_bw()

# New facet label names for dispersal distances
disp.labs <- paste(c(0.5, 0.75, 1, 2),"m",sep="")
names(disp.labs) <- c(0.5, 0.75, 1, 2)
# New facet label names for outcrossing
out.labs <- paste(c(0.005, 0.02, 0.04, 0.08)*100, "%", sep="")
names(out.labs) <- c(0.005, 0.02, 0.04, 0.08)

sim_di_by_year %>% 
  ggplot(aes( x = harvest, y = mean, colour=as.factor(density), linetype=as.factor(dormancy) )) +
  geom_line() +
  theme_bw()+ 
  labs(
    x = "Harvest year",
    y = "Prob. same DGG",
    linetype="Dormancy",
    color ="Density"
  ) +
  theme(
    axis.text.x =element_text(angle = 45, hjust = 1)
  ) +
  facet_grid(
    outcrossing ~ dispersal,
    labeller = labeller(outcrossing = out.labs, dispersal = disp.labs)
  )

ggsave(
  filename = "simulations/structured/sim_di_by_year.png",
  device = "png", width = 16.9, units = 'cm'
  )

