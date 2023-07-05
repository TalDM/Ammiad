library('tidyverse')
library('simmiad')

source("data_processing/import_field_occupancy_data.R")
source('data_processing/ammiad_spatial_clustering.R')
source('data_processing/ammiad_temporal_stability.R')

# Import simulation data
sim_habitat_by_year <- readRDS("simulations/structured/output/sim_habitat_by_year.rds")

d <- obs_geno %>%
  split(.$Year)

outcrossing_values <- unique(sim_habitat_by_year$outcrossing)
dispersal_values   <- unique(sim_habitat_by_year$dispersal)
observed_fst <- vector('list', length(outcrossing_values) * length(dispersal_values))

counter <- 1
for( o in outcrossing_values){
  for( disp in dispersal_values){
    observed_fst[[counter]] <- data.frame(
      outcrossing = o,
      dispersal = disp,
      harvest = as.numeric(names(d)),
      mean = sapply(d, function(x) habitat_clustering (x$DGG, x$Habitat) )
    )
    counter <- counter + 1
  }
}
observed_fst <- do.call('rbind', observed_fst)


# New facet label names for dispersal distances
disp.labs <- paste(c(0.5, 0.75, 1, 2),"m",sep="")
names(disp.labs) <- c(0.5, 0.75, 1, 2)
# New facet label names for outcrossing
out.labs <- paste(c(0.005, 0.02, 0.04, 0.08)*100, "%", sep="")
names(out.labs) <- c(0.005, 0.02, 0.04, 0.08)

sim_habitat_by_year %>% 
  # filter(dormancy == 0.3) %>% 
  ggplot(aes( x = harvest, y = mean, colour=as.factor(density), linetype=as.factor(dormancy) )) +
  geom_line(
    data = observed_fst, aes( x = harvest, y = mean), colour='grey', linetype = 1
  ) +
  # geom_point(
  #   data = observed_fst, aes( x = harvest, y = mean), colour='grey'
  # ) +
  geom_line() +
  theme_bw()+ 
  labs(
    x = "Harvest year",
    y = expression("F"[st]),
    linetype="Dormancy",
    color ="Density"
  ) +
  lims(
    x= c(1980, 2020)
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_grid(
    outcrossing ~ dispersal,
    labeller = labeller(outcrossing = out.labs, dispersal = disp.labs)
  )

ggsave(
  filename = "simulations/structured/structured_simulations.pdf",
  device = "pdf", width = 16.9, units = 'cm'
  )

