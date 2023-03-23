library('tidyverse')
library('simmiad')

source("data_processing/import_field_occupancy_data.R")
source('data_processing/ammiad_spatial_clustering.R')
source('data_processing/ammiad_temporal_stability.R')

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

simfiles <- Sys.glob("simulations/structured/output/**/di_by_year.csv")
summarise_di_by_year <- function(file){
  m <- read_csv(file, col_names = FALSE, show_col_types = FALSE)
  out <- data.frame(
    transect = "sim",
    start_year = strsplit(file, "_")[[1]][3],
    harvest = 1983 + (1:ncol(m)),
    mean = colMeans(m),
    lower = apply(m, 2, quantile, 0.025),
    upper = apply(m, 2, quantile, 0.975)
  )
  
  out
}
sim_di_by_year <- lapply(simfiles, summarise_di_by_year) %>% 
  do.call(what = 'rbind')

obs_di_by_year %>% 
  ggplot(aes( x = harvest, y = mean, colour = transect, group = transect )) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    data = sim_di_by_year %>% filter(start_year == '1984'), 
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
  theme_bw()
