#' Calculate the probability of observing matching genotypes at sampling points 
#' at increasing distances in the same harvest year, averaged over all pairs of 
#' sampling points at a given distance. Sampling points are grouped into 5m
#' bins.
#' 
#' Then, repeat this for each year separately, and keep the probability of 
#' finding identical genotypes for sampling points within 5m only.
#' 
#' This assumes you have imported and formatted the data using 
#' source("data_processing/import_field_occupancy_data.R").
#' 
#' Returns a tibble `obs_di` ('di' is short for 'distance identity') giving 
#' the spatial bin (in metres) from `cut`, the probability that two sampling 
#' points are occupied by the same DGG (`mean`), number of pairs of sampling
#' that could be compared, an integer value for the bin, and transect.
#' 
#' Tom Ellis


if ( !exists('obs_geno') ){
  stop("Table of observed genotypes not found. Please run `source('data_processing/import_field_occupancy_data.R')` before running this script.")
}

# Empty list with an element for each transect.
obs_di <- vector('list', 4)
names(obs_di) <- LETTERS[1:4]
# Loop over all four transects
for(tr in LETTERS[1:4]){
  d <- obs_geno %>%
    filter(transect==tr)
  # Dataframe giving all pairs of sampling points in the transect
  # Shows how often two sampling points 'matched' (same IGG at both point in a year)
  # This is averaged over years (n=number of years they could be compared)
  di <- distance_identities(
    genotypes = split(d$IGG, d$Year), 
    positions = split(d$dist, d$Year)[[1]]
  ) %>% 
    filter( !is.nan(matches))
  # Bin those match probabilities into 5m bins
  obs_di[[tr]] <- di %>%
    mutate(bin = cut(distances, seq(0,max(distances)+5,5))) %>% 
    group_by(bin) %>% 
    summarise(
      mean=mean(matches),
      n = n(),
      .groups='drop'
    ) %>%
    mutate(
      distance =  (1:nrow(.)) * 5,
      transect = tr
    )
}
# Join the list elements into a single dataframe
obs_di <- do.call('rbind', obs_di)



rm(d, di)


