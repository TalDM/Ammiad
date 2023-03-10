# Calculate the probability of plants of the same DGG in the same sampling 
# point at harvests separated by different number of years for each transect.
#
#' This assumes you have imported and formatted the data using 
#' source("data_processing/import_field_occupancy_data.R").
#' 
#' Returns a tibble `obs_stability` with columns for transect, number of years
#' between harvests being compared, the probability that a sampling point is 
#' occupied by the same DGG in those years (stability), and the number of pairs
#' of sampling point that could be compared.
#' 
#' Tom Ellis

if ( !exists('obs_geno') ){
  stop("Table of observed genotypes not found. Please run source('data_processing/import_field_occupancy_data.R') before running this script.")
}

# Lists of unique sample points in each transect
sampling_pts <- lapply(LETTERS[1:4], function(x) unique(obs_geno$Position[obs_geno$transect == x]))
names(sampling_pts) <- LETTERS[1:4]
sampling_pts <- lapply(sampling_pts, na.exclude)

# Empty list to store transition data for each transect.
stability_all <- vector('list', 4)
names(stability_all) <- LETTERS[1:4]

# Loop over transects to retrieve transitions
for(t in LETTERS[1:4]){
  # empty vector to store transitions for this transect
  this_trans <- vector("list", (length(years) *(length(years)-1)) / 2) 
  # loop over the seven transitions between eight years.
  # For each pair, return the labels of each year, prob. of observing the same IGG in
  # both years, and the number of pairs that could be compared.
  counter <- 1
  for(i in 1:length(years)){
    for (j in 2:length(years)) {
      if(i < j){
        # Get genotype occupancies for year i, and index which are also in sampling_pts
        year_one <- obs_geno %>% 
          filter(transect == t & Year == years[i])
        ix <- match(sampling_pts[[t]],
                    year_one %>% select(Position) %>% pull()
        )
        # Get genotype occupancies for year j, and index which are also in sampling_pts
        year_two <- obs_geno %>% 
          filter(transect == t & Year == years[j])
        jx <- match(sampling_pts[[t]],
                    year_two %>% select(Position) %>% pull()
        )
        # Calculate stability for this transition, and the number of valid comparisons.
        this_trans[[counter]] <- c(
          year_one = years[i],
          year_two = years[j],
          stability = transect_stability(year_one$IGG[ix], year_two$IGG[jx]),
          n = (!is.na(ix) & !is.na(jx)) %>% sum()
        )
        counter <- counter +1
      }
    }
  }
  # Bind list, and add column for transect.
  # Add columns for the transect and number of years being compared
  this_trans <- do.call('rbind', this_trans) %>%
    as.data.frame() %>% 
    add_column(transect = t, .before = 'stability') %>% 
    mutate(
      n_years = year_two - year_one
    )
  # add_column(transition = years[i+1], .before = 'stability')
  # Send to bigger list
  stability_all[[t]] <- this_trans
}
# Bind list and add st error
stability_all <- do.call("rbind", stability_all)
stability_all$sterr <- sqrt(stability_all$stability * (1-stability_all$stability)) / (stability_all$n-1)

obs_stability <- stability_all %>% 
  group_by(transect) %>% 
  mutate(bins = cut(n_years, seq(1,40,4))) %>% 
  group_by(transect, n_years) %>% 
  summarise(
    stability = mean(stability),
    n = n(),
    .groups='drop'
  )

rm(this_trans, year_one, year_two, sampling_pts, stability_all)
