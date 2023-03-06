# Script to import data on which samples were found in which sampling position
# in each year.
library("tidyverse")
# Import data one which genotypes are found at each sampling point in each year.
years <- c(1984, 1988, 1992, 1996, 2002, 2014, 2016, 2018, 2020)

# Import data on genotype occupancy.
# Filter to give unique ID, sampling point, unique genotype, year, transect and distance
obs_geno <- read_csv(
  "data/output.csv",
  col_types = cols()
) %>%
  dplyr::select(Sample, Position, IGG, Year, Habitat, dist) %>%
  # Zavitan is the ref genome control; TTD is plants from other populations.
  filter(Position != "Zavitan" & Position != "TTD")

# Ensure you have the same number of entries for each transect-year
obs_geno <- obs_geno %>%
  # Column for transect labels
  mutate(
    transect = ifelse(
      test = grepl("West", Position),
      yes = "A West",
      no = substr(Position, 1 ,1)
    )
  ) %>%
  # Join the observed data to a synthetic dataframe that has all possible years and
  # sampling points. If any are missing in reality, these will appear as NA.
  right_join(
    expand.grid(Position= unique(obs_geno$Position), Year=as.character(years)),
    by = c("Position", "Year")
  ) %>%
  # Fill in the missing data for sampling points that were empty in certain years
  group_by(Position) %>%
  fill(dist, transect, Habitat, .direction = 'downup') %>%
  ungroup() %>%
  filter(!is.na(dist)) %>%
  # reorder
  arrange(Year, transect, dist)