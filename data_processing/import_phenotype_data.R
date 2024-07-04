#' Import data on phenotypes for plants grown in a common-garden nethouse
#' experiment
library('tidyverse')

# function to convert string dates into integer days since 1st January 2020.
str2date <- function(dates) {
  year <- as.numeric(substr(dates, 1,4)) - 2020
  month <- as.numeric(substr(dates, 5,6))
  day <- as.numeric(substr(dates, 7,8))
  
  n_days <- numeric(length(month))
  n_days[month == 2] <- 28
  n_days[month %in% c(4,6,9,11)] <- 30
  n_days[month %in% c(1,3,5,7,8,10,12)] <- 31
  
  year*365 + month*n_days + day
}

# Data file with eleven vegetative phenotypes
pheno <- read_csv("data/raw_data/vegetative_means.csv", show_col_types = FALSE)

#' DGG calling was updated during the experiment with a lower threshold for
#' identity. Based on the new calls, three pairs of old DGGs were called as a
#' single DGG, meaning they are overrepresented. Remove one set of replicates
#' and remove the underscores.
pheno <- pheno %>%
  filter( !grepl("_2", .$dgg), dgg != "zavitan") %>%
  mutate(
    dgg = str_replace(dgg, "_[1,2]$", "")
  )

# Format those phenotypes
pheno <- pheno %>%
  # Change dates to usable integers and coleptile colour to binary
  mutate(
    germination_date = ifelse( grepl("trans", germination_date), NA, germination_date ), # Remove nonsense entries
    # Tidy dates
    germination_date = str2date(germination_date),
    anthesis_date = str2date(anthesis_date),
    # Binarise coleptile colour
    coleptile_colour = ifelse(coleptile_colour == "P", 1, 0),
    block = as.character(block),
    # improve normality by removing rightward skew
    main_flag_leaf_width = log(main_flag_leaf_width)
  )

# Import data on seed traits
seeds   <- read_csv(
  "data/raw_data/mean_median_var_seed_size.csv", show_col_types = FALSE
) %>% 
  filter( !grepl("ZAVITAN", X1), !grepl("N2", X1) ) %>% 
  mutate(
    block = substr(X1, 6,6),
    dgg = substr(X1, 7,10)
  ) %>% 
  rename(seed_area = Area...2, seed_length = Length...3, seed_width = Width...4) %>% 
  dplyr::select(block, dgg, seed_area, seed_length, seed_width)

# Join phenotype and block files.
pheno <- pheno %>% 
  left_join(seeds, by =c("block", "dgg"))

pheno %>% 
  write_csv("data/dataset_S2.csv")

rm(seeds)