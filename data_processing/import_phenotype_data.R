#' Import data on phenotypes for plants grown in a common-garden nethouse
#' experiment
library('tidyverse')
# Data file with eleven vegetative phenotypes
pheno <- read_csv("data/pheno_20_21_data.csv", show_col_types = FALSE)
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
# Format those phenotypes
pheno <- pheno %>% 
  # Change dates to usable integers and coleptile colour to binary
  mutate(
    germination_date = ifelse( grepl("trans", germination_date), NA, germination_date ), # Remove nonsense entries
    # Tidy dates
    germination_date = str2date(germination_date),
    anthesis_date = str2date(anthesis_date),
    # Binarise coleptile colour
    coleptile_color = ifelse(coleptile_color == "P", 1, 0),
    block = as.character(block),
    # improve normality by removing rightward skew
    main_flag_leaf_width = log(main_flag_leaf_width)
  )

# Import data on seed traits
seeds   <- read_csv(
  "data/mean_median_var_seed_size.csv", show_col_types = FALSE
) %>% 
  filter( !grepl("ZAVITAN", X1), !grepl("N2", X1) ) %>% 
  mutate(
    block = substr(X1, 6,6),
    new_dgg = substr(X1, 7,10)
  ) %>% 
  rename(seed_area = Area...2, seed_length = Length...3, seed_width = Width...4) %>% 
  dplyr::select(block, new_dgg, seed_area, seed_length, seed_width)


# Blocks are defined as two rows of plants on the same irrigation line.
# These lines are parallel to each other.
# To allow for variation perpendicular to that, let's also add an extra 'row'
# level.
# Import a map of the experiment
blocks <- read_csv("data/20210110_blocks_common_garden_nethouse.csv")
# Flatten the map, and add eight strata of six plants per block.
blocks <- blocks %>% 
  mutate(
    row = rep(1:8, each = 6)
  ) %>% 
  pivot_longer(block1...1: block6...12, values_to = "new_dgg", names_to = 'block') %>% 
  mutate(
    block = substr(block, 6,6)
  )

# Join phenotype and block files.
pheno <- pheno %>% 
  left_join(seeds, by =c("block", "new_dgg")) %>% 
  left_join(blocks, by =c("block", "new_dgg"))

rm(blocks, seeds)