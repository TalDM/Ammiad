#' Predict phenotypes in the field
#' 
#' This script uses the fitted models to generate posterior draws for plants in
#' the real transects, and estimates the variance explained by habitat and year
#' for each trait. It then rotates the vector of phenotypes
#' 
#'

library(lme4)
library(tidyverse)
source("data_analysis/habitat_permutations/permutation_functions.R")

# Linear models were fitted in this script and saved. Rerunning the script
# imports them.
# Import the saved models, or return an error if they don't exist
if(length(Sys.glob("data_analysis/phenotypes/brms_fits/*rds")) == 14){
  source("data_analysis/phenotypes/001_linear_models.R")
} else{
  stop("Model fits not found. Please run data_analysis/phenotypes/001_linear_models.R first.")
}

# Create a data frame to pass to posterior_predict.
# After removing duplicates, we are left with 801 samples
# Need to merge with DGG names from pheno, because the field data uses old names.
new_data <- pheno %>% 
  dplyr::select(old_dgg, new_dgg) %>% 
  distinct() %>% 
  right_join(
    obs_geno, by = c("old_dgg" = "IGG")
  ) %>%
  filter( !is.na(new_dgg) ) %>%
  arrange(Year, transect, dist) %>% 
  mutate(block = sample(1:4, 1)) # we need to pretend there is a block so `predict` will work

ndraws = 1000 # Number of posterior draws for each plant.
# empty list for results for each trait
sim_data <- vector('list', length(model_fits))
names(sim_data) <- names(model_fits)


# For each trait, generate a ndraws x 801 matrix of posterior draws.
# Predictions are on the link scale, which means we can use gaussian mixed models
# for everything later
for(m in names(model_fits) ){
  mod <- model_fits[[m]]
  # Simulate the data
  sim_data[[m]] <- posterior_predict(
    mod,
    newdata = new_data,
    ndraws = ndraws,
    summary = FALSE,
    scale = "linear",
    re_formula = ~ (1 | new_dgg), # Don't simulate block effects
    # Allows for genotypes in the field that were not in the common garden
    allow_new_levels = TRUE,
    # Phenotypes for new genotypes are drawn from the Gaussian dist of other genotypes
    sample_new_levels = "gaussian"
  )
  # Centre and scale each set of draws.
  sim_data[[m]] <- apply(sim_data[[m]], 1, scale)
}


#' Estimate variance explained by Habitat for each trait separately
#' Permute habitat positions within each transect, within each year, and reestimate variance explained.

#' Empty list to store model output for each trait.
#' Model output is a dataframe of var explained for Year, Habitat and residuals
pred_sd <- null_sd <- vector('list', length(model_fits))
names(pred_sd) <- names(null_sd) <- names(model_fits)

for(m in names(model_fits)){
  # Empty arrays to store results for each model
  # Row for each iteration, and columns for Year, Habitat and residuals
  pred_sd[[m]] <- array(NA, dim = c(ndraws, 3))
  null_sd[[m]] <- array(NA, dim = c(ndraws, 3))
  
  # pb <- txtProgressBar(min = 0, max = ndraws, style = 3)
  for( i in 1:ndraws){
    # setTxtProgressBar(pb,i)
    # Bind data from a single iteration to `new_data`
    glm_data <- cbind( new_data, z = sim_data[[m]][,i])
    # Add an extra column which permutes phenotypes within each transect within each year
    glm_data$z_perm <- glm_data %>%
      split(paste(.$Year)) %>% 
      sapply(., function(x) rotate_vector(x$z, sample(1:nrow(x), size = 1))) %>% 
      # sapply(. , function(x) x[sample(1:nrow(x), size = nrow(x), replace = FALSE), 'z'] ) %>%
      unlist()
    # Fit GLMs for simulated and permuted data
    suppressMessages({
      pred_fit <- lmer( z ~      (1 | Habitat) + (1 | Year), data = glm_data)
      null_fit <- lmer( z_perm ~ (1 | Habitat) + (1 | Year), data = glm_data)
    })
    # Pull out variance components
    pred_sd[[m]][i,] <- as.data.frame(VarCorr(pred_fit))$sdcor
    null_sd[[m]][i,] <- as.data.frame(VarCorr(null_fit))$sdcor
  }
  # close(pb)
}

# Normalise by row sums to get proportion variance explained.
pred_sd <- lapply(pred_sd, function(x) x / rowSums(x) )
null_sd <- lapply(null_sd, function(x) x / rowSums(x) )

# Save the output for later
saveRDS(pred_sd, file = "data_analysis/phenotypes/pred_sd.rds")
saveRDS(null_sd, file = "data_analysis/phenotypes/null_sd.rds")


