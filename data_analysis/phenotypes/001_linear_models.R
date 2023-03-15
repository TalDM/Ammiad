#' Script to fit linear models in BRMS to estimate genetic values for each of
#' 14 traits estimated in a nethouse experiment. 
#' 
#' Plants were grown a six randomised 'blocks' 96 genotypes (95 DGGs plus the 
#' Zavitan reference genotype) randomised over two rows per block on the same 
#' irrigation line.

source("data_processing/import_field_occupancy_data.R")
source("data_processing/import_phenotype_data.R")

library('brms')

# There are 14 traits with different error distributions
# Continuous variables that we can use a Gaussian link function for
gaussian_vars <- c(
  "germination_date", "anthesis_date", "main_length", "main_flag_leaf_width", "weight_10_seeds",
  "seed_area", "seed_length", "seed_width")
# Poisson variables
poisson_vars <- c("main_number_of_spikes", "number_of_tillers")
# Variables on an ordinal scale
ordinal_vars <- c("erect_1_5", "pigment_1_3", "spike_angle_1_3")
# Centre continuous variables and scale by standard deviations
pheno <- pheno %>% 
  mutate(
    across(
      all_of(gaussian_vars),
      ~scale(., center = T, scale = T)
    )
  )

# Folder to save the models
dir.create("data_analysis/phenotypes/brms_fits", showWarnings = FALSE)
# Set up the CPU
chains <- 4
cores <- 4
iter <- 6000


# Empty list to hold model fits
model_fits <- vector('list', length = length(c(gaussian_vars, poisson_vars, ordinal_vars)) + 1)
names(model_fits) <- c(gaussian_vars, poisson_vars, ordinal_vars, "coleptile_color")

# Fit models for normally-distributed variables.
for(var in gaussian_vars){
  cat("Fitting the BRMS model for ", var, "... ", sep="")
  t0 <- Sys.time()
  # Regression formula
  frmula <- formula(paste(var, "~ 1 + (1 | block) + (1 | new_dgg)"))
  # Path to save the model output
  model_path <- paste0("data_analysis/phenotypes/brms_fits/", var, "_fit.rds")
  # Fit the model
  model_fits[[var]] <- brm(
    formula = frmula, 
    data = pheno,
    family = gaussian(),
    chains = chains, cores = cores, iter = iter,
    control = list(adapt_delta = 0.97, max_treedepth = 12), seed=8,
    file= model_path
  )
  t1 <- Sys.time()
  cat("Completed in", round(as.numeric(t1-t0)/60, 3), "minutes.\n\n")
}

# Fit models for poisson variables
for(var in poisson_vars){
  cat("Fitting the BRMS model for", var, ".\n")
  t0 <- Sys.time()
  # Regression formula
  frmula <- formula(paste(var, "~ 1 + (1 | block) + (1 | new_dgg)"))
  # Path to save the model output
  model_path <- paste0("data_analysis/phenotypes/brms_fits/", var, "_fit.rds")
  # Fit the model
  model_fits[[var]] <- brm(
    formula = frmula,
    data = pheno,
    family = poisson,
    chains = chains, cores = 16, iter = iter,
    file= model_path,
    control = list(adapt_delta = 0.97, max_treedepth = 12), seed=8
  )
  t1 <- Sys.time()
  cat("Completed in", round(as.numeric(t1-t0)/60, 3), "minutes.\n\n")
}


# Fit the models for ordinal variables
for(var in ordinal_vars){
  cat("Fitting the BRMS model for", var, ".\n")
  t0 <- Sys.time()
  # Regression formula
  frmula <- formula(paste(var, "~ 1 + (1 | block) + (1 | new_dgg)"))
  # Path to save the model output
  model_path <- paste0("data_analysis/phenotypes/brms_fits/", var, "_fit.rds")
  # Fit the model
  model_fits[[var]] <- brm(
    formula = frmula,
    data = pheno,
    family=cumulative("logit"),
    chains = chains, cores = 16, iter = iter,
    control = list(adapt_delta = 0.97, max_treedepth = 12), seed=8,
    file= model_path
  )
  t1 <- Sys.time()
  cat("Completed in", round(as.numeric(t1-t0)/60, 3), "minutes.\n\n")
}

# Fit the model for the single binary phenotype, coleptile color
var <- "coleptile_color"
t0 <- Sys.time()
# Regression formula
frmula <- formula(paste(var, "~ 1 + (1 | block) + (1 | new_dgg)"))
# Path to save the model output
model_path <- paste0("data_analysis/phenotypes/brms_fits/", var, "_fit.rds")
# Fit the model
cat("Fitting the BRMS model for coleptile colour.")
model_fits[[var]] <- brm(
  formula = frmula,
  data = pheno,
  family=bernoulli,
  chains = chains, cores = 16, iter = iter,
  control = list(adapt_delta = 0.97, max_treedepth = 12), seed=8,
  file= model_path
)
t1 <- Sys.time()
cat("Completed in", round(as.numeric(t1-t0)/60, 3), "minutes.\n\n")
