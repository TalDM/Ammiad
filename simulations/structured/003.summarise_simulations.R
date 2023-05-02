#' Tom Ellis, 15th April 2021
#' 
#' Import simulation results, and average over replicates.
#' Save the output as and RDS file for quick re-use
library("tidyverse")

# Where the simulation results are
simdir <- Sys.glob("simulations/structured/output/*")
# Function to import results and parameter values, and return a table with a row
# for each generation giving input parameters, plus the probability of identical
# genotypes at adjacent sampling points averaged over repicates, with 95% CIs.
summarise_di_by_year <- function(indir){
  m <- read_csv(
    paste0(indir, "/di_by_year.csv"), 
    col_names = FALSE, show_col_types = FALSE
    )
  p <- read_csv(
    paste0(indir, "/parameters.csv"), 
    col_names = c('param', 'value'), show_col_types = FALSE
  )
  out <- data.frame(
    transect = "sim",
    start_year = strsplit(file, "_")[[1]][3],
    harvest = 1983 + (1:ncol(m)),
    dispersal   = p$value[p$param == "mean_dispersal_distance"],
    outcrossing = p$value[p$param == "outcrossing_rate"],
    density     = p$value[p$param == "density"],
    dormancy    = p$value[p$param == "dormancy"],
    mean = colMeans(m, na.rm=TRUE),
    lower = apply(m, 2, quantile, 0.025, na.rm=TRUE),
    upper = apply(m, 2, quantile, 0.975, na.rm=TRUE)
  )
  
  out
}

sim_di_by_year <- lapply(simdir, summarise_di_by_year) %>% 
  do.call(what = 'rbind')

saveRDS(sim_di_by_year, "simulations/structured/output/simmiad_summary.Rds")
