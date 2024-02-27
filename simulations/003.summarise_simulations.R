#' Tom Ellis, 15th April 2021
#' 
#' Import simulation results, and average over replicates.
#' Save the output as and RDS file for quick re-use
library("tidyverse")

# Where the simulation results are
simdir <- Sys.glob("simulations/output/*/")
# Ignore results for mean dispersal = 10cm because they yield nonsense
simdir <- grep("d010", simdir, invert = TRUE, value=TRUE)

items_to_ignore <- grep("d[0-9]+_o[0-9]+", simdir, invert = TRUE, value = TRUE)
if( length(items_to_ignore) >0 ){
  cat(length(items_to_ignore), "item(s) in the target folder don't seem to simulation results and will be ignored:\n")
  items_to_ignore
  cat("\n\n")
}
# Exclude entries that don't match a pattern
simdir <- grep("d[0-9]+_o[0-9]+", simdir, value = TRUE)

# Function to import results and parameter values, and return a table with a row
# for each generation giving input parameters, plus the probability of identical
# genotypes at adjacent sampling points averaged over repicates, with 95% CIs.
summarise_by_year <- function(indir, file_name = 'di_by_year.csv'){
  path <- paste0(indir, "/", file_name)
  
  if( !file.exists(path) ){
    warning("File not found:", path, "\n")
    return(NULL)
  } else {
  m <- read_csv(
    path, 
    col_names = FALSE, show_col_types = FALSE
  )
  p <- read_csv(
    paste0(indir, "/parameters.csv"), 
    col_names = c('param', 'value'), show_col_types = FALSE
  )
  out <- data.frame(
    transect = "sim",
    start_year = strsplit(indir, "_")[[1]][3],
    harvest = 1983 + (1:ncol(m)),
    dispersal   = p$value[p$param == "mean_dispersal_distance"],
    outcrossing = p$value[p$param == "outcrossing_rate"],
    density     = p$value[p$param == "density"],
    dormancy    = p$value[p$param == "dormancy"],
    mean = colMeans(m, na.rm=TRUE),
    lower = apply(m, 2, quantile, 0.025, na.rm=TRUE),
    upper = apply(m, 2, quantile, 0.975, na.rm=TRUE)
  )
  
  return(out)
  }
}

cat("Summarising di_by_year.csv.\n")
sim_di_by_year <- lapply(simdir, summarise_by_year, file_name = 'di_by_year.csv') %>% 
  do.call(what = 'rbind')

cat("Summarising clustering_by_habitat.csv\n")
sim_habitat_by_year <- lapply(simdir, summarise_by_year, file_name = "clustering_by_habitat.csv") %>% 
  do.call(what = 'rbind')

cat("Saving to disk.\n")
saveRDS(sim_di_by_year,      "simulations/output/sim_di_by_year.rds")
saveRDS(sim_habitat_by_year, "simulations/output/sim_habitat_by_year.rds")




if( file.exists('simulations/output/sim_di_by_year.rds') ){
  cat(
    "File sim_di_by_year.rds created and has size",
    file.size('simulations/output/sim_di_by_year.rds'),
    "\n"
    )
} else {
  cat("File sim_di_by_year.rds not found.\n")
}

if( file.exists('simulations/output/sim_habitat_by_year.rds') ){
  cat(
    "File sim_habitat_by_year.rds created and has size",
    file.size('simulations/output/sim_habitat_by_year.rds'),
    '\n'
    )
} else {
  cat("File sim_habitat_by_year.rds not found.")
}
