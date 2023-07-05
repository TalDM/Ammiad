#' Get broad-sense heritability estimates from the BRMS models using QGglmm.
#' 
#' The code below creates a function to first determine the link function, then
#' prepares the necessary variance components and passes them to QGparams. The
#' function tells QGparams whether to calculate heritability on the link or data
#' scale. This is then plotted.
#' 
#' Note that we only present heritability on the scale of the link function.
#' This is because (1) for ordinal traits QGglmm returns a separate heritability
#' estimate for each level, but (2) for other traits there was little or no
#' change in the heritability between the data and link-function scales (except
#' coleptile colour, which was wildly different). As such we felt presenting 
#' the data-scale heritabilities increased the complexity of the information 
#' substantially without adding much useful information. However, I left the 
#' code to do this here in case anyone is interested, but commented it out.
#' 
#' Tom Ellis

library(QGglmm)

# Linear models were fitted in this script and saved. Rerunning the script
# imports them.
# Import the saved models, or return an error if they don't exist
# If you haven't run the models yet, that will take some time...
if(length(Sys.glob("data_analysis/phenotypes/brms_fits/*rds")) == 14){
  source("data_analysis/phenotypes/001_linear_models.R")
} else{
  stop("Model fits not found. Please run data_analysis/phenotypes/002_linear_models.R first.")
}

# Function to set up and run the heritability calculation
qgparams_brms <- function(brms_model, genotype = "dgg", linkfun=NULL, link_scale_only = FALSE){
  # If linkfun is not given explicitly, try and work it out from the model itself.
  if( is.null(linkfun) ) {
    family <- brms_model$family$family
    link   <- brms_model$family$link
    # Check the link function from the BRMS model.
    if( family == "gaussian" ) {
      linkfun <- "Gaussian"
    } else if( (family == "poisson") & (link == "log") ) {
      linkfun <- "Poisson.log"
    } else if( family == "bernoulli" ) {
      if( link == "logit" ) linkfun <- "binom1.logit"
      if( link == "probit" ) linkfun <- "binom1.probit"
    } else if( family == "binomial" ){
      if( link == "logit"  ) linkfun <- "binomN.logit"
      if( link == "probit" ) linkfun <- "binomN.probit"
    } else if( family == "cumulative" ){
      linkfun <- "ordinal"
    } else {
      stop("Cannot determine link function. See `?QGparams` for details of how to specify this.")
    }
  }
  # Dataframe of the intercept, genetic variance and phenotypic variance
  vc <- VarCorr(brms_model, summary = FALSE)
  mod_params <- data.frame(
    vg = vc[[genotype]]$sd,
    vp = rowSums(sapply(vc, function(x) x$sd)),
    mu = fixef(brms_model, summary = FALSE)
  )
  colnames(mod_params) <- c("vg", "vp", "mu")
  mod_params$h2 <- mod_params$vg / mod_params$vp
  # If heritability is only needed one the link scale, return mod_params
  if( link_scale_only ){
    return( mod_params )
  } else if ( linkfun != "ordinal") {
    var_params <- do.call("rbind", apply(mod_params, 1, function(row){
      QGparams(
        mu = row[["mu"]],
        var.a = row[["vg"]],
        var.p = row[["vp"]],
        model = linkfun,
        verbose = FALSE
      )
    }))
    return( var_params )
  } else{
    return(
      data.frame(mean.obs=NA, var.obs=NA, var.a.obs=NA, h2.obs=NA)
    )
  }
}

# Estimate heritability for each trait.
# The line commented out quantifies heritability on the data scale.
h2 <- list(
  # data = lapply(model_fits, qgparams_brms), # Heritability on the data scale: not run.
  link = lapply(model_fits, qgparams_brms, link_scale_only = TRUE)
  
)

# Plot heritability on the link scale only.
plot_heritability <- data.frame(
  trait = names(model_fits),
  mean = sapply(h2$link, function(x) mean(x$h2)),
  lower = sapply(h2$link, function(x) quantile(x$h2, 0.025)),
  upper = sapply(h2$link, function(x) quantile(x$h2, 0.975))
) %>%
  arrange(trait) %>% 
  mutate(trait_nice  = c(
    "Anthesis date", "Coleptile colour", "Erect growth habit", "Germination date",
    "Main flag-leaf width", "Main tiller length", "Main-tiller spikelet number",
    "Tiller number", "Plant-base pigmentation", "Seed area", "Seed length",
    "Seed width", "Main-tiller spike angle", "Seed weight"
  )) %>% 
  ggplot( aes(y=mean, x = trait_nice)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0) +
  labs(
    y = "Broad-sense heritability"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# ggsave(
#   device = 'jpg', width = 16.9, height = 12, units = 'cm',
#   filename = "data_analysis/phenotypes/figures/heritability.jpg",
# )

# Plot heritability on both the link and data scales.
# Not run, because we thought this was confusing and unnecessary.

# pd <- position_dodge(0.5)
# 
# rbind(
#   data.frame(
#     trait = names(model_fits),
#     scale = "link",
#     mean = sapply(h2$link, function(x) mean(x$h2)),
#     lower = sapply(h2$link, function(x) quantile(x$h2, 0.025)),
#     upper = sapply(h2$link, function(x) quantile(x$h2, 0.975))
#   ),
#   data.frame(
#     trait = names(model_fits),
#     scale = "data",
#     mean = sapply(h2$data, function(x) mean(x$h2)),
#     lower = sapply(h2$data, function(x) quantile(x$h2, 0.025, na.rm=T)),
#     upper = sapply(h2$data, function(x) quantile(x$h2, 0.975, na.rm=T))
#   )
# ) %>%
#   ggplot( aes(x=trait, y = mean, group=scale, colour=scale)) +
#   geom_point(position = pd) +
#   geom_errorbar(aes(ymin=lower, ymax=upper), width=0, position = pd) +
#   labs(
#     y = "Heritability"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )

