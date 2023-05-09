#' 

library("tidyverse")
library("ggpubr")

source("data_analysis/phenotypes/001_linear_models.R")

pred_sd <- readRDS("data_analysis/phenotypes/pred_sd.rds")
null_sd <- readRDS("data_analysis/phenotypes/null_sd.rds")

pd <- position_dodge(.5)

#' A helper function to get the means and CIs for variance explained by year (col=1)
#' or habitat (col=2) from the posterior predictions.
prepare_predictions_for_plotting <- function(col){
  # Tidy table of means and CIs for simulated and permuted data
  rbind(
    # predicted data
    data.frame(
      trait = names(model_fits),
      Data  = "Predicted",
      mean  = sapply(pred_sd, function(x) mean(x[,col])) * 100,
      lower = sapply(pred_sd, function(x) quantile(x[,col], 0.025)) * 100,
      upper = sapply(pred_sd, function(x) quantile(x[,col], 0.975)) * 100
    ),
    # Rotated
    data.frame(
      trait = names(model_fits),
      Data  = "Rotated",
      mean  = sapply(null_sd, function(x) mean(x[,col])) * 100,
      lower = sapply(null_sd, function(x) quantile(x[,col], 0.025)) * 100,
      upper = sapply(null_sd, function(x) quantile(x[,col], 0.975)) * 100
    )
  ) %>% 
    arrange(Data, trait) %>% 
    mutate(trait  = rep(c(
      "Anthesis date", "Coleptile colour", "Erect growth habit", "Germination date",
      "Main flag-leaf width", "Main tiller length", "Main-tiller spikelet number",
      "Tiller number", "Plant-base pigmentation", "Seed area", "Seed length",
      "Seed width", "Main-tiller spike angle", "Seed weight"
    ), 2)
    ) 
}

plot_pve_year <- prepare_predictions_for_plotting(1) %>% 
  ggplot( aes( y = trait, x = mean, colour = Data )) +
  geom_errorbar( aes(xmin = lower, xmax=upper), position= pd, width=0 ) +
  geom_point(position=pd) +
  labs(
    x = ""
  ) +
  lims(
    x = c(0, 40)
  ) +
  theme_bw() +
  scale_colour_manual(values = c("#4285f4", "#ea4335", "#fbbc05", "#34a853")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

plot_pve_habitat <- prepare_predictions_for_plotting(2) %>% 
  ggplot( aes( y = trait, x = mean, colour = Data )) +
  geom_errorbar( aes(xmin = lower, xmax=upper), position= pd, width=0 ) +
  geom_point(position=pd) +
  labs(
    x = ""
  ) +
  theme_bw() +
  scale_colour_manual(values = c("#4285f4", "#ea4335", "#fbbc05", "#34a853")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

plot_pve <- ggpubr::ggarrange(
  plot_pve_year, plot_pve_habitat,
  common.legend = TRUE, legend = 'top',
  labels = "AUTO",
  hjust = c(-12, -1),
  vjust = 0.3,
  widths = c(1, 0.625)
)
plot_pve <- annotate_figure(
  plot_pve,
  bottom = text_grob("Percent variance explained by habitat")
)

plot_pve

ggsave(plot = plot_pve,
  filename = "data_analysis/phenotypes/figures/phenotypic_differentiation.jpg",
  device = "jpg",
  width = 16.9, height = 10, units="cm"
)
