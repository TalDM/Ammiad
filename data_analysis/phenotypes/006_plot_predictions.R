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
      mean  = sapply(pred_sd, function(x) mean(x[,col])),
      lower = sapply(pred_sd, function(x) quantile(x[,col], 0.025)),
      upper = sapply(pred_sd, function(x) quantile(x[,col], 0.975))
    ),
    # Rotated
    data.frame(
      trait = names(model_fits),
      Data  = "Rotated",
      mean  = sapply(null_sd, function(x) mean(x[,col])),
      lower = sapply(null_sd, function(x) quantile(x[,col], 0.025)),
      upper = sapply(null_sd, function(x) quantile(x[,col], 0.975))
    )
  ) %>% 
    arrange(Data, trait) %>% 
    mutate(trait  = rep(c(
      "Anthesis date", "Coleoptile colour", "Erect growth habit", "Germination date",
      "Main flag-leaf width", "Main tiller length", "Main-tiller spikelet number",
      "Tiller number", "Plant-base pigmentation", "Seed area", "Seed length",
      "Seed width", "Main-tiller spike angle", "Seed weight"
    ), 2)
    ) 
}

# plot_pve_year <- prepare_predictions_for_plotting(1) %>% 
#   ggplot( aes( y = trait, x = mean, colour = Data )) +
#   geom_errorbar( aes(xmin = lower, xmax=upper), position= pd, width=0 ) +
#   geom_point(position=pd) +
#   labs(
#     x = ""
#   ) +
#   lims(
#     x = c(0, 40)
#   ) +
#   theme_bw() +
#   scale_colour_manual(values = c("#4285f4", "#ea4335", "#fbbc05", "#34a853")) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank()
#   )

plot_pve_habitat <- prepare_predictions_for_plotting(2) %>%
  ggplot( aes( x = trait, y = mean, colour = Data )) +
  geom_errorbar( aes(ymin = lower, ymax=upper), position= pd, width=0 ) +
  geom_point(position=pd) +
  labs(
    x = "",
    y = "Variance explained by habitat"
  ) +
  theme_bw() +
  scale_colour_manual(values = c("#4285f4", "#ea4335", "#fbbc05", "#34a853")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'bottom'
  ) +
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
  

plot_pve <- ggpubr::ggarrange(
  plot_heritability, plot_pve_habitat,
  nrow = 2, ncol = 1,
  common.legend = FALSE,
  labels = "AUTO",
  # hjust = c(-12, -1),
  # widths = c(1.3, 1),
  heights = c(0.6, 1)
  )

# plot_pve <- annotate_figure(
#   plot_pve,
#   bottom = text_grob("Percent variance explained by habitat")
# )

plot_pve

ggsave(plot = plot_pve,
  filename = "data_analysis/phenotypes/figures/phenotypic_differentiation.pdf",
  device = "pdf",
  width = 12, height = 20, units="cm"
)
