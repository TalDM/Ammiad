library(ggpubr)

source("data_analysis/phenotypes/006_plot_predictions.R")
source("data_analysis/phenotypes/007_linear_discriminant_analysis.R")

ggarrange(plot_pve_habitat, plot_lda, nrow=2, labels= "AUTO")

ggsave(plot = plot_pve,
       filename = "data_analysis/phenotypes/figures/phenotypic_differentiation.jpg",
       device = "jpg",
       width = 16.9, height = 10, units="cm"
)
