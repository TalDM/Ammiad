

library(tidyverse)
library(ggpubr)

mean_phenotypes <- read_csv(file = "data_analysis/phenotypes/mean_phenotypes_for_lda.csv", col_types = cols())
# Remove phenotypes with very low heritability
mean_phenotypes <- mean_phenotypes %>%
  dplyr::select(-germination_date, -main_number_of_spikes, -number_of_tillers)

lda_ammiad <- lda(habitat ~ ., data = mean_phenotypes)

lda_for_plotting <- tibble(
  habitat = predict(lda_ammiad)$class,
  as.data.frame(predict(lda_ammiad)$x)
) %>% 
  mutate(
    habitat = factor(
      habitat,
      levels = c(
        'Upper_N_facing_slope',
        'Middle_N_facing_slope',
        'Narrow_valley',
        'S_facing_slope',
        'Plateau',
        'E_facing_slope',
        'Valley',
        'Karst'
        )
    ),
    habitat = recode(
      habitat,
      E_facing_slope = "East-facing slope",
      Karst = "Karst",
      Middle_N_facing_slope = "North-facing slope (mid.)",
      Narrow_valley = "Narrow valley",
      Plateau = "Plateau",
      S_facing_slope = "South-facing slope",
      Upper_N_facing_slope = "North-facing slope (upp.)",
      Valley = "Valley"
    )
  )

ld_12 <- lda_for_plotting %>% 
  ggplot(aes(LD1, LD2, colour = habitat )) +
  geom_point() +
  theme_bw() +
  scale_colour_manual(
    values = c("#aa0000", "#ff5555", "#ff9955", "#ffe680",
               "#5fd35f","#2a7fff", "#9955ff", "#f7a8b8")) 

ld_23 <- lda_for_plotting %>% 
  ggplot(aes(LD2, LD3, colour = habitat )) +
  geom_point() +
  theme_bw() +
  scale_colour_manual(
    values = c("#aa0000", "#ff5555", "#ff9955", "#ffe680",
               "#5fd35f","#2a7fff", "#9955ff", "#f7a8b8")) 

plot_lda <- ggarrange(ld_12, ld_23, ncol = 2, common.legend = TRUE)

lda_ammiad$scaling %>% 
  as.data.frame() %>% 
  write_csv("data_analysis/phenotypes/ld_loadings.csv")

