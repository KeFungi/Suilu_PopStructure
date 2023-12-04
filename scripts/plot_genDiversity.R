library(tidyverse)
library(cowplot)
library(Cairo)
source("scripts/geoinfo.R")
wind_size <- 10000
tajima_Ddf <- read_csv(paste0("Fstats/tajimaD_", wind_size, ".csv"))

tajima_Ddf_input <-
  tajima_Ddf %>%
  filter(!pop %in% c("Slu1", "Slu2")) %>%
  mutate(FID=case_when(
    pop %in% c("NAm_MA") ~ "NAm",
    pop %in% c("CEu_Belgium") ~ "CEu",
    pop %in% c("AU_WA", "AU_other") ~ "AU",
    pop %in% c("SAm_Bariloche") ~ "SAm",
    TRUE ~ pop)
  ) %>%
  left_join(introduce_tb) %>%
  mutate(subpop=case_when(
    pop == "NAm_MA" ~ "North America - MA",
    pop == "AU_WA" ~ "Australia - WA",
    pop == "AU_other" ~ "Australia - other",
    pop == "CEu_Belgium" ~ "CEu - Belgium",
    pop == "SAm_Bariloche" ~ "SAm - Bariloche",
    TRUE ~ FullName
  )
  ) %>%
  rename(Population=FullName) %>%
  mutate(subpop=factor(subpop)) %>%
  mutate(subpop=fct_relevel(subpop, c("Australia", "Australia - WA", "Australia - other")))

div_theme <-
  list(scale_color_manual(values=pop_color_full[unique(tajima_Ddf_input$Population)], drop=TRUE, name="Population"),
       theme_classic(),
       theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       xlab("")
  )

m_pi <-
  tajima_Ddf_input %>%
  filter(subpop == "Central Europe") %>%
  pull(pi) %>%
  median()

m_the <-
  tajima_Ddf_input %>%
  filter(subpop == "Central Europe") %>%
  pull(Wt) %>%
  median()

pi_plot <-
  ggplot(tajima_Ddf_input) +
  geom_boxplot(aes(x=subpop, y=pi, color=Population), outlier.shape = NA, notch=TRUE) +
  geom_hline(yintercept = m_pi, color="red", lty=2) +
  theme_classic() +
  div_theme +
  facet_grid(cols=vars(introduced), scales="free_x", space="free") +
  ylab("Nucleotide Diversity (π)") +
  xlab("")

theta_plot <-
  ggplot(tajima_Ddf_input) +
  geom_boxplot(aes(x=subpop, y=Wt, color=Population), outlier.shape = NA, notch=TRUE) +
  geom_hline(yintercept = m_the, color="red", lty=2) +
  theme_classic() +
  div_theme +
  ylab("Watterson's θ") +
  facet_grid(cols=vars(introduced), scales="free_x", space="free")

D_plot <-
  ggplot(tajima_Ddf_input) +
  geom_boxplot(aes(x=subpop, y=D, color=Population), outlier.shape = NA, notch=TRUE) +
  geom_hline(yintercept = 0, lty=2) +
  theme_classic() +
  div_theme +
  ylab("Tajima's D") +
  facet_grid(cols=vars(introduced), scales="free_x", space="free")

legend_plot <-
  get_legend(theta_plot + theme(legend.key.size = unit(c(0.1, 0.2), "inch")))

final_plot <-
  plot_grid(
    pi_plot,
    theta_plot,
    D_plot,
    ncol=1,
    greedy=FALSE)

ggsave("plots/genDiversity.pdf", final_plot, width=6, height=10, device = cairo_pdf)

ggplot(tajima_Ddf_input) +
  geom_boxplot(aes(x=pop, y=D, color=pop), outlier.shape = NA) +
  geom_hline(yintercept = 0, color="black", lty=2) +
  theme_classic() +
  #ylim(-2, 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
