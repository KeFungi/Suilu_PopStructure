library(cowplot)
library(adegenet)
source("scripts/geoinfo.R")

mapWorld <- borders("world", colour="gray80", fill="gray80")

PCA_theme_subpop <-
  list(scale_color_manual(values=unname(pop_color), drop=TRUE),
       theme_classic(),
       coord_fixed()
  )

map_theme <-
  list(scale_color_manual(values=unname(pop_color), drop=TRUE),
       coord_map(),
       xlab("longitude"),
       ylab("latitude")
  )

pca_list <- readRDS("pca/pca_list.rds")

CEu_pca_plot <-
  pca_list$CEu %>%
  left_join(geoinfo_tb) %>%
  mutate(Locality=fct_infreq(Locality)) %>%
  arrange(Locality) %>%
  ggplot(aes(x=PC1, y=PC2, color=Locality, group=IID)) +
  geom_point() +
  PCA_theme_subpop +
  xlab(paste0("PC1 (", round(attr(pca_list$CEu, "per1"), 2), "%)")) +
  ylab(paste0("PC2 (", round(attr(pca_list$CEu, "per2"), 2), "%)"))

CEu_point <-
  pca_list$CEu %>%
  left_join(geoinfo_tb) %>%
  mutate(Locality=fct_infreq(Locality)) %>%
  arrange(Locality) %>%
  select(Longitude, Latitude, Locality) %>%
  unique() %>%
  ggplot() +
  mapWorld +
  geom_jitter(aes(x=Longitude, y=Latitude, color=Locality)) +
  xlim(-10, 100) +
  ylim(30, 75) +
  map_theme

SAm_pca_plot <-
  pca_list$SAm %>%
  left_join(geoinfo_tb) %>%
  ggplot(aes(x=PC1, y=PC2, color=Locality, group=IID)) +
  geom_point() +
  PCA_theme_subpop +
  xlab(paste0("PC1 (", round(attr(pca_list$SAm, "per1"), 2), "%)")) +
  ylab(paste0("PC2 (", round(attr(pca_list$SAm, "per2"), 2), "%)"))


SAm_point <-
  pca_list$SAm %>%
  left_join(geoinfo_tb) %>%
  select(Longitude, Latitude, Locality) %>%
  unique() %>%
  ggplot() +
  mapWorld +
  geom_point(aes(x=Longitude, y=Latitude, color=Locality)) +
  xlim(-80, -60) +
  ylim(-55, -35) +
  map_theme

NZ_pca_plot <-
  pca_list$NZ %>%
  left_join(geoinfo_tb) %>%
  ggplot(aes(x=PC1, y=PC2, color=Locality, group=IID)) +
  geom_point() +
  PCA_theme_subpop +
  xlab(paste0("PC1 (", round(attr(pca_list$NZ, "per1"), 2), "%)")) +
  ylab(paste0("PC2 (", round(attr(pca_list$NZ, "per2"), 2), "%)"))


NZ_point <-
  pca_list$NZ %>%
  left_join(geoinfo_tb) %>%
  select(Longitude, Latitude, IID, Locality) %>%
  unique() %>%
  ggplot() +
  mapWorld +
  geom_jitter(aes(x=Longitude, y=Latitude, color=Locality)) +
  xlim(165, 180) +
  ylim(-47.5, -34) +
  map_theme

AU_pca_plot <-
  pca_list$AU %>%
  left_join(geoinfo_tb) %>%
  ggplot(aes(x=PC1, y=PC2, color=Locality, group=IID)) +
  geom_point() +
  PCA_theme_subpop +
  xlab(paste0("PC1 (", round(attr(pca_list$AU, "per1"), 2), "%)")) +
  ylab(paste0("PC2 (", round(attr(pca_list$AU, "per2"), 2), "%)"))


AU_point <-
  pca_list$AU %>%
  left_join(geoinfo_tb) %>%
  select(Longitude, Latitude, Locality) %>%
  unique() %>%
  ggplot() +
  mapWorld +
  geom_jitter(aes(x=Longitude, y=Latitude, color=Locality)) +
  xlim(110, 160) +
  ylim(-42.5, -10) +
  map_theme

NAm_pca_plot <-
  pca_list$NAm %>%
  left_join(geoinfo_tb) %>%
  mutate(Locality=fct_infreq(Locality)) %>%
  arrange(Locality) %>%
  ggplot(aes(x=PC1, y=PC2, color=Locality, group=IID)) +
  geom_point() +
  PCA_theme_subpop +
  xlab(paste0("PC1 (", round(attr(pca_list$NAm, "per1"), 2), "%)")) +
  ylab(paste0("PC2 (", round(attr(pca_list$NAm, "per2"), 2), "%)"))


NAm_point <-
  pca_list$NAm %>%
  left_join(geoinfo_tb) %>%
  select(Longitude, Latitude, Locality) %>%
  unique() %>%
  arrange(Locality) %>%
  ggplot() +
  mapWorld +
  geom_jitter(aes(x=Longitude, y=Latitude, color=Locality)) +
  xlim(-160, -30) +
  ylim(10, 75) +
  map_theme

gtests <- readRDS("Fstats/gtests_within999x10000.rds")
ibd_plots <- readRDS("ibd/ibd_plots.rds")

plot_list <-
  list(CEu_point+theme(legend.position = "none"), get_legend(CEu_point), CEu_pca_plot+theme(legend.position = "none"), ~plot(gtests$CEu, main=NULL, xlab="G", sub=paste0("p=", gtests$CEu$pvalue)), ibd_plots$CEu+labs(title = ""),
       NZ_point+theme(legend.position = "none"), get_legend(NZ_point), NZ_pca_plot+theme(legend.position = "none"), ~plot(gtests$NZ, main=NULL, xlab="G", sub=paste0("p=", gtests$NZ$pvalue)), ibd_plots$NZ+labs(title = ""),
       AU_point+theme(legend.position = "none"), get_legend(AU_point), AU_pca_plot+theme(legend.position = "none"), ~plot(gtests$AU, main=NULL, xlab="G", sub=paste0("p=", gtests$AU$pvalue)), ibd_plots$AU+labs(title = ""),
       NAm_point+theme(legend.position = "none"), get_legend(NAm_point), NAm_pca_plot+theme(legend.position = "none"), ~plot(gtests$NAm, main=NULL, xlab="G", sub=paste0("p=", gtests$NAm$pvalue)), ibd_plots$NAm+labs(title = ""),
       SAm_point+theme(legend.position = "none"), get_legend(SAm_point), SAm_pca_plot+theme(legend.position = "none"), ~plot(gtests$SAm, main=NULL, xlab="G", sub=paste0("p=", gtests$SAm$pvalue)), ibd_plots$SAm+labs(title = "")
       )

final_plot <-
  plot_grid(plotlist = plot_list, ncol = 5, rel_widths=c(1,0.5, 1, 1, 1), greedy = FALSE, labels=c("A", "", "B", "C", "D"))

ggsave2("plots/substructures.pdf", final_plot, width = 16, height=16)
ggsave2("plots/AUnoWA_sub.pdf", ibd_plots$AU_noWA , width = 4, height=4)

