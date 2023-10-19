library(tidyverse)
library(plotly)
library(cowplot)
library(magrittr)
library(adegenet)
library(hierfstat)
library(KRIS)

source("scripts/geoinfo.R")
source("scripts/gstat.randtest.R")

PCA_theme <-
  list(scale_fill_manual(values=pop_color_full[-1], drop=TRUE, name="Population"),
       theme_classic(),
       coord_fixed(clip = "off")
  )

genlight_pca <-
  function(genlight) {
    pca <-
      dudi.pca(tab(genlight), scale = FALSE, scannf=FALSE, nf=2)

    per <-
      100 * pca$eig/sum(pca$eig)

    pca_tb <-
      pca$li %>%
      rownames_to_column("IID") %>%
      rename(PC1=Axis1, PC2=Axis2) %>%
      bind_cols(tibble(FID=genlight$pop))

    attr(pca_tb, "per1") <- per[1]
    attr(pca_tb, "per2") <- per[2]
    return(pca_tb)
  }

s208ad_bed <-
  read.PLINK("plink/s208_admixture.raw")

s189ad_bed <-
  s208ad_bed[!s208ad_bed$pop %in% c("Slu1", "Slu2"),]

s189ad_pop <-
  seppop(s189ad_bed)

s208ad_pca <-
  genlight_pca(s208ad_bed)

s208ad_pca_plot <-
  s208ad_pca %>%
  left_join(introduce_tb) %>%
  mutate(Australia=
           ifelse(IID %in% filter(geoinfo_tb, FID=="AU", Locality=="WA")$IID,
                  'Western Australia (WA)',
                  'other Australia'
           )
  ) %>%
  arrange(match(FID, c("AF", "NAm"))) %>%
  arrange(-row_number()) %>%
  ggplot(aes(x=PC1, y=PC2, fill=FullName, color=Australia)) +
  geom_point(shape=21)+
  PCA_theme +
  scale_color_manual(values=c("Western Australia (WA)"="red"), na.value = "transparent") +
  xlab(paste0("PC1 (", round(attr(s208ad_pca, "per1"), 2), "%)")) +
  ylab(paste0("PC2 (", round(attr(s208ad_pca, "per2"), 2), "%)"))

s189ad_pca <-
  genlight_pca(s189ad_bed)

s189ad_pca_plot <-
  s189ad_pca %>%
  left_join(introduce_tb) %>%
  mutate(Australia=
           ifelse(IID %in% filter(geoinfo_tb, FID=="AU", Locality=="WA")$IID,
                  'Western Australia (WA)',
                  'other Australia'
           )
  ) %>%
  sample_n(nrow(.)) %>%
  arrange(match(FID, "AF")) %>%
  arrange(-row_number()) %>%
  ggplot(aes(x=PC1, y=PC2, fill=FullName, color=Australia)) +
  geom_point(shape=21) +
  PCA_theme +
  scale_color_manual(values=c("Western Australia (WA)"="red"), na.value = "transparent") +
  xlab(paste0("PC1 (", round(attr(s189ad_pca, "per1"), 2), "%)")) +
  ylab(paste0("PC2 (", round(attr(s189ad_pca, "per2"), 2), "%)"))

pca1_legend <-
  get_legend(s189ad_pca_plot + theme(legend.key.size = unit(c(0.1, 0.1), "inch")))

PCA_plot1 <-
  plot_grid(s208ad_pca_plot+theme(legend.position = "none"),
            ggplot()+theme_nothing(),
            plot_grid(pca1_legend),
            s189ad_pca_plot+theme(legend.position = "none"),
            ggplot()+theme_nothing(),
            ggplot()+theme_nothing(),
            scale = c(1, 0.9, 1,
                      1, 0.9, 1),
            rel_widths = c(1, 0.4, 0.1,
                           1, 0.4, 0.1),
            labels=c("A", "", "",
                     "B", "", ""),
            nrow=2, ncol=3,
            greedy=FALSE)

ggsave("plots/PCA.pdf", PCA_plot1, height = 6, width = 6)

applying_pop <-
  which(table(s189ad_bed$pop) > 5) %>%
  names() %>%
  set_names(., .)

pca_list <-
  map(applying_pop, ~genlight_pca(s189ad_pop[[.x]]))

saveRDS(pca_list, "pca/pca_list.rds")
