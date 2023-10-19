library(tidyverse)
library(magrittr)
library(adegenet)
library(hierfstat)

source("scripts/geoinfo.R")
source("scripts/gstat.randtest.R")

set.seed(12345)
nsim <- 999
nsite <- commandArgs(trailingOnly=TRUE)[1]

s208ad_bed <-
  read.PLINK("plink/s208_admixture.raw")

s189ad_bed <-
  s208ad_bed[!s208ad_bed$pop %in% c("Slu1", "Slu2"),]

keep_subpops <-
  tibble(IID=s189ad_bed$ind.names) %>%
  left_join(select(geoinfo_tb, IID, FID, Locality)) %>%
  group_by(FID, Locality) %>%
  summarise(n=n()) %>%
  filter(n>=5) %>%
  select(FID, Locality)

keep_ind <-
  tibble(IID=s189ad_bed$ind.names) %>%
  left_join(select(geoinfo_tb, IID, FID, Locality)) %>%
  right_join(keep_subpops) %>%
  pull(IID)

keep_ind_bed <-
  s189ad_bed[s189ad_bed$ind.names %in% keep_ind,]

applying_pop <-
  which(table(keep_ind_bed$pop) > 5) %>%
  names() %>%
  set_names(., .)

within_g_bed <-
  seppop(keep_ind_bed)

for (pop in applying_pop) {
  within_g_bed[[pop]]$pop <-
    geoinfo_tb$Locality[match(within_g_bed[[pop]]$ind.names, geoinfo_tb$IID)] %>%
    as.factor()

  #within_g_bed[[pop]] <-
  #  genlight_dropInv(within_g_bed[[pop]])
}

applying_genlights <-
  c(within_g_bed, list(s208=s208ad_bed, s189=s189ad_bed))

gtests_within <-
  map(applying_genlights, ~gstat_wrapper(.x, nsim, nsite))

gtests_within_value  <-
  map_df(gtests_within, ~get_randtest_value(.x), .id = "pop") %>%
  set_colnames(c("pop", "G", "pvalue"))

saveRDS(gtests_within, paste0("Fstats/gtests_within", nsim, "x", nsite, ".rds"))
write_csv(gtests_within_value,
          paste0("Fstats/wGstats_", nsim, "x", nsite, ".csv")
          )
