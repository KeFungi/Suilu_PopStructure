library(adegenet)
library(poppr)
library(dartR)
library(tidyverse)

gl <-
  read.PLINK("plink/s213_admixture.raw")

gi <-
  gl2gi(gl)

gif <- gi[!gi$pop %in% c("AF"),]

genpop <-
  genind2genpop(gif)

pop_dist <-
  dist.genpop(genpop, method=3)

pop_dist %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  write_csv("Fstats/reynold.csv")
