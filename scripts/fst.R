library(tidyverse)
library(KRIS)

fst_2pop <-
  function(plink_bed, pop1, pop2){
    pop1_pos <- which(plink_bed$ind.info$FamID==pop1)
    pop2_pos <- which(plink_bed$ind.info$FamID==pop2)
    fst.hudson(plink_bed$snp, pop1_pos, pop2_pos)
  }

Dxy_2pop <-
  function(plink_bed, pop1, pop2){
    pop1_pos <- which(plink_bed$ind.info$FamID==pop1)
    pop2_pos <- which(plink_bed$ind.info$FamID==pop2)
    Dxy(plink_bed$snp, pop1_pos, pop2_pos)
  }

fis_pop <-
  function(plink_bed, pop1){
    pop1_pos <- which(plink_bed$ind.info$FamID==pop1)
    fis.hudson(plink_bed$snp, pop1_pos)
  }

fis.hudson <-
  function (X, idx) {
    prestep.fis.one.marker <- function(alleles, idx) {
      G = alleles[idx]
      no.AA = length(which(G == 0))
      no.AB = length(which(G == 1))
      no.BB = length(which(G == 2))
      n = (no.AA + no.AB + no.BB)
      p1 = (no.AA * 2 + no.AB)/(n*2)
      Fis = 1 - (no.AB/n)/(2*(p1-p1^2))
      return(Fis)
    }
    set.fis = apply(X, 2, prestep.fis.one.marker, idx = idx)
    fis = mean(set.fis, na.rm = T)
    return(fis)
  }

Dxy <-
  function (X, idx.p1, idx.p2) {
    prestep.dxy.one.marker <-
      function(alleles, idx.p1, idx.p2) {
        G = alleles[idx.p1]
        no.AA = length(which(G == 0))
        no.AB = length(which(G == 1))
        no.BB = length(which(G == 2))
        n1 = (no.AA + no.AB + no.BB) * 2
        p.A = (no.AA * 2 + no.AB)/n1
        p1 = p.A
        G = alleles[idx.p2]
        no.AA = length(which(G == 0))
        no.AB = length(which(G == 1))
        no.BB = length(which(G == 2))
        n2 = (no.AA + no.AB + no.BB) * 2
        p.A = (no.AA * 2 + no.AB)/n2
        p2 = p.A
        Dxy = p1*(1-p2) + p2*(1-p1)
        return(Dxy)
      }
  set.Dxy = apply(X, 2, prestep.dxy.one.marker, idx.p1 = idx.p1,
                  idx.p2 = idx.p2)
  m_Dxy = mean(set.Dxy, na.rm = T)
  return(m_Dxy)
  }


fst_wrapper<-
  function(plink_bed){
    pop_min5 <-
      plink_bed$ind.info %>%
      group_by(FamID) %>%
      summarise(n=n()) %>%
      filter(n>=5) %>%
      pull(FamID) %>%
      expand.grid(., .)

    fst_vec <-
      map2_dbl(pop_min5$Var1, pop_min5$Var2, ~fst_2pop(plink_bed, .x, .y))

    bind_cols(pop_min5, tibble(fst=fst_vec)) %>%
      mutate(fst=ifelse(Var1==Var2, NA, fst)) %>%
      pivot_wider(names_from = "Var1", values_from = "fst")
  }

dxy_wrapper<-
  function(plink_bed){
    pop_min5 <-
      plink_bed$ind.info %>%
      group_by(FamID) %>%
      summarise(n=n()) %>%
      filter(n>=5) %>%
      pull(FamID) %>%
      expand.grid(., .)

    dxy_vec <-
      map2_dbl(pop_min5$Var1, pop_min5$Var2, ~Dxy_2pop(plink_bed, .x, .y))

    bind_cols(pop_min5, tibble(Dxy=dxy_vec)) %>%
      mutate(Dxy=ifelse(Var1==Var2, NA, Dxy)) %>%
      pivot_wider(names_from = "Var1", values_from = "Dxy")
  }


fis_wrapper <-
  function(plink_bed){
    pop_min5 <-
      plink_bed$ind.info %>%
      group_by(FamID) %>%
      summarise(n=n()) %>%
      filter(n>=5) %>%
      pull(FamID)

    fis_vec <-
      map_dbl(pop_min5, ~fis_pop(plink_bed, .x))

    tibble(FamID=pop_min5) %>%
      bind_cols(tibble(Fis=fis_vec))
  }


s213ad_bed <-
  read.bed("plink/s213_admixture.bed",
           "plink/s213_admixture.bim",
           "plink/s213_admixture.fam",
           only.snp = FALSE)

fst_s213 <-
  fst_wrapper(s213ad_bed)

fis_s213 <-
  fis_wrapper(s213ad_bed)

dxy_s213 <-
  dxy_wrapper(s213ad_bed)

write_csv(fst_s213, "Fstats/s213_fst.csv")
write_csv(fis_s213, "Fstats/s213_fis.csv")
write_csv(dxy_s213, "Fstats/s213_Dxy.csv")

