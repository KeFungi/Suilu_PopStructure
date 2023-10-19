library(tidyverse)
library(adegenet)

source("scripts/geoinfo.R")

harmonic <- function(n) {
  sum(1 / (1:n))
}

harmonic_sq <- function(n) {
  sum(1 / (1:n)^2)
}

tajimaD.win <-
  function(focalpop_matrix){
    n <- rowSums(!is.na(focalpop_matrix), na.rm = TRUE)
    n_s <- round(mean(n), 0)

    if (n_s == 0) {list(pi=0, Wt=0, D=0)} else {

      n1 <- rowSums(focalpop_matrix, na.rm = TRUE)
      pi <- sum(((n1)*(2*n-n1))/((2*n-1)*n), na.rm=TRUE)

      S <-
        (rowSums(focalpop_matrix, na.rm = TRUE)>=1 & rowSums(focalpop_matrix, na.rm = TRUE)<=(2*n-1)) %>%
        sum(na.rm=TRUE)

      a1 <- map_dbl(2*n_s-1, ~harmonic(.x))
      Wt <- S/a1

      a2 <- map_dbl(2*n_s-1, ~harmonic_sq(.x))
      b1 <- (2*n_s+1)/(3*(2*n_s-1))
      b2 <- 2*(4*n_s^2+2*n_s+3)/(9*2*n_s*(2*n_s-1))
      c1 <- b1 - 1/a1
      c2 <- b2 - (2*n_s+2)/(a1*2*n_s) + a2/a1^2
      e1 <- c1/a1
      e2 <- c2/(a1^2+a2)

      D <- (pi- Wt)/(e1*S+e2*S*(S-1))^0.5

      return(list(pi=pi, Wt=Wt, D=D))
    }
  }

cut_rows <-
  function(matrix, size){
    nrow <- dim(matrix)[1]
    nblocks <- ifelse(nrow %% size ==0, nrow/size, nrow%/%size + 1)
    starts <- size*(1L:nblocks) - (size-1)
    ends <- size*(1L:nblocks)
    ends[nblocks] <- nrow
    out <- list()
    for (i in 1L:nblocks){
      out <- c(out, list(matrix[starts[i]:ends[i],]))
    }
    return(out)
  }

wind_size <- as.numeric(commandArgs(trailingOnly=TRUE)[1])

s208con_pop <-
  readRDS("s208con_pop.rds")

CEu_Belgium_ids <-
  geoinfo_tb %>%
  filter(FID=="CEu", Locality=="Belgium") %>%
  filter(!IID %in% c("UH-Slu-Ew1", "UH-Slu-Ew2", "UH-Slu-Ew7")) %>%
  pull(IID)

NAm_MA_ids <-
  geoinfo_tb %>%
  filter(FID=="NAm", Locality=="MA") %>%
  pull(IID)

AU_other_ids <-
  geoinfo_tb %>%
  filter(FID=="AU", Locality!="WA") %>%
  pull(IID)

AU_WA_ids <-
  geoinfo_tb %>%
  filter(FID=="AU", Locality=="WA") %>%
  pull(IID)

SAm_Bariloche_ids <-
  geoinfo_tb %>%
  filter(FID=="SAm", Locality=="Bariloche") %>%
  pull(IID)

s208con_pop[["CEu_Belgium"]] <-
  s208con_pop[["CEu"]][s208con_pop[["CEu"]]$ind.names %in% CEu_Belgium_ids,]

s208con_pop[["NAm_MA"]] <-
  s208con_pop[["NAm"]][s208con_pop[["NAm"]]$ind.names %in% NAm_MA_ids,]

s208con_pop[["SAm_Bariloche"]] <-
  s208con_pop[["SAm"]][s208con_pop[["SAm"]]$ind.names %in% SAm_Bariloche_ids,]

s208con_pop[["AU_other"]] <-
  s208con_pop[["AU"]][s208con_pop[["AU"]]$ind.names %in% AU_other_ids,]

s208con_pop[["AU_WA"]] <-
  s208con_pop[["AU"]][s208con_pop[["AU"]]$ind.names %in% AU_WA_ids,]

applying_pop <-
  c("CEu", "CEu_Belgium", "NAm", "NAm_MA", "SAm", "AU", "AU_other", "AU_WA", "NZ", "SAm_Bariloche")

tajima_Ddf <- tibble()
for (pop in applying_pop){
  message(paste0("start ", pop))
  site_matrix <-
    t(as.matrix(s208con_pop[[pop]]))

  d <-
    site_matrix %>%
    cut_rows(wind_size) %>%
    map_df(tajimaD.win) %>%
    mutate(pop=pop)

  tajima_Ddf <- rbind(tajima_Ddf, d)
}

tajima_Ddf <- rbind(tajima_Ddf, d)

write_csv(tajima_Ddf, paste0("Fstats/tajimaD_", wind_size, ".csv"))
