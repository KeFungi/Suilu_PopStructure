require(hierfstat)
require(adegenet)
require(tidyverse)
require(cowplot)

gstat.randtest <- function(X,
                           method=c("global","within","between"),
                           nsim=499, subpop){
  met <- tolower(method[1])

  obs <- g.stats.glob(X)$g.stats

  pop <- X[,1]
  X <- X[,-1]

  sim <- vector(mode="numeric", length=nsim)

  if(met=="global"){

    sim <- sapply(1:nsim, function(i) g.stats.glob.mat(cbind(sample(pop),X))$g.stats)

  } else if(met=="within"){

    if(length(sup.pop) != length(pop)) stop("pop and sup.pop do not have the same length.")
    sim <- sapply(1:nsim, function(i) g.stats.glob.mat(cbind(pop,X[samp.within(sup.pop),]))$g.stats)

  } else if(met=="between"){

    if(length(sub.pop) != length(pop)) stop("pop and sub.pop do not have the same length.")
    sim <- sapply(1:nsim, function(i) g.stats.glob.mat(cbind(pop,X[samp.between(sub.pop),]))$g.stats)

  } else {
    stop("Unknown method requested.")
  }

  prevcall <- match.call()

  res <- as.randtest(sim=sim, obs=obs, call=prevcall)

  return(res)

}

genlight_dropInv <-
  function(genlight){
    genlight[,!(colSums(as.matrix(genlight)) %in% c(0, 2*dim(genlight)[1]))]
  }

g.stats.glob.mat <-
  function(heir, diploid=TRUE) {
    popl <- unique(heir[,1])
    names(popl) <- popl

    if (diploid) {
      X <- apply(heir[,-1], MARGIN = 2, genot2al)
      pop <-  rep(heir[,1], 2)
    }
    else {
      pop <- data[, 1]
      X <- data[,-1]
    }
    rows <- map(popl, ~pop==.x)
    cols <- list(A1=(X==1), A2=(X==2))
    A1_pops <- list()
    A2_pops <- list()
    pop_sums <- list()
    for (popi in popl){
      popi_A1 <-
        apply(cols$A1, MARGIN = 2, function(x){x&rows[[popi]]}) %>%
        apply(MARGIN = 2, sum, na.rm=TRUE)
      popi_A2 <-
        apply(cols$A2, MARGIN = 2, function(x){x&rows[[popi]]}) %>%
        apply(MARGIN = 2, sum, na.rm=TRUE)
      popi_A12 <- popi_A1+popi_A2
      A1_pops <- c(A1_pops, list(popi_A1))
      A2_pops <- c(A2_pops, list(popi_A2))
      pop_sums <- c(pop_sums, list(popi_A12))
    }

    nt <- colSums(!is.na(X), na.rm=TRUE)
    A1_sums <- colSums(cols$A1, na.rm=TRUE)
    A2_sums <- colSums(cols$A2, na.rm=TRUE)
    pop_A1exp <- map(pop_sums, ~.x*(A1_sums/nt))
    pop_A2exp <- map(pop_sums, ~.x*(A2_sums/nt))
    g_A1 <- map2(A1_pops, pop_A1exp, ~.x*log(.x/.y))
    g_A1 <- map(g_A1, ~ifelse(is.na(.x), 0, .x))
    g_A2 <- map2(A2_pops, pop_A2exp, ~.x*log(.x/.y))
    g_A2 <- map(g_A2, ~ifelse(is.na(.x), 0, .x))
    g.stats.l <- 2 * (reduce(g_A1,`+`) + reduce(g_A2, `+`))
    g.stats <- sum(g.stats.l)
    return(list(g.stats.l = g.stats.l, g.stats = g.stats))
  }

g.stats.glob <-
  function(data, diploid=TRUE) {
    nl <- dim(data)[2] - 1
    g.stats.l <- vector(length = nl)
    g.stats <- 0
    for (i in 1:nl) {
      if (diploid) {
        y <- genot2al(data[, (i + 1)])
        x <- rep(data[, 1], 2)
      }
      else {
        x <- data[, 1]
        y <- data[, i + 1]
      }
      obs <- table(x, y)
      nt <- sum(obs)
      s.r <- apply(obs, 1, sum)
      s.c <- apply(obs, 2, sum)
      expe <- s.r %*% t(s.c)/nt
      g.stats.l[i] <- 2 * sum(obs * log(obs/expe), na.rm = TRUE)
    }
    g.stats <- sum(g.stats.l)
    list(g.stats.l = g.stats.l, g.stats = g.stats)
  }

genlight2hierfst <-
  function(genlight, nsubsite=10000){
    genlight_mat <-
    as.matrix(genlight)
    nsite <- dim(genlight_mat)[2]
    if (!is.na(nsubsite)) {
      genlight_mat <- genlight_mat[,sample(1:nsite, size = nsubsite)]
    }
    genlight_mat[genlight_mat==0] <- 11
    genlight_mat[genlight_mat==1] <- 12
    genlight_mat[genlight_mat==2] <- 22

    X <-
      cbind(data.frame(pop=genlight$pop),
            as.data.frame(genlight_mat))
    return(X)
  }

gstat_wrapper <-
  function(genlight, nsim, nsite){
    hfst <-
      genlight2hierfst(genlight, nsite)
    gtest <-
      gstat.randtest(hfst, nsim = nsim)
    return(gtest)
  }

get_randtest_value <-
  function(randtest){
    return(list(obs=randtest$obs, pvalue=randtest$pvalue))
  }

get_randtest_plot <-
  function(randtest){
    return(as_grob(~plot(randtest)))
  }
