library(plotly)
library(cowplot)
library(magrittr)
library(adegenet)
library(hierfstat)
library(geosphere)
library(poppr)
library(dartR)

source("scripts/geoinfo.R")
source("scripts/gstat.randtest.R")

genetic.distance <-
  function(genlight) {
    provesti.dist(as.matrix(genlight))
  }

s208ad_bed <-
  read.PLINK("plink/s208_admixture.raw")

s189ad_bed <-
  s208ad_bed[!s208ad_bed$pop %in% c("Slu1", "Slu2"),]

s189ad_pop <-
  seppop(s189ad_bed)

applying_pop <-
  which(table(s189ad_bed$pop) > 5) %>%
  names() %>%
  set_names(., .)

ibd_plots <- list()

for (pop in applying_pop){
  set.seed(12345)
  subgen <-
    s189ad_pop[[pop]] %>%
    .[,sample(1:dim(.)[2], size=10000)]

  subgen$pop <-
    geoinfo_tb$Locality[match(subgen$ind.names, geoinfo_tb$IID)] %>%
    as.factor()

  gen.dist <- dist(subgen)

  geo.dist <-
    geoinfo_tb[match(subgen$ind.names, geoinfo_tb$IID),] %>%
    select(IID, Longitude, Latitude) %>%
    column_to_rownames("IID") %>%
    as.matrix() %>%
    distm() %>%
    as.dist()

  ibd <- mantel.randtest(gen.dist,geo.dist)
  #plot(geo.dist, gen.dist)
  dist_lm <- lm(as.vector(gen.dist) ~ as.vector(geo.dist))
  #abline(dist_lm, col="red", lty=2)
  #title(introduce_tb$FullName[match(pop, introduce_tb$FID)],
  #      paste0("mantel p=",ibd$pvalue)
  #      )

  m.A <-
    dist(subgen) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("IID") %>%
    pivot_longer(-IID, values_to = "gen")

  m.B <-
    geoinfo_tb[match(subgen$ind.names, geoinfo_tb$IID),] %>%
    select(IID, Longitude, Latitude) %>%
    column_to_rownames("IID") %>%
    as.matrix() %>%
    distm() %>%
    set_colnames(subgen$ind.names) %>%
    set_rownames(subgen$ind.names) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("IID") %>%
    pivot_longer(-IID, values_to = "geo")

  p <-
    left_join(m.A, m.B) %>%
    left_join(geoinfo_tb) %>%
    left_join(geoinfo_tb, by=c("name"="IID")) %>%
    mutate(IID=ordered(IID, levels = unique(IID)),
           name=ordered(name, levels = unique(IID))
    ) %>%
    filter(IID>name) %>%
    ggplot() +
    geom_jitter(aes(x=geo, y=gen)) +
    geom_abline(slope = dist_lm$coefficients[2],
                intercept = dist_lm$coefficients[1],
                color="red"
    ) + xlab("geographic distance") +
    ylab("genetic distance") +
    labs(title = introduce_tb$FullName[match(pop,introduce_tb$FID )],
         caption = paste0("matel p=", ibd$pvalue)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5))
  ibd_plots[[pop]] <- p
}

# AU
pop <- "AU"
set.seed(12345)
subgen <-
  s189ad_pop[[pop]] %>%
  .[,sample(1:dim(.)[2], size=10000)]

subgen$pop <-
  geoinfo_tb$Locality[match(subgen$ind.names, geoinfo_tb$IID)] %>%
  as.factor()

gen.dist <- dist(subgen)

geo.dist <-
  geoinfo_tb[match(subgen$ind.names, geoinfo_tb$IID),] %>%
  select(IID, Longitude, Latitude) %>%
  column_to_rownames("IID") %>%
  as.matrix() %>%
  distm() %>%
  as.dist()

ibd <- mantel.randtest(gen.dist, geo.dist)
#plot(geo.dist, gen.dist)
dist_lm <- lm(as.vector(gen.dist) ~ as.vector(geo.dist))
#abline(dist_lm, col="red", lty=2)
#title(introduce_tb$FullName[match(pop, introduce_tb$FID)],
#      paste0("mantel p=",ibd$pvalue)
#      )

m.A <-
  gen.dist %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("IID") %>%
  pivot_longer(-IID, values_to = "gen")

m.B <-
  geo.dist %>%
  as.matrix() %>%
  set_colnames(subgen$ind.names) %>%
  set_rownames(subgen$ind.names) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("IID") %>%
  pivot_longer(-IID, values_to = "geo")


p <-
  left_join(m.A, m.B) %>%
  left_join(geoinfo_tb) %>%
  left_join(geoinfo_tb, by=c("name"="IID")) %>%
  mutate(IID=ordered(IID, levels = unique(IID)),
         name=ordered(name, levels = unique(IID))
  ) %>%
  filter(IID>name) %>%
  mutate(pop=case_when(Locality.x=="WA"&Locality.y=="WA" ~ "within WA",
                       Locality.x=="WA"|Locality.y=="WA" ~ "to WA",
                       TRUE ~ "others")
  ) %>%
  mutate(pop=factor(pop, levels = c("within WA", "to WA", "others"))) %>%
  arrange(pop) %>%
  ggplot() +
  geom_jitter(aes(x=geo, y=gen, shape=pop, color=pop)) +
  scale_shape_manual(name="pairs",
                     values=c("within WA"=1,
                              "to WA" = 2,
                              "others"=4)) +
  scale_color_manual(name="pairs",
                     values=c("within WA"="blue",
                              "to WA"="blue",
                              "others"="black")) +
  geom_abline(slope = dist_lm$coefficients[2],
              intercept = dist_lm$coefficients[1],
              color="red"
  ) + xlab("geographic distance") +
  ylab("genetic distance") +
  labs(title = introduce_tb$FullName[match(pop,introduce_tb$FID )],
       caption = paste0("matel p=", ibd$pvalue)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

ibd_plots[[pop]] <- p

#AU(no WA)
pop <- "AU_noWA"
set.seed(12345)
subgen <-
  s189ad_pop[["AU"]] %>%
  .[,sample(1:dim(.)[2], size=10000)]

subgen$pop <-
  geoinfo_tb$Locality[match(subgen$ind.names, geoinfo_tb$IID)] %>%
  as.factor()

subgen <-
  subgen[subgen$pop!="WA",]

gen.dist <- dist(subgen)

geo.dist <-
  geoinfo_tb[match(subgen$ind.names, geoinfo_tb$IID),] %>%
  select(IID, Longitude, Latitude) %>%
  column_to_rownames("IID") %>%
  as.matrix() %>%
  distm() %>%
  as.dist()

ibd <- mantel.randtest(gen.dist, geo.dist)
#plot(geo.dist, gen.dist)
dist_lm <- lm(as.vector(gen.dist) ~ as.vector(geo.dist))
#abline(dist_lm, col="red", lty=2)
#title(introduce_tb$FullName[match(pop, introduce_tb$FID)],
#      paste0("mantel p=",ibd$pvalue)
#      )

m.A <-
  gen.dist %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("IID") %>%
  pivot_longer(-IID, values_to = "gen")

m.B <-
  geo.dist %>%
  as.matrix() %>%
  set_colnames(subgen$ind.names) %>%
  set_rownames(subgen$ind.names) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("IID") %>%
  pivot_longer(-IID, values_to = "geo")


p <-
  left_join(m.A, m.B) %>%
  left_join(geoinfo_tb) %>%
  left_join(geoinfo_tb, by=c("name"="IID")) %>%
  mutate(IID=ordered(IID, levels = unique(IID)),
         name=ordered(name, levels = unique(IID))
  ) %>%
  filter(IID>name) %>%
  mutate(pop=case_when(Locality.x=="WA"&Locality.y=="WA" ~ "within WA",
                       Locality.x=="WA"|Locality.y=="WA" ~ "to WA",
                       TRUE ~ "others")
  ) %>%
  mutate(pop=factor(pop, levels = c("within WA", "to WA", "others"))) %>%
  arrange(pop) %>%
  ggplot() +
  geom_jitter(aes(x=geo, y=gen, shape=pop, color=pop)) +
  scale_shape_manual(name="pairs",
                     values=c("within WA"=1,
                              "to WA" = 2,
                              "others"=4)) +
  scale_color_manual(name="pairs",
                     values=c("within WA"="blue",
                              "to WA"="blue",
                              "others"="black")) +
  geom_abline(slope = dist_lm$coefficients[2],
              intercept = dist_lm$coefficients[1],
              color="red"
  ) + xlab("geographic distance") +
  ylab("genetic distance") +
  labs(title = "Australia (no WA)",
       caption = paste0("matel p=", ibd$pvalue)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

ibd_plots[[pop]] <- p

# SAm
pop <- "SAm"
set.seed(12345)
subgen <-
  s189ad_pop[[pop]] %>%
  .[,sample(1:dim(.)[2], size=10000)]

subgen$pop <-
  geoinfo_tb$Locality[match(subgen$ind.names, geoinfo_tb$IID)] %>%
  as.factor()

gen.dist <- dist(subgen)

geo.dist <-
  geoinfo_tb[match(subgen$ind.names, geoinfo_tb$IID),] %>%
  select(IID, Longitude, Latitude) %>%
  column_to_rownames("IID") %>%
  as.matrix() %>%
  distm() %>%
  as.dist()

ibd <- mantel.randtest(gen.dist, geo.dist)
#plot(geo.dist, gen.dist)
dist_lm <- lm(as.vector(gen.dist) ~ as.vector(geo.dist))
#abline(dist_lm, col="red", lty=2)
#title(introduce_tb$FullName[match(pop, introduce_tb$FID)],
#      paste0("mantel p=",ibd$pvalue)
#      )

m.A <-
  gen.dist %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("IID") %>%
  pivot_longer(-IID, values_to = "gen")

m.B <-
  geo.dist %>%
  as.matrix() %>%
  set_colnames(subgen$ind.names) %>%
  set_rownames(subgen$ind.names) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("IID") %>%
  pivot_longer(-IID, values_to = "geo")


p <-
  left_join(m.A, m.B) %>%
  left_join(locality_tb) %>%
  left_join(locality_tb, by=c("name"="IID")) %>%
  mutate(IID=ordered(IID, levels = unique(IID)),
         name=ordered(name, levels = unique(IID))
  ) %>%
  filter(IID>name) %>%
  mutate(pop=case_when(Country.x=="Argentina"&Country.y=="Argentina" ~ "within Argentina",
                       Country.x=="Chile"&Country.y=="Chile" ~ "within Chile",
                       TRUE ~ "between")
  ) %>%
  mutate(pop=factor(pop, levels = c("within Argentina", "within Chile", "between"))) %>%
  arrange(pop) %>%
  ggplot() +
  geom_jitter(aes(x=geo, y=gen, shape=pop, color=pop)) +
  scale_shape_manual(name="pairs",
                     values=c("within Argentina"=1,
                              "within Chile" = 2,
                              "between"=4)) +
  scale_color_manual(name="pairs",
                     values=c("within Argentina"="blue",
                              "within Chile"="blue",
                              "between"="black")) +
  geom_abline(slope = dist_lm$coefficients[2],
              intercept = dist_lm$coefficients[1],
              color="red"
  ) + xlab("geographic distance") +
  ylab("genetic distance") +
  labs(title = introduce_tb$FullName[match(pop,introduce_tb$FID )],
       caption = paste0("matel p=", ibd$pvalue)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

ibd_plots[[pop]] <- p

saveRDS(ibd_plots, "ibd/ibd_plots.rds")
