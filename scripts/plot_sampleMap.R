library(tidyverse)
library(geosphere)
#library(DT)
library(maps)
#library(scatterpie)

source("scripts/geoinfo.R")

clades <-
  c("S. luteus, Clade Central Europe (introduced)",
    "S. luteus, Clade Central Europe (native)",
    "S. luteus, Clade Northern Europe",
    "S. luteus, Clade Asia",
    "S. brunnescens"
)

colors <-
  colour("muted")(9)

color_plate <-
  colors[c(5,2,7,3,9)]

names(color_plate) <-
  clades

keep_id <-
  read_table("plink/s214.fam", col_names = FALSE) %>%
  pull(X2)

raw_tb <-
  read_csv("metadata/parsed_geo.csv") %>%
  filter(IID %in% keep_id)

sum_country_tb <-
  raw_tb %>%
  group_by(FID) %>%
  summarise(n=n(), Country=paste(unique(Country), collapse = ", ")) %>%
  left_join(introduce_tb) %>%
  select(`Population/Clade`=FullName,
         Introduced=introduced,
         Country,
         n) %>%
  ungroup() %>%
  arrange(Introduced) %>%
  arrange(match(`Population/Clade`, c("S. brunnescens", "Clade Asia", "Clade Northern Europe")))

write_csv(sum_country_tb, "tables/sum_country_tb.csv")

mapWorld <- borders("world", colour="gray80", fill="gray80")
ggplot(raw_tb) +
  mapWorld +
  geom_point(aes(x=Longitude, y=Latitude), color="red")

# cluster
```{r}
k=13

GPS_points <-
  select(raw_tb, Longitude, Latitude)

mdist <-
  distm(GPS_points)

hc <- hclust(as.dist(mdist))

GPS_points_tb <-
  GPS_points %>%
  cbind(select(raw_tb, FID)) %>%
  mutate(group=cutree(hc, k=k)) %>%
  mutate(group=as.factor(group))

centroid_tb <-
  GPS_points_tb %>%
  group_by(group) %>%
  summarise(mean_long=mean(Longitude), mean_lat=mean(Latitude), .groups="keep")

ggplot(centroid_tb) +
  mapWorld +
  geom_point(aes(x=mean_long, y=mean_lat), color="red")

# all samples
cluster_sum <-
  GPS_points_tb %>%
  group_by(group, FID) %>%
  summarise(n(), .groups="keep") %>%
  ungroup() %>%
  rename(n=`n()`)

cluster_sum_wide <-
  cluster_sum %>%
  left_join(centroid_tb)

ggplot(cluster_sum_wide) +
  mapWorld +
  geom_point(aes(x=mean_long,
                 y=mean_lat,
                 size=n
  ), alpha=0.8, color="red") +
  geom_text(aes(mean_long, y=mean_lat, label=group))

set.seed(0)
sum_plot <-
  cluster_sum_wide %>%
  mutate(`Clade`=
           case_when(
             FID=="SAm" ~ "S. luteus, Clade Central Europe (introduced)",
             FID=="AF" ~ "S. luteus, Clade Central Europe (introduced)",
             FID=="NAm" ~ "S. luteus, Clade Central Europe (introduced)",
             FID=="AU" ~ "S. luteus, Clade Central Europe (introduced)",
             FID=="NZ" ~ "S. luteus, Clade Central Europe (introduced)",
             FID=="Sbo" ~ "S. brunnescens",
             FID=="Slu2" ~ "S. luteus, Clade Northern Europe",
             FID=="Slu1" ~ "S. luteus, Clade Asia",
             FID=="CEu" ~ "S. luteus, Clade Central Europe (native)"
           )
  ) %>%
  mutate(Clade=ordered(Clade, levels = clades)) %>%
  arrange(Clade) %>%
  ggplot() +
  mapWorld +
  geom_point(aes(x=mean_long,
                 y=mean_lat,
                 size=n,
                 color=`Clade`
  ),
  alpha=0.8,
  width=5,
  height=5) +
  scale_color_manual(values=color_plate) +
  scale_size_area(breaks=c(1, 10, 20, 40)) +
  coord_quickmap() +
  theme() +
  xlab("Longitude") +
  ylab("Latitude")

ggsave("plots/sample_map.pdf", sum_plot, width = 6.5*1.2, height = 5.72*1.2)
