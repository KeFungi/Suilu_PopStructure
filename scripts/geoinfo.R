require(tidyverse)
require(khroma)

fam_tb <-
  read_delim("metadata/pop.fam.txt", col_names = FALSE, delim=" ") %>%
  select(IID=X4, FID=X3)

locality_tb <-
  read_csv("metadata/parsed_geo.csv")

introduce_tb <-
  tribble(~ "FID", ~"introduced", ~"FullName",
          "Sbo", "Native", "S. brunnescens",
          "Slu1", "Native", "Clade Asia",
          "Slu2", "Native", "Clade Northern Europe",
          "CEu", "Native", "Central Europe",
          "AU", "Introduced", "Australia",
          "NZ", "Introduced", "New Zealand",
          "SAm", "Introduced", "South America",
          "NAm", "Introduced", "North America",
          "AF", "Introduced", "Africa"
          ) %>%
  mutate(introduced=factor(introduced, levels = c("Native", "Introduced")))%>%
  arrange(introduced)

geoinfo_tb <-
  tibble() %>%
  bind_rows(
    select(filter(locality_tb, FID=="SAm"),
           IID, FID, Longitude, Latitude,
           Locality=Locality)
    ) %>%
  bind_rows(
    select(filter(locality_tb, FID=="NAm"),
           IID, FID, Longitude, Latitude,
           Locality=StateProvince)
  ) %>%
  bind_rows(
    select(filter(locality_tb, FID=="NZ"),
           IID, FID, Longitude, Latitude,
           Locality=Locality)
  ) %>%
  bind_rows(
    select(filter(locality_tb, FID=="AU"),
           IID, FID, Longitude, Latitude,
           Locality=StateProvince)
  ) %>%
  bind_rows(
    select(filter(locality_tb, FID=="CEu"),
           IID, FID, Longitude, Latitude,
           Locality=Country)
  ) %>%
  bind_rows(
    select(filter(locality_tb, FID=="Slu1"),
           IID, FID, Longitude, Latitude,
           Locality=Country)
  ) %>%
  bind_rows(
    select(filter(locality_tb, FID=="Slu2"),
           IID, FID, Longitude, Latitude,
           Locality=Country)
  ) %>%
  bind_rows(
    select(filter(locality_tb, FID=="Sbo"),
           IID, FID, Longitude, Latitude,
           Locality=Country)
  ) %>%
  bind_rows(
    select(filter(locality_tb, FID=="AF"),
           IID, FID, Longitude, Latitude,
           Locality=Country)
  )

muted <- colour("muted")

pop_names <-
  introduce_tb$FID

pop_names_full <-
  introduce_tb$FullName

pop_color <-
  muted(9) %>%
  set_names(pop_names)

pop_color_full <-
  muted(9) %>%
  set_names(pop_names_full)
