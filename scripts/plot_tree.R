library(treeio)
library(ggtree)
library(aplot)
library(tanggle)
library(phangorn)

source("scripts/geoinfo.R")

clades <-
  c("S. brunnescens",
    "Clade Northern Europe",
    "Clade Asia",
    "Clade Central Europe",
    "Central Europe",
    "North America",
    "South America",
    "Australia",
    "Western Australia",
    "New Zealand",
    "Africa")

colors <-
  colour("muted")(9)

color_plate <-
  c(colors[c(9, 7, 3)], "#000000", colors[c(4, 5, 6, 2, 2, 8,1)])

names(color_plate) <-
  clades

tree_id <-
  read_csv("metadata/S. luteus resequencing metadata_main.csv")

phy_tree <-
  read.newick("phylogeny/RAxML_result.s214-CAT") %>%
  phytools::reroot(node.number = 105, position=0.05)

phylo_data <-
  phy_tree %>%
  as_tibble() %>%
  left_join(select(tree_id, GenomeCode, City, Country), by=c("label"="GenomeCode")) %>%
  left_join(select(geoinfo_tb, IID, FID),  by=c("label"="IID")) %>%
  left_join(introduce_tb) %>%
  mutate(Clade=factor(FullName,
                      levels =clades
  )
  ) %>%
  mutate(long_label=
           ifelse(is.na(label), NA,
             paste0(label,
                    "-",
                    City,
                    "(",
                    Country,
                    ")"
             )
           )
  ) %>%
  mutate(long_label=
           ifelse(label=="Sb2", "S. brevipes", long_label),
         branch.length=
           ifelse(label!="Sb2"|is.na(label), branch.length, 0.2)
  ) %>%
  as.treedata()

id_clade_Eu <-
  phylo_data %>%
  as_tibble() %>%
  filter(FullName %in% clades[5:10]) %>%
  pull(label)

id_clade_WA <-
  pull(filter(geoinfo_tb, Locality=="WA", FID=="AU"), IID)

tree_plot <-
  ggtree(phylo_data) +
  geom_tiplab(aes(label=long_label, color=Clade), size=1) +
  geom_treescale(width=0.05) +
  scale_color_manual(values=color_plate) +
  theme(legend.position = "none") #+
  #xlim(0, 0.75)

location_dots <-
  locality_tb %>%
  left_join(introduce_tb) %>%
  select(IID, FullName) %>%
  rbind(tibble(IID="Sb2", FullName=NA))

location_dots_expand <-
  rbind(location_dots,
        tibble(IID=location_dots$IID[location_dots$FullName %in% c("Central Europe",
                                  "Africa",
                                  "South America",
                                  "Australia",
                                  "New Zealand",
                                  "North America")],
         FullName="Clade Central Europe")
  ) %>%
  rbind(tibble(IID=id_clade_WA,
               FullName="Western Australia")
  ) %>%
  mutate(Clade=factor(FullName,
                      levels = clades

  )
  ) %>%
  arrange(Clade)

dot_plot <-
  location_dots_expand %>%
  ggplot(aes(x=Clade, y=IID, color=Clade)) +
  geom_point(size=0.5) +
  scale_color_manual(values=color_plate) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  scale_x_discrete(position = "top",
                   limits=clades) +
  theme(axis.text.y =element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, size=14, hjust = 0),
        axis.title.x =element_blank(),
        legend.key.size = unit(c(16, 4), "points"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        panel.grid.major = element_line(color="grey80", size = 0.2)
        )

insert_left(dot_plot+theme(legend.position = "none"), tree_plot, 4.5/1.5) %>%
  ggsave("plots/dot_tree.pdf", ., width = 12, height = 12)

dummay1 <-
  ggplot(tibble(A=1:11, Clade=clades)) +
  geom_point(aes(x=A, y=A, color=Clade)) +
  scale_color_manual(values=color_plate[1:4])

dummay2 <-
  ggplot(tibble(A=1:11, Population=clades)) +
  geom_point(aes(x=A, y=A, color=Population)) +
  scale_color_manual(values=color_plate[5:10])

ggsave("plots/dot_tree_legend1.pdf", dummay1)
ggsave("plots/dot_tree_legend2.pdf", dummay2)

Nnet <- phangorn::read.nexus.networx("/Users/keyh/Suilu_thesis/splitstree/s189_admixture.splitstree.nex")

nt_plot <-
  ggsplitnet(Nnet, size=0.1) +
  geom_tiplab2(aes(label=label), size=2)

ggsave("plots/nt_plot.pdf", nt_plot, width = 8, height = 8)
Nnet %>%
  as_tibble() %>%
  as.phylo() %>%
as.treedata() %>%
  as.phylo %>%
  ggsplitnet()
ggplot(Nnet, aes(x, y))  + geom_splitnet() #+ theme_tree()+
  #scale_color_manual(values=rainbow(15)) +
  geom_tiplab2(aes(label=label, color=label))+
  theme(legend.position="none")

ggplot(Nnet, aes(x, y))  + geom_splitnet() + theme_tree() +

  #scale_color_manual(values=rainbow(2)) +

  geom_tiplab2(aes(label=label, color=label))+

  theme(legend.position="none")



Nnet %>%
  as_tibble() %>%
  as.phylo() %>%
  ggtree()
  as.networx() %>%
ggsplitnet() #+
  #geom_tiplab2()
#ggtree(layout="daylight", aes(color=Clade)) +
  geom_tiplab(aes(label=long_label, color=label), size=1)
