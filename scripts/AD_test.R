library(tidyverse)
library(cowplot)

source("scripts/geoinfo.R")

AF1_freq <-
  function(A1, A2) {
      return(A1/(A1+A2))
  }

AF_freq_plot <-
  function(file, color="black", title=NA) {
    read_delim(file, delim=" ", col_names = c("GT", "AD1","AD2", "DP")) %>%
      filter(GT!="./.", GT!="0/0", GT!="1/1", GT!="0|0", GT!="1|1") %>%
      mutate(AF1=AF1_freq(AD1, AD2)) %>%
      mutate(AF=map2(AF1, 1-AF1, ~c(.x, .y))) %>%
      mutate(seed=sample(c(1, 2), size=n(), replace=TRUE)) %>%
      mutate(AF=map2_dbl(AF, seed, ~.x[.y])) %>%
      ggplot(.) +
      scale_x_continuous(breaks = c(0, 0.5, 1.0), limits = c(0, 1)) +
      stat_density(aes(x=AF), fill=color) +
      xlab(ifelse(is.na(title), file, title)) +
      ylab("") +
      theme(axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(size=6),
            title = element_text(size=10)
      )
  }

IIDs <- read_tsv("metadata/AD_list.txt", col_names = c("IID", "group"))

hap_ID <-
  read_table("metadata/hap.fam.txt", col_names = c("FID", "IID")) %>%
  mutate(group="hap")

mixed_ID <-
  read_table("metadata/mixed.fam.txt", col_names = c("FID", "IID")) %>%
  mutate(group="mixed")

other_ID <-
  read_table("plink/s275.fam", col_names = c("FID", "IID")) %>%
  filter(!(IID %in% hap_ID$IID | IID %in% mixed_ID$IID)) %>%
  mutate(group="other")

ind_groups <-
  bind_rows(hap_ID, mixed_ID, other_ID)

MIX_plot <-
  map(ind_groups$IID[ind_groups$group=="mixed"], ~AF_freq_plot(paste0("AF/", .x, ".AD.txt"), color = "red", title=.x))

normal_plots <-
  map(ind_groups$IID[ind_groups$group=="other"], ~AF_freq_plot(paste0("AF/", .x, ".AD.txt"), color = "blue", title=.x))

hap_plots <-
  map(ind_groups$IID[ind_groups$group=="hap"], ~AF_freq_plot(paste0("AF/", .x, ".AD.txt"), color = "black", title=.x))

plotlist <-
  c(normal_plots, hap_plots, MIX_plot) %>%
  map(~.x+theme(title = element_text(size=7)))

dummy_plot <-
  tibble(y=1:3, Category=factor(c("Normal", "Haploid", "Mixed"))) %>%
  ggplot(aes(x=y, fill=Category)) +
  geom_bar() +
  scale_fill_manual(
    values = c(
      "Normal"="blue",
      "Haploid"="black",
      "Mixed"="red"
    )
  ) +
  theme(legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.1, "inches")
  )

dummy_legend <-
  get_legend(dummy_plot)

n <-
  dim(ind_groups)[1]

n_subplot <-
  (n %/% (6*8)) + 1

ggblank <-
  ggplot() +
  theme_nothing()

plotlist[(n+1):(n_subplot*6*8)] <-
  list(ggblank)

plotlist <-
  append(plotlist, list(dummy_legend), 0)

range_start <-
  seq(1, n, 6*8)

range_end <-
  seq(6*8, n_subplot*6*8, 6*8)

AD_plots <-
  map(
    1:n_subplot,
    ~plot_grid(
      plotlist = plotlist[range_start[.x]:range_end[.x]],
      ncol = 6
      )
  )

walk(1:n_subplot,
  ~ggsave(paste0("plots/AD_density", .x, ".pdf"), AD_plots[[.x]], width=6, height = 8)
)
