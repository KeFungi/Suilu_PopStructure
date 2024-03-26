library(tidyverse)
library(cowplot)
library(Cairo)
library(multcompView)
library(exactRankTests)
library(rcompanion)
library(magrittr)
source("scripts/geoinfo.R")

wilcox_wrapper <-
  function(.data, pre_col, res_col){
    pre_col <- rlang::sym(pre_col)
    res_col <- rlang::sym(res_col)

    .data <- filter(.data, !is.na(!!pre_col), !is.na(!!res_col))
    x <- pull(.data, !!pre_col) %>% as.factor()
    y <- pull(.data, !!res_col)

    w_results <- wilcox.exact(y ~ x)

    tibble(var=as.character(res_col), W=w_results$statistic, p=w_results$p.value)
  }



box_group <-
  function(.intb, group, refgroup, var, color_scale=NULL, p.adjust.method="holm"){
    .tb <-
      select(.intb , all_of(c(group, var))) %>%
      set_colnames(c("group", "var")) %>%
      mutate(group=fct_relevel(group, refgroup))

    range_data <- max(pull(.tb, var)) - min(pull(.tb, var))

    groups <-
      .tb$group %>%
      levels()

    top_values <-
      groups %>%
      map(~filter(.tb, group==.x)$var) %>%
      map(~tibble(med=median(.x), q100=quantile(.x, 1), q75=quantile(.x, 0.75), q25=quantile(.x, 0.25))) %>%
      reduce(bind_rows) %>%
      mutate(IQR=q75-q25) %>%
      mutate(top=q75+1.5*IQR)

    pairs <-
      expand_grid(g1=refgroup, g2=groups[!groups %in% refgroup]) %>%
      rowwise() %>%
      mutate(vec=list(c(g1, g2)))

    p <- c()
    for (pair in pairs$vec){

      fit <-
        filter(.tb, group %in% pair) %>%
        mutate(group=fct_drop(group)) %>%
        wilcox_wrapper("group", "var")

      p <-
        pull(fit, p) %>%
        c(p, .)
    }

    p.adj <-
      bind_cols(select(pairs, -vec), tibble(p=p.adjust(p, method = p.adjust.method)))

    p.matrix <-
      mutate(p.adj, ng1=g2, ng2=g1) %>%
      select(g1=ng1, g2=ng2, p) %>%
      bind_rows(p.adj) %>%
      full_join(filter(expand_grid(g1=groups, g2=groups))) %>%
      rowwise() %>%
      mutate(p=ifelse(is.na(p), 1, p)) %>%
      pivot_wider(names_from = "g1", values_from = "p") %>%
      arrange(factor(g2, groups)) %>%
      select(g2, any_of(groups)) %>%
      column_to_rownames("g2") %>%
      as.matrix()

    cld_letters <-
      multcompLetters(p.matrix)$Letters

    cld_dt <-
      tibble(group=names(cld_letters), cld=cld_letters)
    return(cld_dt)
  }

wind_size <- 10000
tajima_Ddf <- read_csv(paste0("Fstats/tajimaD_", wind_size, ".csv"))

tajima_Ddf_input <-
  tajima_Ddf %>%
  filter(!pop %in% c("Slu1", "Slu2")) %>%
  mutate(FID=case_when(
    pop %in% c("NAm_MA") ~ "NAm",
    pop %in% c("CEu_Belgium") ~ "CEu",
    pop %in% c("AU_WA", "AU_other") ~ "AU",
    pop %in% c("SAm_Bariloche") ~ "SAm",
    TRUE ~ pop)
  ) %>%
  left_join(introduce_tb) %>%
  mutate(subpop=case_when(
    pop == "NAm_MA" ~ "North America - MA",
    pop == "AU_WA" ~ "Australia - WA",
    pop == "AU_other" ~ "Australia - other",
    pop == "CEu_Belgium" ~ "CEu - Belgium",
    pop == "SAm_Bariloche" ~ "SAm - Bariloche",
    TRUE ~ FullName
  )
  ) %>%
  rename(Population=FullName) %>%
  mutate(subpop=factor(subpop)) %>%
  mutate(subpop=fct_relevel(subpop, c("Australia", "Australia - WA", "Australia - other")))

div_theme <-
  list(scale_color_manual(values=pop_color_full[unique(tajima_Ddf_input$Population)], drop=TRUE, name="Population"),
       theme_classic(),
       theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       xlab("")
  )

m_pi <-
  tajima_Ddf_input %>%
  filter(subpop == "Central Europe") %>%
  pull(pi) %>%
  median()

m_the <-
  tajima_Ddf_input %>%
  filter(subpop == "Central Europe") %>%
  pull(Wt) %>%
  median()

pi_letters_tb <-
  box_group(.intb = tajima_Ddf_input, group = "subpop",
            refgroup = c("Central Europe", "CEu - Belgium"),
            var = "pi") %>%
  rename(subpop=group) %>%
  left_join(unique(select(tajima_Ddf_input, introduced, Population, subpop)))

D_letters_tb <-
  box_group(.intb = tajima_Ddf_input, group = "subpop",
            refgroup = c("Central Europe", "CEu - Belgium"),
            var = "Wt") %>%
  rename(subpop=group) %>%
  left_join(unique(select(tajima_Ddf_input, introduced, Population, subpop)))

theta_letters_tb <-
  box_group(.intb = tajima_Ddf_input, group = "subpop",
            refgroup = c("Central Europe", "CEu - Belgium"),
            var = "D") %>%
  rename(subpop=group) %>%
  left_join(unique(select(tajima_Ddf_input, introduced, Population, subpop)))


pi_plot <-
  ggplot(tajima_Ddf_input) +
  geom_boxplot(aes(x=subpop, y=pi, color=Population), outlier.shape = NA, notch=TRUE) +
  geom_text(data=pi_letters_tb, aes(x=subpop, y=150, label=cld, color=Population)) +
  geom_hline(yintercept = m_pi, color="red", lty=2) +
  theme_classic() +
  div_theme +
  facet_grid(cols=vars(introduced), scales="free_x", space="free") +
  ylab("Nucleotide Diversity (π)") +
  xlab("")

theta_plot <-
  ggplot(tajima_Ddf_input) +
  geom_boxplot(aes(x=subpop, y=Wt, color=Population), outlier.shape = NA, notch=TRUE) +
  geom_text(data=theta_letters_tb, aes(x=subpop, y=150, label=cld, color=Population)) +
  geom_hline(yintercept = m_the, color="red", lty=2) +
  theme_classic() +
  div_theme +
  ylab("Watterson's θ") +
  facet_grid(cols=vars(introduced), scales="free_x", space="free")

D_plot <-
  ggplot(tajima_Ddf_input) +
  geom_boxplot(aes(x=subpop, y=D, color=Population), outlier.shape = NA, notch=TRUE) +
  geom_text(data=D_letters_tb, aes(x=subpop, y=-1.2, label=cld, color=Population)) +
  geom_hline(yintercept = 0, lty=2, color="red") +
  theme_classic() +
  div_theme +
  ylab("Tajima's D") +
  facet_grid(cols=vars(introduced), scales="free_x", space="free")

legend_plot <-
  get_legend(theta_plot + theme(legend.key.size = unit(c(0.1, 0.2), "inch")))

final_plot <-
  plot_grid(
    pi_plot,
    theta_plot,
    D_plot,
    ncol=1,
    greedy=FALSE)

ggsave("plots/genDiversity.pdf", final_plot, width=6, height=10, device = cairo_pdf)

ggplot(tajima_Ddf_input) +
  geom_boxplot(aes(x=pop, y=D, color=pop), outlier.shape = NA) +
  geom_hline(yintercept = 0, color="black", lty=2) +
  theme_classic() +
  #ylim(-2, 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
