library(dplyr)
library(tidyverse)
library(flextable)

#### SETUP ####
Tobfuscate <- readRDS(here::here("Objects/Talent obfuscate.rds"))

Aobfuscate <- readRDS(here::here("Objects/Ashland obfuscate.rds"))

t1 <- readRDS(here::here('Objects/Talent sum table part 1 (spatial centrality).RDS'))%>%
  dplyr::select(
    -c(Total, week, dg_cn, btwn, closeness, pgrnk, rev_pgrnk)
  )

t1 <- t1 %>%
  dplyr::rename(
    btwn = "btwn_in",
    closeness = "closeness_in",
    pgrnk = "pgrnk_in",
    dg_cn = "deg_in"
  )

t2 <- readRDS(here::here('Objects/Talent sum table part 2 (genomic centrality).RDS'))

t3 <- readRDS(here::here('Objects/Talent sum table part 3 (dispersal).RDS'))

tal_sum <- list(t1, t2, t3) %>%
  reduce(full_join, by = "Site") %>%
  arrange(Site) %>%
  left_join(Tobfuscate, by = c("Site" = "Tsite"))

tal_long <- tal_sum %>%
  dplyr::select(-c(Site)) %>%
  dplyr::mutate(Tnum = as.numeric(stringr::str_extract(Tlabel, "\\d+"))) %>%
  dplyr::arrange(Tnum) %>%
  tidyr::pivot_longer(
    cols = c(
      dg_cn,
      btwn,
      closeness,
      pgrnk,
      dg_cn0,
      btwn0,
      closeness0,
      pgrnk0,
      dg_cnA,
      btwnA,
      closenessA,
      pgrnkA
    ),
    names_to = c("metric", "network"),
    names_pattern = "^(dg_cn|btwn|closeness|pgrnk)(0|A)?$",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    network = dplyr::case_when(
      network == "0" ~ "Genetic Null",
      network == "A" ~ "Genetic Alternative",
      TRUE ~ "Spatial network"
    ),
    metric = dplyr::recode(
      metric,
      dg_cn = "Degree",
      btwn = "Betweenness",
      closeness = "Closeness",
      pgrnk = "PageRank"
    )
  ) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value) %>%
  dplyr::select(-Tnum)

tal_long %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  relocate(Tlabel) %>%
  flextable() %>%
  set_header_labels(Tlabel = "Site", network = "Network") %>%
  merge_v(j = "Tlabel") %>%
  flextable::compose(part = "header",
                     j = "FS",
                     value = as_paragraph(">\u03C3", as_sub("FS"))) %>%
  flextable::compose(part = "header",
                     j = "GG",
                     value = as_paragraph(">\u03C3", as_sub("GG"))) %>%
  flextable::compose(part = "header",
                     j = "PO",
                     value = as_paragraph(">\u03C3", as_sub("PO"))) %>%
  theme_vanilla() %>%
  flextable::autofit()

a1 <- readRDS(here::here('Objects/Ashland sum table part 1 (spatial centrality).RDS'))
a2 <- readRDS(here::here('Objects/Ashland sum table part 2 (genomic centrality).RDS'))
a3 <- readRDS(here::here('Objects/Ashland sum table part 3 (dispersal).RDS'))

ash_sum <- list(a1, a2, a3) %>%
  reduce(full_join, by = "Site") %>%
  arrange(Site) %>%
  left_join(Aobfuscate, by = c("Site" = "Asite")) %>%
  dplyr::mutate(Alabel = ifelse(is.na(Alabel), "A14", Alabel),
                label = ifelse(is.na(label), "A14", label))

ash_long <- ash_sum %>%
  dplyr::select(-c(Site, label)) %>%
  dplyr::mutate(Tnum = as.numeric(stringr::str_extract(Alabel, "\\d+"))) %>%
  dplyr::arrange(Tnum) %>%
  tidyr::pivot_longer(
    cols = c(
      dg_cn,
      btwn,
      closeness,
      pgrnk,
      dg_cn0,
      btwn0,
      closeness0,
      pgrnk0,
      dg_cnA,
      btwnA,
      closenessA,
      pgrnkA
    ),
    names_to = c("metric", "network"),
    names_pattern = "^(dg_cn|btwn|closeness|pgrnk)(0|A)?$",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    network = dplyr::case_when(
      network == "0" ~ "Genetic Null",
      network == "A" ~ "Genetic Alternative",
      TRUE ~ "Spatial network"
    ),
    metric = dplyr::recode(
      metric,
      dg_cn = "Degree",
      btwn = "Betweenness",
      closeness = "Closeness",
      pgrnk = "PageRank"
    )
  ) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value) %>%
  dplyr::select(-Tnum)

ash_long %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  relocate(Alabel) %>%
  flextable() %>%
  set_header_labels(Alabel = "Site", network = "Network") %>%
  merge_v(j = "Alabel") %>%
  flextable::compose(part = "header",
                     j = "FS",
                     value = as_paragraph(">\u03C3", as_sub("FS"))) %>%
  flextable::compose(part = "header",
                     j = "GG",
                     value = as_paragraph(">\u03C3", as_sub("GG"))) %>%
  flextable::compose(part = "header",
                     j = "PO",
                     value = as_paragraph(">\u03C3", as_sub("PO"))) %>%
  theme_vanilla() %>%
  autofit()


# correlation graph #####
ggplot() +
  ggplot2::geom_point(data = tal_long, aes(x = Tlabel, y = Degree, color = network))

ggplot() +
  ggplot2::geom_point(data = ash_long, aes(x = Alabel, y = Degree, color = network))

### correlation coefficient (spearman) #####
ash_cor <- ash_sum %>%
  dplyr::mutate(
    dg_cn = dg_cn / max(dg_cn, na.rm = TRUE),
    dg_cnA = dg_cnA / max(dg_cnA, na.rm = TRUE),
    dg_cn0 = dg_cn0 / max(dg_cn0, na.rm = TRUE),
    btwn = btwn / max(btwn, na.rm = TRUE),
    btwnA = btwnA / max(btwnA, na.rm = TRUE),
    btwn0 = btwn0 / max(btwn0, na.rm = TRUE),
    closeness = closeness / max(closeness, na.rm = TRUE),
    closenessA = closenessA / max(closenessA, na.rm = TRUE),
    closeness0 = closeness0 / max(closeness0, na.rm = TRUE),
    pgrnk = pgrnk / max(pgrnk, na.rm = TRUE),
    pgrnkA = pgrnkA / max(pgrnkA, na.rm = TRUE),
    pgrnk0 = pgrnk0 / max(pgrnk0, na.rm = TRUE),
    city = "Ashland",
    Alabel = label
  ) %>%
  dplyr::select(
    label,
    dg_cn,
    dg_cn0,
    dg_cnA,
    btwn,
    btwn0,
    btwnA,
    closeness,
    closeness0,
    closenessA,
    pgrnk,
    pgrnk0,
    pgrnkA,
    city
  )

tal_cor <- tal_sum %>%
  dplyr::mutate(
    dg_cn = dg_cn / max(dg_cn, na.rm = TRUE),
    dg_cnA = dg_cnA / max(dg_cnA, na.rm = TRUE),
    dg_cn0 = dg_cn0 / max(dg_cn0, na.rm = TRUE),
    btwn = btwn / max(btwn, na.rm = TRUE),
    btwnA = btwnA / max(btwnA, na.rm = TRUE),
    btwn0 = btwn0 / max(btwn0, na.rm = TRUE),
    closeness = closeness / max(closeness, na.rm = TRUE),
    closenessA = closenessA / max(closenessA, na.rm = TRUE),
    closeness0 = closeness0 / max(closeness0, na.rm = TRUE),
    pgrnk = pgrnk / max(pgrnk, na.rm = TRUE),
    pgrnkA = pgrnkA / max(pgrnkA, na.rm = TRUE),
    pgrnk0 = pgrnk0 / max(pgrnk0, na.rm = TRUE),
    city = "Talent",
    label = Tlabel
  ) %>%
  dplyr::select(
    label,
    dg_cn,
    dg_cn0,
    dg_cnA,
    btwn,
    btwn0,
    btwnA,
    closeness,
    closeness0,
    closenessA,
    pgrnk,
    pgrnk0,
    pgrnkA,
    city
  )


ash_cor_l <- pivot_longer(
  ash_cor,
  cols = c(
    dg_cn,
    dg_cn0,
    dg_cnA,
    btwn,
    btwn0,
    btwnA,
    closeness,
    closeness0,
    closenessA,
    pgrnk,
    pgrnk0,
    pgrnkA,
  ),
  values_to = "value"
) %>%
  mutate(
    type = case_when(
      endsWith(name, "A") ~ "Gene alternative",
      endsWith(name, "0") ~ "Gene null",
      TRUE ~ "Spatial"
    ),
    metric = sub("(0|A)$", "", name)
  )

tal_cor_l <- pivot_longer(
  tal_cor,
  cols = c(
    dg_cn,
    dg_cn0,
    dg_cnA,
    btwn,
    btwn0,
    btwnA,
    closeness,
    closeness0,
    closenessA,
    pgrnk,
    pgrnk0,
    pgrnkA,
  ),
  values_to = "value"
) %>%
  mutate(
    type = case_when(
      endsWith(name, "A") ~ "Gene alternative",
      endsWith(name, "0") ~ "Gene null",
      TRUE ~ "Spatial"
    ),
    metric = sub("(0|A)$", "", name)
  )

tal_cor_l <- tal_cor_l %>%
  dplyr::mutate(
    label = ifelse(
      is.na(label), "T29", label
    )
  )

comb <- rbind(tal_cor_l, ash_cor_l) %>%
  dplyr::mutate(metric = ifelse(
    metric == "btwn",
    "Betweenness",
    ifelse(
      metric == "closeness",
      "Closeness",
      ifelse(
        metric == "dg_cn",
        "Degree",
        ifelse(metric == "pgrnk", "Pagerank", metric)
      )
    )
  ))

comb <- comb[str_order(comb$label, numeric = TRUE), , drop = FALSE]

# COMPARING NETWORK METRICS #####
ggplot(comb, aes(
  x = fct_inorder(label),
  y = value,
  color = metric,
  shape = metric
)) + # HAVE TO USE FCT_INORDER OR X AXIS IS IN A RANDOM ORDER!
  geom_point(size = 2, alpha = 0.5) +
  scale_shape_manual(values = c(
    Degree = 16,
    Betweenness = 17,
    Closeness = 15,
    Pagerank = 3
  )) +
  facet_wrap( ~ city + type, scales = "free_x") +
  labs(
    x = "Site",
    y = "Value (normalized to 1)",
    shape = "Metric",
    color = "Metric"
  ) +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    size = 5
  ))

get_cors_table <- function(df, type) {
  metrics <- c("dg_cn", "btwn", "closeness", "pgrnk")
  
  map_df(metrics, function(m) {
    mat <- df %>%
      dplyr::select(all_of(c(m, paste0(m, "0"), paste0(m, "A")))) %>%
      psych::corr.test(method = type)
    
    data.frame(
      Metric = m,
      `Spatial vs 0` = round(mat$r[1, 2], digits = 2),
      `S0p` = round(mat$p[1, 2], digits = 3),
      `Spatial vs A` = round(mat$r[1, 3], digits = 2),
      `SAp` = round(mat$p[1, 3], digits = 3),
      `0 vs A` = round(mat$r[2, 3], digits = 2),
      `0 vs A` = round(mat$p[2, 3], digits = 3)
    )
  })
}
get_cors_table(ash_cor, "pearson")
get_cors_table(tal_cor, "pearson")
