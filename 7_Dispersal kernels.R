library(kindisperse)
library(ape)
library(igraph)
library(timeordered)
library(dplyr)
library(here)
library(flextable)
library(MMWRweek)
source(here::here("From_to_functions.R"))
source(here::here("Creating_networks.R"))

ash <- readRDS(here::here("Objects/kinship/ash full graph.RDS")) %>%
  gene_dist() %>%
  dplyr::mutate(
    Dist = as.numeric(Dist),
    date_diff = as.numeric(date_diff),
    epi_wk = MMWRweek(To_Date)$MMWRweek,
    TopRel = ifelse(
      # a guesstimate
      date_diff <= 13 & Type == "First",
      "FS",
      ifelse(
        between(date_diff, 14, 28) & Type == "First",
        "PO",
        ifelse(
          date_diff >= 28 & Type == "Second",
          "GG",
          "other"
        )
      )
    ),
    Type = ifelse(
      TopRel == "Second" & date_diff < 28, "other", Type
    )
  )

# null network (undirected)
ag0 <- gene_graph_wk_it_fun(ash, max(ash$epi_wk), direct = FALSE)
agA <- gene_graph_wk_it_fun(ash, max(ash$epi_wk), alternative = TRUE)

ash_di <- ash %>%
  dplyr::select(From, To, TopRel, Dist) %>%
  dplyr::rename(
    id1 = "From",
    id2 = "To",
    kinship = "TopRel",
    distance = "Dist"
  ) %>%
  dplyr::mutate(distance = as.numeric(distance)) %>%
  st_drop_geometry()  %>%
  dplyr::filter(kinship != "other")

FS <- df_to_kinpair(ash_di, kinship = "FS")
GP <- df_to_kinpair(ash_di, kinship = "GG")
PO <- df_to_kinpair(ash_di, kinship = "PO")

dyad_l <- list(FS, GP, PO)

relationships <- c("full", "grand", "parent")
names(dyad_l) <- relationships
comp <- c(2, 1, 1)

sigma_l <- vector(mode = "list", length = length(relationships))
names(sigma_l) <- relationships

for (i in 1:length(dyad_l)) {
  set.seed(123)
  sigma_l[[i]] <- axpermute(dyad_l[[i]], composite = comp[i])
}

### ASHLAND t.2 #######
kerns <- do.call(rbind, sigma_l) %>% as.data.frame()
kerns

fs_mean <- sigma_l$full[['mean']]
gg_mean <- sigma_l$grand[['mean']]
po_mean <- sigma_l$parent[['mean']]

kern_df <- data.frame("FS" = fs_mean, "GG" = gg_mean, "PO" = po_mean)
kern_df

fs_edge <- E(ag0)[TopRel == "FS" & Dist > fs_mean]
gg_edge <- E(ag0)[TopRel == "GG" & Dist > gg_mean]
po_edge <- E(ag0)[TopRel == "PO" & Dist > po_mean]

el_edge <- c(fs_edge, gg_edge, po_edge)

ash_di[ash_di$kinship == "FS" & ash_di$distance > fs_mean, ]
ash_di[ash_di$kinship == "GG" & ash_di$distance > gg_mean, ]
ash_di[ash_di$kinship == "PO" & ash_di$distance > po_mean, ]

par(mfrow = c(1, 1))
pn_fun_kin2(ash, ag0, NULL, lab = FALSE)

pn_fun_kin3(
  ash,
  subgraph_from_edges(ag0, E(ag0)[TopRel == "FS"], delete.vertices = TRUE),
  NULL,
  legend = FALSE,
  lab = FALSE,
  buffer_m = fs_mean,
  type = "FS"
)

fs_edge <- E(agA)[TopRel == "FS" & Dist > fs_mean]
gg_edge <- E(agA)[TopRel == "GG" & Dist > gg_mean]
po_edge <- E(agA)[TopRel == "PO" & Dist > po_mean]

el_edge <- c(fs_edge, gg_edge, po_edge)

### build last part of table for node summary #####
ash_di <- ash %>%
  dplyr::select(From, To, From_site, To_site, TopRel, Dist) %>%
  dplyr::rename(
    id1 = From,
    id2 = To,
    kinship = TopRel,
    distance = Dist
  ) %>%
  dplyr::mutate(distance = as.numeric(distance)) %>%
  st_drop_geometry() %>%
  dplyr::filter(kinship != "other")

kernel_table <- ash_di %>%
  dplyr::mutate(
    threshold = dplyr::case_when(
      kinship == "FS" ~ kern_df$FS,
      kinship == "GG" ~ kern_df$GG,
      kinship == "PO" ~ kern_df$PO
    ),
    exceed = distance > threshold
  ) %>%
  dplyr::filter(exceed) %>%
  tidyr::pivot_longer(
    cols = c(From_site, To_site),
    names_to = NULL,
    values_to = "Site"
  ) %>%
  dplyr::distinct(id1, id2, kinship, Site) %>%
  dplyr::count(Site, kinship) %>%
  tidyr::pivot_wider(
    names_from = kinship,
    values_from = n,
    values_fill = 0
  )

saveRDS(kernel_table, here::here('Objects/Ashland sum table part 3 (dispersal).RDS'))

### talent ####
tal <- readRDS(here::here("Objects/kinship/tal full graph.RDS")) %>%
  gene_dist() %>%
  dplyr::mutate(
    Dist = as.numeric(Dist),
    date_diff = as.numeric(date_diff),
    epi_wk = MMWRweek(To_Date)$MMWRweek,
    TopRel = ifelse(
      # a guesstimate
      date_diff <= 13 & Type == "First",
      "FS",
      ifelse(
        between(date_diff, 14, 28) & Type == "First",
        "PO",
        ifelse(
          date_diff >= 28 & Type == "Second",
          "GG",
          "other"
        )
      )
    ),
    Type = ifelse(
      TopRel == "Second" & date_diff < 28, "other", Type
    )
  )

# null network (undirected)
tg0 <- gene_graph_wk_it_fun(tal, max(tal$epi_wk), direct = FALSE)
tgA <- gene_graph_wk_it_fun(tal, max(tal$epi_wk), alternative = TRUE)

tal_di <- tal %>%
  dplyr::filter(!(TopRel %in% c("2nd", "MZ"))) %>%
  dplyr::select(From, To, TopRel, Dist) %>%
  dplyr::rename(
    id1 = "From",
    id2 = "To",
    kinship = "TopRel",
    distance = "Dist"
  ) %>%
  dplyr::mutate(distance = as.numeric(distance)) %>%
  st_drop_geometry()  %>%
  dplyr::filter(kinship != "other")

FS <- df_to_kinpair(tal_di, kinship = "FS")
GP <- df_to_kinpair(tal_di, kinship = "GG")
PO <- df_to_kinpair(tal_di, kinship = "PO")

dyad_l <- list(FS, GP, PO)

relationships <- c("full", "grand", "parent")
names(dyad_l) <- relationships
comp <- c(2, 1, 1)

sigma_l <- vector(mode = "list", length = length(relationships))
names(sigma_l) <- relationships

for (i in 1:length(dyad_l)) {
  set.seed(123)
  sigma_l[[i]] <- axpermute(dyad_l[[i]], composite = comp[i])
}

kerns <- do.call(rbind, sigma_l) %>% as.data.frame()
kerns

fs_mean <- sigma_l$full[['mean']]
gg_mean <- sigma_l$grand[['mean']]
po_mean <- sigma_l$parent[['mean']]

kern_df <- data.frame("FS" = fs_mean, "GG" = gg_mean, "PO" = po_mean)

fs_edge <- E(tg0)[TopRel == "FS" & Dist > fs_mean]
gg_edge <- E(tg0)[TopRel == "GG" & Dist > gg_mean]
po_edge <- E(tg0)[TopRel == "PO" & Dist > po_mean]

el_edge <- c(fs_edge, gg_edge, po_edge)

par(mfrow = c(1, 1))
pn_fun_kin2(tal, tg0, NULL, lab = FALSE)

fs_edge <- E(tgA)[TopRel == "FS" & Dist > fs_mean]
gg_edge <- E(tgA)[TopRel == "GG" & Dist > gg_mean]
po_edge <- E(tgA)[TopRel == "PO" & Dist > po_mean]

el_edge <- c(fs_edge, gg_edge, po_edge)

par(mfrow = c(1, 1))
pn_fun_kin3(
  ash,
  subgraph.edges(tg0, E(tg0)[TopRel == "FS"], delete.vertices = TRUE),
  NULL,
  legend = FALSE,
  lab = FALSE,
  buffer_m = fs_mean,
  type = "FS"
)

tal_di <- tal %>%
  dplyr::select(From, To, From_site, To_site, TopRel, Dist) %>%
  dplyr::rename(
    id1 = From,
    id2 = To,
    kinship = TopRel,
    distance = Dist
  ) %>%
  dplyr::mutate(distance = as.numeric(distance)) %>%
  st_drop_geometry() %>%
  dplyr::filter(kinship != "other")

kernel_table <- tal_di %>%
  dplyr::mutate(
    threshold = dplyr::case_when(
      kinship == "FS" ~ kern_df$FS,
      kinship == "GG" ~ kern_df$GG,
      kinship == "PO" ~ kern_df$PO
    ),
    exceed = distance > threshold
  ) %>%
  dplyr::filter(exceed) %>%
  tidyr::pivot_longer(
    cols = c(From_site, To_site),
    names_to = NULL,
    values_to = "Site"
  ) %>%
  dplyr::distinct(id1, id2, kinship, Site) %>%
  dplyr::count(Site, kinship) %>%
  tidyr::pivot_wider(
    names_from = kinship,
    values_from = n,
    values_fill = 0
  )
kernel_table

saveRDS(kernel_table, here::here('Objects/Talent sum table part 3 (dispersal).RDS'))
