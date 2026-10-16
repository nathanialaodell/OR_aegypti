library(ape)
library(adegenet)
library(genetics)
library(ade4)
library(RColorBrewer)
library(hierfstat)
library(pegas)
library(vegan)
library(poppr)
library(genetics)
library(hierfstat)
library(tidyr)
library(igraph)
library(stringr)
library(dplyr)
library(patchwork)
library(ggtext)
source(here::here("From_to_functions.R"))
source(here::here("Creating_networks.R"))
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")
library(StAMPP)
library(dartR)

# for plot layout (Talent only right now) ########
or_dat_f_c <- readRDS(here::here("Objects/Talent Trap Counts.rds")) %>%
  dplyr::select(Site, Date, epi_wk, Total)

# talent shapefile
or_sf_cnct <- read_sf(here::here("Objects/Talent connectivity NA omitted-copy.shp")) 

or_sf_cnct <- or_sf_cnct[!(or_sf_cnct$From == or_sf_cnct$To), ]

Tsite <- unique(c(or_sf_cnct$From, or_sf_cnct$To))

Tlabel <- paste("T", 1:length(Tsite), sep = "")

Tobfuscate <- data.frame(
  Tsite,
  Tlabel
)

saveRDS(Tobfuscate, here::here("Objects/Talent obfuscate.rds"))

ash <- read_sf(here::here("Objects/Ashland connectivity NA omitted.shp"))%>%
  dplyr::filter(Dist <= 700)
ash <- ash[!(ash$From == ash$To), ]
Asite <- unique(c(ash$From, ash$To))

Alabel <- paste("A", 1:length(Asite), sep = "")

Aobfuscate <- data.frame(
  Asite,
  Alabel
)

######## tube codes ######
ash_sf <- readRDS(here::here("Objects/kinship/ash full graph.RDS"))

tal_sf <- readRDS(here::here("Objects/kinship/tal full graph.RDS"))

mydataNA = read.genepop(
  here::here("Objects/kinship/Oregon_ALL_2026_GP.gen"),
  ncode = 3L,
  quiet = FALSE
)

ash_codes <- read.csv(here::here('Objects/kinship/Ashland_codes for predigree.csv')) %>%
  dplyr::mutate(Tube = paste(Tube, Sex, sep = "-")) %>%
  dplyr::select(Tube, Trap.location)

id_reg <- "[A-Z]{2,}[0-9]{0,}_[0-9]{1,}[A-Z]{0,}$"
tal_codes <- read.csv(here::here(
  'Objects/kinship/Talent_codes for predigree_copy.csv'
)) %>%
  dplyr::mutate(Tube = str_extract(Tube, id_reg)) %>%
  dplyr::select(Tube, Trap.location)

## ashland ######
# grouping by site
ash <- mydataNA[1:50]
names <- as_tibble(indNames(ash))
names <- left_join(names, ash_codes, by = c("value" = "Tube")) %>%
  dplyr::select(-value)
sites <- as.factor(names$Trap.location)

strata(ash) <- names
nameStrata(ash) <- ~ Site
setPop(ash) <- ~ Site
site_counts <- table(ash@pop)

ash2 <- dartR::gi2gl(ash, parallel = FALSE, verbose = NULL)

A.dist <- stamppFst(ash2, nboots = 99, percent = 95)
A.dist2 <- A.dist[[1]]
A.dist2[A.dist2 < 0] <- 1e-6
A.dist2[is.nan(A.dist2)] <- NA

A.dist[[1]][!lower.tri(A.dist[[1]])] <- NA
A.dist[[1]][is.nan(A.dist[[1]])] <- NA
A.dist[[1]][A.dist[[1]] < 0] <- 1e-6
A.dist[[1]][A.dist[[2]] > 0.05] <- NA

unique(A.dist[[1]])

ash_fst <- melt(A.dist[[1]], varnames = c("from", "to"), value.name = "fst")
ash_fst <- na.omit(ash_fst)
ash_fst <- ash_fst[as.character(ash_fst$from) < as.character(ash_fst$to), ]
rownames(ash_fst) <- 1:nrow(ash_fst)
ash_fst

ash_fst2 <- na.omit(melt(A.dist2, varnames = c("from", "to"), value.name = "fst"))
ash_fst2 <- ash_fst[as.character(ash_fst$from) < as.character(ash_fst$to), ]
rownames(ash_fst2) <- 1:nrow(ash_fst2)

ag0 <- graph_from_data_frame(
  ash_fst,
  directed = FALSE
)

E(ag0)$weight <- ash_fst$value

V(ag0)$label <- Aobfuscate$Alabel[match(V(ag0)$name, Aobfuscate$Asite)]
V(ag0)$label[12] <- "A14" # missing site in spatial data 

ag1 <- ag0

msn_graph <- igraph::mst(ag0)
E(msn_graph)$color <- "#D81B60"

adding_el <- ash_sf %>%
  dplyr::select(From_site, To_site, date_diff, Type) %>%
  dplyr::mutate(date_diff = as.numeric(date_diff))
adding_el <- adding_el %>%
  dplyr::filter(date_diff >= 14 & Type == "First") %>%
  distinct()

adding_el <- anti_join(
  adding_el %>%
    dplyr::select(-c(date_diff, Type)),
  igraph::as_data_frame(msn_graph) %>% dplyr::select(-weight),
  by = c("From_site" = "from", "To_site" = "to")
) %>%
  dplyr::group_by(From_site, To_site) %>%
  dplyr::mutate(
    weight = A.dist2[From_site, To_site][1]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    from = "From_site",
    to = "To_site"
  ) %>%
  dplyr::mutate(
    color = "#004D40",
    weight = ifelse(weight == 1e-6, NA, weight)
  ) %>%
  na.omit()

msn_df <- igraph::as_data_frame(msn_graph) %>%
  dplyr::mutate(
    color = "#D81B60"
  ) %>%
  dplyr::select(
    -value
  )

msn_graph1 <- graph_from_data_frame(rbind(msn_df, adding_el), directed = FALSE)
V(msn_graph1)$label <- Aobfuscate$Alabel[match(V(msn_graph1)$name, Aobfuscate$Asite)]
V(msn_graph1)$label[10] <- "A14"


V(msn_graph)$Lat <- ifelse(is.na(ash_sf$From_Y[match(V(msn_graph)$name, ash_sf$From_site)]),
                           ash_sf$To_Y[match(V(msn_graph)$name, ash_sf$To_site)],
                           ash_sf$From_Y[match(V(msn_graph)$name, ash_sf$From_site)])
V(msn_graph)$Long <- ifelse(is.na(ash_sf$From_X[match(V(msn_graph)$name, ash_sf$From_site)]),
                            ash_sf$To_X[match(V(msn_graph)$name, ash_sf$To_site)],
                            ash_sf$From_X[match(V(msn_graph)$name, ash_sf$From_site)])

lo <- matrix(c(V(msn_graph)$Long, V(msn_graph)$Lat), ncol = 2) %>%
  na.omit()
bbox <- st_bbox(read_sf(here::here("Objects/Ashland connectivity NA omitted.shp")))

V(msn_graph1)$Lat <- ifelse(is.na(ash_sf$From_Y[match(V(msn_graph1)$name, ash_sf$From_site)]),
                            ash_sf$To_Y[match(V(msn_graph1)$name, ash_sf$To_site)],
                            ash_sf$From_Y[match(V(msn_graph1)$name, ash_sf$From_site)])
V(msn_graph1)$Long <- ifelse(is.na(ash_sf$From_X[match(V(msn_graph1)$name, ash_sf$From_site)]),
                             ash_sf$To_X[match(V(msn_graph1)$name, ash_sf$To_site)],
                             ash_sf$From_X[match(V(msn_graph1)$name, ash_sf$From_site)])


### for summary table (08/25/2026) ####
ash_sum0 <- data.frame(
  Site = V(msn_graph)$name,
  dg_cn0 = as.numeric(igraph::degree(msn_graph, mode = "out")),
  btwn0 = as.numeric(igraph::betweenness(msn_graph, directed = FALSE)),
  closeness0 = as.numeric(igraph::closeness(msn_graph, mode = "all")),
  pgrnk0 = as.numeric(page_rank(msn_graph, directed = FALSE)$vector)
) %>%
  rbind(
    data.frame(
      Site = Site, # PII... obfuscated
      dg_cn0 = 0,
      btwn0 = 0,
      closeness0 = 0,
      pgrnk0 = 0
    )
  ) %>%
  arrange(desc(Site))

ash_sum1 <- data.frame(
  Site = V(msn_graph1)$name,
  dg_cnA = as.numeric(igraph::degree(msn_graph1, mode = "out")),
  btwnA = as.numeric(igraph::betweenness(msn_graph1, directed = FALSE)),
  closenessA = as.numeric(igraph::closeness(msn_graph1, mode = "all")),
  pgrnkA = as.numeric(page_rank(msn_graph1, directed = FALSE)$vector)
) %>%
  arrange(desc(Site)) %>%
  dplyr::select(
    -Site
  )

saveRDS(cbind(ash_sum0, ash_sum1), here::here('Objects/Ashland sum table part 2 (genomic centrality).RDS'))

######

lo2 <- matrix(c(V(msn_graph1)$Long, V(msn_graph1)$Lat), ncol = 2) %>%
  na.omit()

# for the legend..
shared_range <- range(c(E(msn_graph)$weight, E(msn_graph1)$weight), na.rm = TRUE)

p1 <- ggraph(msn_graph, layout = "manual", x = lo[,1], y = lo[,2]) +
  geom_edge_arc(aes(width = weight, alpha = weight, edge_colour = color), strength = 0.15) +
  scale_edge_width(range = c(1,0.3), breaks = c(0.04, 0.08, 0.12, 0.16, 0.20)) +
  scale_edge_alpha(range = c(1,0.3), breaks = c(0.04, 0.08, 0.12, 0.16, 0.20)) +
  scale_edge_colour_identity()+
  geom_node_point(size = 4, color = "#FFC107", show.legend = FALSE, shape = 20) +
  scale_size_continuous(range = c(2,7)) + 
  geom_node_text(
    aes(label = label), repel = TRUE, size = 5,
    max.overlaps = Inf, box.padding = 0.6,
    segment.size = 0.2, segment.color = "grey50"
  ) +
  theme_void()

p2 <- ggraph(msn_graph1, layout = "manual", x = lo2[,1], y = lo2[,2]) +
  geom_edge_arc(aes(width = weight, alpha = weight, edge_colour = color), strength = 0.15) +
  scale_edge_width(range = c(1,0.3), breaks = c(0.04, 0.08, 0.12, 0.16, 0.20)) +
  scale_edge_alpha(range = c(1,0.3), breaks = c(0.04, 0.08, 0.12, 0.16, 0.20)) +
  scale_edge_colour_identity() +
  geom_node_point(size = 4, color = "#FFC107", show.legend = FALSE, shape = 20) +
  scale_size_continuous(range = c(2,7)) + 
  geom_node_text(
    aes(label = label), repel = TRUE, size = 5,
    max.overlaps = Inf, box.padding = 0.6,
    segment.size = 0.2, segment.color = "grey50"
  ) +
  theme_void()

p1 + p2 +
  plot_layout(guides = "collect")


### talent #####

tal <- mydataNA[51:87]
names <- as_tibble(indNames(tal))
names <- left_join(names, tal_codes, by = c("value" = "Tube")) %>%
  dplyr::select(-value)
sites <- as.factor(names$Trap.location)

strata(tal) <- names
nameStrata(tal) <- ~ Site
setPop(tal) <- ~ Site
site_counts <- table(tal@pop)


tal2 <- dartR::gi2gl(tal, parallel = FALSE, verbose = NULL)

T.dist <- stamppFst(tal2, nboots = 99, percent = 95)
T.dist2 <- T.dist[[1]]
T.dist2[T.dist2 < 0] <- 1e-6
T.dist2[is.nan(T.dist2)] <- NA

T.dist[[1]][!lower.tri(T.dist[[1]])] <- NA
T.dist[[1]][is.nan(T.dist[[1]])] <- NA
T.dist[[1]][T.dist[[1]] < 0] <- 1e-6
T.dist[[1]][T.dist[[2]] > 0.05] <- NA

tal_fst <- melt(T.dist[[1]], varnames = c("from", "to"), value.name = "fst")
tal_fst <- na.omit(tal_fst)
tal_fst <- tal_fst[as.character(tal_fst$from) < as.character(tal_fst$to), ]
rownames(tal_fst) <- 1:nrow(tal_fst)
tal_fst

tal_fst2 <- na.omit(melt(T.dist2, varnames = c("from", "to"), value.name = "fst"))
tal_fst2 <- tal_fst[as.character(tal_fst$from) < as.character(tal_fst$to), ]
rownames(tal_fst2) <- 1:nrow(tal_fst2)

tg0 <- graph_from_data_frame(
  tal_fst,
  directed = FALSE
)

E(tg0)$weight <- tal_fst$value
V(tg0)$label <- Tobfuscate$Tlabel[match(V(tg0)$name, Tobfuscate$Tsite)]
tg1 <- tg0

msn_graph <- igraph::mst(tg0)
E(msn_graph)$color <- "#D81B60"

adding_el <- tal_sf %>%
  dplyr::select(From_site, To_site, date_diff, Type) %>%
  dplyr::mutate(date_diff = as.numeric(date_diff))
adding_el <- adding_el %>%
  dplyr::filter(date_diff >= 14 & Type == "First") %>%
  distinct()

adding_el <- anti_join(
  adding_el %>%
    dplyr::select(-c(date_diff, Type)),
  igraph::as_data_frame(msn_graph) %>% dplyr::select(-weight),
  by = c("From_site" = "from", "To_site" = "to")
) %>%
  dplyr::group_by(From_site, To_site) %>%
  dplyr::mutate(
    weight = T.dist2[From_site, To_site]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    from = "From_site",
    to = "To_site"
  ) %>%
  dplyr::mutate(
    color = "#004D40",
    weight = ifelse(weight == 1e-6, NA, weight)
  ) %>%
  na.omit()

msn_df <- igraph::as_data_frame(msn_graph) %>%
  dplyr::mutate(
    color = "#D81B60"
  ) %>%
  dplyr::select(
    -value
  )

msn_graph1 <- graph_from_data_frame(rbind(msn_df, adding_el), directed = FALSE)
V(msn_graph1)$label <- Tobfuscate$Tlabel[match(V(msn_graph1)$name, Tobfuscate$Tsite)]

V(msn_graph)$Lat <- ifelse(is.na(or_sf_cnct$From_Y[match(V(msn_graph)$name, or_sf_cnct$From)]),
                           or_sf_cnct$To_Y[match(V(msn_graph)$name, or_sf_cnct$To)],
                           or_sf_cnct$From_Y[match(V(msn_graph)$name, or_sf_cnct$From)])
V(msn_graph)$Long <- ifelse(is.na(or_sf_cnct$From_X[match(V(msn_graph)$name, or_sf_cnct$From)]),
                            or_sf_cnct$To_X[match(V(msn_graph)$name, or_sf_cnct$To)],
                            or_sf_cnct$From_X[match(V(msn_graph)$name, or_sf_cnct$From)])


### for summary table (08/25/2026) ####
tal_sum0 <- data.frame(
  Site = V(msn_graph)$name,
  dg_cn0 = as.numeric(igraph::degree(msn_graph, mode = "all")),
  btwn0 = as.numeric(igraph::betweenness(msn_graph, directed = FALSE)),
  closeness0 = as.numeric(igraph::closeness(msn_graph, mode = "all")),
  pgrnk0 = as.numeric(page_rank(msn_graph, directed = FALSE)$vector)
) %>%
  arrange(desc(Site))

tal_sum1 <- data.frame(
  Site = V(msn_graph1)$name,
  dg_cnA = as.numeric(igraph::degree(msn_graph1, mode = "all")),
  btwnA = as.numeric(igraph::betweenness(msn_graph1, directed = FALSE)),
  closenessA = as.numeric(igraph::closeness(msn_graph1, mode = "all")),
  pgrnkA = as.numeric(page_rank(msn_graph1, directed = FALSE)$vector)
) %>%
  arrange(desc(Site)) %>%
  dplyr::select(
    -Site
  )

saveRDS(cbind(tal_sum0, tal_sum1), here::here('Objects/Talent sum table part 2 (genomic centrality).RDS'))

lo <- matrix(c(V(msn_graph)$Long, V(msn_graph)$Lat), ncol = 2) %>%
  na.omit()
bbox <- st_bbox(or_sf_cnct)

V(msn_graph1)$Lat <- ifelse(is.na(or_sf_cnct$From_Y[match(V(msn_graph1)$name, or_sf_cnct$From)]),
                            or_sf_cnct$To_Y[match(V(msn_graph1)$name, or_sf_cnct$To)],
                            or_sf_cnct$From_Y[match(V(msn_graph1)$name, or_sf_cnct$From)])
V(msn_graph1)$Long <- ifelse(is.na(or_sf_cnct$From_X[match(V(msn_graph1)$name, or_sf_cnct$From)]),
                             or_sf_cnct$To_X[match(V(msn_graph1)$name, or_sf_cnct$To)],
                             or_sf_cnct$From_X[match(V(msn_graph1)$name, or_sf_cnct$From)])


lo2 <- matrix(c(V(msn_graph1)$Long, V(msn_graph1)$Lat), ncol = 2) %>%
  na.omit()
shared_range <- range(c(E(msn_graph)$weight, E(msn_graph1)$weight), na.rm = TRUE)

# both plotted #####
p3 <- ggraph(msn_graph, layout = "manual", x = lo[,1], y = lo[,2]) +
  geom_edge_arc(aes(width = weight, alpha = weight, edge_colour = color), strength = 0.15) +
  scale_edge_width(range = c(1,0.3), breaks = c(0.04, 0.08, 0.12, 0.16, 0.20)) +
  scale_edge_alpha(range = c(1,0.3), breaks = c(0.04, 0.08, 0.12, 0.16, 0.20)) +
  scale_edge_colour_identity()+
  geom_node_point(size = 4, color = "#FFC107", show.legend = FALSE, shape = 20) +
  scale_size_continuous(range = c(2,7)) + 
  geom_node_text(
    aes(label = label), repel = TRUE, size = 5,
    max.overlaps = Inf, box.padding = 0.6,
    segment.size = 0.2, segment.color = "grey50"
  ) +
  theme_void()

p4 <- ggraph(msn_graph1, layout = "manual", x = lo2[,1], y = lo2[,2]) +
  geom_edge_arc(aes(width = weight, alpha = weight, edge_colour = color), strength = 0.15) +
  scale_edge_width(range = c(1,0.3), breaks = c(0.04, 0.08, 0.12, 0.16, 0.20)) +
  scale_edge_alpha(range = c(1,0.3), breaks = c(0.04, 0.08, 0.12, 0.16, 0.20)) +
  scale_edge_colour_identity() +
  geom_node_point(size = 4, color = "#FFC107", show.legend = FALSE, shape = 20) +
  scale_size_continuous(range = c(2,7)) + 
  geom_node_text(
    aes(label = label), repel = TRUE, size = 5,
    max.overlaps = Inf, box.padding = 0.6,
    segment.size = 0.2, segment.color = "grey50"
  ) +
  theme_void()

(p1 + p2)/(p3+p4) +
  plot_layout(guides = "collect")
