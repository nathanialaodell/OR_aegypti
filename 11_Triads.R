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

# set up############
urel <- c("First", "Second")
pnw_love <- PNWColors::pnw_palette(name = "Sunset2",
                                   n = length(urel),
                                   type = "continuous")
pal <- setNames(pnw_love, urel)

collect_triads <- function(vids, isoclass) {
  if (isoclass == target_class) {
    found_triplets <- c(found_triplets, list(vids))
  }
  return(FALSE) 
}

or_dat_f_c <- readRDS(here::here("Objects/Talent Trap Counts.rds")) %>%
  dplyr::select(Site, Date, epi_wk, Total)

or_sf_cnct <- read_sf(here::here("Objects/Talent connectivity NA omitted.shp"))

or_sf_cnct <- or_sf_cnct[!(or_sf_cnct$From == or_sf_cnct$To), ]

Tobfuscate <- readRDS(here::here("Objects/Talent obfuscate.rds"))
Aobfuscate <- readRDS(here::here("Objects/Ashland obfuscate.rds"))

or_sf_cnct <- or_sf_cnct %>%
  left_join(Tobfuscate, by = c("From" = "Tsite")) %>%
  dplyr::rename(From_label = "Tlabel") %>%
  left_join(Tobfuscate, by = c("To" = "Tsite")) %>%
  dplyr::rename(To_label = "Tlabel")

ash_sf <- readRDS(here::here("Objects/kinship/ash full graph.RDS")) %>%
  dplyr::mutate(
    date_diff = as.numeric(date_diff)
  ) %>%
  dplyr::filter(
    (date_diff >= 14 & Type == "First") | (date_diff >= 28 & Type == "Second")
  )

#### TRIAD ANALYSIS #####
ag1 <- igraph::graph_from_data_frame(ash_sf, directed = TRUE)
pattern <- make_graph(edges = c(1,2,  2,3,  1,3), directed = TRUE)
trio <- subgraph_isomorphisms(pattern, ag1, method = "vf2")

target_030T <- do.call(rbind, lapply(trio, function(m) {
  node_names <- V(ag1)$name[m]
  data.frame(v1 = node_names[1], v2 = node_names[2], v3 = node_names[3])
}))

triad_vertices <- unique(unlist(target_030T[, c("v1","v2","v3")]))
ag1_030T <- induced_subgraph(ag1, vids = triad_vertices)

# A -> C connections that are 2nd degree KING
ac_edges <- do.call(rbind, lapply(seq_len(nrow(target_030T)), function(i) {
  trio <- as.character(unlist(target_030T[i, c("v1", "v2", "v3")]))
  sub <- induced_subgraph(ag1, trio)
  
  outd <- igraph::degree(sub, mode = "out")
  ind  <- igraph::degree(sub, mode = "in")
  
  A <- V(sub)$name[outd == 2] # sends to both others
  C <- V(sub)$name[ind == 2]  # receives from both others
  
  # pull the actual A->C edge back out of the full graph
  eid <- get_edge_ids(ag1, c(A, C))
  
  data.frame(A = A, C = C, Type = E(ag1)$Type[eid], stringsAsFactors = FALSE)
}))

table(ac_edges$Type)
sum(ac_edges$Type == "Second", na.rm = TRUE)

V(ag1_030T)$Lat <- ifelse(is.na(ash_sf$From_Y[match(V(ag1_030T)$name, ash_sf$From_Ind)]),
                          ash_sf$To_Y[match(V(ag1_030T)$name, ash_sf$To_Ind)],
                          ash_sf$From_Y[match(V(ag1_030T)$name, ash_sf$From_Ind)])
V(ag1_030T)$Long <- ifelse(is.na(ash_sf$From_X[match(V(ag1_030T)$name, ash_sf$From_Ind)]),
                           ash_sf$To_X[match(V(ag1_030T)$name, ash_sf$To_Ind)],
                           ash_sf$From_X[match(V(ag1_030T)$name, ash_sf$From_Ind)])
E(ag1_030T)$color <- pal[E(ag1_030T)$Type]
V(ag1_030T)$label <- Aobfuscate$Alabel[match(str_extract(V(ag1_030T)$name, "^.*?(?=_[0-9]{1,2}[A-Z])"), Aobfuscate$Asite)]
V(ag1_030T)$label[8] <- "A14"
V(ag1_030T)$label[9] <- "A14"

n_triads <- nrow(target_030T)
n_col <- ceiling(sqrt(n_triads))
n_row <- ceiling(n_triads / n_col)

##### FIGURE CREATION (ASHLAND) #####
set.seed(123)

jitter_amount <- 0.0005
p_list <- list()
for (i in seq_len(n_triads)) {
  
  trio <- unlist(target_030T[i, c("v1", "v2", "v3")], use.names = FALSE)
  trio <- as.character(trio)
  
  sub  <- induced_subgraph(ag1_030T, trio)
  
  lo <- matrix(c(V(sub)$Long, V(sub)$Lat), ncol = 2)
  
  lo <- lo + matrix(
    runif(length(lo), -jitter_amount, jitter_amount),
    ncol = 2
  )
  
  p_list[[i]] <- ggraph(sub, layout = "manual", x = lo[, 1], y = lo[, 2]) +
    geom_edge_link(aes(edge_color = I(color)), # have to use I() otherwise can't read hex codes
                   arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
                   end_cap = circle(2, "mm")
    ) +
    geom_node_point(
      color = "orange",
      shape = 16, 
      size = 1
    ) +
    geom_node_text(
      aes(label = paste0(label, "_", str_extract(name, "[0-9]{1,2}[A-Z]{1}"))), repel = TRUE, size = 3.75,
      max.overlaps = Inf, box.padding = 1,
      segment.size = 0.2, segment.color = "black"
    ) +
    theme_void()
}

wrap_plots(p_list)

ABCs <- c(1,8,11,13:14, 16, 20, 21, 23)

p_list <- list()
for (i in ABCs) {
  
  trio <- unlist(target_030T[i, c("v1", "v2", "v3")], use.names = FALSE)
  trio <- as.character(trio)
  
  sub  <- induced_subgraph(ag1_030T, trio)
  
  lo <- matrix(c(V(sub)$Long, V(sub)$Lat), ncol = 2)
  
  lo <- lo + matrix(
    runif(length(lo), -jitter_amount, jitter_amount),
    ncol = 2
  )
  
  p_list[[i]] <- ggraph(sub, layout = "manual", x = lo[, 1], y = lo[, 2]) +
    geom_edge_link(aes(edge_color = I(color)), # have to use I() otherwise can't read hex codes
      arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
      end_cap = circle(2, "mm")
    ) +
    geom_node_point(
      color = "orange",
      shape = 16, 
      size = 1
    ) +
    geom_node_text(
      aes(label = paste0(label, "_", str_extract(name, "[0-9]{1,2}[A-Z]{1}"))), repel = TRUE, size = 3.75,
      max.overlaps = Inf, box.padding = 2.5,
      segment.size = 0.2, segment.color = "black"
    ) +
    theme_void()
}
p_list <- Filter(Negate(is.null), p_list)

wrap_plots(p_list, ncol = 3)

############ talent ##########
tal_sf <- readRDS(here::here("Objects/kinship/tal full graph.RDS")) %>%
  dplyr::mutate(
    date_diff = as.numeric(date_diff)
  ) %>%
  dplyr::filter(
    (date_diff >= 14 & Type == "First") | (date_diff >= 28 & Type == "Second")
  )

tg1 <- igraph::graph_from_data_frame(tal_sf, directed = TRUE)
trio <- subgraph_isomorphisms(pattern, tg1, method = "vf2")

target_030T <- do.call(rbind, lapply(trio, function(m) {
  node_names <- V(tg1)$name[m]
  data.frame(v1 = node_names[1], v2 = node_names[2], v3 = node_names[3])
}))

triad_vertices <- unique(unlist(target_030T[, c("v1","v2","v3")]))
tg1_030T <- induced_subgraph(tg1, vids = triad_vertices)

# A -> C connections that are 2nd degree KING
ac_edges <- do.call(rbind, lapply(seq_len(nrow(target_030T)), function(i) {
  trio <- as.character(unlist(target_030T[i, c("v1", "v2", "v3")]))
  sub <- induced_subgraph(tg1, trio)
  
  outd <- degree(sub, mode = "out")
  ind  <- degree(sub, mode = "in")
  
  A <- V(sub)$name[outd == 2]   # sends to both others
  C <- V(sub)$name[ind == 2]    # receives from both others
  
  # pull the actual A->C edge back out of the full graph
  eid <- get_edge_ids(tg1, c(A, C))
  
  data.frame(A = A, C = C, Type = E(tg1)$Type[eid], stringsAsFactors = FALSE)
}))

V(tg1_030T)$Lat <- ifelse(is.na(tal_sf$From_Y[match(V(tg1_030T)$name, tal_sf$From_Ind)]),
                          tal_sf$To_Y[match(V(tg1_030T)$name, tal_sf$To_Ind)],
                          tal_sf$From_Y[match(V(tg1_030T)$name, tal_sf$From_Ind)])
V(tg1_030T)$Long <- ifelse(is.na(tal_sf$From_X[match(V(tg1_030T)$name, tal_sf$From_Ind)]),
                           tal_sf$To_X[match(V(tg1_030T)$name, tal_sf$To_Ind)],
                           tal_sf$From_X[match(V(tg1_030T)$name, tal_sf$From_Ind)])
E(tg1_030T)$color <- pal[E(tg1_030T)$Type]

lo2 <- matrix(c(V(tg1_030T)$Long, V(tg1_030T)$Lat), ncol = 2) %>%
  na.omit()
bbox <- st_bbox(or_sf_cnct)

igraph::plot.igraph(
  tg1_030T,
  layout = lo2,
  rescale = FALSE,
  asp = 0,
  vertex.color = "orange",
  vertex.size = 0.025,
  vertex.shape = "circle",
  vertex.frame.color = "black",
  vertex.label.color = "black",
  vertex.label.cex = 0.65,
  edge.arrow.size = 0.75,
  vertex.label = V(tg1_030T)$name,
  edge.label.color = "black",
  xlim = c(bbox[1], bbox[3]),
  ylim = c(bbox[2], bbox[4])
)
