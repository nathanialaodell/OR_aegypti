library(sf)
library(igraph)
library(timeordered)
library(dplyr)
library(pals)
library(TNC)
source(here::here('Creating_networks.R'))

tal <- read_sf(here::here("Objects/Talent connectivity NA omitted.shp"))%>%
  dplyr::filter(Dist <= 700)
tal <- tal[!(tal$From == tal$To), ]

Tobfuscate <- readRDS(here::here("Objects/Talent obfuscate.rds"))

# coding convention seen in Figure 3
tal <- tal %>%
  left_join(Tobfuscate, by = c("From" = "Tsite")) %>%
  dplyr::rename(From_label = "Tlabel") %>%
  left_join(Tobfuscate, by = c("To" = "Tsite")) %>%
  dplyr::rename(To_label = "Tlabel")

ash <- read_sf(here::here("Objects/Ashland connectivity NA omitted.shp"))%>%
  dplyr::filter(Dist <= 700)
ash <- ash[!(ash$From == ash$To), ]
Asite <- unique(c(ash$From, ash$To))

Alabel <- paste("A", 1:length(Asite), sep = "")

Aobfuscate <- data.frame(
  Asite,
  Alabel
)
Aobfuscate
saveRDS(Aobfuscate, here::here("Objects/Ashland obfuscate.rds"))

ash <- ash %>%
  left_join(Aobfuscate, by = c("From" = "Asite")) %>%
  dplyr::rename(From_label = "Alabel") %>%
  left_join(Aobfuscate, by = c("To" = "Asite")) %>%
  dplyr::rename(To_label = "Alabel")

# get time-ordered networks
g_33 <- graph_wk_with_nn_fun(tal, 33)
g_34 <- graph_wk_with_nn_fun(tal, 34)
g_35 <- graph_wk_with_nn_fun(tal, 35)
g_36 <- graph_wk_with_nn_fun(tal, 36)
g_37 <- graph_wk_with_nn_fun(tal, 37)

# and a static
git_37 <- graph_wk_it_fun(tal, 37)

# there are two sites that don't appear until the final week
# they need to be added to the previous time-ordered network for figures
# and temporal network centrality measures
names = V(git_37)$name[27:28]
Lats = V(git_37)$Lat[27:28]
Longs = V(git_37)$Long[27:28]

g_33 <- add_vertices(
  g_33,
  nv = 2,
  name = names,
  Lat = Lats,
  Long = Longs
)
g_34 <- add_vertices(
  g_34,
  nv = 2,
  name = names,
  Lat = Lats,
  Long = Longs
)
g_35 <- add_vertices(
  g_35,
  nv = 2,
  name = names,
  Lat = Lats,
  Long = Longs
)
g_36 <- add_vertices(
  g_36,
  nv = 2,
  name = names,
  Lat = Lats,
  Long = Longs
)

ash_wks <- unique(ash$epi_wk)
length(ash_wks)

ga_1 <- graph_wk_with_nn_fun(ash, ash_wks[1])
ga_2 <- graph_wk_with_nn_fun(ash, ash_wks[2])
ga_3 <- graph_wk_with_nn_fun(ash, ash_wks[3])
ga <- graph_wk_it_fun(ash, ash_wks[3])

names_36 <- V(ga_3)[!(V(ga_3)$name %in% V(ga_1)$name)]$name
names_37 <- V(ga_3)[!(V(ga_3)$name %in% V(ga_2)$name)]$name
Lats_36 <- V(ga_3)[!(V(ga_3)$name %in% V(ga_1)$name)]$Lat
Lats_37 <- V(ga_3)[!(V(ga_3)$name %in% V(ga_2)$name)]$Lat
Long_36 <- V(ga_3)[!(V(ga_3)$name %in% V(ga_1)$name)]$Long
Long_37 <- V(ga_3)[!(V(ga_3)$name %in% V(ga_2)$name)]$Long

ga_1 <- add_vertices(
  ga_1,
  nv = length(names_36),
  name = names_36,
  Lat = Lats_36,
  Long = Long_36
)

ga_2 <- add_vertices(
  ga_2,
  nv = length(names_37),
  name = names_37,
  Lat = Lats_37,
  Long = Long_37
)

# for the spread analysis, it's easiest to just use the `timeordered` package syntax
or_sf_cnct_time <- tal %>%
  mutate(
    VertexFrom = From_label,
    VertexTo = To_label,
    TimeStart = epi_wk - 1, # technically epi_wk is a week ahead. Go back 1
    TimeStop = epi_wk
  ) %>%
  distinct() %>%
  dplyr::select(
    VertexFrom, VertexTo, TimeStart, TimeStop, Dist
  ) %>%
  st_drop_geometry()
g <- generatetonetwork(or_sf_cnct_time)

ash_sf_cnct_time <- ash %>%
  mutate(
    VertexFrom = From_label,
    VertexTo = To_label,
    TimeStart = epi_wk - 1,
    TimeStop = epi_wk
  ) %>%
  distinct() %>%
  dplyr::select(
    VertexFrom, VertexTo, TimeStart, TimeStop, Dist
  ) %>%
  st_drop_geometry()
gat <- generatetonetwork(ash_sf_cnct_time)

time_ordered <- list(
  "or_sf_cnct_time" = g,
  "ash_sf_cnct_time" = gat
)

time_df <- list(
  "or_sf_cnct_time" = or_sf_cnct_time,
  "ash_sf_cnct_time" = ash_sf_cnct_time
)

# weight g by distance
time_ordered <- lapply(names(time_df), function(graph_name){
  x <- time_df[[graph_name]]
  g_obj <- time_ordered[[graph_name]]
  
  weights <- x %>%
    dplyr::select(VertexFrom, VertexTo, Dist)
  
  connections <- data.frame(
    From = E(g_obj)$VertexFrom,
    To = E(g_obj)$VertexTo,
    stringsAsFactors = FALSE
  )
  
  matched_idx <- match(
    paste(connections$From, connections$To, sep = "_"),
    paste(weights$VertexFrom, weights$VertexTo, sep = "_")
  )
  
  edge_weights <- ifelse(
    !is.na(matched_idx),
    weights$Dist[matched_idx],
    0
  )
  
  E(g_obj)$weight <- ifelse(edge_weights == 0, 0, 1/edge_weights)
  
  g_obj  
})
names(time_ordered) <- names(time_df)
list2env(time_ordered, envir = .GlobalEnv)


#### SUPPLEMENTAL FIGURE 5 #####
layout.matrix <- matrix(c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10),
                        nrow = 2, ncol = 5)

layout.matrix

layout(
  mat = layout.matrix,
  heights = c(1, 1),
  widths = c(1.25, 1.25, 1.25, 1.25, 1.5)
)

pn_nn_fun(tal, g_33, git_37, "A)") # 32-33
pn_nn_fun(tal, g_34, git_37, "B)") # 33-34
pn_nn_fun(tal, g_35, git_37, "C)") # 34-35
pn_nn_fun(tal, g_36, git_37, "D)") # 35-36
pn_nn_fun(tal, g_37, git_37, "E)") # 36-37
pn_fun(tal, git_37, "F)", lab = FALSE)

pn_nn_fun(ash, ga_1, ga, "G)") # 35-36
pn_nn_fun(ash, ga_2, ga,  "H)") # 36-37
pn_nn_fun(ash, ga_3, ga, "I)") # 38-39
pn_fun(ash, ga, "J)", lab = FALSE)

# panel 4
or_dat_f_c <- readRDS(here::here("Objects/Talent Trap Counts.rds"))
or_dat_f_sf <- or_dat_f_c %>%
  dplyr::group_by(Site, Lat, Long) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  dplyr::select(Site) %>%
  arrange(desc(Site))
or_dat_f_sf <- left_join(or_dat_f_sf, Tobfuscate, by = c("Site" = "Tsite"))
or_dat_f_sf[is.na(or_dat_f_sf$Tlabel),][1,3] <- "T1"
or_dat_f_sf[is.na(or_dat_f_sf$Tlabel),][1,3] <- "T26"
or_dat_f_sf[is.na(or_dat_f_sf$Tlabel),][1,3] <- "T7"

### FIGURES #####
# summary descriptive node table ####
a.36 <- as_adj_list(ga_1)
a.37 <- as_adj_list(ga_2)
a.39 <- as_adj_list(ga_3)
As <- list(a.36, a.37, a.39)

t.33 <- as_adj_list(g_33)
t.34 <- as_adj_list(g_34)
t.35 <- as_adj_list(g_35)
t.36 <- as_adj_list(g_36)
t.37 <- as_adj_list(g_37)
Ts <- list(t.33, t.34, t.35, t.36, t.37)

talent_sum <- data.frame(
  Site = V(git_37)$name,
  dg_cn = as.numeric(igraph::degree(git_37, mode = "out") / (((
    length(V(git_37)$name) - 1
  ) * 5))),
  btwn = as.numeric(tbc(Ts, type = "L", directed = TRUE)),
  closeness = ifelse(is.na(as.numeric(
    tcc(Ts, type = "L", directed = TRUE)
  )), 0, as.numeric(tcc(
    Ts, type = "L", directed = TRUE
  ))),
  pgrnk = as.numeric(page_rank(git_37, directed = TRUE)$vector)
) %>%
  arrange(desc(Site)) %>%
  left_join(Tobfuscate, by = c("Site" = "Tsite")) %>%
  dplyr::rename(
    label = "Tlabel"
  )
saveRDS(talent_sum, here::here('Objects/Talent sum table part 1 (spatial centrality).RDS'))

ash_sum <- data.frame(
  Site = V(ga)$name,
  dg_cn = as.numeric(igraph::degree(ga, mode = "out") / (((
    length(V(ga)$name) - 1
  ) * 5))),
  btwn = as.numeric(tbc(As, type = "L", directed = TRUE)),
  closeness = ifelse(is.na(as.numeric(
    tcc(As, type = "L", directed = TRUE)
  )), 0, as.numeric(tcc(
    As, type = "L", directed = TRUE
  ))),
  pgrnk = as.numeric(page_rank(ga, directed = TRUE)$vector)
) %>%
  arrange(desc(Site)) %>%
  left_join(Aobfuscate, by = c("Site" = "Asite")) %>%
  dplyr::rename(
    label = "Alabel"
  )
saveRDS(ash_sum, here::here('Objects/Ashland sum table part 1 (spatial centrality).RDS'))
