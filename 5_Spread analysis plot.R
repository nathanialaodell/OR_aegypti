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

ash <- ash %>%
  left_join(Aobfuscate, by = c("From" = "Asite")) %>%
  dplyr::rename(From_label = "Alabel") %>%
  left_join(Aobfuscate, by = c("To" = "Asite")) %>%
  dplyr::rename(To_label = "Alabel")

g_33 <- graph_wk_with_nn_fun(tal, 33)
g_34 <- graph_wk_with_nn_fun(tal, 34)
g_35 <- graph_wk_with_nn_fun(tal, 35)
g_36 <- graph_wk_with_nn_fun(tal, 36)
g_37 <- graph_wk_with_nn_fun(tal, 37)
git_37 <- graph_wk_it_fun(tal, 37)

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

or_sf_cnct_time <- tal %>%
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

# spread analysis
set.seed(123)
sa <- spreadanalysis(g, seq(1, 6, by = 1), 100) # number of random events
boxplot(sa[, -1], xlab = "Time delay", ylab = "Fraction reached")
b <- transformspreadbyindividual(sa)
# numeric order for sites just for flare
ord <- order(as.numeric(gsub("\\D", "", colnames(b))))
b <- b[, ord]

# panel 4
or_dat_f_c <- readRDS(here::here("Objects/Talent Trap Counts.rds"))
or_dat_f_sf <- or_dat_f_c %>%
  dplyr::group_by(Site, Lat, Long) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  dplyr::select(Site) %>%
  dplyr::filter(Site != "Wagner Creek") %>%
  arrange(desc(Site))
or_dat_f_sf <- left_join(or_dat_f_sf, Tobfuscate, by = c("Site" = "Tsite"))
or_dat_f_sf[is.na(or_dat_f_sf$Tlabel),][1,3] <- "T1"
or_dat_f_sf[is.na(or_dat_f_sf$Tlabel),][1,3] <- "T26"
or_dat_f_sf[is.na(or_dat_f_sf$Tlabel),][1,3] <- "T7"

#### FIGURES ##########
# have to finangle patchwork to make base plots (spread analysis) play nice
# with ggplot
cols <- rep("grey", ncol(b))
c_indices <- head(order(b[6, ], decreasing = TRUE), 10)

cols[c_indices] <- cols25(n = 10)
