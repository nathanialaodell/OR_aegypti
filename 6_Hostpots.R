library(tidyverse)
library(sf)
library(tidygraph)
library(igraph)
library(ggraph)
library(gridExtra)
library(stringr)
library(tigris)
library(purrr)
library(nngeo)
library(reshape)
library(ggpubr)
library(timeordered)
library(stplanr)
library(flextable)
library(readxl)
library(MMWRweek)
library(hrbrthemes)
library(EpiContactTrace)
library(patchwork)
library(ggrepel)
library(ggplotify)
source(here::here('Creating_networks.R'))
source(here::here('5_spread analysis plot.R'))
# spread analysis
p5 <- as.ggplot(~{
  plot(ts(b, delta = 1), plot.type = "single", col = cols,
       xlab = "Time step (week)", ylab = "Fraction of sites reached")
  legend(colnames(b[, c_indices]), lwd = 1, cex = 0.75,
         col = cols25(n = 10), bg = "white", x = "topleft")
})
p5

# cleaned up trap data
or_dat_f_c <- readRDS(here::here("Objects/Talent Trap Counts.rds")) %>%
  dplyr::select(Site, Date, epi_wk, Total)

# talent shapefile
or_sf_cnct <- read_sf(here::here("Objects/Talent connectivity NA omitted.shp")) %>%
  dplyr::filter(Dist <= 700)

or_sf_cnct <- or_sf_cnct[!(or_sf_cnct$From == or_sf_cnct$To), ]

site <- unique(c(or_sf_cnct$From, or_sf_cnct$To))

label <- paste("T", 1:length(site), sep = "")

obfuscate <- data.frame(
  site,
  label
)
obfuscate

first.aeg <- or_sf_cnct %>%
  dplyr::group_by(From)%>%
  dplyr::summarise(week = min(epi_wk))

names = V(git_37)$name[27:28]
Lats = V(git_37)$Lat[27:28]
Longs = V(git_37)$Long[27:28]

g_33 <- add_vertices(
  graph_wk_with_nn_fun(or_sf_cnct, 33),
  nv = 2,
  name = names,
  Lat = Lats,
  Long = Longs
)
g_34 <- add_vertices(
  graph_wk_with_nn_fun(or_sf_cnct, 34),
  nv = 2,
  name = names,
  Lat = Lats,
  Long = Longs
)
g_35 <- add_vertices(
  graph_wk_with_nn_fun(or_sf_cnct, 35),
  nv = 2,
  name = names,
  Lat = Lats,
  Long = Longs
)
g_36 <- add_vertices(
  graph_wk_with_nn_fun(or_sf_cnct, 36),
  nv = 2,
  name = names,
  Lat = Lats,
  Long = Longs
)
g_37 <- graph_wk_with_nn_fun(or_sf_cnct, 37)
git_37 <- graph_wk_it_fun(or_sf_cnct, 37)

abundance <- or_dat_f_c %>%
  dplyr::filter(epi_wk > 30) %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
  arrange(desc(Site))

a.33 <- as_adj_list(g_33)
a.34 <- as_adj_list(g_34)
a.35 <- as_adj_list(g_35)
a.36 <- as_adj_list(g_36)
a.37 <- as_adj_list(g_37)
As <- list(a.33, a.34, a.35, a.36, a.37)


# getting per-site metrics
net.met <- data.frame(
  Site = V(git_37)$name,
  dg_cn = as.numeric(igraph::degree(git_37, mode = "out") / (((
    length(V(git_37)$name) - 1
  ) * 5))),
  btwn = as.numeric(tbc(As, type = "L", directed = TRUE)),
  closeness = ifelse(is.na(as.numeric(
    tcc(As, type = "L", directed = TRUE)
  )), 0, as.numeric(tcc(
    As, type = "L", directed = TRUE
  ))),
  pgrnk = as.numeric(page_rank(git_37, directed = TRUE)$vector),
  rev_pgrnk = rev(as.numeric(page_rank(git_37, directed = TRUE)$vector))
) %>%
  arrange(desc(Site))

net.met <- left_join(net.met, abundance, by = "Site")
# making sure order is the same
first.aeg <- first.aeg[match(net.met$Site, first.aeg$From),]
first.aeg$week <- ifelse(is.na(first.aeg$week), 0, first.aeg$week)

net.met$week <- ifelse(first.aeg$week != 0, first.aeg$week - 31, 7)


or_dat_f_c <- readRDS("Objects/Talent Trap Counts.rds") %>%
  dplyr::rename(To = Site) %>%
  dplyr::select(To, Date, epi_wk, Total)

or_tab <- left_join(or_sf_cnct, or_dat_f_c, by = c("To", "epi_wk", "Total")) %>%
  dplyr::rename(
    source = From,
    destination = To,
    t = Date,
    n = Total
  ) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  dplyr::select(source, destination, t, n) %>%
  na.omit()

###### RISK INDEX #######

# degree centrality index is the centrality at site k/max centrality in all
# of Talent (i.e., spatial extent)
# same train of thought for other indices
net.met.in <- net.met %>%
  dplyr::mutate(
    btwn_in = (btwn/max(btwn))* (Total/max(Total)),
    closeness_in = (closeness/max(closeness))* (Total/max(Total)),
    pgrnk_in = ((pgrnk/max(pgrnk)))* (Total/max(Total)),
    deg_in = ((dg_cn/max(dg_cn))) * (Total/max(Total))
  )
saveRDS(net.met.in, here::here('Objects/Talent sum table part 1 (spatial centrality).RDS'))

net.met.r <- net.met.in %>%
  dplyr::mutate(
    btwn_in = btwn_in/max(btwn_in),
    closeness_in = closeness_in/max(closeness_in),
    pgrnk_in = pgrnk_in/max(pgrnk_in),
    deg_in = deg_in/max(deg_in)
  )


### plotting indices ####
or_dat_f_c <- readRDS(here::here("Objects/Talent Trap Counts.rds"))
or_dat_f_sf <- or_dat_f_c %>%
  dplyr::group_by(Site, Lat, Long) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  dplyr::select(Site) %>%
  dplyr::filter(Site != "Wagner Creek") %>%
  arrange(desc(Site))


or_dat_f_sf <- left_join(or_dat_f_sf, net.met.r, by = "Site")
or_dat_f_sf <- left_join(or_dat_f_sf, obfuscate, by = c("Site" = "site")) %>%
  dplyr::mutate(
    shape = ifelse(
      label %in% colnames(b[, c_indices]), "18", "20"
    )
  )

set.seed(123)
or_dat_f_sf <- or_dat_f_sf %>%
  st_jitter(amount = 0.0003)

bbox <- st_bbox(or_sf_cnct)

big_street.or <- roads(state = "OR", county = "Jackson")
big_street.or <- big_street.or %>%
  st_transform(crs = st_crs(bbox))


p1 <- ggplot(or_dat_f_sf) +
  geom_sf(
    data = big_street.or$geometry,
    inherit.aes = FALSE,
    color = "darkgrey"
  ) +
  geom_sf(
    aes(color = deg_in),
    size = 3,
    alpha = .75) +
  stat_sf_coordinates(
    geom = "text_repel",
    aes(label = label),
    box.padding = 0.6
  ) +
  scale_color_gradient(high = "#F4D03F", low = "#16A085") +
  coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
  labs(color = "Index value") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p2 <- ggplot(or_dat_f_sf) +
  geom_sf(
    data = big_street.or$geometry,
    inherit.aes = FALSE,
    color = "darkgrey"
  ) +
  geom_sf(
    aes(color = btwn_in, shape = shape),
    size = 4,
    alpha = .75) +
  scale_color_gradient(high = "#F4D03F", low = "#16A085") +
  coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
  scale_shape_manual(values = c("18" = 18, "20" = 20)) +
  labs(color = "Index value") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p3 <- ggplot(or_dat_f_sf) +
  geom_sf(
    data = big_street.or$geometry,
    inherit.aes = FALSE,
    color = "darkgrey"
  ) +
  geom_sf(
    aes(color = closeness_in, shape = shape),
    size = 4,
    alpha = .75
  ) +
  scale_color_gradient(high = "#F4D03F", low = "#16A085") +
  coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
  scale_shape_manual(values = c("18" = 18, "20" = 20)) +
  labs(color = "Index value") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p4 <- ggplot(or_dat_f_sf) +
  geom_sf(
    data = big_street.or$geometry,
    inherit.aes = FALSE,
    color = "darkgrey"
  ) +
  geom_sf(
    aes(color = pgrnk_in, shape = shape),
    size = 4,
    alpha = .75) +
  scale_color_gradient(high = "#F4D03F", low = "#16A085") +
  coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
  scale_shape_manual(values = c("18" = 18, "20" = 20)) +
  labs(color = "Index value") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

top_row <- ggpubr::ggarrange(
  p1, p2, p3, p4,
  ncol = 4, nrow = 1,
  labels = c("A)", "B)", "C)", "D)"),
  common.legend = TRUE,
  legend = "bottom"
)

git_37 <- graph_wk_it_fun(or_sf_cnct, 37)

top_row / (p5  + ~pn_fun(or_sf_cnct, git_37, title = NULL, spat_legend = TRUE, lab = FALSE, temp = TRUE))

pn_fun(or_sf_cnct, git_37, title = NULL, spat_legend = TRUE, lab = FALSE, temp = TRUE)
