library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(dplyr)
library(marmap) 
library(tigris)
library(cowplot)

us_states <- ne_states(country = "united states of america", returnclass = "sf") %>%
  st_transform(
    crs = 4326
  )
oregon <- us_states %>% filter(name == "Oregon")
oregon_counties <- counties(state = "OR", cb = TRUE, resolution = "20m") %>%
  st_transform(crs = 4326)

talent <- read_sf(here::here("Objects/Talent connectivity NA omitted.shp")) %>%
  dplyr::filter(Dist <= 700) %>%
  dplyr::mutate(
    city = "Talent"
  ) %>%
  st_drop_geometry()

ash <- read_sf(here::here("Objects/Ashland connectivity NA omitted.shp")) %>%
  dplyr::filter(Dist <= 700) %>%
  dplyr::mutate(
    city = "Ashland"
  ) %>%
  st_drop_geometry()

sites <- rbind(talent, ash)

from_coords <- sites %>% select(X = From_X, Y = From_Y, city)
to_coords   <- sites %>% select(X = To_X, Y = To_Y, city)

coords <- bind_rows(from_coords, to_coords) %>%
  dplyr::distinct(X, Y, city) %>%
  dplyr::filter(!is.na(X) & !is.na(Y)) %>% 
  sf::st_as_sf(coords = c("X", "Y"), crs = 4326) 


bbox <- st_bbox(coords)
pad <- 5
t_bbox <- st_bbox(coords %>% dplyr::filter(city == "Talent")) %>%
  st_as_sfc()
a_bbox <- st_bbox(coords %>% dplyr::filter(city == "Ashland")) %>%
  st_as_sfc()

xlims <- c(bbox["xmin"] - pad, bbox["xmax"] + pad)
ylims <- c(bbox["ymin"] - pad, bbox["ymax"] + pad)
plot_boundary <- st_as_sfc(st_bbox(c(xmin = xlims[[1]], xmax = xlims[[2]], 
                                     ymin = ylims[[1]], ymax = ylims[[2]]), crs = 4326))

main <- ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0) +
  geom_sf(data = t_bbox, fill = "#5d74a9", color = NA, alpha = 0.5, size = 2) +
  geom_sf(data = a_bbox, fill = "#a45851", color = NA, alpha = 0.5, size = 2) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_orienteering()) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(0, 0, 0, 0), # have to add the below to remove white bordering (annoying...)
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )
main  


study_area_box <- st_as_sfc(bbox, crs = 4326)

inset <- ggplot() +
  geom_sf(data = oregon, fill = "white", color = "black", linewidth = 0.4) +
  geom_sf(data = study_area_box, fill = NA, color = "red", linewidth = 0.3) +
  geom_sf(data = oregon_counties, fill = NA, color = "black", linewidth = 0.15) +
  theme_void()

final_map <- ggdraw() +
  draw_plot(main) +
  draw_plot(inset, x = 0.64, y = 0.65, width = 0.2, height = 0.2)

final_map
