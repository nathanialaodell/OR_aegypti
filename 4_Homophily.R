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
library(osmdata)
library(hrbrthemes)
library(flextable)

or_dat_f_c <- readRDS(here::here("Objects/Talent Trap Counts.rds")) %>%
  dplyr::rename(
    To = Site
  ) %>%
  dplyr::select(
    To, Date, epi_wk, Total
  )

attributes <- readRDS(here::here("Objects/Talent attributes.rds"))%>%
  dplyr::select(
    Site, Date, epi_wk, Total, road_density, mn.hum, tdmean, tmax, tmean, ppt
  ) %>%
  dplyr::mutate(
    road_density = as.numeric(road_density)
  )

# talent shapefile
or_sf_cnct <- read_sf(here::here("Objects/Talent connectivity NA omitted.shp")) %>%
  dplyr::filter(Dist <= 700)
or_sf_cnct <- or_sf_cnct[!(or_sf_cnct$From == or_sf_cnct$To),]

# ashland load
ash <- read_sf(here::here("Objects/Ashland connectivity NA omitted.shp"))%>%
  dplyr::filter(Dist <= 700)
ash <- ash[!(ash$From == ash$To), ]
ash_attributes <- readRDS(here::here("Objects/Ashland Clim attributes.rds"))%>%
  dplyr::select(
    Site, Date, epi_wk, F_Total, road_density, mn.hum, tdmean, tmax, tmean, ppt
  ) %>%
  dplyr::mutate(
    road_density = as.numeric(road_density)
  ) %>%
  dplyr::rename(
    Total = "F_Total"
  )

# since we already have this object set up with From To connections,
# it should be easy to do a pearson coefficient test
tal_m <- left_join(or_sf_cnct, attributes %>%
                      dplyr::rename(From = Site) %>%
                      dplyr::mutate(epi_wk = epi_wk + 1), by = c("From", "epi_wk"))

tal_m <- left_join(tal_m, attributes %>%
                      dplyr::rename(To = Site) %>%
                      dplyr::mutate(epi_wk = epi_wk), by = c("To", "epi_wk"))

tal_m <- tal_m %>%
  rename_with(~ gsub(".x$", "_From", .x), ends_with(".x"))

tal_m <- tal_m %>%
  rename_with(~ gsub(".y$", "_To", .x), ends_with(".y"))

ash_m <- left_join(ash, ash_attributes %>%
                     dplyr::rename(From = Site) %>%
                     dplyr::mutate(epi_wk = epi_wk + 1), by = c("From", "epi_wk"))

ash_m <- left_join(ash_m, ash_attributes %>%
                     dplyr::rename(To = Site) %>%
                     dplyr::mutate(epi_wk = epi_wk), by = c("To", "epi_wk"))

ash_m <- ash_m %>%
  rename_with(~ gsub(".x$", "_From", .x), ends_with(".x"))

ash_m <- ash_m %>%
  rename_with(~ gsub(".y$", "_To", .x), ends_with(".y"))

cor.df <- rbind(ash_m, tal_m)

names <- c(
  "Road Density",
  "Relative Humidity",
  "Mean Dew Point (C°)",
  "Max Temperature (C°)",
  "Mean Temperature (C°)",
  "Precipitation (cm)"
)

from_cols <- names(cor.df)[c(15:20)]
to_cols   <- names(cor.df)[c(23:28)]

# initialize vectors to store results
cor_values <- numeric(length(from_cols))
p_values <- numeric(length(from_cols))

clean_names <- gsub("_From$", "", from_cols)

# loop through column pairs
for (i in seq_along(from_cols)) {
  test <- cor.test(cor.df[[from_cols[i]]], cor.df[[to_cols[i]]])
  cor_values[i] <- test$estimate
  p_values[i] <- test$p.value
}
cor.df.2 <- setNames(as.data.frame(as.list(cor_values)), clean_names)
conf <- setNames(as.data.frame(as.list(p_values)), clean_names)


conf.df <- conf %>%
  pivot_longer(
    cols = 1:6, names_to = "Attribute", values_to = "Pval"
  )

cor.df.2 <- cor.df.2 %>%
  pivot_longer(cols = 1:6, names_to = "Attribute", values_to = "Correlation")

cor.df.2$Attribute <- names
cor.df.2$Pval <- conf.df$Pval

flextable(cor.df.2) %>%
                theme_alafoli()%>%
                autofit() %>%
                set_caption("Network Homophily") %>%
                colformat_double(j = "Correlation", digits = 3) %>%
                colformat_double(j = "Pval", digits = 5) %>%
                set_header_labels(Correlation = "Pearson coefficient",
                                  Pval = "P-value") %>%
                set_formatter("Pval" = function(x) { formatC(x, format = "e", digits = 2) }) %>%
                footnote(j = 1, part = "header",
                         value = as_paragraph("All environmental attributes except relative humidity were extracted at 100m buffers."),
                         ref_symbols = "*") %>%
                footnote(
                  j = 3, i = c(2:6), part = "body",
                  value = as_paragraph("Significant at α = 0.05."), 
                  ref_symbols = "**"
                )

cor.test(tal_m$ppt_From, tal_m$ppt_To, method = "pearson")
cor.test(tal_m$tdmean_From, tal_m$tdmean_To, method = "pearson")
cor.test(tal_m$tmax_From, tal_m$tmax_To, method = "pearson")
cor.test(tal_m$tmean_From, tal_m$tmean_To, method = "pearson")
