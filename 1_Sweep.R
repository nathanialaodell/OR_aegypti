library(tidyverse)
library(MMWRweek)
library(sf)
library(tidygeocoder)
library(maps)
library(tidygraph)
library(igraph)
library(ggraph)
library(gridExtra)
library(stringr)
library(tigris)
library(purrr)
library(nngeo)
library(readxl)
library(here)
source(here::here("From_to_functions.R"))

#--------
# talent
#---------

dat <- readxl::read_excel(here::here('Objects/talent_traps_OR.xlsx'))
talent <- dat[-1,] # no site info in first row, so remove

#### IGNORE COLUMN NAME MESSAGE (see below) #####

talent_f <- talent[, -c(6, 9, 12)]

names(talent_f) <- c(
  "Site", "Lat", "Long", "Trap_Type", "Detection_1",
  "Tot_1", "Detection_2", "Tot_2",
  "Detection_3", "Tot_3"
)

talent_f <- talent_f %>%
  mutate(
    Lat = as.numeric(Lat),
    Long = as.numeric(Long),
    Tot_1 = as.numeric(Tot_1),
    Tot_2 = as.numeric(Tot_2),
    Tot_3 = as.numeric(Tot_3),
    Detection_1 = as.Date(Detection_1, format = "%m/%d/%Y"),
    Detection_2 = as.Date(Detection_2, format = "%m/%d/%Y"),
    Detection_3 = as.Date(Detection_3, format = "%m/%d/%Y"),
    Site = gsub(".*-", "", Site), # remove the 6 digit prefix from sites
  )

talent_f$Site <- str_trim(talent_f$Site)

talent_f <- talent_f %>%
  pivot_longer(
    cols = starts_with("Detection_"),
    names_to = "Detection",
    values_to = "Date"
  ) %>%
  mutate(
    Total = case_when(
      str_detect(Detection, "1") ~ Tot_1,
      str_detect(Detection, "2") ~ Tot_2,
      str_detect(Detection, "3") ~ Tot_3
    ),
    City = "Talent",
    State = "Oregon"
  ) %>%
  dplyr::select(
    Site, Lat, Long, City,
    State, Trap_Type, Date, Total
  )
talent_f$Total <- ifelse(is.na(talent_f$Total), 0, talent_f$Total)

talent_f <- talent_f %>%
  na.omit()

# grab epidemiological week and append
talent_f$epi_wk <- MMWRweek(talent_f$Date)$MMWRweek

talent_f <- talent_f %>%
  dplyr::filter(epi_wk > 30) %>%
  dplyr::group_by(Site, Lat, Long, epi_wk) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE))

talent_f_sf <- talent_f %>%
  dplyr::select(Site, Lat, Long, epi_wk, Total) %>%
  dplyr::group_by(Site, Lat, Long) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
  dplyr::rename(Address = Site) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

t_p <- talent_f %>%
  mutate(Obs = ifelse(Total > 0, 1, 0)) %>%
  arrange(epi_wk, Site)

# presence only
t_p_1 <- t_p %>%
  filter(Obs == 1) %>%
  arrange(epi_wk, Site) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326, remove = FALSE)

# ls2 <- list()
# Adr <- levels(as.factor(t_p_1$Site))
# for (i in 1:length(Adr)) {
#   ls2[[i]] <- t_p_1 %>%
#     filter(Site == Adr[i]) %>%
#     dplyr::arrange(epi_wk)
# }

### date of first detection
t_p_2 <- t_p_1 %>% dplyr::arrange(epi_wk) %>% 
  dplyr::arrange(epi_wk)

dts.nw <- c(unique(t_p_2$epi_wk))[-1]
dts.all <- unique(t_p_1$epi_wk)

t_cnct <- rbind(
  new.fun(from = t_p_1, to = t_p_2, x = 1, y = 1), # week 1
  new.fun(from = t_p_1, to = t_p_2, x = 2, y = 2), # 2
  new.fun(from = t_p_1, to = t_p_2, x = 3, y = 3), # so on....
  new.fun(from = t_p_1, to = t_p_2, x = 4, y = 4),
  new.fun(from = t_p_1, to = t_p_2, x = 5, y = 5)
) %>%
  na.omit()

write_sf(t_cnct, here::here("Objects/Talent connectivity NA omitted.shp"))

#---------
# Ashland
#---------

ash <- read.csv(here::here("Objects/Ashland Traps Raw.csv")) %>%
  filter(
    F_Total != 0 # only want female collections (removes one address!)
  ) %>%
  dplyr::rename(Total = F_Total) %>%
  mutate(Date = as.Date(Date, "%m/%d/%Y")) %>%
  dplyr::select(Address, City, State, Total, Date)

ash$epi_wk <- MMWRweek(ash$Date)$MMWRweek

ash_geo <- ash %>%
  geocode(
    street = Address,
    city = City,
    state = State,
    method = 'census',
    full_results = FALSE
  ) %>%
  dplyr::select(Address, Total, epi_wk, lat, long) %>%
  dplyr::rename(Site = Address, Lat = lat, Long = long)

ash_geo <- ash_geo %>%
  dplyr::group_by(Site, Lat, Long, epi_wk) %>%
  dplyr::summarise(Total = sum(Total)) %>%
  dplyr::arrange(epi_wk) %>%
  ungroup()

write.csv(ash_geo, here::here("Objects/Ashland Cleaned.csv"))

ash_sf <- ash_geo %>%
  dplyr::rename(Address = Site) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

ash_pres <- ash_geo %>%
  dplyr::mutate(Obs = ifelse(Total > 0, 1, 0)) %>%
  arrange(epi_wk, Site)

ash_pres_sf <- ash_pres %>%
  dplyr::select(Site, Lat, Long, epi_wk, Obs) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

ash_pres_sf$Lat <- ash_pres$Lat
ash_pres_sf$Long <- ash_pres$Long

Lat = ash_pres$Lat
Long = ash_pres$Long
ash_pres <- ash_pres %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

ash_pres$Lat = Lat
ash_pres$Long = Long

ash_pres.p2 <- ash_pres %>% dplyr::arrange(epi_wk) %>%
dplyr::arrange(epi_wk)

dts.nw <- c(unique(ash_pres.p2$epi_wk))[-1]
dts.all <- unique(ash_pres$epi_wk)

ash_sf_cnct <- rbind(
  new.fun(from = ash_pres, to = ash_pres.p2, x = 1, y = 1),
  new.fun(from = ash_pres, to = ash_pres.p2, x = 2, y = 2),
  new.fun(from = ash_pres, to = ash_pres.p2, x = 3, y = 3)
) %>%
  na.omit()

write_sf(ash_sf_cnct, here::here("Objects/Ashland connectivity NA omitted.shp"))
