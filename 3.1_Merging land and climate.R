library(sf)
library(stars)
library(tidyverse)
library(ggplot2)
library(sp)
library(tmap)
library(spdep)
library(raster)
library(fasterize)
library(snow)
library(RColorBrewer)
library(rasterVis)    
library(colorspace)
library(readr)
library(gridExtra)
library(readxl)
library(MMWRweek)
library(maps)
library(plotrix)
library(terra)
library(readxl)
library(tidygeocoder)
library(writexl)
library(phylin)
library(spatstat)
library(tidycensus)
library(prism)
library(reshape)
library(stars)
library(osmdata)

or.sf <- readRDS(here::here("Objects/Talent Trap Counts.rds")) %>%
  dplyr::select(
    Site, Date, epi_wk, Total
  )

or.ln <- readRDS(here::here("Objects/Talent with landcover.rds"))
  
or.clim <- readRDS(here::here("Objects/Talent with climate.rds")) %>%
  dplyr::mutate(
    epi_wk = MMWRweek(date)$MMWRweek
  )

or.clim <- or.clim %>%
  dplyr::group_by(Site, epi_wk) %>%
  dplyr::summarise(ppt = mean(ppt, na.rm = TRUE),
            tdmean = mean(tdmean, na.rm = TRUE),
            tmax = mean(tmax, na.rm = TRUE),
            tmean = mean(tmean, na.rm = TRUE))

or.rd <- readRDS(here::here("Objects/Talent with roads.rds"))[,-2]

or.hum <- readRDS(here::here("Objects/Talent with rel humidity.rds"))

or.sf <- left_join(or.sf, or.ln, by = "Site")

or.sf <- left_join(or.sf, or.rd, by = "Site")

or.sf <- left_join(or.sf, or.hum, by = c("Site", "epi_wk"))

or.sf <- left_join(or.sf, or.br, by = c("Site", "epi_wk"))

or.sf <- left_join(or.sf, or.clim%>%
                     dplyr::select(Site, tdmean, tmax, tmean, ppt, epi_wk), 
                   by = c("Site", "epi_wk")
)

write_rds(or.sf, here::here("Objects/Talent attributes.rds"))