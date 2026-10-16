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
library(geosphere)
library(openmeteo)
library(here)

or_dat_f_c <- readRDS(here::here("Objects/Talent Trap Counts.rds")) %>%
  arrange(desc(Site))

# talent shapefile
or.sf <- or_dat_f_c %>%
  dplyr::group_by(Site, Lat, Long) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  arrange(desc(Site))

or.sf <- left_join(or.sf, or_dat_f_c %>%
                     group_by(Site, Lat, Long) %>%
                     dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
                     dplyr::select(Site, Lat, Long), by = "Site")


or.sf <- st_transform(or.sf, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(st_geometry(or.sf), axes=T)

state_map <- map('county', fill = TRUE, plot = FALSE) %>%
  st_as_sf()

wrd.pp.20 <- raster(here::here('Objects/gpw_v4_population_density_rev11_2020_30_sec_2020.tif'))

or_sf1 <- st_transform(or.sf, crs(wrd.pp.20))
state_map2 <- st_transform(state_map, crs(wrd.pp.20))
e <- extent(or_sf1)

or.pop.20 <- crop(wrd.pp.20, e)

# getting population density
ns.line = st_sfc(st_linestring(rbind(c(st_bbox(or.sf)$xmin, st_bbox(or.sf)$ymin), 
                                     c(st_bbox(or.sf)$xmin, st_bbox(or.sf)$ymax))), 
                 crs = crs(or.sf))
st_length(ns.line)/1000

ew.line = st_sfc(st_linestring(rbind(c(st_bbox(or.sf)$xmin, st_bbox(or.sf)$ymin), 
                                     c(st_bbox(or.sf)$xmax, st_bbox(or.sf)$ymin))), 
                 crs = crs(or.sf))
st_length(ew.line)/1000  

r <- raster(e, nrow=86, ncol=56, vals=1,crs=crs(or.sf))

dens.fun <- function(data) {
  pop.mean.20 <- raster::extract(x=data, y=or.sf)
  pop.mean.20.2a <- unlist(lapply(pop.mean.20, mean))
  
  or_sf1.x <- or.sf
  or_sf1.x$Pop.Dens <- pop.mean.20.2a
  return(or_sf1.x)
}

dens.20 <- dens.fun(or.pop.20)

##### CLIMATE VARIABLES
or.df <- or_dat_f_c %>%
  dplyr::group_by(Site, Lat, Long) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
  arrange(desc(Site))

# set temporary wd (with PRISM files)
setwd('~/Extracted 800 m PRISM')

bils <- paste(getwd(),"/",list.files(pattern = "\\.bil$"), sep = "")

climate_data <- stack(
  bils
)

# extract project coordinates from raster stack
climate_data@crs@projargs <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
climate_crs <- climate_data@crs@projargs

coordinates(or.df) <- c("Long", "Lat")
proj4string(or.df) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

climate_or <- data.frame(
  coordinates(or.df),
  or.df$Site,
  extract(climate_data, or.df)
)

### cleaning up columns
climate_or <- climate_or %>%
  gather(date, value, 4:ncol(climate_or))

# The column header includes the date and data type, but also some other metadata that we don't need
# Here, I remove the extra info from the column header
climate_or$date <- gsub('prism_', '', climate_or$date) %>% 
  gsub('us_30s_', '', .)

# Split header into type (precipitation or temperature), year, and month
climate_or <- separate(climate_or, 'date', 
                         into = c('type', 'date'), 
                         sep = '_')
climate_or$date <- ymd(as.numeric(climate_or$date))
climate_or$date <- as.Date(climate_or$date, format = "%Y/%m/%d")

# Reshape data-- make a separate column for temperature and precipitation
climate_or <- unique(climate_or)
climate_or <- climate_or %>% 
  spread(type, value) %>%
  dplyr::rename(Site = or.df.Site)

## save
write_rds(climate_or, here::here("Objects/Talent with climate.rds"))

or.sf <- readRDS(here::here("Objects/Talent with climate.rds"))


## land cover
rasterdf <- function(x, aggregate = 1) {
  resampleFactor <- aggregate        
  inputRaster <- x    
  inCols <- ncol(inputRaster)
  inRows <- nrow(inputRaster)
  # Compute numbers of columns and rows in the resampled raster
  resampledRaster <- rast(ncol=(inCols / resampleFactor), 
                          nrow=(inRows / resampleFactor),
                          crs = crs(inputRaster))
  # Match to the extent of the original raster
  ext(resampledRaster) <- ext(inputRaster)
  # Resample data on the new raster
  y <- resample(inputRaster,resampledRaster,method='near')
  # Extract cell coordinates into a data frame
  coords <- xyFromCell(y, seq_len(ncell(y)))
  # Extract layer names
  dat <- stack(values(y, dataframe = TRUE))
  # Add names - 'value' for data, 'variable' for different
  # layer names in a multilayer raster
  names(dat) <- c('value', 'variable')
  dat <- cbind(coords, dat)
  dat
}

nlcd <- raster(
  here::here("Objects/Talent_NLCD_2024.tiff")
)

rasterdf(nlcd)

or.df.1 <- spTransform(or.df, crs(nlcd))

nlcd.crp <- crop(nlcd, or.df.1)

plot(nlcd.crp)
plot(or.df.1, pch = 1, col = "red",
     add = TRUE)

landcover <- extract(nlcd.crp, or.df.1, buffer = 100) # give prop of each class in a 100m buffer

landcover_proportions <- lapply(landcover, function(x) {
  counts_x <- table(x)
  proportions_x <- prop.table(counts_x)
})
landcover_proportions

or.lnd <- data.frame(
  Site = rep(or.df.1$Site, lapply(landcover_proportions, length)),
  cover = names(unlist(landcover_proportions)),
  percent = unlist(landcover_proportions)
)

or.lnd$type <- ifelse( # ugly....
  or.lnd$cover == "21", "Developed 1",
  ifelse(
    or.lnd$cover == "22", "Developed 2",
    ifelse(
      or.lnd$cover == "23", "Developed 3",
      ifelse(
        or.lnd$cover == "24", "Developed 4",
        ifelse(
          or.lnd$cover == "81", "Argiculture 1",
          "Agriculture 2"
        )
      )
    )
  )
)

or.lnd.piv <- or.lnd %>%
  dplyr::select(Site, type, percent) %>%
  pivot_wider( names_from = type, 
                          values_from = percent)

or.lnd.piv[is.na(or.lnd.piv)] <- 0

write_rds(or.lnd.piv, "Talent with landcover.rds")

# road density
calculate_road_density <- function(data, roads, dist) {
  
  # generate centroids for each unique address
  # this allows us to avoid calculating a road density for every single 
  # observation (i.e., saves us some time)
  data_unique <- data %>%
    group_by(Site) %>% # name place column
    slice(1) %>% 
    ungroup()
  
  centroids <- st_centroid(data_unique)
  
  # create buffer around each centroid
  buffer <- st_buffer(centroids, dist = dist)
  
  # find the roads that intersect the buffer
  roads_intersect <- st_intersection(roads, buffer)
  
  # calculate the length of the intersecting roads
  roads_intersect$length <- st_length(roads_intersect)
  
  # aggregate the total road length by address
  roads_aggregate <- aggregate(roads_intersect$length, 
                               by = list(roads_intersect$Site), FUN = sum)
  colnames(roads_aggregate) <- c('Site', 'length')
  
  # calculate the area of the buffer
  buffer$area <- st_area(buffer) / 10000 # convert to km^2
  
  # merge the buffer data with the road lengths
  buffers_merge <- merge(buffer, roads_aggregate, by = 'Site', all = TRUE)
  
  # fill NAs with 0 for lengths
  buffers_merge$length[is.na(buffers_merge$length)] <- 0
  
  # convert length to kilometers
  buffers_merge$length <- buffers_merge$length / 1000
  
  # calculate road density (km of road per km^2)
  buffers_merge$road_density <- buffers_merge$length / buffers_merge$area
  
  # select columns from buffers_merge
  Site  <- buffers_merge$Site
  area <- buffers_merge$area
  length <- buffers_merge$length
  road_density <- buffers_merge$road_density
  
  # create a new data frame that has measurements
  new_df <- data.frame(Site, road_density) 
  
  tx.rd <- left_join(data, new_df, by = "Site")
  
  # returning that appended table and we are done!
  return(tx.rd)
}

streets <- getbb("Talent Oregon")%>%
  opq()%>%
  add_osm_feature(key = "highway", 
                  value = c("residential", "living_street",
                            "unclassified",
                            "service", "footway","secondary", "tertiary", 
                            "secondary_link", "tertiary_link",
                            "motorway", "primary", "motorway_link", "primary_link"
                  )) %>%
  osmdata_sf()

or.rd <- calculate_road_density(or.sf, 
                       small_streets$osm_lines%>%
                         dplyr::select(geometry), # just need the geometries 
                       100)

write_rds(or.rd, "Talent with roads.rds")



## 08/22/2025; post joining the above variables
# using humidity instead of dew point
or_sf_cnct.1 <- readRDS("Talent NA with env.rds")

or.df <- readRDS("Talent Trap Counts.rds") %>%
  arrange(desc(Site))

coordinates(or.df) <- c("Lat", "Long")

hum.ls <- list()

for (i in 1:length(unique(or.df$Site))){
  hum.ls[[i]] <- data.frame(Site = unique(or.df$Site)[i],
                            weather_history(
    location = unique(coordinates(or.df))[i,],
    start = min(or.df$Date) - 7, # so we can appropriately average week 30
    end = max(or.df$Date) + 7, # same as above, but wk 37
    hourly = "relative_humidity_2m"
  )
  )
  
  hum.ls[[i]] <- hum.ls[[i]] %>%
    dplyr::mutate(
      Date = as.Date(datetime)
    ) %>%
    dplyr::group_by(Site, Date) %>%
    dplyr::summarise(
      hourly_relative_humidity_2m = mean(hourly_relative_humidity_2m)
    ) %>%
    dplyr::rename(
      mn.humidity = hourly_relative_humidity_2m
    )
}

hum.df <- bind_rows(hum.ls) %>%
  dplyr::mutate(
    epi_wk = MMWRweek(Date)$MMWRweek
  ) %>%
  group_by(
    Site, epi_wk
  ) %>%
  dplyr::summarise(
    mn.hum = mean(mn.humidity)
  )

write_rds(hum.df, "Talent with rel humidity.rds")

# getting bearing of connections (given daily)
or.df <- readRDS("Talent Trap Counts.rds") %>%
  arrange(desc(Site))

coordinates(or.df) <- c("Lat", "Long")

bear.ls <- list()
for (i in 1:length(unique(or.df$Site))){
  bear.ls[[i]] <- data.frame(Site = unique(or.df$Site)[i],
                            weather_history(
                              location = unique(coordinates(or.df))[i,],
                              start = min(or.df$Date) - 7, # so we can appropriately average week 30
                              end = max(or.df$Date) + 7, # same as above, but wk 37
                              daily = "wind_direction_10m_dominant"
                            )
  )
  
  bear.ls[[i]] <- bear.ls[[i]] %>%
    dplyr::mutate(
      Date = as.Date(date)
    ) %>%
    dplyr::group_by(Site, Date) %>%
    dplyr::summarise(
      daily_wind_direction_10m_dominant = mean(daily_wind_direction_10m_dominant)
    ) %>%
    dplyr::rename(
      mn.bear = daily_wind_direction_10m_dominant
    )
}

bear.df <- bind_rows(bear.ls) %>%
  dplyr::mutate(
    epi_wk = MMWRweek(Date)$MMWRweek
  )

write_rds(bear.df, "Talent with wind bearing.rds")

# wind speed
speed.ls <- list()
for (i in 1:length(unique(or.df$Site))){
  speed.ls[[i]] <- data.frame(Site = unique(or.df$Site)[i],
                             weather_history(
                               location = unique(coordinates(or.df))[i,],
                               start = min(or.df$Date) - 7, # so we can appropriately average week 30
                               end = max(or.df$Date) + 7, # same as above, but wk 37
                               daily = "windspeed_10m_max"
                             )
  )
  
  speed.ls[[i]] <- speed.ls[[i]] %>%
    dplyr::mutate(
      Date = as.Date(date)
    ) %>%
    dplyr::group_by(Site, Date) %>%
    dplyr::summarise(
      daily_windspeed_10m_max = mean(daily_windspeed_10m_max)
    ) %>%
    dplyr::rename(
      mn.speed = daily_windspeed_10m_max
    )
}

speed.df <- bind_rows(speed.ls) %>%
  dplyr::mutate(
    epi_wk = MMWRweek(Date)$MMWRweek
  )

wind_data <- inner_join(speed.df, bear.df, by = c("Site", "Date", "epi_wk"))

# calculate the U (east) and V (north) components
wind_data <- wind_data %>%
  mutate(
    u = mn.speed * sin(mn.bear * pi / 180),
    v = mn.speed * cos(mn.bear * pi / 180)
  )

# compute weekly means of u and v
weekly_wind <- wind_data %>%
  group_by(Site, epi_wk) %>%
  summarise(
    mean_u = mean(u, na.rm = TRUE),
    mean_v = mean(v, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # convert back to polar coordinates
  mutate(
    mean_speed = sqrt(mean_u^2 + mean_v^2),
    mean_direction = (atan2(mean_u, mean_v) * 180 / pi) %% 360
  )
head(weekly_wind)

write_rds(weekly_wind, "Talent with wind bearing.rds")


### ASHLAND ####
ash <- read.csv("Ashland Traps Raw.csv")%>%
  mutate(
    Date = as.Date(Date, "%m/%d/%Y"),
    City = "Ashland",
    State = "Oregon"
  )

ash$epi_wk <- MMWRweek(ash$Date)$MMWRweek

ash <- ash %>%
  tidygeocoder::geocode(
    street = Address,
    city = City,
    state = State,
    method = 'census',
    full_results = FALSE
  ) %>%
  dplyr::select(
    Address,
    F_Total,
    epi_wk,
    lat,
    long
  ) %>%
  dplyr::rename(
    Site = Address,
    Lat = lat,
    Long = long,
    Total = F_Total
  )

# ashland shapefile
ash.sf <- ash %>%
  dplyr::group_by(Site, Lat, Long) %>%
  dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  arrange(desc(Site))

ash.sf <- left_join(ash.sf, ash %>%
                      group_by(Site, Lat, Long) %>%
                      dplyr::summarise(Total = sum(Total, na.rm = TRUE)) %>%
                      dplyr::select(Site, Lat, Long), by = "Site")


ash.sf <- st_transform(ash.sf, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

##### CLIMATE VARIABLES
# set temporary wd (with PRISM files)
setwd('~/Objects')

bils <- paste(getwd(),"/",list.files(pattern = "\\.bil$"), sep = "")

climate_data <- stack(
  bils
)

# extract project coordinates from raster stack
climate_data@crs@projargs <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
climate_crs <- climate_data@crs@projargs

coordinates(ash) <- c("Long", "Lat")
proj4string(ash) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

climate_or <- data.frame(
  coordinates(ash),
  ash$Site,
  extract(climate_data, ash)
)

### cleaning up columns
climate_or <- climate_or %>%
  gather(date, value, 4:ncol(climate_or))

# The column header includes the date and data type, but also some other metadata that we don't need
# Here, I remove the extra info from the column header
climate_or$date <- gsub('prism_', '', climate_or$date) %>% 
  gsub('us_30s_', '', .)

# Split header into type (precipitation or temperature), year, and month
climate_or <- separate(climate_or, 'date', 
                       into = c('type', 'date'), 
                       sep = '_')
climate_or$date <- ymd(as.numeric(climate_or$date))
climate_or$date <- as.Date(climate_or$date, format = "%Y/%m/%d")

# Reshape data-- make a separate column for temperature and precipitation
climate_or <- unique(climate_or)
climate_or <- climate_or %>% 
  spread(type, value) %>%
  dplyr::rename(Site = ash.Site)

climate_or <- climate_or %>%
  dplyr::mutate(
    epi_wk = MMWRweek(date)$MMWRweek
  )

## save
write_rds(climate_or, here::here("Objects/Ashland with climate.rds"))

climate_or <- climate_or %>%
  dplyr::group_by(Site, epi_wk) %>%
  dplyr::summarise(ppt = mean(ppt, na.rm = TRUE),
                   tdmean = mean(tdmean, na.rm = TRUE),
                   tmax = mean(tmax, na.rm = TRUE),
                   tmean = mean(tmean, na.rm = TRUE))

ash <- read.csv(here::here("Objects/Ashland Traps Raw.csv"))%>%
  mutate(
    Date = as.Date(Date, "%m/%d/%Y")
  ) %>%
  dplyr::rename(
    Site = Address
  )

ash$epi_wk <- MMWRweek(ash$Date)$MMWRweek

or.sf <- left_join(ash, climate_or%>%
                     dplyr::select(Site, tdmean, tmax, tmean, ppt, epi_wk), 
                   by = c("Site", "epi_wk")
)

write_rds(or.sf, here::here("Objects/Ashland Clim attributes.rds"))

ash_att <- readRDS(here::here("Objects/Ashland Clim attributes.rds"))


# road density
calculate_road_density <- function(data, roads, dist) {
  
  # generate centroids for each unique address
  # this allows us to avoid calculating a road density for every single 
  # observation (i.e., saves us some time)
  data_unique <- data %>%
    group_by(Site) %>% # name place column
    slice(1) %>% 
    ungroup()
  
  centroids <- st_centroid(data_unique)
  
  # create buffer around each centroid
  buffer <- st_buffer(centroids, dist = dist)
  
  # find the roads that intersect the buffer
  roads_intersect <- st_intersection(roads, buffer)
  
  # calculate the length of the intersecting roads
  roads_intersect$length <- st_length(roads_intersect)
  
  # aggregate the total road length by address
  roads_aggregate <- aggregate(roads_intersect$length, 
                               by = list(roads_intersect$Site), FUN = sum)
  colnames(roads_aggregate) <- c('Site', 'length')
  
  # calculate the area of the buffer
  buffer$area <- st_area(buffer) / 10000 # convert to km^2
  
  # merge the buffer data with the road lengths
  buffers_merge <- merge(buffer, roads_aggregate, by = 'Site', all = TRUE)
  
  # fill NAs with 0 for lengths
  buffers_merge$length[is.na(buffers_merge$length)] <- 0
  
  # convert length to kilometers
  buffers_merge$length <- buffers_merge$length / 1000
  
  # calculate road density (km of road per km^2)
  buffers_merge$road_density <- buffers_merge$length / buffers_merge$area
  
  # select columns from buffers_merge
  Site  <- buffers_merge$Site
  area <- buffers_merge$area
  length <- buffers_merge$length
  road_density <- buffers_merge$road_density
  
  # create a new data frame that has measurements
  new_df <- data.frame(Site, road_density) 
  
  tx.rd <- left_join(data, new_df, by = "Site")
  
  # returning that appended table and we are done!
  return(tx.rd)
}

streets <- getbb("Ashland Oregon")%>%
  opq()%>%
  add_osm_feature(key = "highway", 
                  value = c("residential", "living_street",
                            "unclassified",
                            "service", "footway","secondary", "tertiary", 
                            "secondary_link", "tertiary_link",
                            "motorway", "primary", "motorway_link", "primary_link"
                  )) %>%
  osmdata_sf()

or.rd <- calculate_road_density(ash.sf, 
                                streets$osm_lines%>%
                                  dplyr::select(geometry), # just need the geometries 
                                100) %>%
  dplyr::select(
    Site, road_density
  )

ash_att <- left_join(
  ash_att, or.rd, by = "Site"
)

write_rds(ash_att, "Ashland Clim attributes.rds")