nn.fun <- function(x, y, z) {
  # create From and To DFs
  dt.from <- from %>%
    filter(epi_wk == dts.all[x] |
             epi_wk == dts.all[y]) %>%
    dplyr::group_by(Site) %>%
    summarise(Obs2 = sum(Obs))
  dt.to <- from %>% filter(epi_wk == dts.nw[z])
  
  # nearest neighbor (all possible nearest neighbors)
  nn <- st_nn(dt.from, dt.to, k = nrow(dt.to), returnDist = FALSE)
  nhbr <- st_connect(dt.from, dt.to, k = nrow(dt.to), ids = nn)
  
  ### nearest neighbor but in df format
  nn2 <- st_nn(dt.from, dt.to, k = nrow(dt.to), returnDist = TRUE) %>%
    set_names("To", "Dist") %>%
    map_df(., unlist)
  
  ## converting sf line joins into an sf with distance between points
  nhbr2 <- nhbr %>%
    st_as_sf()
  nhbr2$Dist <- nn2$Dist
  return(nhbr2)
}

new.fun <- function(from, to, x, y) {
  dt.from <- from %>%
    dplyr::filter(epi_wk <= dts.all[x]) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(Obs2 = sum(Obs))
  dt.to <- to %>%
    dplyr::filter(epi_wk == dts.nw[y])
  
  nn <- st_nn(dt.from, dt.to, k = nrow(dt.to), returnDist = FALSE)
  nhbr <- st_connect(dt.from, dt.to, k = nrow(dt.to), ids = nn)
  
  nn2_raw <- st_nn(dt.from, dt.to, k = nrow(dt.to), returnDist = TRUE)
  nn2 <- data.frame(
    From_idx = rep(seq_along(nn2_raw[[1]]), lengths(nn2_raw[[1]])),
    To       = unlist(nn2_raw[[1]]),
    Dist     = unlist(nn2_raw[[2]])
  )
  
  nhbr2 <- nhbr %>% st_as_sf()
  nhbr2$From_idx <- nn2$From_idx
  nhbr2$To       <- nn2$To
  nhbr2$Dist     <- nn2$Dist
  
  Adr2 <- levels(as.factor(dt.from$Site))
  nhbr2$From   <- Adr2[nhbr2$From_idx]
  nhbr2$From_idx <- NULL
  nhbr2$epi_wk <- rep(dts.nw[y], nrow(nhbr2))
  
  dt.from2 <- from %>%
    dplyr::filter(epi_wk == dts.all[x]) %>%
    st_drop_geometry() %>%
    dplyr::select(Site, Total, Lat, Long) %>%
    dplyr::rename(
      Prior.A = Total,
      From = Site,
      From.Y = Lat,
      From.X = Long
    )
  nhbr3 <- left_join(nhbr2, dt.from2, by = "From")
  
  To_idx <- 1:nrow(dt.to)
  dt.to$To <- To_idx
  dt.to <- dt.to %>%
    st_drop_geometry() %>%
    dplyr::select(Site, Total, To, Lat, Long) %>%
    dplyr::rename(To.Y = Lat, To.X = Long)
  nhbr4 <- left_join(nhbr3, dt.to, by = "To") %>%
    dplyr::select(-To) %>%
    dplyr::rename(To = Site)
  
  nhbr4$Brng.Rel.N <- nhbr3 %>% stplanr::line_bearing()
  
  nhbr5 <- nhbr4 %>%
    dplyr::select(From,
                  To,
                  epi_wk,
                  Dist,
                  Brng.Rel.N,
                  Prior.A,
                  Total,
                  x,
                  To.Y,
                  To.X,
                  From.Y,
                  From.X)
  return(nhbr5)
}

gene_dist <- function(data) {
  dt.from <- data %>%
    st_as_sf(
      coords = c("From_X", "From_Y"),
      remove = FALSE,
      crs = 4326
    )
  
  dt.to <- data %>%
    st_as_sf(coords = c("To_X", "To_Y"),
             remove = FALSE,
             crs = 4326)
  
  Dist <- st_distance(
    st_transform(dt.from, 4087),
    st_transform(dt.to, 4087),
    which = "Euclidean",
    by_element = TRUE
  )
  
  lines <- st_sfc(lapply(seq_len(nrow(data)), function(i) {
    st_linestring(rbind(st_coordinates(dt.from[i, ]), st_coordinates(dt.to[i, ])))
  }), crs = 4326)
  
  done <- st_sf(data, geometry = lines)
  done$Brng.Rel.N <- stplanr::line_bearing(done)
  done$Dist <- Dist
  
  return(done)
}
