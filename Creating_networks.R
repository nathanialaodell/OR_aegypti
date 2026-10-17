library(tidyverse)
library(dplyr)
library(MMWRweek)
library(plotrix)
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
library(reshape)
library(ggpubr)
library(osmdata)
library(stplanr)
library(igraph)
library(readxl)
library(timeordered) # for timeordered networks
library(EpiContactTrace)
library(TNC)
library(dplyr)
# devtools::install_github("jakelawlor/PNWColors")
library(PNWColors)

graph_wk_fun <- function(data, wk) {
  graph <- graph_from_data_frame(
    data %>%
      filter(epi_wk == wk) %>%
      dplyr::select(From, To) %>% # keep it simple; no attributes (distance) for now
      st_drop_geometry(),
    directed = T
  )
  return(graph)
}

get_nearest_neighbors <- function(data, wk) {
  data %>%
    filter(epi_wk == wk) %>%
    dplyr::group_by(From) %>%
    slice_min(order_by = Dist,
              n = 1,
              with_ties = FALSE) %>%
    mutate(type = "nearest") %>%
    dplyr::select(From, To, type) %>%
    ungroup() %>%
    st_drop_geometry()
}

graph_wk_with_nn_fun <- function(data, wk) {
  or_sf_cnct <- data
  fil.dat <- or_sf_cnct %>%
    filter(epi_wk == wk) %>%
    dplyr::select(From, To, Dist) %>%
    st_drop_geometry() %>%
    distinct()
  
  nn.dat <- get_nearest_neighbors(or_sf_cnct, wk)
  
  # combine filtered original and nearest edges
  all_edges <- left_join(fil.dat, nn.dat, by = c("From", "To"))
  all_edges$type <- ifelse(is.na(all_edges$type), "original", all_edges$type)
  
  # build the graph
  graph <- graph_from_data_frame(d = all_edges, directed = TRUE)
  
  # add missing vertices (From sites not present as nodes)
  rest <- setdiff(or_sf_cnct$From, V(graph)$name)
  graph <- add_vertices(graph, length(rest), attr = list(name = rest))
  
  # assign coordinates to vertices
  V(graph)$Lat <- ifelse(is.na(or_sf_cnct$From_Y[match(V(graph)$name, or_sf_cnct$From)]),
                         or_sf_cnct$To_Y[match(V(graph)$name, or_sf_cnct$To)],
                         or_sf_cnct$From_Y[match(V(graph)$name, or_sf_cnct$From)])
  
  V(graph)$Long <- ifelse(is.na(or_sf_cnct$From_X[match(V(graph)$name, or_sf_cnct$From)]),
                          or_sf_cnct$To_X[match(V(graph)$name, or_sf_cnct$To)],
                          or_sf_cnct$From_X[match(V(graph)$name, or_sf_cnct$From)])
  
  # assign edge attributes
  E(graph)$weight <- 1 / all_edges$Dist
  E(graph)$type <- all_edges$type
  
  return(graph)
}

pn_nn_fun <- function(data, graph, full_graph, title) {
  or_sf_cnct <- data
  
  to_add <- setdiff(V(full_graph)$name, V(graph)$name)
  
  atts <- vertex_attr(full_graph, index = to_add)
  atts$name <- NULL
  
  graph <- add_vertices(graph,
                        nv = length(to_add),
                        name = to_add,
                        attr = atts)
  
  lo <- matrix(c(V(graph)$Long, V(graph)$Lat), ncol = 2) %>%
    na.omit()
  
  bbox <- st_bbox(or_sf_cnct)
  
  edge_colors <- ifelse(E(graph)$type == "nearest", "red", "grey")
  edge_lty <- ifelse(E(graph)$type == "nearest", "dashed", "solid")
  
  # 08/17/26 changing shape of node based on whether it is receiving or sending
  # connection
  
  # TRIANGLE CREATION FROM https://r.igraph.org/reference/shapes.html
  mytriangle <- function(coords, v = NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    
    symbols(
      x = coords[, 1],
      y = coords[, 2],
      bg = vertex.color,
      stars = cbind(vertex.size, vertex.size, vertex.size),
      add = TRUE,
      inches = FALSE
    )
  }
  # clips as a circle
  add_shape("triangle",
            clip = shapes("circle")$clip,
            plot = mytriangle)
  
  in_deg  <- degree(graph, mode = "in")
  out_deg <- degree(graph, mode = "out")
  V(graph)$shape <- ifelse(in_deg > out_deg,
                           "square",
                           ifelse(in_deg == 0 &
                                    out_deg == 0, "triangle", "circle"))
  V(graph)$color <- ifelse(in_deg > out_deg,
                           "firebrick",
                           ifelse(in_deg == 0 &
                                    out_deg == 0, "grey", "dodgerblue"))
  
  plot(
    graph,
    layout = lo,
    rescale = FALSE,
    asp = 0,
    vertex.frame.color = "black",
    vertex.label.color = "black",
    vertex.label.cex = 0.65,
    vertex.size = 0.025,
    vertex.label = NA,
    edge.arrow.size = 0.35,
    edge.color = edge_colors,
    edge.lty = edge_lty,
    xlim = c(bbox[1], bbox[3]),
    ylim = c(bbox[2], bbox[4])
  )
  title(title, adj = 0)
}

pn_fun <- function(data,
                   graph,
                   title,
                   kin_legend = FALSE,
                   spat_legend = FALSE,
                   lab = TRUE,
                   temp = FALSE) {
  if (!lab) {
    label <- NA
  } else{
    label <- V(graph)$name
  }
  
  mystar <- function(coords, v = NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    norays <- params("vertex", "norays")
    if (length(norays) != 1 && !is.null(v)) {
      norays <- norays[v]
    }
    
    mapply(
      coords[, 1],
      coords[, 2],
      vertex.color,
      vertex.size,
      norays,
      FUN = function(x, y, bg, size, nor) {
        symbols(
          x = x,
          y = y,
          bg = bg,
          stars = matrix(c(size, size / 2), nrow = 1, ncol = nor * 2),
          add = TRUE,
          inches = FALSE
        )
      }
    )
  }
  # no clipping, edges will be below the vertices anyway
  add_shape(
    "star",
    clip = shape_noclip,
    plot = mystar,
    parameters = list(vertex.norays = 5)
  )
  
  if (temp) {
    min_edge_val <- sapply(V(graph), function(v) {
      inc_edges <- incident(graph, v, mode = "out")
      ifelse(!is.infinite(min(E(graph)[inc_edges]$epi_wk)), min(E(graph)[inc_edges]$epi_wk) - 1, 37)
    })
    
    weeks <- c(32, unique(E(graph)$epi_wk))
    
    pnw_love <- PNWColors::pnw_palette(name = "Sunset2",
                                       n = length(weeks),
                                       type = "continuous")
    pal <- cut(min_edge_val, breaks = length(weeks), labels = FALSE)
    V(graph)$color <- pnw_love[pal]
    
  }
  
  or_sf_cnct <- data
  lo <- matrix(c(V(graph)$Long, V(graph)$Lat), ncol = 2) %>%
    na.omit()
  bbox <- st_bbox(or_sf_cnct)
  plot.igraph(
    graph,
    layout = lo,
    rescale = FALSE,
    asp = 0,
    vertex.shape = "star",
    vertex.frame.color = "black",
    vertex.label.color = "black",
    vertex.label.cex = 0.65,
    vertex.size = 0.03,
    vertex.label = label,
    edge.arrow.size = 0.35,
    vertex.label.dist = .003,
    xlim = c(bbox[1], bbox[3]),
    ylim = c(bbox[2], bbox[4])
  )
  title(title, adj = 0)
  
  if (kin_legend) {
    urel <- c("PO", "FS", "HS", "GP", "FA", "2nd")
    pnw_love <- PNWColors::pnw_palette(name = "Bay",
                                       n = length(urel),
                                       type = "continuous")
    pal <- setNames(pnw_love, urel)
    
    legend(
      "topright",
      legend = urel,
      col = pal,
      lty = 1,
      lwd = 4,
      title = "Relationship",
      bty = "n",
      cex = 0.65
    )
  }
  
  if (spat_legend) {
    pal <- setNames(pnw_love, weeks)
    
    legend(
      "topright",
      legend = weeks,
      col = pal,
      lty = 1,
      lwd = 4,
      title = "First detection (epidemiological week)",
      bty = "n",
      cex = 0.65
    )
  }
  
}

pn_fun_kin2 <- function(data,
                        graph,
                        title,
                        legend = TRUE,
                        lab = TRUE) {
  if (!lab) {
    label <- NA
  } else{
    label <- V(graph)$name
  }
  or_sf_cnct <- data
  lo <- matrix(c(V(graph)$Long, V(graph)$Lat), ncol = 2) %>%
    na.omit()
  bbox <- st_bbox(or_sf_cnct)
  plot.igraph(
    graph,
    layout = lo,
    rescale = FALSE,
    asp = 0,
    vertex.color = "orange",
    vertex.shape = "circle",
    vertex.frame.color = "black",
    vertex.label.color = "black",
    vertex.label.cex = 0.65,
    vertex.size = 0.01,
    vertex.label = label,
    edge.arrow.size = 0.35,
    edge.width = 3,
    vertex.label.dist = .003,
    xlim = c(bbox[1], bbox[3]),
    ylim = c(bbox[2], bbox[4]),
    main = title
  )
  
  if (legend) {
    urel <- c("PO", "FS", "GG")
    pnw_love <- PNWColors::pnw_palette(name = "Sunset2",
                                       n = length(urel),
                                       type = "continuous")
    pal <- setNames(pnw_love, urel)
    
    legend(
      "topright",
      legend = urel,
      col = pal,
      lty = 1,
      lwd = 4,
      title = "Relationship",
      bty = "n",
      cex = 0.65
    )
  }
}

pn_fun_kin3 <- function(data,
                        graph,
                        title,
                        legend = TRUE,
                        lab = TRUE,
                        buffer_m = sigma,
                        type = c("PO", "GG", "FS")) {
  if (!lab) {
    label <- NA
  } else {
    label <- V(graph)$name
  }
  
  urel <- c("PO", "FS", "GG")
  pnw_love <- PNWColors::pnw_palette(name = "Sunset2",
                                     n = length(urel),
                                     type = "continuous")
  pal <- setNames(pnw_love, urel)
  
  buffer_border <- pal[type]
  buffer_col <- scales::alpha(pal[type], 0.1)
  
  
  or_sf_cnct <- data
  
  coords_df <- data.frame(Long = V(graph)$Long, Lat = V(graph)$Lat) %>%
    na.omit()
  
  lo <- as.matrix(coords_df)
  
  pts_sf <- st_as_sf(coords_df, coords = c("Long", "Lat"), crs = 4326) %>%
    distinct()
  buff_sf <- st_buffer(pts_sf, dist = buffer_m)
  
  bbox <- st_bbox(buff_sf)
  
  
  plot(
    st_geometry(buff_sf),
    col = buffer_col,
    border = buffer_border,
    xlim = c(bbox[1], bbox[3]),
    ylim = c(bbox[2], bbox[4]),
    main = title,
    asp = 0
  )
  
  
  plot.igraph(
    graph,
    layout = lo,
    rescale = FALSE,
    add = TRUE,
    asp = 0,
    vertex.color = "orange",
    vertex.shape = "circle",
    vertex.frame.color = "black",
    vertex.label.color = "black",
    vertex.label.cex = 0.65,
    vertex.size = 0.01,
    vertex.label = label,
    edge.arrow.size = 0.35,
    edge.color = "black",
    vertex.label.dist = .003,
    xlim = c(bbox[1], bbox[3]),
    ylim = c(bbox[2], bbox[4])
  )
  
  if (legend) {
    urel <- c("PO", "FS", "GG")
    pnw_love <- PNWColors::pnw_palette(name = "Spring",
                                       n = length(urel),
                                       type = "continuous")
    pal <- setNames(pnw_love, urel)
    
    legend(
      "topright",
      legend = urel,
      col = pal,
      lty = 1,
      lwd = 4,
      title = "Relationship",
      bty = "n",
      cex = 0.65
    )
  }
}

# iterative graphing function
graph_wk_it_fun <- function(data, wk) {
  or_sf_cnct <- data
  fil.dat <- or_sf_cnct %>%
    filter(epi_wk <= wk) %>% # iteratively add weeks; for the purposes of paths
    dplyr::select(From, To, Dist, epi_wk) %>% # distance as attribute
    st_drop_geometry()
  fil.dat <- fil.dat[!duplicated(fil.dat), ]
  
  graph <- graph_from_data_frame(fil.dat, directed = T)
  
  # ensure all the sites are present (for plotting)
  rest <- setdiff(or_sf_cnct$From, V(graph)$name)
  
  graph <- add_vertices(graph, length(rest), attr = list(name = rest))
  
  # lat/long attributes for each vertex in a graph
  V(graph)$Lat <- ifelse(is.na(or_sf_cnct$From_Y[match(V(graph)$name, or_sf_cnct$From)]),
                         or_sf_cnct$To_Y[match(V(graph)$name, or_sf_cnct$To)],
                         or_sf_cnct$From_Y[match(V(graph)$name, or_sf_cnct$From)])
  V(graph)$Long <- ifelse(is.na(or_sf_cnct$From_X[match(V(graph)$name, or_sf_cnct$From)]),
                          or_sf_cnct$To_X[match(V(graph)$name, or_sf_cnct$To)],
                          or_sf_cnct$From_X[match(V(graph)$name, or_sf_cnct$From)])
  
  # assign Dist as edge weight
  E(graph)$weight <- 1 / fil.dat$Dist
  
  return(graph)
}

gene_graph_wk_it_fun <- function(data,
                                 wk,
                                 direct = TRUE,
                                 alternative = FALSE) {
  or_sf_cnct <- data
  fil.dat <- or_sf_cnct %>%
    filter(epi_wk <= wk) %>%
    dplyr::select(From, To, Dist, epi_wk, date_diff, TopRel, Kin2) %>%
    st_drop_geometry()
  fil.dat <- fil.dat[!duplicated(fil.dat), ]
  
  graph <- graph_from_data_frame(fil.dat, directed = direct)
  
  rest <- setdiff(or_sf_cnct$From, V(graph)$name)
  graph <- add_vertices(graph, length(rest), attr = list(name = rest))
  
  V(graph)$Lat <- ifelse(is.na(or_sf_cnct$From_Y[match(V(graph)$name, or_sf_cnct$From)]),
                         or_sf_cnct$To_Y[match(V(graph)$name, or_sf_cnct$To)],
                         or_sf_cnct$From_Y[match(V(graph)$name, or_sf_cnct$From)])
  V(graph)$Long <- ifelse(is.na(or_sf_cnct$From_X[match(V(graph)$name, or_sf_cnct$From)]),
                          or_sf_cnct$To_X[match(V(graph)$name, or_sf_cnct$To)],
                          or_sf_cnct$From_X[match(V(graph)$name, or_sf_cnct$From)])
  
  E(graph)$weight   <- 1 / fil.dat$Kin2
  E(graph)$kin      <- fil.dat$TopRel
  E(graph)$date_diff <- as.numeric(fil.dat$date_diff)
  E(graph)$TopRel <- fil.dat$TopRel
  E(graph)$Dist <- as.numeric(fil.dat$Dist)
  
  urel <- c("PO", "FS", "GG")
  pnw_love <- PNWColors::pnw_palette(name = "Sunset2",
                                     n = length(urel),
                                     type = "continuous")
  pal <- setNames(pnw_love, urel)
  
  E(graph)$color <- pal[E(graph)$kin]
  width <- ifelse(E(graph)$date_diff > 0, 1 / E(graph)$date_diff, 1.25)
  E(graph)$width <- scales::rescale(width, to = c(0.45, 2))
  
  if (alternative) {
    el <- as_edgelist(graph)
    
    diff_idx <- which(E(graph)$date_diff >= 14)
    
    rev_attrs <- lapply(edge_attr(graph), function(x)
      x[diff_idx])
    
    bidirect <- unname(t(el[diff_idx, c(2, 1), drop = FALSE]))
    
    graph <- add_edges(graph, bidirect, attr = rev_attrs)
    
    graph <- igraph::simplify(
      graph,
      edge.attr.comb = list(
        weight = "first",
        kin = "first",
        TopRel = "first",
        Dist = "first",
        date_diff = "first",
        color = "first",
        width = "first",
        "ignore"
      )
    )
  }
  
  return(graph)
}

edge_list_fun <- function(graph, kern_type) {
  kern <- switch(
    kern_type,
    "FS" = kerns$FS,
    "HS" = kerns$HS,
    "GG" = kerns$GG,
  )
  
  E(graph)[Dist > kern]
}
