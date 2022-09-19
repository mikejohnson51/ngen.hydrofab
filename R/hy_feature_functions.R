#' Apply Nexus Topology
#' This function enforces the nexus-->flowpath topology and adds a catchment and flowpath
#' edgelist to the network_list object. Additonally, nexus locations are identified and
#' added as well.
#' @param network_list  a list containing flowpath and catchment `sf` objects
#' @param nexus_prefix  character prefix for nexus IDs
#' @param terminal_nexus_prefix character prefix for terminal nexus IDs
#' @param catchment_prefix character prefix for catchment IDs
#' @param waterbody_prefix character prefix for catchment IDs
#' @param term_cut terminal IDs begin above a defined threshold
#' @return list
#' @importFrom sf read_sf
#' @importFrom dplyr select mutate left_join everything distinct
#' @export


apply_nexus_topology = function(gpkg,
                                nexus_prefix = "nex-",
                                terminal_nexus_prefix = "tnex-",
                                catchment_prefix = "cat-",
                                waterbody_prefix = "wb-",
                                term_cut = 1e9,
                                verbose = TRUE){

  hyaggregate_log("INFO", "\n--- Applying HY_feature topology ---\n", verbose)

  network_list = read_hydrofabric(gpkg, verbose = verbose)

  # Slope is unitless (e.g. m/m)
  network_list$flowpaths = add_slope(network_list$flowpaths)

  network_list$catchment_edge_list <- get_catchment_edges_terms(flowpaths = network_list$flowpaths,
                                                      nexus_prefix = nexus_prefix,
                                                      terminal_nexus_prefix = terminal_nexus_prefix,
                                                      catchment_prefix = catchment_prefix,
                                                      term_cut = term_cut )

  network_list$flowpath_edge_list  <- get_catchment_edges_terms(network_list$flowpaths,
                                                      nexus_prefix = nexus_prefix,
                                                      terminal_nexus_prefix = terminal_nexus_prefix,
                                                      catchment_prefix = waterbody_prefix,
                                                      term_cut = term_cut)

  network_list$crosswalk = read_sf(gpkg, "lookup_table") %>%
    select(id = aggregated_ID,
           NHDPlusV2_COMID, NHDPlusV2_COMID_part,
           reconciled_ID, mainstem,
           POI_ID, POI_TYPE, POI_VALUE) %>%
    mutate(id = paste0(waterbody_prefix, id)) %>%
    left_join(network_list$flowpath_edge_list, by = "id") %>%
    select(id, toid, everything()) %>%
    distinct()

  network_list$nexus =  get_nexus(fline = network_list$flowpaths,
                                  term_cut = term_cut,
                                  nexus_prefix = nexus_prefix,
                                  terminal_nexus_prefix = terminal_nexus_prefix) %>%
    left_join(network_list$catchment_edge_list, by = "id")

  hyaggregate_log("INFO", glue("Created {nrow(network_list$nexus)} nexus locations"), verbose)

  network_list$catchments <- get_catchment_data(catchment = network_list$catchments,
                                                catchment_edge_list = network_list$catchment_edge_list,
                                                catchment_prefix = catchment_prefix)


  network_list$flowpaths  <- get_flowpath_data( fline = network_list$flowpaths,
                                                catchment_edge_list = network_list$catchment_edge_list,
                                                waterbody_prefix = waterbody_prefix,
                                                catchment_prefix = catchment_prefix)

  filter(network_list$flowpaths, id == 10404) %>%
    select(id, toid)

  network_list

}



#' Get catchment edge list
#' @description get a edge list for catchments
#' @param flowpaths  sf data.frame containing hyRefactor or hyAggregate output.
#' @param nexus_prefix  character prefix for nexus IDs
#' @param terminal_nexus_prefix character prefix for terminal nexus IDs
#' @param catchment_prefix character prefix for catchment IDs
#' @param term_cut terminal IDs begin above a defined threshold
#' @return data.frame
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate left_join bind_rows `%>%`
#' @importFrom nhdplusTools rename_geometry

get_catchment_edges_terms = function(flowpaths,
                                     nexus_prefix = "nex-",
                                     terminal_nexus_prefix = "tnex-",
                                     catchment_prefix = "cat-",
                                     term_cut = 1e9) {

  fline = select(st_drop_geometry(flowpaths), id, toid)

  filter(fline, id == 10404)

  fline = flush_prefix(fline, c("id", "toid"))

  obj1 = fline %>%
    mutate(id = paste0(catchment_prefix, .data$id),
           toid = paste0(
             ifelse(.data$toid > term_cut, terminal_nexus_prefix, nexus_prefix), .data$toid
           ))

  filter(obj1, id == "cat-10404")

  obj2 =  data.frame(id = unique(fline$toid)) %>%
    left_join(mutate(select(fline, id), toid = id), by = "id") %>%
    mutate(toid = ifelse(is.na(.data$toid), 0, .data$toid)) %>%
    mutate(id =  paste0(
      ifelse(.data$id > term_cut, terminal_nexus_prefix, nexus_prefix),
      .data$id
    ),
    toid = paste0(catchment_prefix, .data$toid))

  bind_rows(obj1, obj2)
}


assign_nex_ids = function(fline, term_cut = 1e9) {

  term_node = filter(fline, toid == 0 | is.na(toid)) %>%
    mutate(toid = term_cut + 1:n())


  no_term = filter(fline, !id %in% term_node$id)

  bind_rows(term_node, no_term) %>%
    rename_geometry("geometry")

}


#' @title get nexuses
#' @title get nexuses
#' @param fline sf data.frame NHDPlus Flowlines or hyRefactor output.
#' @param nexus_prefix character prefix for nexus IDs
#' @importFrom sf st_coordinates st_as_sf st_crs
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter ungroup select n row_number rename
#' @export
#'

get_nexus <- function(fline, term_cut = 1e9,
                      nexus_prefix = "nex-",
                      terminal_nexus_prefix = "tnx-") {

  nexus <- fline %>%
    st_cast("MULTILINESTRING") %>%
    st_coordinates() %>%
    as.data.frame()

  if("L2" %in% names(nexus)) {
    nexus <- rename(nexus, GG = .data$L2)
  } else {
    nexus <- rename(nexus, GG = .data$L1)
  }

  fline <- check_nexus(fline)

  nexus <- nexus %>%
    group_by(.data$GG) %>%
    filter(row_number() == n()) %>%
    ungroup() %>%
    select(.data$X, .data$Y) %>%
    st_as_sf(coords = c("X", "Y"), crs = st_crs(fline))

  print(nrow(nexus))

  nexus$id <- fline$to_nID
  nexus$type <- ifelse(is.na(fline$poi_id), "infered", "poi")

  nexus$id = ifelse(nexus$id >= term_cut, paste0(terminal_nexus_prefix, nexus$id), paste0(nexus_prefix, nexus$id))

  if(length(unique(nexus$id)) < nrow(nexus)) {
    nexus <- group_by(nexus, .data$id) %>%
      filter(row_number() == 1) %>%
      ungroup()
  }

  write_sf(nexus, data/test_nex.gpkg)
  return(nexus)
}

check_nexus <- function(fline) {

  fline$from_nID <- fline$id

  fline <- left_join(fline,
                     select(st_drop_geometry(fline), .data$id, to_nID = .data$from_nID),
                     by = c("toid" = "id"))

  fline$to_nID[is.na(fline$to_nID)] <- fline$toid[is.na(fline$to_nID)]

  fline

}


#' Get waterbody edge list
#' @description get a edge list for waterbodies
#' @param flowpaths  sf data.frame containing hyRefactor or hyAggregate output.
#' @param wb_prefix  character prefix for waterbody IDs
#' @param terminal_wb_prefix character prefix for terminal waterbody IDs
#' @param term_cut terminal IDs begin above a defined threshold
#' @return data.frame
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate

get_waterbody_edges_terms = function(flowpaths,
                                     wb_prefix = "wb-",
                                     terminal_wb_prefix = "twb-",
                                     term_cut = 1e9) {

  fline = select(st_drop_geometry(flowpaths), id, toid)

  fline = flush_prefix(fline, "id")
  fline = flush_prefix(fline, "toid")

  fline %>% select(.data$id, .data$toid) %>%
    mutate(id = paste0(
      ifelse(.data$id > term_cut, terminal_wb_prefix, wb_prefix),
      .data$id
    ),
    toid = paste0(
      ifelse(.data$toid > term_cut, terminal_wb_prefix, wb_prefix),
      .data$toid
    ))
}


#' Get Catchment Data
#' @description get a edge list for waterbodies
#' @param catchment  sf data.frame containing hyRefactor or hyAggregate output.
#' @param wb_prefix  character prefix for waterbody IDs
#' @param terminal_wb_prefix character prefix for terminal waterbody IDs
#' @param cutoff terminal IDs begin above a defined threshold
#' @return data.frame
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate

get_catchment_data = function(catchment,
                              catchment_edge_list,
                              catchment_prefix = "cat-") {
  catchment %>%
    mutate(id = paste0(catchment_prefix, .data$id),
           toid = NULL) %>%
    left_join(catchment_edge_list,  by = "id")
}


get_flowpath_data = function(fline,
                             catchment_edge_list,
                             waterbody_prefix = "wb-",
                             catchment_prefix = "cat-") {

  if ("main_id" %in% names(fline)) {
    fline = rename(fline, levelpathid = main_id)
  }

  if (!"slope" %in% names(fline)) {
    fline = add_slope(flowpaths = fline)
  }

  select(
      fline,
      id = .data$id,
      lengthkm = .data$lengthkm,
      slope_percent = .data$slope,
      main_id = .data$levelpathid,
      member_comid = .data$member_comid,
      tot_drainage_areasqkm = tot_drainage_areasqkm,
      order = order
    ) %>%
    mutate(id = paste0(waterbody_prefix, .data$id),
           slope_percent = 100 * .data$slope_percent) %>%
    mutate(realized_catchment = gsub(waterbody_prefix, catchment_prefix, id)) %>%
    left_join(catchment_edge_list, by = c("realized_catchment" = "id"))

}


#' Merge VPU NextGen files into single CONUS file
#'
#' @param gpkg a set of geopackage paths
#' @param outfile a file.path to write to
#' @return outfile path
#' @export
#' @importFrom  sf read_sf write_sf
#' @importFrom dplyr bind_rows

national_merge = function(gpkg, outfile){

  # catchments
  lapply(gpkg, read_sf, "divides") %>%
    bind_rows() %>%
    write_sf(outfile, "divides")

  # flowpaths
  lapply(gpkg, read_sf, "flowpaths") %>%
    bind_rows() %>%
    write_sf(outfile, "flowpaths")

   # nexus
  lapply(gpkg, read_sf, "nexus") %>%
    bind_rows() %>%
    write_sf(outfile, "nexus")

  # flowpath_edge_list
  lapply(gpkg, read_sf, "flowpath_edge_list") %>%
    bind_rows() %>%
    write_sf(outfile, "flowpath_edge_list")

  # flowpath_attributes
  lapply(gpkg, read_sf, "flowpath_attributes") %>%
    bind_rows() %>%
    write_sf(outfile, "flowpath_attributes")

  # crosswalk
  lapply(gpkg, read_sf, "crosswalk") %>%
    bind_rows() %>%
    write_sf(outfile, "crosswalk")

  # cfe_noahowp_attributes
  lapply(gpkg, read_sf, "cfe_noahowp_attributes") %>%
    bind_rows() %>%
    write_sf(outfile, "cfe_noahowp_attributes")

  outfile
}


realign_topology = function(network_list, term_cut = 1e9){

  tmp = network_list$flowpaths

  term_fl      <- filter(tmp, .data$toid == 0)
  term_fl$toid <- (nrow(tmp) + term_cut + 1:nrow(term_fl))

  tmp     = bind_rows(term_fl, filter(tmp, !.data$id %in% term_fl$id))
  nrow(tmp)

  ends = tmp  %>%
    nhdplusTools::rename_geometry('geometry') %>%
    mutate(geometry = nhdplusTools::get_node(., "end")$geometry)

  nrow(ends)


  starts_ends = bind_rows(get_node(tmp, "start"), get_node(tmp, "end"))

  emap     = st_intersects(ends, starts_ends)
  tmp$type = ifelse(lengths(emap) > 1, "nex", "jun")
  tmp$type = ifelse(tmp$toid > term_cut, "term", tmp$type)

  ends2 = left_join(st_drop_geometry(select(ends, id)), st_drop_geometry(select(tmp, id, type)), by = "id")

  tmap = st_intersects(ends, tmp)

  df = data.frame(
    id            = rep(ends2$id, times = lengths(tmap)),
    type          = rep(ends2$type, times = lengths(tmap)),
    touches       = tmp$id[unlist(tmap)],
    touches_toID  = tmp$toid[unlist(tmap)]) %>%
    filter(.data$id != .data$touches) %>%
    left_join(st_drop_geometry(select(tmp, .data$id, .data$toid)), by = "id") %>%
    mutate(real_toID = ifelse(.data$type == "jun", .data$touches_toID, .data$toid)) %>%
    select(.data$id, toid = .data$real_toID) %>%
    distinct() %>%
    left_join(select(tmp, -.data$toid), by = "id")  %>%
    bind_rows(term_fl) %>%
    st_as_sf() %>%
    select(-.data$type)

}

