add_nonnetwork_nexus_location  = function(divides,
                                          coastal_nexus_prefix = "cnx-",
                                          internal_nexus_prefix = "inx-"){

  coastal_nex = suppressWarnings({ st_point_on_surface(filter(divides, type == "coastal") ) }) %>%
    mutate(id = paste0(coastal_nexus_prefix, id),
           type = "coastal",
           toid = "NA") %>%
    select(id, toid, type) %>%
    rename_geometry("geometry")

  internal_nex = suppressWarnings({ st_point_on_surface(filter(divides, type == "internal") ) }) %>%
    mutate(id = paste0(internal_nexus_prefix, id),
           type = "internal",
           toid = "NA") %>%
    select(id, toid, type) %>%
    rename_geometry("geometry")

  return(bind_rows(coastal_nex, internal_nex))

}


#' Apply Nexus Topology
#' This function enforces the nexus-->flowpath topology and adds a catchment and flowpath
#' edgelist to the network_list object. Additionally, nexus locations are identified and
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
                                terminal_nexus_prefix = "tnx-",
                                coastal_nexus_prefix = "cnx-",
                                internal_nexus_prefix = "inx-",
                                catchment_prefix = "cat-",
                                waterbody_prefix = "wb-",
                                verbose = TRUE){

  hyaggregate_log("INFO", "\n--- Applying HY_feature topology ---\n", verbose)

  network_list = read_hydrofabric(gpkg, verbose = verbose)

  ngen_flows   = realign_topology(network_list,
                                  nexus_prefix = "nex-",
                                  terminal_nexus_prefix = "tnx-",
                                  coastal_nexus_prefix = "cnx-",
                                  internal_nexus_prefix = "inx-",
                                  catchment_prefix = "cat-",
                                  waterbody_prefix = "wb-")

  ngen_flows$catchment_edge_list = add_prefix(ngen_flows$topo , hf_prefix = "cat-", nexus_prefix = "nex-")

  ngen_flows$flowpath_edge_list  = add_prefix(ngen_flows$topo , hf_prefix = "wb-", nexus_prefix = "nex-")

  ngen_flows$topo = NULL

  ngen_flows$lookup_table = read_sf(gpkg, "lookup_table") %>%
    select(id = aggregated_flowpath_ID,
           NHDPlusV2_COMID, NHDPlusV2_COMID_part,
           reconciled_ID, mainstem,
           POI_ID, POI_TYPE, POI_VALUE) %>%
    mutate(id = paste0(waterbody_prefix, id)) %>%
    left_join(ngen_flows$flowpath_edge_list, by = "id") %>%
    select(id, toid, everything()) %>%
    distinct()

  ngen_flows

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


realign_topology = function(network_list,
                            nexus_prefix = NULL,
                            terminal_nexus_prefix = NULL,
                            coastal_nexus_prefix = NULL,
                            internal_nexus_prefix = NULL,
                            catchment_prefix = NULL,
                            waterbody_prefix = NULL){

  # Isolate the flow network
  iso = select(network_list$flowpaths, id, toid, hydroseq, poi_id)

  # Cast flow network to end nodes
  ends = iso  %>%
    rename_geometry('geometry') %>%
    mutate(geometry = get_node(., "end")$geometry) %>%
    left_join(st_drop_geometry(select(iso, id)), by = "id")

  # Get all start and end node geometries
  starts_ends = bind_rows(get_node(iso, "start"), get_node(iso, "end"))

  # Find the locations where the end points interestect with starting/ending points
  emap     = st_intersects(ends, starts_ends)

  # If more then one intersection occurs its a nexus,
  #  otherwise it is a junction
  ends$type = ifelse(lengths(emap) > 1, "nex", "jun")

  # Now, intersect the typed ends with the isolated flow network
  tmap = st_intersects(ends, iso)

  # Build a data.frame that stored the following:
    # 1. ID - the ID of the end node
    # 2. type - the type of the end node
    # 3. touches - the flowlines the end node touches
    # 4. hs - the hydrosequence of the end node
    # 5. tocuhes_toID - the toID of the flowline touched by the endnode

  # The data.frame is grouped by the ID (1) and the total entries are tabulated.
  # The data.frame is then joined to the isolated flownetwork to append the topologic toID

  df = data.frame(
    id            = rep(ends$id, times = lengths(tmap)),
    type          = rep(ends$type, times = lengths(tmap)),
    hs            = rep(ends$hydroseq, times = lengths(tmap)),
    touches       = iso$id[unlist(tmap)],
    touches_toID  = iso$toid[unlist(tmap)]) %>%
    group_by(id) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    left_join(st_drop_geometry(select(iso, id, toid)), by = "id")

  ### --- TERMINALS --- ###
  # Terminals are those where an end point only touches itself
  # ID/toID topology persists

  terminals = filter(df, n == 1 & .data$id == .data$touches) %>%
    select(id, toid  = toid) %>%
    mutate(type = "terminal")

  ### --- NEXUS 1 --- ###
  # easy nexuses are those typed as nex where the id and touches ID are not equal
  # ID/toID topology persists

  nexus    = filter(df, type == "nex") %>%
    filter(.data$id != .data$touches) %>%
    select(id, toid  = toid) %>%
    distinct() %>%
    mutate(type = "nexus")

  # Terminals and easy nexuses make the first set of nexus locations
  tmp1 = bind_rows(terminals, nexus)

  ### --- Junctions --- ###
  # When flowlines have junctions involved (e.g. more then one incoming flowline)
  # We need to ensure the most Upstream is selected as the nexus AND
  # That the others are topologically pushed downstream

  # First, all junctions are found
  juns    = filter(df, type == "jun") %>%
    filter(.data$id != .data$touches)

  # The "top" junctions are those that do no touches IDs in tmp1
  # AND that have the largest hydrosequnece
  # Here we shift the toID to the flowline it touches.

  top_juns = juns %>%
    filter(!touches %in% c(tmp1$toid)) %>%
    group_by(touches) %>%
    slice_max(hs) %>%
    ungroup() %>%
    select(id, toid  = touches) %>%
    distinct() %>%
    mutate(type = "nexus")

  # The "inner" junctions are those that are not top junctions
  # Here we shift the toID to the toID of the flowline it touches

  inner_juns = juns %>%
    filter(!id %in% top_juns$id) %>%
    select(id, toid  = touches_toID) %>%
    distinct() %>%
    mutate(type = "junction")

  # The complete nexus topo network is the combination of the first set,
  # The top junctions and the inner junctions
  # Collectively, these define the fl --> nex network topology
  topo = bind_rows(tmp1, top_juns, inner_juns)  %>%
    mutate(topo_type = "network")

  # We'll use the fl-->nex topo to modify the input flow network toIDs
  # Additionally we will add the Nextgen required prefixes.
  fl =  left_join(select(network_list$flowpaths, -.data$toid),topo, by = "id")  %>%
    st_as_sf() %>%
    rename_geometry('geometry') %>%
    mutate(id = paste0(waterbody_prefix, id),
           toid = paste0(ifelse(type == "terminal", terminal_nexus_prefix, nexus_prefix), toid),
           realized_catchment = gsub(waterbody_prefix, catchment_prefix, id)) %>%
    rename(main_id = levelpathid)


  divide =  left_join(select(network_list$catchments, -.data$toid),
                  select(topo, -type), by = "id")  %>%
    st_as_sf() %>%
    rename_geometry('geometry') %>%
    mutate(id = paste0(catchment_prefix, id),
           toid = paste0(ifelse(type == "terminal", terminal_nexus_prefix, nexus_prefix), toid)) %>%
    select(id, toid, type, areasqkm)

  # The nexuses defined so far are those part of the dendritic network,
  # We also want to add those that are coastal or inland.
  # We need to make POINT locations for
  # all nexus and terminal, coastal and inland divides

  nex = filter(fl, type %in% c("nexus",  "terminal")) %>%
    mutate(geometry = get_node(., "end")$geometry,
           id = toid) %>%
    select(id, toid, poi_id, type) %>%
    flush_prefix(col = c("id", "toid")) %>%
    distinct() %>%
    mutate(id = paste0(ifelse(type == "terminal", terminal_nexus_prefix, nexus_prefix), id ),
           toid = paste0(waterbody_prefix, toid)) %>%
    bind_rows(
      add_nonnetwork_nexus_location(
        divide,
        coastal_nexus_prefix = coastal_nexus_prefix,
        internal_nexus_prefix = internal_nexus_prefix
      )
    )

  # We then add the nex --> fl topo to the existing fl --> nex topo
  topo = suppressWarnings({
        nex %>%
          flush_prefix(col = c('id', 'toid')) %>%
          select(id, toid, type) %>%
          st_drop_geometry() %>%
          mutate(topo_type = "nexus") %>%
          bind_rows(topo) %>%
          distinct()
  })


  return(list(flowpaths = select(fl, -type, -topo_type),
              divides   = divide,
              nexus     = nex,
              topo      = topo))

}


add_prefix = function(topo = topo, hf_prefix = "cat-", nexus_prefix = "nex-"){

  t1 = filter(topo, topo_type == "network") %>%
    mutate(id = paste0(hf_prefix, id),
           toid = paste0(nexus_prefix, toid))

  t2 = filter(topo, topo_type == "nexus") %>%
    mutate(id = paste0(nexus_prefix, id),
           toid = paste0(hf_prefix, toid))

  bind_rows(t1,t2) %>%
    select(id, toid) %>%
    arrange(id)

}
