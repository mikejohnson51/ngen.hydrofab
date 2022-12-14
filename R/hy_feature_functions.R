add_nonnetwork_nexus_location  = function(divides,
                                          coastal_nexus_prefix = "cnx-",
                                          internal_nexus_prefix = "inx-",
                                          waterbody_prefix = "wb-"){

  coastal_nex = suppressWarnings({ st_point_on_surface(filter(divides, type == "coastal") ) }) %>%
    flush_prefix(c('id', 'toid')) %>%
    mutate(id = paste0(coastal_nexus_prefix, id),
           type = "coastal",
           toid = paste0(waterbody_prefix, 0)) %>%
    select(id, toid, type) %>%
    rename_geometry("geometry")

  internal_nex = suppressWarnings({ st_point_on_surface(filter(divides, type == "internal") ) }) %>%
    flush_prefix(c('id', 'toid')) %>%
    mutate(id = paste0(internal_nexus_prefix, id),
           type = "internal",
           toid = paste0(waterbody_prefix, 0)) %>%
    select(id, toid, type) %>%
    rename_geometry("geometry")

  return(bind_rows(coastal_nex, internal_nex))

}


#' Apply Nexus Topology
#' This function enforces the nexus-->flowpath topology and adds nexus locations,
#' a catchment edge list, a flowpath edge list, and a lookup_table to the
#' network_list object.
#' @param network_list          list containing flowpath and catchment `sf` objects
#' @param nexus_prefix          character prefix for nexus IDs
#' @param terminal_nexus_prefix character prefix for terminal nexus IDs
#' @param coastal_nexus_prefix  character prefix for coastal nexus IDs
#' @param internal_nexus_prefix character prefix for internal nexus IDs
#' @param catchment_prefix      character prefix for catchment IDs
#' @param waterbody_prefix      character prefix for catchment IDs
#' @param export_gpkg           file path to write new data. If NULL list object is returned
#' @return list or file path
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
                                verbose = TRUE,
                                export_gpkg = NULL){

  hyaggregate_log("INFO", "\n--- Applying HY_feature topology ---\n", verbose)

  network_list = read_hydrofabric(gpkg, verbose = verbose)

  ngen_flows   = realign_topology(network_list,
                                  nexus_prefix = nexus_prefix,
                                  terminal_nexus_prefix = terminal_nexus_prefix,
                                  coastal_nexus_prefix = coastal_nexus_prefix,
                                  internal_nexus_prefix = internal_nexus_prefix,
                                  catchment_prefix = catchment_prefix,
                                  waterbody_prefix = waterbody_prefix)

  ngen_flows$flowpaths = ngen_flows$flowpaths %>%
    select(id, toid, mainstem = main_id, lengthkm, areasqkm,
           tot_drainage_areasqkm, order, hydroseq, divide_id = realized_catchment)

  ngen_flows$divides = ngen_flows$divides %>%
    select(divide_id = id, toid, divide_type = type, areasqkm) %>%
    left_join(st_drop_geometry(select(ngen_flows$flowpaths, id,divide_id)))

  ngen_flows$POIs = read_sf(gpkg, "mapped_POIs") %>%
    select(poi_id) %>%
    inner_join(st_drop_geometry(select(ngen_flows$nexus, id, poi_id)), by = "poi_id") %>%
    rename_geometry("geometry")

  ngen_flows$network = add_prefix(ngen_flows$topo ,
                                   hf_prefix = waterbody_prefix,
                                   nexus_prefix = nexus_prefix,
                                   terminal_nexus_prefix = terminal_nexus_prefix,
                                   coastal_nexus_prefix = coastal_nexus_prefix,
                                   internal_nexus_prefix = internal_nexus_prefix) %>%
    left_join(st_drop_geometry(select(ngen_flows$flowpaths, id, divide_id, lengthkm, areasqkm, tot_drainage_areasqkm, mainstem)), by = "id") %>%
    left_join(st_drop_geometry(select(ngen_flows$POIs, toid = id, poi_id)), by = "toid")

  ngen_flows$topo = NULL

  ngen_flows$lookup_table = read_sf(gpkg, "lookup_table") %>%
    select(id = aggregated_flowpath_ID,
           hf_id = NHDPlusV2_COMID,
           hf_id_part = NHDPlusV2_COMID_part,
           mainstem,
           poi_id = POI_ID,
           poi_type = POI_TYPE,
           poi_value = POI_VALUE) %>%
    mutate(hf_source = "NHDPlusV2_COMID",
      id = paste0(waterbody_prefix, id)) %>%
    left_join(select(ngen_flows$network, id, toid, divide_id), by = "id") %>%
    select(id, toid, everything()) %>%
    distinct()

  if(!is.null(export_gpkg)){
    write_hydrofabric(ngen_flows, export_gpkg, enforce_dm = TRUE)
    return(export_gpkg)
  } else {
    return(ngen_flows)
  }
}


#' Merge VPU NextGen files into single CONUS file
#' @param gpkg a set of geopackage paths
#' @param outfile a file.path to write to
#' @param verbose emit messaging?
#' @return outfile path
#' @export
#' @importFrom  sf read_sf write_sf
#' @importFrom dplyr bind_rows

national_merge = function(gpkg, outfile, verbose = TRUE){

  all = st_layers(files[1])$name

  for(i in 1:length(all)){
    hyaggregate_log("INFO", glue("Merging: {all[i]}"), verbose)
    lapply(gpkg, read_sf, all[i]) %>%
      bind_rows() %>%
      write_sf(outfile, all[i])
  }

  outfile
}


#' Realign Topology to a nexus network
#' @param network_list          list containing flowpath and catchment `sf` objects
#' @param nexus_prefix          character prefix for nexus IDs
#' @param terminal_nexus_prefix character prefix for terminal nexus IDs
#' @param coastal_nexus_prefix  character prefix for coastal nexus IDs
#' @param internal_nexus_prefix character prefix for internal nexus IDs
#' @param catchment_prefix      character prefix for catchment IDs
#' @param waterbody_prefix      character prefix for catchment IDs
#' @return list
#' @export
#' @importFrom  dplyr select mutate left_join bind_rows group_by mutate ungroup filter distinct bind_rows slice_max rename case_when
#' @importFrom  nhdplusTools rename_geometry get_node
#' @importFrom  sf st_drop_geometry st_intersects st_as_sf

realign_topology = function(network_list,
                            nexus_prefix = NULL,
                            terminal_nexus_prefix = NULL,
                            coastal_nexus_prefix = NULL,
                            internal_nexus_prefix = NULL,
                            catchment_prefix = NULL,
                            waterbody_prefix = NULL){

  # Isolate the flow network
  iso = select(network_list$flowpaths,
               id, toid, hydroseq, poi_id)

  # Cast flow network to end nodes, these are the outlets of the
  ends = iso  %>%
    rename_geometry('geometry') %>%
    mutate(geometry = get_node(., "end")$geometry) %>%
    left_join(st_drop_geometry(select(iso, id)), by = "id")

  # Get all start and end node geometries
  starts_ends = bind_rows(get_node(iso, "start"),
                          get_node(iso, "end"))

  # Find the locations where the end points intersect with starting/ending points
  emap     = st_intersects(ends, starts_ends)

  # If more then one intersection occurs its a nexus,
  #  otherwise it is a junction
  ends$type = ifelse(lengths(emap) > 1, "nex", "jun")

  # Now, intersect the typed ends with the isolated flow network
  tmap = st_intersects(ends, iso)

  # Build a data.frame that stores the following:
    # 1. ID           - the ID of the end node
    # 2. type         - the type of the end node
    # 3. touches      - the flowlines the end node touches
    # 4. hs           - the hydrosequence of the end node
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
  # nexuses are those typed as nex where the id and touches ID are not equal
  # ID/toID topology persists

  nexus    = filter(df, type == "nex") %>%
    filter(.data$id != .data$touches) %>%
    select(id, toid  = toid) %>%
    distinct() %>%
    mutate(type = "nexus")

  # Terminals and easy nexuses make the first set of nexus locations
  term_nex = bind_rows(terminals, nexus)

  ### --- Junctions --- ###
  # When flowlines have junctions involved (e.g. more then one incoming flowline)
  # we need to ensure the most upstream is selected as the nexus AND
  # that the others are topologically pushed downstream

  # First, all junctions are found
  juns    = filter(df, type == "jun") %>%
    filter(.data$id != .data$touches)

  # The "top" junctions are those that do not touch any IDs in tmp1
  # AND that have the largest hydrosequnece
  # Here we shift the toID to the flowline it touches.

  top_juns = juns %>%
    filter(!touches %in% c(term_nex$toid)) %>%
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
  topo = bind_rows(term_nex, top_juns, inner_juns)  %>%
    mutate(topo_type = "network")

  ## NEW!!! ##
    xx =  select(topo, toid, type) %>%
      group_by(toid) %>%
      mutate(type = ifelse(any(type == 'terminal'), "terminal", type)) %>%
      slice(1) %>%
      ungroup()

    topo = topo %>%
      mutate(type = NULL) %>%
      left_join(xx, by = "toid") %>%
      filter(id != toid)

  # We'll use the fl-->nex topo to modify the input flow network toIDs
  # Additionally we will add the Nextgen required prefixes.
  fl =  left_join(select(network_list$flowpaths, -.data$toid), topo, by = "id")  %>%
    st_as_sf() %>%
    rename_geometry('geometry') %>%
    mutate(id = paste0(waterbody_prefix, id),
           toid = paste0(ifelse(type == "terminal", terminal_nexus_prefix, nexus_prefix), toid),
           realized_catchment = gsub(waterbody_prefix, catchment_prefix, id)) %>%
    rename(main_id = levelpathid)

  divide =  left_join(select(network_list$catchments, -.data$toid), select(topo, id, toid, net_type = type), by = "id")  %>%
    st_as_sf() %>%
    rename_geometry('geometry') %>%
    mutate(toid = ifelse(type %in% c("coastal", "internal"), id, toid),
           net_type = ifelse(is.na(net_type), type, net_type),
           type = ifelse(net_type == "terminal", "terminal", type)) %>%
    mutate(nex_pre = case_when(
      type == "terminal" ~ terminal_nexus_prefix,
      type == "coastal"  ~ coastal_nexus_prefix,
      type == "internal" ~ internal_nexus_prefix,
      TRUE ~ nexus_prefix
    )) %>%
    mutate(id = paste0(catchment_prefix, id),
           toid = paste0(nex_pre, toid),
           topo_type = NULL,
           nex_pre = NULL) %>%
    select(id, toid, type, areasqkm)

  # The nexuses defined so far are those part of the dendritic network,
  # We also want to add those that are coastal or inland.
  # We need to make POINT locations for
  # all nexus and terminal, coastal and inland divides

  nex = filter(fl, type %in% c("nexus",  "terminal")) %>%
    group_by(toid) %>%
    slice_max(hydroseq) %>%
    ungroup() %>%
    mutate(geometry = get_node(., "end")$geometry,
           id = toid) %>%
    select(id, toid, poi_id, type) %>%
    flush_prefix(col = c("id", "toid")) %>%
    distinct() %>%
    mutate(nex_pre = case_when(
      type == "terminal" ~ terminal_nexus_prefix,
      type == "coastal"  ~ coastal_nexus_prefix,
      type == "internal" ~ internal_nexus_prefix,
      TRUE ~ nexus_prefix
    )) %>%
    mutate(id = paste0(nex_pre, id),
           toid = paste0(waterbody_prefix, toid),
           nex_pre = NULL) %>%
    mutate(toid = ifelse(type == "terminal", paste0(waterbody_prefix, 0), toid)) %>%
    bind_rows(
      add_nonnetwork_nexus_location(
        divide,
        coastal_nexus_prefix = coastal_nexus_prefix,
        internal_nexus_prefix = internal_nexus_prefix,
        waterbody_prefix = waterbody_prefix
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


  if(sum(!divide$toid %in% nex$id) != 0){
    stop('Divides flow to nexus locations that do not exisit!')
  }

  if(sum(!fl$toid %in% nex$id) != 0 ){
    stop('Flowpaths flow to nexus locations that do not exisit!')
  }

  if(sum(duplicated(divide$id)) != 0 ){
    stop('Divides have duplicated IDs!!')
  }

  if(sum(duplicated(fl$id)) != 0 ){
    stop('Divides flow to nexus locations that do not exisit!')
  }


  return(list(flowpaths = select(fl, -type, -topo_type),
              divides   = divide,
              nexus     = nex,
              topo      = topo))

}


add_prefix = function(topo, hf_prefix = "cat-",
                      nexus_prefix = "nex-",
                      terminal_nexus_prefix = "tnx-",
                      coastal_nexus_prefix  = "cnx-",
                      internal_nexus_prefix = "inx-"){

  t1 = filter(topo, topo_type == "network") %>%
    mutate(nex_pre = case_when(
      type == "terminal" ~ terminal_nexus_prefix,
      type == "coastal"  ~ coastal_nexus_prefix,
      type == "internal" ~ internal_nexus_prefix,
      TRUE ~ nexus_prefix
    )) %>%
    mutate(id = paste0(hf_prefix, id),
           toid = paste0(nex_pre, toid),
           nex_pre = NULL)

  t2 = filter(topo, topo_type == "nexus") %>%
    mutate(nex_pre = case_when(
      type == "terminal" ~ terminal_nexus_prefix,
      type == "coastal"  ~ coastal_nexus_prefix,
      type == "internal" ~ internal_nexus_prefix,
      TRUE ~ nexus_prefix
    )) %>%
    mutate(id    = paste0(nex_pre,id),
           toid = paste0(hf_prefix, toid),
           nex_pre = NULL)

  bind_rows(t1,t2)   %>%
    select(id, toid, type)

}
