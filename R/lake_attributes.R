#' Add Lake Parameters
#' @param gpkg a hydrofabric gpkg
#' @param lake_vars a set of lake variables to extract for lake_path
#' @param lake_path a NetCDF path to an NWM LAKEPARAM file
#' @return gpkg filepath
#' @export

add_lake_attributes  = function(gpkg,
                           lake_vars = c("Dam_Length",
                                         "ifd",
                                         "LkArea", "LkMxE",
                                         "OrificeA", "OrificeC", "OrificeE",
                                         "time",
                                         "WeirC", "WeirE", "WeirL"),
                           lake_path = NULL){

  if(layer_exists(gpkg, "hydrolocations")){

  if(is.null(lake_path)){
    stop("lake_path cannot be NULL")
  }

  all = c("lake_id", lake_vars)

  nc = open.nc(lake_path)

  out = lapply(1:length(all), FUN = function(x) { var.get.nc(nc, all[x]) } )

  out = data.frame(do.call("cbind", out))

  names(out) = all

  net = read_sf(gpkg, "network") %>%
    select(id, toid) %>%
    distinct()

    read_sf(gpkg, "hydrolocations") %>%
      filter(hl_reference %in% c("WBIn", 'WBOut')) %>%
      mutate(hl_link = as.numeric(hl_link)) %>%
      left_join(out, by = c("hl_link" = "lake_id")) %>%
      left_join(net, by = "id") %>%
      select(id, toid, hl_id, hl_reference, hl_link, hl_uri, all_of(lake_vars)) %>%
      write_sf(gpkg, "lakes")

  } else {
    message("Can only write LAKES when hydrolocation layer exisits in gpkg.")
  }

  gpkg

}
