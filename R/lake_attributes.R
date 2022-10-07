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

  if(is.null(lake_path)){
    stop("lake_path cannot be NULL")
  }

  all = c("lake_id", lake_vars)

  nc = open.nc(lake_path)

  out = lapply(1:length(all), FUN = function(x) { var.get.nc(nc, all[x]) } )

  out = data.frame(do.call("cbind", out))

  names(out) = all

  read_sf(gpkg, "lookup_table") %>%
    filter(POI_TYPE %in% c("WBIn", 'WBOut')) %>%
    mutate(POI_VALUE = as.numeric(POI_VALUE)) %>%
    left_join(out, by = c("POI_VALUE" = "lake_id")) %>%
    select(id, toid, lake_vars) %>%
    write_sf(gpkg, "lake_attributes")

  gpkg

}
