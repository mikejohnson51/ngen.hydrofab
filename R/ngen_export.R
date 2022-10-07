#' Write GeoJSON
#' @param x sf object
#' @param file file path
#' @export
#' @importFrom sf write_sf st_make_valid st_transform
#'
write_geojson <- function(x, file) {
  names(x) <- tolower(names(x))
  unlink(file)
  write_sf(st_make_valid(st_transform(x, 4326)), file,
           layer_options = c("ID_FIELD=id", "ID_TYPE=String"))
}

#' @title Write NextGen files from GPKG
#' @param gpkg path to geopackage
#' @param dir directory to write data to, if NULL, then its written at the same level as the gpkg
#' @param catchment_name name of catchment layer
#' @return NULL
#' @export
#' @importFrom logger log_info log_success
#' @importFrom sf read_sf st_drop_geometry
#' @importFrom jsonlite write_json
#' @importFrom dplyr select mutate left_join group_by arrange ungroup
#' @importFrom nhdplusTools get_vaa
#' @importFrom tidyr unnest_longer
#' @importFrom hydrofab layer_exists


write_ngen_dir = function(gpkg,
                          dir = NULL,
                          catchment_name = "divides",
                          export_shapefiles = FALSE,
                          verbose = TRUE,
                          overwrite = FALSE){

  if(!layer_exists(gpkg, catchment_name)){
    stop("Need ", catchment_name, " in gpkg", call. = FALSE)
  }

  if(!layer_exists(gpkg, "nexus")){
    stop("Need ", "nexus", " in gpkg", call. = FALSE)
  }

  if(!layer_exists(gpkg, "flowpath_edge_list")){
    stop("Need ", "flowpath_edge_list", " in gpkg", call. = FALSE)
  }

  if(!layer_exists(gpkg, "flowpath_attributes")){
    stop("Need ", "flowpath_attributes", " in gpkg", call. = FALSE)
  }

  if(is.null(dir)){
    dir = file.path(dirname(gpkg),strsplit(basename(gpkg), "\\.")[[1]][1])
  }

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  hyaggregate_log("INFO", glue("Writing to: {dir}"))

  if(any(!file.exists(file.path(dir, "catchment_data.geojson")), overwrite)){
    write_geojson(
      read_sf(gpkg, catchment_name),
      file.path(dir, "catchment_data.geojson")
    )
    hyaggregate_log("SUCCESS", glue("Completed: {file.path(dir, 'catchment_data.geojson')}"), verbose)
  }

  if(any(!file.exists(file.path(dir, "nexus_data.geojson")), overwrite)){
    write_geojson(
      read_sf(gpkg, "nexus"),
      file.path(dir, "nexus_data.geojson")
    )
    hyaggregate_log("SUCCESS", glue("Completed: {file.path(dir, 'nexus_data.geojson')}"), verbose)
  }

  if(any(!file.exists(file.path(dir, "flowpath_edge_list.json")), overwrite)){

    write_json(
      read_sf(gpkg, "flowpath_edge_list"),
      file.path(dir, "flowpath_edge_list.json"),
      pretty = TRUE
    )

    hyaggregate_log("SUCCESS", glue("Completed: {file.path(dir, 'flowpath_edge_list.json')}"), verbose)

  }


  if(any(!file.exists(file.path(dir, "flowpath_attributes.json")), overwrite)){

    wb_feilds =  read_sf(gpkg, "flowpath_attributes")

    wb_feilds2 <- split(select(wb_feilds, -id), seq(nrow(wb_feilds)))

    names(wb_feilds2) = wb_feilds$id

    write_json(
      wb_feilds2,
      file.path(dir, "flowpath_attributes.json"),
      pretty = TRUE
    )

    hyaggregate_log("SUCCESS", glue("Completed: {file.path(dir, 'flowpath_params.json')}"), verbose)
  }

  if(any(!file.exists(file.path(dir, "crosswalk")), overwrite) &
     layer_exists(gpkg,"crosswalk")){
    write.csv(
      read_sf(gpkg, "crosswalk"),
      file.path(dir, "crosswalk.csv"),
      row.names = FALSE
    )

    hyaggregate_log("SUCCESS", glue("Completed: {file.path(dir, 'crosswalk.csv')}"), verbose)
  }


  if(any(!file.exists(file.path(dir, "aorc_weights")), overwrite) &
     layer_exists(gpkg,"aorc_weights")){
    write.csv(
      read_sf(gpkg, "aorc_weights"),
      file.path(dir, "aorc_weight_grids.csv"),
      row.names = FALSE
    )

    hyaggregate_log("SUCCESS", glue("Completed: {file.path(dir, 'aorc_weights.csv')}"), verbose)
  }

  if(any(!file.exists(file.path(dir, "cfe_noahowp_attributes")), overwrite) &
     layer_exists(gpkg,"cfe_noahowp_attributes")){
    write.csv(
      read_sf(gpkg, "cfe_noahowp_attributes"),
      file.path(dir, "cfe_noahowp_attributes.csv"),
      row.names = FALSE
    )

    hyaggregate_log("SUCCESS", glue("Completed: {file.path(dir, 'cfe_noahowp_attributes.csv')}"), verbose)
  }

  if(export_shapefiles){ write_shapefile_dir(gpkg, dir = dir) }

  return(dir)
}


#' @title Write NextGen geopackage as directory of shapefiles
#' @param gpkg path to geopackage
#' @param dir directory path to create a 'shp' folder for output
#' @param verbose should messages be emmited?
#' @return NULL
#' @export
#' @importFrom logger log_info
#' @importFrom sf gdal_utils

#' @importFrom tidyr unnest_longer

write_shapefile_dir = function(gpkg, dir, verbose = TRUE){

  outpath = file.path(dir, "shps")
  dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

  hyaggregate_log("INFO", glue("Writing shapefiles to: {outpath}"), verbose)
  log_info("Writing shapefiles to: {outpath}")

  gdal_utils(
    util = "vectortranslate",
    source = gpkg,
    destination = outpath, # output format must be specified for GDAL < 2.3
    options = c("-f", "ESRI Shapefile", "-overwrite")
)
}


