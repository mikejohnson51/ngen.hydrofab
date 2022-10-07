#' Add grid mapping to geopackage
#' @param gpkg a geopackage file path
#' @param layer_name gpkg layer name with polygonal geometries
#' @param template a filepath or SpatRast object
#' @param add_to_gpkg should the weight grid be added to the gpkg?
#' @param grid_name layer name of exported layer when add_to_gpkg = TRUE
#' @param log should a log be kept
#' @return file.path or data.frame
#' @export
#' @importFrom logger log_appender appender_file appender_console
#' @importFrom zonal weight_grid
#' @importFrom sf write_sf

add_grid_mapping = function(gpkg = NULL,
                            layer_name = "aggregate_divides",
                            template = '/Users/mjohnson/Downloads/AORC-OWP_2012063021z.nc4',
                            grid_name = NULL,
                            add_to_gpkg = TRUE,
                            log = TRUE){

  if(!is.logical(log)){
    log_appender(appender_file(log))
    verbose = TRUE
  } else {
    log_appender(appender_console)
    verbose = log
  }

  hy = read_hydrofabric(gpkg)

  out = weight_grid(template,
                    geom = hy$catchments,
                    ID = "id",
                    progress = verbose)

  if(add_to_gpkg){
    if(is.null(grid_name)){ stop("To write this file to a gpkg, a `grid_name` must be provided ...") }
    write_sf(out, gpkg, grid_name)
  } else{
    return(out)
  }

  hyaggregate_log("INFO", glue('Build AORC weight grid from {template}'),  verbose)

  log_appender(appender_console)

}

