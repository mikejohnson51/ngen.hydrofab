#' Add CFE and NOAH-OWP attributes
#' @param gpkg a geopackage with aggregation units
#' @param catchment_name the layer name of the aggregation units
#' @param flowline_name the layer name of the flowpath units
#' @param single_layer should only the top layer of a multilayer parameter be processed?
#' @param dir directory of expected nwm_grid_file (see hyAggregate::get_nwm_grids())
#' @param precision the precision of the computations
#' @return NULL
#' @export
#' @importFrom sf read_sf st_drop_geometry st_is_empty st_layers
#' @importFrom zonal weight_grid execute_zonal
#' @importFrom dplyr select mutate filter bind_cols rename inner_join group_by summarize across everything right_join
#' @importFrom tidyr unnest_longer
#' @importFrom RNetCDF open.nc var.get.nc close.nc
#' @importFrom stats weighted.mean complete.cases setNames
#' @importFrom logger log_info log_success
#' @importFrom terra sources
#' @importFrom hydrofab read_hydrofabric hyaggregate_log

add_cfe_noahowp_attributes = function(gpkg = NULL,
                                 nwm_dir = NULL,
                                 catchment_name = NULL,
                                 flowline_name  = NULL,
                                 precision = 9,
                                 outfile = NULL,
                                 verbose = TRUE,
                                 log = TRUE) {

  if(!is.logical(log)){
    log_appender(appender_file(log))
    verbose = TRUE
  } else {
    log_appender(appender_console)
    verbose = log
  }

  if(is.null(nwm_dir)){
    stop("nwm_dir cannot be NULL")
  }

  lyr = "cfe_noahowp_attributes"

  .SD <- . <- .data <-  NULL

  if(is.null(gpkg)){
    stop('gpkg cannot be missing', call. = FALSE)
  }

  hf = read_hydrofabric(gpkg, catchment_name, flowline_name)

  hyaggregate_log("INFO", "Getting NWM grid data", verbose)

  f = list.files(nwm_dir, full.names = TRUE)
  soils = correct_nwm_spatial(grep('soilproperties_CONUS_FullRouting', f, value = TRUE))
  wrf = correct_nwm_spatial(grep('wrfinput_CONUS', f, value = T))

  hyaggregate_log("INFO", glue("Building weighting grid from: {terra::sources(soils)[1]}"),  verbose)

  nwm_w_1000m = weight_grid(soils, hf$catchments,  ID = "divide_id", progress = FALSE)

  soils_exe = list()

  ### soil_properties
  soil_mode_var = c("bexp")
  wrf_mode_var = c("IVGTYP", "ISLTYP")
  soil_gm_var   = c("dksat", "psisat")
  soil_mean_var = c("slope", "smcmax", "smcwlt", "refkdt", 'cwpvt', 'vcmx25', 'mp', 'mfsno', "quartz")

  soils_exe[[1]] = execute_zonal(
    data = soils[[grepl(paste(soil_mode_var, collapse = "|"), names(soils))]],
    w = nwm_w_1000m,
    ID = 'divide_id',
    fun = "mode"
  )

  hyaggregate_log("INFO", glue('Getting mode: {paste(soil_mode_var, collapse = ", ")}'),  verbose)

  soils_exe[[2]] = execute_zonal(
    data = wrf[[grepl(paste(wrf_mode_var, collapse = "|"), names(wrf))]],
    w = nwm_w_1000m,
    ID = 'divide_id',
    drop = 'divide_id',
    fun = "mode"
  )

  hyaggregate_log("INFO", glue('Getting mode: {paste(wrf_mode_var, collapse = ", ")}'),  verbose)

  soils_exe[[3]] = execute_zonal(
    data = soils[[grepl(paste(soil_gm_var, collapse = "|"), names(soils))]],
    w = nwm_w_1000m,
    ID = 'divide_id',
    drop = 'divide_id',
    fun = zonal:::geometric_mean
  )

  hyaggregate_log("INFO", glue('Getting geometric_mean: {paste(soil_gm_var, collapse = ", ")}'),  verbose)

  soils_exe[[4]] = execute_zonal(
    data = soils[[grepl(paste(soil_mean_var, collapse = "|"), names(soils))]],
    w = nwm_w_1000m,
    ID = 'divide_id',
    drop = 'divide_id',
    fun = "mean"
  )

  hyaggregate_log("INFO", glue('Getting mean: {paste(soil_mean_var, collapse = ", ")}'),  verbose)

  exe <-  do.call('cbind', soils_exe)
  names(exe) = gsub("_Time=1", "", names(exe))

  ####

    crosswalk = read_sf(gpkg, "network") %>%
      select(divide_id, comid  = hf_id) %>%
      distinct() %>%
      filter(complete.cases(.))

    gwnc = open.nc(grep("GWBUCKPARM_CONUS_FullRouting.nc", f, value = TRUE))

    on.exit(close.nc(gwnc))

    vars      = c("Area_sqkm", "ComID", "Coeff",  "Zmax")
    vars_mode = c("Area_sqkm", "ComID", "Expon")

    gwparams_means = suppressMessages({
      lapply(vars, function(x)
        var.get.nc(gwnc, x)) %>%
        bind_cols() %>%
        setNames(vars) %>%
        rename(comid = ComID) %>%
        mutate(comid = as.integer(comid)) %>%
        inner_join(select(crosswalk, divide_id, comid), by = 'comid') %>%
        filter(complete.cases(.)) %>%
        filter(!duplicated(.)) %>%
        group_by(divide_id) %>%
        summarize(across(everything(), ~ round(
          weighted.mean(.x, w = Area_sqkm, na.rm = TRUE),
          precision
        ))) %>%
        select(-comid, -Area_sqkm)
    })

    hyaggregate_log("INFO", glue('Getting weighted mean: {paste(vars, collapse = ", ")}'),  verbose)

    getmode = function(x){
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }

    gwparams_mode = suppressMessages({
      lapply(vars_mode, function(x)
        var.get.nc(gwnc, x)) %>%
        bind_cols() %>%
        setNames(vars_mode) %>%
        rename(comid = ComID) %>%
        inner_join(select(crosswalk, divide_id, comid), by = 'comid') %>%
        filter(complete.cases(.)) %>%
        filter(!duplicated(.)) %>%
        group_by(divide_id) %>%
        summarize(Expon = getmode(floor(Expon)))
    })

    hyaggregate_log("INFO", glue('Getting mode: {paste(vars_mode, collapse = ", ")}'),  verbose)

    traits = left_join(gwparams_means, gwparams_mode, by = "divide_id") %>%
      setNames(c('divide_id', paste0('gw_', names(.)[-1]))) %>%
      full_join(exe, by = 'divide_id')

  log_appender(appender_console)

  #write_sf(traits, export_gpkg, lyr)
  arrow::write_parquet(traits, outfile)

  outfile
}


