#' Find Latest Version of NWM on NCEP
#' @return character
#' @export

latest_nwm_version = function(){
  ncep = 'https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/'
  ver = grep("nwm", readLines(ncep), value = TRUE)
  gsub('^.*href=\"\\s*|\\s*/.*$', '', ver)
}


#' hyAggregate data direcory
#' @param dir if not supplies will default to `get("ngen_dat_dir", envir = ngen_env)`
#' @return character
#' @export

ngen_data_dir  = function (dir = NULL){
  if (is.null(dir)) {
    return(get("ngen_dat_dir", envir = ngen_env))
  }
  else {
    assign("ngen_dat_dir", dir, envir = ngen_env)
    return(invisible(get("ngen_dat_dir", envir = ngen_env)))
  }
}

#' Get Routelink Path
#' @param dir if not supplies will default to `get("ngen_dat_dir", envir = ngen_env)`
#' @param build if TRUE, and the file does not exist, should it be built?
#' @return character
#' @export

get_routelink_path = function(dir = ngen_data_dir(), build = TRUE){

  ver = latest_nwm_version()

  local_netcdf = file.path(dir, paste0("RouteLink_", gsub("\\.", "_", ver), ".nc"))
  local_fst    = gsub(".nc", ".fst", local_netcdf)

  if(!file.exists(local_fst)){
    if(build){

      rl = paste0('https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/', ver, '/parm/domain/RouteLink_CONUS.nc')

      logger::log_info("Downloading: \n{rl} \n to \n{local_netcdf}")

      dir.create(dir, recursive = TRUE, showWarnings = FALSE)

      httr::GET(rl, httr::write_disk(local_netcdf, overwrite = TRUE), httr::progress())

      logger::log_info("Extracting NetCDF variables.")

      rl_vars = c("link", "from", "to", "alt", "order", "Qi", "MusK", "MusX", "Length", "n", "So", "ChSlp", "BtmWdth", "time", "Kchan", "nCC", "TopWdthCC", "TopWdth")

      nc = open.nc(local_netcdf)
      on.exit(close.nc(nc))

      df = data.frame(do.call(cbind, lapply(rl_vars, function(x) var.get.nc(nc, x))))
      names(df) = rl_vars

      df2 = data.frame(do.call(cbind, lapply(c("link", "NHDWaterbodyComID", "gages"), function(x) var.get.nc(nc, x))))
      names(df2) = c("link", "NHDWaterbodyComID", "gages")

      df  = df2 %>%
        mutate(link = as.numeric(link)) %>%
        right_join(df, by = "link")  %>%
        rename(comid = link) %>%
        mutate(gages = trimws(.data$gages),
               gages = ifelse(.data$gages == "", NA, .data$gages),
               NHDWaterbodyComID = ifelse(NHDWaterbodyComID == -9999, NA, NHDWaterbodyComID),
               NHDWaterbodyComID = as.numeric(.data$NHDWaterbodyComID))

      logger::log_info("Converting NetCDF to fst")
      fst::write.fst(df,  local_fst)
      unlink(local_netcdf)
    } else {
      stop("File does not exisit, set `build = TRUE`", call. = FALSE)
    }
  }
    return(local_fst)
}

get_routelink_names = function(path){

  if(is.null(path)){
    path <- get_routelink_path()
  }

  fst::metadata_fst(path)[["columnNames"]]
}


#' Routlink Attribute Subset
#' @param atts character The variable names you would like, always includes comid
#' @param path	character path where the file should be saved. Default is a persistent system data as retrieved by hyAggregate_data_dir. Also see: get_routlink_path
#' @param build if TRUE, and the file does not exist, should it be built?
#' @return character
#' @export

get_routelink = function(atts = NULL, path = get_routelink_path(), build = TRUE){

  available_names = get_routelink_names(path)

  if (is.null(atts)) {
    atts <- available_names
  } else {
    bad_atts = atts[!atts %in% available_names]
    atts = atts[atts %in% available_names]
    if (length(bad_atts) > 0) {
      message(paste(bad_atts, collapse = ", "), " not in routelink data. Ignoring...")
    }
  }

  return(fst::read_fst(path, c("comid", atts[atts != "comid"])))
}


#' Add Length Mapping from VAA
#' @param flowpaths an sf object
#' @return data.frame
#' @export
#' @importFrom dplyr select mutate filter left_join right_join arrange group_by summarize
#' @importFrom sf st_drop_geometry st_as_sf st_cast st_length
#' @importFrom tidyr unnest

build_length_map = function (flowpaths, length_table) {

  select(st_drop_geometry(flowpaths), .data$id, comid = .data$member_comid) %>%
    mutate(comid = strsplit(.data$comid, ",")) %>%
    unnest(cols = .data$comid) %>%
    mutate(comid = as.numeric(gsub("\\..*","", comid))) %>%
    left_join(length_table, by = "comid") %>%
    group_by(id) %>%
    mutate(totLength = sum(lengthkm),
           perLength = round(lengthkm / totLength, 3),
           totLength = NULL, lengthkm = NULL) %>%
    ungroup()

}

#' Add Slope to Flowpaths
#' @param flowpaths sf object (LINESTRING)
#' @return sf object
#' @export
#' @importFrom nhdplusTools get_vaa
#' @importFrom dplyr select mutate right_join group_by summarize
#' @importFrom sf st_drop_geometry st_as_sf
#' @importFrom tidyr unnest
#' @importFrom stats weighted.mean

add_slope = function(flowpaths) {

  flowpaths$slope = NULL

  suppressMessages({
    build_length_map(flowpaths, length_table = get_vaa(c("lengthkm", "slope"), updated_network = FALSE)) %>%
    # To cacluate the true unitless (m/m) slope provided in the
    # NHDplusFlowlineVAA table the units must be divided by 1000 (m/km)
    mutate(slope = slope / 1000) %>%
    group_by(id) %>%
    summarize(slope = round(weighted.mean(.data$slope, w = .data$perLength, na.rm = TRUE), 10)) %>%
    right_join(flowpaths, by = "id") %>%
    st_as_sf()
  })

}


#' Length Average Routelink Variables
#' @param flowpaths sf LINESTRING
#' @param rl_vars routelink variables
#' @param rl_path routelink path (see get_routelink_path())
#' @return data.frame
#' @export
#' @importFrom nhdplusTools get_vaa
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate rename right_join rename right_join across everything summarize `%>%` bind_cols
#' @importFrom stats weighted.mean
#' @importFrom tidyr unnest
#' @importFrom RNetCDF open.nc close.nc var.get.nc

add_flowpath_attributes   = function (gpkg,
                                     rl_vars = c("comid", "Qi", "MusK", "MusX", "n", "So", "ChSlp", "BtmWdth",
                                                 "time", "Kchan", "nCC", "TopWdthCC", "TopWdth"),
                                     rl_path = get_routelink_path()) {




  flowpaths = read_sf(gpkg, "flowpaths") %>%
    #add_slope(flowpaths) %>%
    mutate(length_m = as.numeric(st_length(.)))

  net_map =  build_length_map(flowpaths, length_table = suppressMessages({
    get_vaa(c("lengthkm"),  updated_network = TRUE)}))

  if (!"Length" %in% rl_vars) { rl_vars = c("Length", rl_vars) }

  df = get_routelink(atts = rl_vars, path = rl_path)

  df = df %>%
    right_join(net_map, by = "comid") %>%
    group_by(.data$id) %>%
    summarize(across(everything(), ~ round(
      weighted.mean(x = .,
                    w = .data$perLength, na.rm = TRUE), 8))) %>%
    select(-.data$comid, -.data$Length, -.data$perLength) %>%
    left_join(select(st_drop_geometry(flowpaths), id, length_m), by = 'id')

  df2 = get_routelink(c("comid", "gages", "NHDWaterbodyComID"), path = rl_path) %>%
    right_join(net_map, by = 'comid') %>%
    group_by(.data$id) %>%
    summarize(rl_gages = paste(.data$gages[!is.na(.data$gages)], collapse = ","),
              rl_NHDWaterbodyComID = paste(unique(.data$NHDWaterbodyComID[!is.na(.data$NHDWaterbodyComID)]), collapse = ",")) %>%
    left_join(df, by = "id") %>%
    mutate(rl_gages = ifelse(.data$rl_gages == "", NA, .data$rl_gages),
           rl_NHDWaterbodyComID = ifelse(.data$rl_NHDWaterbodyComID == "", NA, .data$rl_NHDWaterbodyComID)) %>%
    mutate(rl_gages = as.character(.data$rl_gages),
           rl_NHDWaterbodyComID = as.character(.data$rl_NHDWaterbodyComID))

  write_sf(df2, gpkg, "flowpath_attributes")

  gpkg
}
