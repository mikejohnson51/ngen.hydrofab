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
        mutate(gages = trimws(gages),
               gages = ifelse(gages == "", NA, gages),
               NHDWaterbodyComID = ifelse(NHDWaterbodyComID == -9999, NA, NHDWaterbodyComID),
               NHDWaterbodyComID = as.numeric(NHDWaterbodyComID))

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

  select(st_drop_geometry(flowpaths), id, comid = member_comid) %>%
    mutate(comid = strsplit(comid, ",")) %>%
    unnest(cols = comid) %>%
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
    summarize(slope = round(weighted.mean(slope, w = perLength, na.rm = TRUE), 10)) %>%
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
#' @importFrom dplyr select mutate rename right_join rename right_join across everything summarize `%>%` bind_cols collect
#' @importFrom stats weighted.mean
#' @importFrom glue glue
#' @importFrom arrow open_dataset

add_flowpath_attributes   = function (gpkg,
                                      rl_vars = c("hf_id", "rl_Qi_m3s", "rl_MusK", "rl_MusX", "rl_n",
                                                  "rl_So", "rl_ChSlp", "rl_BtmWdth_m",
                                                  "rl_Kchan_mmhr",
                                                  "rl_nCC", "rl_TopWdthCC_m", "rl_TopWdth_m"),
                                      hf_version = "2.2",
                                      source = "s3://lynker-spatial/hydrofabric",
                                      add_to_gpkg = TRUE) {

  net = read_sf(gpkg, "network") %>%
    select(id, hf_id, full_length = lengthkm) %>%
    filter(complete.cases(.))

  net_map =  open_dataset(glue("{source}/v{hf_version}/reference/conus_network")) %>%
    select(hf_id = id, lengthkm) %>%
    inner_join(net, by = "hf_id") %>%
    collect() %>%
    group_by(id) %>%
    mutate(totLength = sum(lengthkm),
           perLength = lengthkm / totLength) %>%
    ungroup() %>%
    select(hf_id, id, perLength)

  if (!"rl_Length_m" %in% rl_vars) { rl_vars = c("rl_Length_m", rl_vars) }

  df = open_dataset(glue('{source}/v{hf_version}/reference/conus_routelink')) %>%
    select(any_of(rl_vars)) %>%
    filter(hf_id %in% net_map$hf_id) %>%
    collect() %>%
    right_join(net_map, by = "hf_id") %>%
    group_by(id) %>%
    summarize(across(everything(), ~ round(
      weighted.mean(x = ., w = perLength, na.rm = TRUE), 8))) %>%
    select(-hf_id, -rl_Length_m, -perLength) %>%
    left_join(select(net, id, lengthkm = full_length), by = 'id') %>%
    mutate(length_m = lengthkm * 1000, lengthkm = NULL) %>%
    distinct()

  # df2 = arrow::open_dataset('/Users/mjohnson/hydrofabric/v2.2/reference/routelink_ls') %>%
  #  select(c("hf_id", #"rl_gages",
  #           "rl_NHDWaterbodyComID")) %>%
  #   filter(hf_id %in% net_map$hf_id) %>%
  #   collect() %>%
  #   right_join(net_map, by = "hf_id")  %>%
  #   group_by(id) %>%
  #   summarize(#rl_gages = paste(gages[!is.na(gages)], collapse = ","),
  #             rl_NHDWaterbodyComID = paste(unique(NHDWaterbodyComID[NHDWaterbodyComID != -9999]), collapse = ",")) %>%
  #   left_join(df, by = "id") %>%
  #   mutate(#rl_gages = ifelse(rl_gages == "", NA, rl_gages),
  #          rl_NHDWaterbodyComID = ifelse(rl_NHDWaterbodyComID == "", NA, rl_NHDWaterbodyComID)) %>%
  #   mutate(#rl_gages = as.character(rl_gages),
  #          rl_NHDWaterbodyComID = as.character(rl_NHDWaterbodyComID))

  if(add_to_gpkg){
    write_sf(df, gpkg, "flowpath_attributes")
    return(gpkg)
  } else {
    return(df)
  }
}
