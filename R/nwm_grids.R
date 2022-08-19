#' Find Latest Version of NWM on NCEP
#' @return character
#' @export

latest_nwm_version = function(){
  ncep = 'https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/'
  ver = grep("nwm", readLines(ncep), value = TRUE)
  gsub('^.*href=\"\\s*|\\s*/.*$', '', ver)
}

#' @title Check to See if NCO is on the system
#' @description NCO is a fantastic open source tool for working with NetCDF files.
#' It can be downloaded from here: http://nco.sourceforge.net/#Source
#' @return a boolean condition
#' @importFrom sys exec_internal
#' @examples
#' check_nco()
#' @export

check_nco = function(){
  tryCatch({
    x = sys::exec_internal("ncks", "--version")
    grepl("version", rawToChar(x$stderr))
  }, error = function(e){
    FALSE
  }, warning = function(w){
    FALSE
  })
}


#' @title Get NWM hydrologcic grids
#' @description Base Nextgen requires summary data from 17 layers across the soilproperties_CONUS_FullRouting
#' GWBUCKPARM_CONUS_FullRouting, and wrfinput_CONUS NWM NetCDF files. These three files have been combined and the needed layers extracted.
#' If this dataset doesn't exist locally, it can be built using this function (download_cache = FALSE) or downloaded from ScienceBase (download_cache = FALSE).
#' The spatial metadata in these files is incomplete and needs to be corrected. This can be done by passing the path returned by this function when spatial = FALSE to
#' `correct_nwm_spatial` or, spatial can be set to TRUE. Using data_cache = FALSE requires NCO to be installed!
#' @param dir directory where file should exist
#' @param spatial should the returned object be a file path (spatial= FALSE, default) or SpatRast (spatial = TRUE)
#' @param download_cache if the file does not exist locally, should it be downloaded (TRUE) or built from raw data (FALSE)
#' @return a file path (spatial = FALSE) or SpatRast object (spatial = TRUE)
#' @importFrom sbtools item_file_download
#' @importFrom logger log_info
#' @importFrom httr GET write_disk progress
#' @export

get_nwm_grids   =  function(dir = NULL, spatial = FALSE, download_cache = TRUE){

if (is.null(dir)) {
  stop("`dir` cannot be NULL", call. = FALSE)
}

name = "ngen_gridded_data.nc"

outfile = file.path(dir, name)

if(!file.exists(outfile) & download_cache){
  item_file_download(
    sb_id = "629a4246d34ec53d276f446d",
    names = name,
    destinations = outfile
  )
} else if(!file.exists(outfile) & !download_cache){

  if(!check_nco()){
    stop("Check your NCO install, or download: http://nco.sourceforge.net/#Source.", call. = FALSE)
  }

  files   = c('soilproperties_CONUS_FullRouting.nc', 'GWBUCKPARM_CONUS_FullRouting.nc', 'wrfinput_CONUS.nc')
  urls    = paste0('https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/', latest_nwm_version(), '/parm/domain/', files)

  local_files = file.path(dir, basename(urls))

  for(i in 1:length(local_files)){
    if(!file.exists(local_files[i])){
      log_info("Downloading: \n", urls[i], "\nto \n", local_files[i])
      GET(urls[i], write_disk(local_files[i], overwrite = TRUE), progress())
    }
  }

  log_info("Downloaded ", length(local_files), " files with a size of ",
           round(sum(unlist(lapply(local_files, file.size))) / 1e9, 2), " GB"
  )

  # Push 3 into 2
  system(paste("ncks -A", local_files[3], local_files[2]))
  # Push 2 into 1
  system(paste("ncks -A", local_files[2], local_files[1]))

  log_info("Merged all data into 1 file with a size of ",
           round(file.size(local_files[1]) / 1e9, 2), " GB"
  )

  # Copy 1 to final

  file.rename(local_files[1], outfile)

  # Remove files
  unlink(local_files, recursive = TRUE)

  # Extact subsets from final
  varnames  = c(
    # soild properties
    'bexp', "dksat", "psisat", "slope", "smcmax",
    "smcwlt", "refkdt", 'cwpvt', 'vcmx25', 'mp', 'mfsno',
    # wrfinput
    "IVGTYP", "ISLTYP",
    # GW
    "Area_sqkm", "ComID", "Coeff",  "Zmax", "Expon"
  )

  system(paste("ncks -C -O -v",  paste(varnames, collapse = ","), outfile, outfile))

  log_info("Extracted ", length(varnames), " varibles. New file size ",
           round(file.size(outfile) / 1e9, 2), " GB")
}

if(spatial){
  return(correct_nwm_spatial(outfile) )
} else {
  return(outfile)
}
}



#' @title Correct Degenerate NWM files
#' @description Degenerate = having lost the physical, mental, or moral qualities considered normal and desirable; showing evidence of decline. NWM files do not store extents or projections...
#' @param path path to NWM file
#' @param subds subdatasets to extract
#' @return SpatRast object
#' @export
#' @importFrom terra ext crs rast

correct_nwm_spatial = function(path, lyrs = NULL){

  template = list(
    ext = ext(
      -2303999.62876143,
      2304000.37123857,
      -1920000.70008381,
      1919999.29991619
    ),
    crs = 'PROJCS["Sphere_Lambert_Conformal_Conic",GEOGCS["GCS_Sphere",DATUM["D_Sphere",SPHEROID["Sphere",6370000.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-97.0],PARAMETER["standard_parallel_1",30.0],PARAMETER["standard_parallel_2",60.0],PARAMETER["latitude_of_origin",40.000008],UNIT["Meter",1.0]];-35691800 -29075200 126180232.640845;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision'
  )


  data = suppressWarnings({ rast(path, lyrs = lyrs) })
  terra::ext(data) <-  template$ext
  terra::crs(data) <-  template$crs

  data
}



