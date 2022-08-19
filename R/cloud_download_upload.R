#' Authenticated Upload to s3
#' You must authenticate s3 for this working session to use this function. Also,
#' the aws.s3 package is required on your machine and is NOT a hyAggregate dependency.
#' @param path a path to a file or directory
#' @param bucket the s3 bucket to upload to
#' @param prefix an optional file object prefix (think of sub directory)
#' @param verbose should messages be emitted?

upload_to_aws = function(path,
                         bucket = "formulations-dev",
                         prefix = "hf_1.0",
                         verbose = TRUE){

  check_pkg("aws.s3")

  prefix = gsub("/$", "", prefix)

  if(!file_test("-f", path)){
    f = list.files(path,
                   full.names = TRUE,
                   recursive = TRUE)

    o = file.path(prefix, basename(path), gsub("/", "", gsub(path, "", f)))

    fin = lapply(1:length(f), function(i) {
      aws.s3::put_object(file = f[i],
                         object = o[i],
                         bucket = bucket,
                         multipart = TRUE,
                         show_progress = verbose)
    })

  } else {

    o = file.path(prefix, basename(path))

    aws.s3::put_object(file = path,
                       object = o,
                       bucket = bucket,
                       multipart = TRUE,
                       show_progress = verbose)
  }

  hyaggregate_log("SUCCESS", glue("{length(o)} file(s) uplaoded to s3!"), verbose)

}


#' Find Processing Unit (Vector or Raster)
#' @param location sf object
#' @param pu either "vpu" or "rpu" (default is "vpu")
#' @return intersection processing units (sf object)
#' @export
#' @importFrom nhdplusTools get_boundaries
#' @importFrom sf st_transform
#' @importFrom dplyr filter slice_min
#' @importFrom sbtools item_file_download

find_pu = function(location, pu = "vpu"){ get_boundaries(pu)[st_transform(location, 4269),] }

#' Extract File Extension
#' @details returns file extension
#' @param x a file path
#' @param prefix character string to precede extracted extension.  default = "". (Usefull if you want to keep the ".")
#' @return character string
#' @export

.getExtension = function (x, prefix = "") {
  ext <- strsplit(basename(x), split = "\\.")[[1]]
  return(paste0(prefix, ext[length(ext)]))
}

#' Authenticated Upload to ScienceBase
#' @param x A string vector of paths to files to be uploaded
#' @param item a character ScienceBase ID corresponding to the item (default = '629a4246d34ec53d276f446d')
#' @importFrom sbtools item_replace_files

upload_to_sb = function(x, item = "629a4246d34ec53d276f446d") {
  authenticate_sb(Sys.getenv("sb_user", unset=""), Sys.getenv("sb_pass", unset=""))

  item_replace_files(item, x)
}
