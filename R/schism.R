process_schism = function(schism_file, outfile, crs = 4326, overwrite = FALSE){

  if(!file.exists(outfile) | overwrite){
    nc = open.nc(schism_file)
    on.exit(close.nc(nc))

    data.frame(
      X = var.get.nc(nc, "x"),
      Y = var.get.nc(nc, "y"),
      Z = var.get.nc(nc, "depth"),
      elev = var.get.nc(nc, "elev")
    ) %>%
      mutate(#nodeID = 1:n(),
        Z = Z + elev,
        elev = NULL) %>%
      data.table::fwrite(outfile)
  }

  return(outfile)
}


interpolate_schism = function(AOI, schism = gpkg, outfile = NULL, algo = "invdist", res = c(.001,.001)){


  if(grepl("gpkg", schism)){
    meta =  st_layers(schism)
    bb   = st_bbox(st_transform(AOI, meta$crs[[1]]))

    g = glue("gdal_grid -a {algo} -clipsrc {bb$xmin} {bb$ymin} {bb$xmax} {bb$ymax} -txe {bb$xmin} {bb$xmax} -tye {bb$ymin} {bb$ymax} -tr {res[1]} {res[1]} -of GTiff -ot Float64 -l {meta$name} -zfield Z {schism} {outfile}")

  } else if(grepl("fgb", schism)){

    bb   = st_bbox(AOI)
    g = glue("gdal_grid -a {algo} -clipsrc {bb$xmin} {bb$ymin} {bb$xmax} {bb$ymax} -txe {bb$xmin} {bb$xmax} -tye {bb$ymin} {bb$ymax} -tr {res[1]} {res[1]} -of GTiff -ot Float64 -zfield Z {schism} {outfile}")
  } else {


    bb   = st_bbox(AOI)

    # The OGR CSV driver returns all attribute columns
    # with a type of string if no field type information
    # file (with .csvt extension) is available.

    vrt  = file.path(getwd(), gsub("csv", 'vrt', schism))
    csvt = file.path(getwd(), gsub("csv", 'csvt', schism))

    csvt_filecon <- file(csvt,"w")
    writeLines('"Real","Real","Real"',con=csvt_filecon)
    close(csvt_filecon)

    vrt_header = c(
      '<OGRVRTDataSource>',
      '<OGRVRTLayer name="schism_nodes">',
      glue('<SrcDataSource>{file.path(getwd(), schism)}</SrcDataSource>'),
      #'<SrcDataSource>/Volumes/Transcend/depths.csv</SrcDataSource>',
      '<GeometryType>wkbPoint</GeometryType>',
      '<LayerSRS>EPSG:4326</LayerSRS>',
      '<GeometryField encoding="PointFromColumns" x="X" y="Y" z="Z"/>',
      '<Field name="Z" src="Z" type="Real"/>',
      '<Field name="X" src="X" type="Real"/>',
      '<Field name="Y" src="Y" type="Real"/>',
      '</OGRVRTLayer>',
      '</OGRVRTDataSource>'
    )


    vrt_filecon <- file(vrt,"w")
    writeLines(vrt_header,con=vrt_filecon)
    close(vrt_filecon)

    g = glue("gdal_grid -a {algo} -clipsrc {bb$xmin} {bb$ymin} {bb$xmax} {bb$ymax} -txe {bb$xmin} {bb$xmax} -tye {bb$ymin} {bb$ymax} -tr {res[1]} {res[1]} -of GTiff -ot Float64 -l schism_nodes -zfield Z {vrt} {outfile}")

  }

  system(g)
  return(outfile)

}
