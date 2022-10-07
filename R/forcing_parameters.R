#' Generate Forcing Parameters
#' @description The NextGen forcing workstream team has been working on regridding the GFS meteorological forcings data into NextGen catchments.
#' There is a preliminary step that we must take to “downscale” the regridded data into the NextGen catchments based on previous NCAR formulas
#' that were integrated into the WRFHydro Forcings Engine for the National Water Model (NWM). The following characteristics are needed to downscale GFS data for each NextGen catchment:
#' Mean elevation (meters); Catchment Latitude centroid; Catchment Longitude Centroid; Mean Catchment Slope (meters/kilometer); Circular mean azimuth (degrees);
#' Catchment area (km2)
#' @param gpkg hydrofabric geopackage
#' @param dem digital elevation model (DEM) SpatRaster to process. Units should be meters
#' @param verbose emit messages?
#' @return data.table
#' @export
#' @importFrom terra crop project vect crs terrain
#' @importFrom zonal weight_grid execute_zonal
#' @importFrom dplyr mutate select right_join left_join
#' @importFrom sf st_centroid st_coordinates st_drop_geometry
#' @importFrom hydrofab add_areasqkm

add_forcing_attributes = function(gpkg, dem, verbose = TRUE){

  geom = read_hydrofabric(gpkg, realization = "catchments")[[1]]

  DEM_m = crop(dem, project(vect(geom), crs(dem)))

  slope_m_km = 1000 * terrain(DEM_m, v = "slope", unit = "radians")

  aspect = terrain(DEM_m, v = "aspect", unit = "degree")

  hyaggregate_log("INFO", "Building Weight Grid", verbose)
  w = weight_grid(DEM_m, geom, ID = "id", progress = verbose)

  hyaggregate_log("INFO", "Building summaries", verbose)
  summary1 = execute_zonal(c(DEM_m, slope_m_km),
                           w = w,
                           fun = "mean",
                           ID = "id")


  names(summary1) = c("id", "elevation", "slope_m_km")


  summary2 = execute_zonal(aspect,
                           w = w,
                           fun = zonal:::circular_mean,
                           ID = "id")

  out = suppressWarnings({
    geom %>%
      mutate(areasqkm = add_areasqkm(.)) %>%
      select(id, areasqkm) %>%
      st_centroid() %>%
      st_transform(4326) %>%
      mutate(cetroid_lon = st_coordinates(.)[,1], centroid_lat = st_coordinates(.)[,2]) %>%
      st_drop_geometry( ) %>%
      right_join(left_join(summary1, summary2, by = "id"), by = "id")
  })


  write_sf(out, gpkg, "forcing_attributes", overwrite = TRUE)

  gpkg
}
