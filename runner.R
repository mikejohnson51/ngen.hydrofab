pacman::p_load(hydrofabric)
devtools::load_all()
sf_use_s2(FALSE)

## -- STEP 1: VPU base NextGen Artifacts --- #

folder = '04_global_uniform'

base = glue('/Volumes/MyBook/nextgen/global_uniform')

nwm_dir = '/Volumes/Transcend/nwmCONUS-v216'

process = data.frame(files = list.files(base, full.names = TRUE, pattern = ".+gpkg$")) %>%
  mutate(VPU = gsub(".gpkg", "", gsub(paste0("uniform", "_"), "", basename(files))),
         hf_outfile =  glue('/Volumes/MyBook/v20/v21/gpkg/nextgen_{VPU}.gpkg'))

glue("{dirname(base)}/05_nextgen/nextgen_{VPU}.gpkg")


# Produce Nextgen Fabrics -------------------------------------------------

for(i in 1:nrow(process)){

  if(!file.exists(process$hf_outfile[i])){
    dir.create(dirname(process$hf_outfile[i]), showWarnings = FALSE)

     o    = apply_nexus_topology(gpkg        = process$files[i],
                                 vpu         = process$VPU[i],
                                 export_gpkg = process$hf_outfile[i]) %>%
         add_flowpath_attributes() %>%
         add_lake_attributes(lake_path = '/Volumes/Transcend/nwmCONUS-v216/LAKEPARM_CONUS.nc') %>%
         append_style(layer_names = c("nexus", "hydrolocations", "flowpaths", "divides", "lakes"))

     bas = read_sf(o, "divides")
     bas = clean_geometry(bas, ID = "divide_id")
     write_sf(bas, o, "divides")
  }
}


national_merge(gpkg = process$hf_outfile, outfile = '/Volumes/MyBook/v20/v21/conus.gpkg')

 # Build National Topology -------------------------------------------------
vaa = get_vaa("hydroseq") %>%
  rename(hf_id = comid, hf_hydroseq = hydroseq)

xx = "/Volumes/MyBook/v20/v21/conus.gpkg" %>%
  read_sf('network') %>%
  left_join(vaa, by = "hf_id")

write_parquet(xx, "/Volumes/MyBook/v20/v21/conus_net.parquet")



# Build Supplementary Data -------------------------------------------------
for(i in 1:nrow(process)){
o    = add_cfe_noahowp_attributes(gpkg        = process$hf_outfile[i],
                                  nwm_dir     = "data",
                                  outfile     = glue('/Volumes/Transcend/ngen/CONUS-hydrofabric/05_nextgen/nextgen_{process$VPU[i]}_cfe_noahowp.parquet'))
message(i)
}

dem = correct_nwm_spatial("data/dem.tif")
forcing_files = list.files('/Users/mjohnson/Downloads/extended_tifs', full.names = TRUE)

for(i in 1:nrow(process)) {

  if (any(
    !layer_exists(process$ext_outfile[i], 'forcing_dem_attributes'),
    !layer_exists(process$ext_outfile[i], 'nwm_weight_grid')
  )) {

    add_forcing_attributes(
      process$hf_outfile[i],
      export_gpkg = process$ext_outfile[i],
      add_grid = 'nwm_weight_grid',
      dem
    )
  }

  hy   = read_hydrofabric(process$hf_outfile[i], realization = "catchments")[[1]]

  for (j in 1:length(forcing_files)) {
    n = paste0(gsub(".tif", "", basename(forcing_files[j])), "_weight_grid")

      if (!layer_exists(process$ext_outfile[i], n)) {
        out = weight_grid(
          forcing_files[j],
          geom = hy,
          ID = "divide_id",
          progress = TRUE
        )

        tmp = terra::rowColFromCell(rast(forcing_files[j]), out$cell) %>%
          as.data.frame() %>%
          setNames(c("row", "col"))

        write_sf(cbind(out, tmp), process$ext_outfile[i], n)
      }
    }

  # if(!layer_exists(process$ext_outfile[i], 'twi_distributions')){
  #
  #   system.time({
  #     o = warp('/vsis3/nextgen-hydrofabric/DEM-products/twi.vrt', hy, 500)
  #
  #     xx =  execute_zonal(rast(o),
  #                         geom = hy,
  #                         fun = zonal:::binned_json,
  #                         ID = "divide_id",
  #                         area_weight = FALSE,
  #                         join = TRUE)
  #     unlink(o)
  #   })
  #
  #   write_sf(xx, process$ext_outfile[i], "twi_distributions")
  # }
}

## -- STEP 2: National Merge --- ##
unlink(glue("{dirname(base)}/05_nextgen/conus.gpkg"))
outfile = nat_merge(gpkg = process$hf_outfile,
                    outfile = glue("{dirname(base)}/05_nextgen/conus.gpkg"))


outfile = nat_merge(gpkg = process$ext_outfile,
                    outfile = glue("{dirname(base)}/05_nextgen/conus_ext.gpkg"))


nat_merge = function(gpkg, outfile, verbose = TRUE){

  all = st_layers(gpkg[1])$name

  all = all[all != "layer_styles"]

  for(i in 1:length(all)){
    hyaggregate_log("INFO", glue("Merging: {all[i]}"), verbose)
    tmp = list()
    for(j in 1:length(gpkg)){
      tmp[[j]] = read_sf(gpkg[j], all[i])
      message("\tread ", j," of ", length(gpkg))
    }

    write_sf(bind_rows(tmp), outfile, all[i])
    rm(tmp)
  }

  outfile
}

geom = read_sf('/Volumes/Transcend/ngen/CONUS-hydrofabric/pre-release/camels/VPU01/gauge_01022500.gpkg',
               "divides")
dem = opendap.catalog::dap('/vsis3/nextgen-hydrofabric/dem.vrt', AOI = geom)

plot(dem)
plot(geom$geom, add=T)


## -- STEP 3: CAMEL EXTRACTS --- ##

vpus = get_boundaries()

c  = read.csv('/Users/mjohnson/github/hydrofabric_attributes/data/camels_compiled.csv') |>
  st_as_sf(coords = c('gauge_lon', 'gauge_lat'), crs = 4269) |>
  mutate(gauge_id = sprintf("%08s", gauge_id)) |>
  st_join(vpus) %>%
  arrange(VPUID)

log_info("{nrow(c)} camels basins to process!")

for(j in 1:nrow(c)){

  cc = c[j,]

  if(is.na(cc$VPUID)){
    new = st_nearest_feature(cc, vpus)
    cc$VPUID = vpus$VPUID[new]
  }

  here = glue("{dirname(base)}/camels/VPU{cc$VPUID}/gauge_{cc$gauge_id}.gpkg")
  dir.create(dirname(here), recursive = TRUE, showWarnings = FALSE)
  gpkg = glue('{dirname(base)}/nextgen/nextgen_{cc$VPUID}.gpkg')

  origin =   read_sf(gpkg, "hydrolocations_lookup") %>%
    filter(hl_link == cc$gauge_id)

  if(nrow(origin) > 0){
    x = subset_hf(gpkg             = gpkg,
                  origin           = origin$id,
                  export_gpkg      = here)
  }


  log_info(here)
}



