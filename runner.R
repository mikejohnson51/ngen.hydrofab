pacman::p_load(glue, hydrofab, archive, sf, nhdplusTools, dplyr)
devtools::load_all()
sf_use_s2(FALSE)

## -- STEP 1: VPU base NextGen Artifacts --- #

#type = 'uniform_national'
type = 'global_uniform'

base = '/Volumes/Transcend/ngen/CONUS-hydrofabric/calibration'

f = list.files(base,
               full.names = TRUE,
               pattern = type)
f = f[grepl(".gpkg$", f)]

dem = correct_nwm_spatial(raster::raster('/Volumes/Transcend/nwmCONUS-v216/geo_em.d01_1km.nc', varname = "HGT_M"))


for(i in 1:length(f)){

  VPU = gsub(".gpkg", "", gsub(paste0(type, "_"), "", basename(f[i])))

  message(VPU)

  outfile = glue("{base}/nextgen_{VPU}.gpkg")

  dir.create(dirname(outfile), showWarnings = FALSE)

  ngen    = apply_nexus_topology(gpkg = f[i])

  outfile = write_hydrofabric(ngen, outfile)
  outfile = add_flowpath_attributes(outfile)
  outfile = add_lake_attributes(outfile, lake_path = '/Volumes/Transcend/nwmCONUS-v216/LAKEPARM_CONUS.nc')

  outfile = add_cfe_noahowp_attributes(gpkg      = outfile,
                        dir                 = '/Volumes/Transcend/nwmCONUS-v216/',
                        add_to_gpkg         = TRUE,
                        add_weights_to_gpkg = FALSE,
                        log                 = TRUE,
                        overwrite = TRUE)

  outfile = add_forcing_attributes(outfile, dem)

  st_layers(outfile)
#
#   dir = write_ngen_dir(outfile, overwrite = TRUE)
#   unlink(paste0(dir, ".zip"))
#   archive_write_dir(paste0(dir, ".zip"), dir)

}


## -- STEP 2: National Merge --- ##

outfile = national_merge(files, glue("{base}/conus.gpkg"))
dir = write_ngen_dir(outfile)
archive_write_dir(paste0(dir, ".zip"), dir)


source("secret/aws_creds.R")
aws.dir = 's3://nextgen-hydrofabric/v1.2/'

system(glue('aws s3 sync {base} {aws.dir}'))


## -- STEP 3: CAMEL EXTRACTS --- ##

gpkg = glue("{base}/nextgen_{VPU}.gpkg") #/Volumes/Transcend/ngen/CONUS-hydrofabric/calibration/nextgen_01.gpkg

c  = read_sf(gpkg, "lookup_table") %>%
  filter(POI_TYPE == "Gages")

dir.create('/Volumes/Transcend/ngen/CONUS-hydrofabric/calibration/subsets')

for(i in 1:nrow(c)){

  here = glue("{base}/subsets/gauge_{c$POI_VALUE[i]}.gpkg")

  x = subset_network(gpkg             = gpkg,
                     origin           = c$toid[i],
                     attribute_layers = c("flowpath_attributes",
                                          "lake_attributes",
                                          "cfe_noahowp_attributes",
                                          "forcing_attributes"),
                     export_gpkg = here)

  message(i, " of ", nrow(c))

}





log_info("{nrow(c)} camels basins to process!")

for(j in 1:nrow(c)){

  cc = c[j,]

  if(is.na(cc$VPUID)){
    new = st_nearest_feature(cc, vpus)
    cc$VPUID = vpus$VPUID[new]
  }

  here = glue("{base}/camels/VPU{cc$VPUID}/gauge_{cc$gauge_id}.gpkg")
  dir.create(dirname(here), recursive = TRUE, showWarnings = FALSE)
  gpkg = glue('{base}/nextgen_{cc$VPUID}.gpkg')

  x = subset_network(gpkg             = gpkg,
                     origin           = find_origin(gpkg, cc, "divides"),
                     export_gpkg      = here,
                     attribute_layers = c("cfe_noahowp_attributes", 'forcing_metadata'),
                     overwrite = FALSE)

  log_info(here)
}






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

  here = glue("{base}/camels/VPU{cc$VPUID}/gauge_{cc$gauge_id}.gpkg")
  dir.create(dirname(here), recursive = TRUE, showWarnings = FALSE)
  gpkg = glue('{base}/nextgen_{cc$VPUID}.gpkg')

  x = subset_network(gpkg             = gpkg,
                     origin           = find_origin(gpkg, cc, "divides"),
                     export_gpkg      = here,
                     attribute_layers = c("cfe_noahowp_attributes", 'forcing_metadata'),
                     overwrite = FALSE)

  log_info(here)
}



