% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfe_noahowp_summary.R
\name{add_cfe_noahowp_attributes}
\alias{add_cfe_noahowp_attributes}
\title{Add CFE and NOAH-OWP attributes}
\usage{
add_cfe_noahowp_attributes(
  gpkg = NULL,
  nwm_dir = NULL,
  catchment_name = NULL,
  flowline_name = NULL,
  precision = 9,
  outfile = NULL,
  verbose = TRUE,
  log = TRUE
)
}
\arguments{
\item{gpkg}{a geopackage with aggregation units}

\item{catchment_name}{the layer name of the aggregation units}

\item{flowline_name}{the layer name of the flowpath units}

\item{precision}{the precision of the computations}

\item{single_layer}{should only the top layer of a multilayer parameter be processed?}

\item{dir}{directory of expected nwm_grid_file (see hyAggregate::get_nwm_grids())}
}
\description{
Add CFE and NOAH-OWP attributes
}
