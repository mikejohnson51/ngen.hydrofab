% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/routing.R
\name{add_flowpath_attributes}
\alias{add_flowpath_attributes}
\title{Length Average Routelink Variables}
\usage{
add_flowpath_attributes(
  gpkg,
  rl_vars = c("comid", "Qi", "MusK", "MusX", "n", "So", "ChSlp", "BtmWdth", "time",
    "Kchan", "nCC", "TopWdthCC", "TopWdth"),
  rl_path = get_routelink_path()
)
}
\arguments{
\item{rl_vars}{routelink variables}

\item{rl_path}{routelink path (see get_routelink_path())}

\item{flowpaths}{sf LINESTRING}
}
\value{
data.frame
}
\description{
Length Average Routelink Variables
}
