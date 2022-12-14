% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realms.R
\name{setVerticality}
\alias{setVerticality}
\title{Add realms to params object.}
\usage{
setVerticality(params, vertical_migration_array, exposure_array = NULL)
}
\arguments{
\item{params}{A mizer params object}

\item{vertical_migration_array}{An array of number of realms x number of
species x number of sizes filled with the time ratio of each species spending
in each realms. Values must be positive and the sum of every realms per
species per size must be one}
}
\description{
set different species specific realms with different temperatures
}
\details{
Vertical integration functions
}
