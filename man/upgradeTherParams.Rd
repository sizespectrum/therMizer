% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/therParams.R
\name{upgradeTherParams}
\alias{upgradeTherParams}
\title{Upgrade to thermizer object}
\usage{
upgradeTherParams(
  params,
  temp_min = NULL,
  temp_max = NULL,
  ocean_temp_array = NULL,
  n_pp_array = NULL,
  vertical_migration_array = NULL,
  exposure_array = NULL,
  aerobic_effect = TRUE,
  metabolism_effect = TRUE
)
}
\arguments{
\item{params}{A mizer params object}

\item{ocean_temp_array}{An array of temperatures}

\item{n_pp_array}{An array of plankton forcing}

\item{vertical_migration_array}{An array of number of realms x number of
species x number of sizes filled with the time ratio of each species spending
in each realms. Values must be positive and the sum of every realms per
species per size must be one}

\item{aerobic_effect}{Boolean value which determines if encounter rate is
affected by temperature. Default is TRUE.}

\item{metabolism_effect}{Boolean value which determines if metabolism rate is
affected by temperature. Default is TRUE.}
}
\description{
Wrapper function making a mizer params object
into something useable by the package
}
