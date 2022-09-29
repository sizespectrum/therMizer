### Function aiming to upgrade a default mizer object to one able to work with the therMizer extension

#' @title Upgrade to thermizer object
#'
#' @description Wrapper function making a mizer params object
#' into something useable by the package
#'
#' @param params A mizer params object
#' @param ocean_temp_array An array of temperatures
#' @param n_pp_array An array of plankton forcing
#' @param vertical_migration_array An array of number of realms x number of species x number of sizes
#' filled with the time ratio of each species spending in each realms.
#' Values must be positive and the sum of every realms per species per size must
#' be one
#'
#' @export

upgradeTherParams <- function(params, temp_min = NULL, temp_max = NULL, ocean_temp_array = NULL,
                              n_pp_array = NULL, vertical_migration_array = NULL, exposure_array = NULL){

  ## temperature parameters
  if(is.null(temp_min)){
    if(is.null(species_params(params)$temp_min)) stop("You need to setup min temperature for your species.")
  } else if(length(temp_min) != length(species_params(params)$species)) { stop("The length of temp_min is not the same as the number of species.")
  } else {species_params(params)$temp_min <- temp_min}


  if(is.null(temp_max)){
    if(is.null(species_params(params)$temp_max)) stop("You need to setup max temperature for your species.")
  } else if(length(temp_max) != length(species_params(params)$species)) { stop("The length of temp_max is not the same as the number of species.")
  } else {species_params(params)$temp_max <- temp_max}

  params <- setEncounterPredScale(params)
  params <- setMetabTher(params)

  ## temperature data arrays
  if(is.null(ocean_temp_array)) stop("You need to specify a temperature array to do the projections.")
  else other_params(params)$ocean_temp <- ocean_temp_array

  if(!is.null(n_pp_array)){
    if(!dim(ocean_temp_array)[1] == dim(n_pp_array)[1])
      stop("The time dimension of ocean_temp_array and n_pp_array must be equal.")

    if(!dim(n_pp_array[2] == length(params@w_full)))
      stop("The size dimension of the n_pp_array must be the same as w_full.")

    params <- setResource(params, resource_dynamics = "plankton_forcing")
    other_params(params)$n_pp_array <- n_pp_array
  }

  if(!is.null(vertical_migration_array)) params <- setVerticality(params, vertical_migration_array, exposure_array = exposure_array)

  ## rate functions
  params <- setRateFunction(params, "Encounter", "therMizerEncounter")
  params <- setRateFunction(params, "PredRate", "therMizerPredRate")
  params <- setRateFunction(params, "EReproAndGrowth", "therMizerEReproAndGrowth")

  ## time dimension
  other_params(params)$t_idx = - as.numeric(dimnames(ocean_temp_array)[[1]][1])

  return(params)

}

#' @title Project thermizer object
#'
#' @description Wrapper function adjusting simulation time and start time
#' for the project function
#'
#' @inheritParams upgradeTherParams
#'
#' @export
#'
therProject <- function(params){
  sim_times <- c(as.numeric(dimnames(other_params(params)$ocean_temp)[[1]][1]),
                 dim(other_params(params)$ocean_temp)[1])

  cat(sprintf("The simulation is set to start in %d and will run for %d years.\n",sim_times[1], sim_times[2]))

  sim <- project(params, t_start = sim_times[1], t_max = sim_times[2]-1)
}
