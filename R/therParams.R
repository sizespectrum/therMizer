### Function aiming to upgrade a default mizer object to one able to work with the therMizer extension


# species_params = data.frame(species = c("speciesA", "speciesB"), w_inf = c(500, 5000), k_vb = c(0.8, 0.3), w_min = c(0.001, 0.001), w_mat = c(5, 50), beta = c(1000,100), sigma = c(3,3))
# species_params$interaction_resource <- c(1,0.5)
# params <- newMultispeciesParams(species_params, no_w = 200, kappa = 0.0001) |>
#   steady(tol = 0.001)
#
# species_params(params)$temp_min <- c(-5, 5)
# species_params(params)$temp_max <- c(10, 20)

# Create temperature array and fill it
# times <- 0:500
# ocean_temp_array <- array(NA, dim = c(length(times), length(realm_names)), dimnames = list(time = times, realm = realm_names))
# temp_inc <- 0
# for (i in 1:501) {
#   ocean_temp_array[i,] <- c(-5 + temp_inc, -5 + temp_inc, -5 + temp_inc, -5 + temp_inc)
#   temp_inc <- temp_inc + 0.1
# }


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
                              n_pp_array = NULL, vertical_migration_array = NULL){

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

  if(is.null(ocean_temp_array)) stop("You need to specify a temperature array to do the projections.")
  else other_params(params)$ocean_temp <- ocean_temp_array

  if(!is.null(n_pp_array)){
    # params <- setResource(params, resource_dynamics = "plankton_forcing")
    other_params(params)$n_pp_array <- n_pp_array
  }

  if(!is.null(vertical_migration_array)) params <- setVerticality(params, vertical_migration_array)

  params <- setRateFunction(params, "Encounter", "therMizerEncounter")
  params <- setRateFunction(params, "PredRate", "therMizerPredRate")
  params <- setRateFunction(params, "EReproAndGrowth", "therMizerEReproAndGrowth")

  return(params)

}
