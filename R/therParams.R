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
#' @param realm_names A character vector for realms names
#'
#' @export

upgradeTherParams <- function(params, temp_min = NULL, temp_max = NULL, ocean_temp_array = NULL, n_pp_array = NULL, realm_names = NULL){

  if(is.null(temp_min)){
    if(is.null(species_params(params)$temp_min)) stop("You need to setup min temperature for your species")
  } else if(length(temp_min) != length(species_params(params)$species)) { stop("The length of temp_min is not the same as the number of species")
  } else {species_params(params)$temp_min <- temp_min}


  if(is.null(temp_max)){
    if(is.null(species_params(params)$temp_max)) stop("You need to setup max temperature for your species")
  } else if(length(temp_max) != length(species_params(params)$species)) { stop("The length of temp_max is not the same as the number of species")
  } else {species_params(params)$temp_max <- temp_max}

  params <- setEncounterPredScale(params)
  params <- setMetabTher(params)

  if(!is.null(ocean_temp_array)) other_params(params)$ocean_temp <- ocean_temp_array
  if(!is.null(plankton_forcing)){
    # params <- setResource(params, resource_dynamics = "plankton_forcing")
    other_params(params)$n_pp_array <- n_pp_array
  }
  realm_names <- c("upper50m","bottom","DVM_day","DVM_night")

  if(!is.null(real_names)) params <- setVerticality(params, real_names)

  params <- setRateFunction(params, "Encounter", "therMizerEncounter")
  params <- setRateFunction(params, "PredRate", "therMizerPredRate")
  params <- setRateFunction(params, "EReproAndGrowth", "therMizerEReproAndGrowth")

  return(params)

}
