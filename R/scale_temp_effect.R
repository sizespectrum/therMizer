### Functions related to creating the scaling parameters

#' @title Set encounter scalar
#'
#' @description Creates the encounterpred_scale parameter which is
#' used for scaling encounter and mortality rates and set the temperature scalar
#' between 0 and 1.
#'
#' @inheritParams scaled_temp_effect
#'
#' @export
#'
setEncounterPredScale <- function(params){
  species_params(params)$encounterpred_scale <- rep(NA, length(species_params(params)$temp_min))
  for (indv in seq(1:length(species_params(params)$temp_min))) {

    # Create a vector of all temperatures each species might encounter
    # Convert from degrees Celsius to Kelvin
    temperature <- seq(species_params(params)$temp_min[indv], species_params(params)$temp_max[indv], by = 0.1) + 273

    # Find the maximum value of the unscaled effect of temperature on encounter and predation rate for each species
    species_params(params)$encounterpred_scale[indv] <-
      max((temperature) * (temperature - (species_params(params)$temp_min[indv] + 273)) *
            ((species_params(params)$temp_max[indv] + 273) - temperature)^(1/2))

  }
  return(params)
}

#' @title Metabolism temperature
#'
#' @description Determine the minimum, maximum, and range of value for the
#' effect of temperature on metabolism.
#'
#' @inheritParams scaled_temp_effect
#'
#' @export
#'

setMetabTher <- function(params){
  min_metab_value <- (exp(25.22 - (0.63/((8.62e-5)*(273 + species_params(params)$temp_min)))))
  max_metab_value <- (exp(25.22 - (0.63/((8.62e-5)*(273 + species_params(params)$temp_max)))))

  species_params(params)$metab_min <- min_metab_value
  species_params(params)$metab_range <- max_metab_value - min_metab_value

  return(params)
}


#' @title Temperature scaling factor
#'
#' @description Calculate the temperature scaling factor for the encounter rate
#' and predation rate.
#'
#' @param params An object of class \linkS4class{MizerParams}.
#' @param t Time
#'
#' @export
#'
scaled_temp_effect <- function(params, t) {
  # print("t")
  # print(t)
  # checking that t is within ocean_temp, defaulting to first value otherwise
  # TODO update with monthly version
  if(!round(t) %in% as.numeric(dimnames(other_params(params)$ocean_temp)[[1]]))
    t = as.numeric(dimnames(other_params(params)$ocean_temp)[[1]])[1] + t

  scaled_temp_effect_realms <- array(NA, dim = c(dim(other_params(params)$vertical_migration)),
                                     dimnames = c(dimnames(other_params(params)$vertical_migration)))

  # Using t+1 to avoid calling ocean_temp[0,] at the first time step
  # Looping through each realm
  nb_realms <- dim(other_params(params)$vertical_migration)[1]
  for (r in seq(1, nb_realms, 1)) {
    index <- which.min(abs(as.numeric(dimnames(other_params(params)$ocean_temp)[[1]])  - t))
    temp_at_t <- other_params(newP)$ocean_temp[index,r] + 273

    # Calculate unscaled temperature effect using a generic polynomial rate equation
    unscaled_temp_effect <-
      temp_at_t * (temp_at_t - (species_params(params)$temp_min + 273)) *
      (species_params(params)$temp_max - temp_at_t + 273)^(1/2)

    # Scale using new parameter
    scaled_temp_effect_r <- unscaled_temp_effect / species_params(params)$encounterpred_scale

    # Set temperature effect to 0 if temperatures are outside thermal tolerance limits
    above_max <- (temp_at_t - 273) > species_params(params)$temp_max
    below_min <- (temp_at_t - 273) < species_params(params)$temp_min

    scaled_temp_effect_r[above_max | below_min] = 0
    scaled_temp_effect_realms[r,,] <- scaled_temp_effect_r *
      other_params(params)$exposure[r,] * other_params(params)$vertical_migration[r,,]
  }

  # print("index")
  # print(index)
  # print("temperature")
  # print(temp_at_t - 273)

  scaled_temp_effect <- colSums(scaled_temp_effect_realms)

  return(scaled_temp_effect)
}
