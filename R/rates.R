
#' @title therMizerEncounter
#'
#' @description calculates the temperature-scaled encounter rate 
#'
#' @param params mizer object
#' @param t time
#'
#' @export

therMizerEncounter <- function(params, t, ...) {

  # Calculate maximum possible encounter rate
  max_encounter <- mizerEncounter(params, t, ...)

  # Apply temperature effect
  return(max_encounter * scaled_temp_effect(params,t))

}


#' @title therMizerPredRate
#'
#' @description calculates the temperature-scaled predation rate
#'
#' @param params mizer object
#' @param t time
#'
#' @export

# Because temperature effect can vary by predator size, we need to insert the temperature scaling into the mizerPredRate function.  This code is the same as that from the mizer github page (24 Aug 22), just with the scaling added in.
therMizerPredRate <- function(params, n, n_pp, n_other, t, feeding_level, ...) {
  no_sp <- dim(params@interaction)[1]
  no_w <- length(params@w)
  no_w_full <- length(params@w_full)

  # If the the user has set a custom pred_kernel we can not use fft.
  # In this case we use the code from mizer version 0.3
  if (!is.null(comment(params@pred_kernel))) {
    n_total_in_size_bins <- sweep(n, 2, params@dw, '*', check.margin = FALSE)
    # The next line is a bottle neck
    pred_rate <- sweep(params@pred_kernel, c(1,2),
                       (1 - feeding_level) * params@search_vol *
                         n_total_in_size_bins,
                       "*", check.margin = FALSE)

    pred_rate <- pred_rate * scaled_temp_effect(params, t)

    # integrate over all predator sizes
    pred_rate <- colSums(aperm(pred_rate, c(2, 1, 3)), dims = 1)
    return(pred_rate)
  }

  # Get indices of w_full that give w
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  # We express the result as a a convolution  involving
  # two objects: Q[i,] and ft_pred_kernel_p[i,].
  # Here Q[i,] is all the integrand of (3.12) except the feeding kernel
  # and theta
  Q <- matrix(0, nrow = no_sp, ncol = no_w_full)
  # We fill the end of each row of Q with the proper values
  Q[, idx_sp] <- sweep( (1 - feeding_level) * params@search_vol * n, 2,
                        params@dw, "*")

  Q[, idx_sp] <- Q[, idx_sp] * scaled_temp_effect(params, t)

  # We do our spectral integration in parallel over the different species
  pred_rate <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_p) *
                                  mvfft(base::t(Q)), inverse = TRUE))) / no_w_full
  # Due to numerical errors we might get negative or very small entries that
  # should be 0
  pred_rate[pred_rate < 1e-18] <- 0

  return(pred_rate * params@ft_mask)
}

#' @title therMizererrepro
#'
#' @description calculates the temperature-scaled energy available for 
#' growth and reproduction 
#' 
#' @param params mizer object
#' @param t time
#'
#' @export

therMizerEReproAndGrowth <- function(params, t, encounter, feeding_level, ...) {

  temp_effect_metab_realms <- array(NA, dim = c(dim(other_params(params)$vertical_migration)), dimnames = c(dimnames(other_params(params)$vertical_migration)))

  # Using t+1 to avoid calling ocean_temp[0,] at the first time step
  # Looping through each realm
  nb_realms <- dim(other_params(params)$exposure)[1]
  for (r in seq(1, nb_realms, 1)) {
    temp_at_t <- other_params(params)$ocean_temp[t + 1,r]

    # Arrhenius equation
    unscaled_temp_effect <- (exp(25.22 - (0.63/((8.62e-5)*(273 + temp_at_t)))))

    # Arrhenius equation scaled to a value between 0 and 1
    temp_effect_metabolism_r <- (unscaled_temp_effect - species_params(params)$metab_min) / species_params(params)$metab_range

    # Set temperature effect to 0 if temperatures are outside thermal tolerance limits
    above_max <- temp_at_t > species_params(params)$temp_max
    below_min <- temp_at_t < species_params(params)$temp_min

    temp_effect_metabolism_r[above_max | below_min] = 0

    temp_effect_metab_realms[r,,] <- temp_effect_metabolism_r*other_params(params)$exposure[r,]*other_params(params)$vertical_migration[r,,]
  }

  temp_effect_metabolism <- colSums(temp_effect_metab_realms)

  # Apply scaled Arrhenius value to metabolism
  sweep((1 - feeding_level) * encounter, 1,
        species_params(params)$alpha, "*", check.margin = FALSE) -
    metab(params)*temp_effect_metabolism

}

#' @title plankton forcing
#'
#' @description allows for time-varying resource
#'
#' @param params mizer object
#' @param t time
#'
#' @export

# Set up new resource forcing "function" - just changes to different time slot in the n_pp_array array
plankton_forcing <- function(params, t, ...) {
  w_cut_off <- 10
  pkt <- 10^(other_params(params)$n_pp_array[t + params@other_params$other$t_idx,])/params@dw_full # converting to density
  pkt[which(as.numeric(names(pkt)) >= w_cut_off)] <- 0

  return(pkt)
}
