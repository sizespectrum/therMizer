### Function aiming to upgrade a default mizer object to one able to work with the therMizer extension

#' @title Upgrade to thermizer object
#'
#' @description Wrapper function making a mizer params object
#' into something useable by the package
#'
#' @param params A mizer params object
#' @param temp_min A vector of numeric values the same length as the number of
#' species in params. Contains the minimum temperature range per species.
#' @param temp_max A vector of numeric values the same length as the number of
#' species in params. Contains the maximum temperature range per species.
#' @param ocean_temp_array A scalar, vector or array of temperature. The first dimension
#' must be time. If it is an array, the second dimension must be realms. Numeric
#' dimnames for time are assumed to be years. Use dates otherwise.
#' @param n_pp_array An array of plankton forcing of dimensions time x sizes
#' @param vertical_migration_array An array of number of realms x number of
#' species x number of sizes filled with the time ratio of each species spending
#' in each realms. Values must be positive and the sum of every realms per
#' species per size must be one. Default is species spending equal time between
#' all realms.
#' @param exposure_array An array of number of realms x species filled with values
#' from 0 to 1. This array scales the effect of temperature for each species per
#' realms. Default is ones.
#' @param aerobic_effect Boolean value which determines if encounter rate is
#' affected by temperature. Default is TRUE.
#' @param metabolism_effect Boolean value which determines if metabolism rate is
#' affected by temperature. Default is TRUE.
#'
#' @export

upgradeTherParams <- function(params, temp_min = NULL, temp_max = NULL,
                              ocean_temp_array = NULL,
                              n_pp_array = NULL,
                              vertical_migration_array = NULL,
                              exposure_array = NULL,
                              aerobic_effect = TRUE,
                              metabolism_effect = TRUE){
  no_sp <- length(species_params(params)$species)

  ## temperature parameters
  if(is.null(temp_min)){
    if(is.null(species_params(params)$temp_min)) stop("You need to setup min temperature for your species.")
  } else if(length(temp_min) != no_sp) {
    stop("The length of temp_min is not the same as the number of species.")
  } else {species_params(params)$temp_min <- temp_min}


  if(is.null(temp_max)){
    if(is.null(species_params(params)$temp_max)) stop("You need to setup max temperature for your species.")
  } else if(length(temp_max) != no_sp) {
    stop("The length of temp_max is not the same as the number of species.")
  } else {species_params(params)$temp_max <- temp_max}

  params <- setEncounterPredScale(params)
  params <- setMetabTher(params)

  ## temperature data array
  if(is.null(ocean_temp_array)) stop("You need to specify a temperature array to do the projections.")

  ## check dimension
  if(is.vector(ocean_temp_array)){
    ## check dimnames
    if(is.null(names(ocean_temp_array))){
      warning("ocean_temp_array has no dates. They are assumed to be successive years staring from 0.")
      names(ocean_temp_array) <- sprintf("%04d", seq(0, length.out = length(ocean_temp_array)))
    }
    date_vec <- names(ocean_temp_array)
    parse_order <- ifelse(nchar(date_vec) <= 4, "%Y",
                          ifelse(nchar(date_vec) == 7, "%Y-%m",
                                 "%Y-%m-%d"))
    date_vec <- parse_date_time(date_vec, orders = parse_order)
    date_vec <- as.numeric(year(date_vec)) + as.numeric(yday(date_vec) - 1) / as.numeric(ifelse(leap_year(date_vec), 366, 365))

    ocean_temp_array <- matrix(ocean_temp_array,
                               nrow = length(ocean_temp_array),
                               ncol = dim(params@species_params)[1],
                               byrow = F,
                               dimnames = list("time" = date_vec, "species" = params@species_params$species))

  } else if(is.array(ocean_temp_array))
  {
    # assuming that the second dimension is realms. If there is an issue with that it will be signaled later
    ## check dimnames
    date_vec <- dimnames(ocean_temp_array)[[1]]

    parse_order <- ifelse(nchar(date_vec) == 4, "%Y",
                          ifelse(nchar(date_vec) == 7, "%Y-%m",
                                 "%Y-%m-%d"))

    date_vec <- parse_date_time(date_vec, orders = parse_order)
    date_vec <- as.numeric(year(date_vec)) + as.numeric(yday(date_vec) - 1) / as.numeric(ifelse(leap_year(date_vec), 366, 365))
    dimnames(ocean_temp_array)[[1]] <- date_vec
  } else stop("The ocean_temp_array is of the wrong format")



  other_params(params)$ocean_temp <- ocean_temp_array

  ## background resource
  if(!is.null(n_pp_array)){
    if(!dim(ocean_temp_array)[1] == dim(n_pp_array)[1])
      stop("The time dimension of ocean_temp_array and n_pp_array must be equal.")

    if(!dim(n_pp_array)[2] == length(params@w_full))
      stop("The size dimension of the n_pp_array must be the same as w_full.")

    if(!identical(dimnames(n_pp_array)[[2]],params@w_full)){
      warning("The dimnames of the second dimension of n_pp_array are not the same
      as the ones in w_full. Since the dimension size is the same, the dimnames
      have been replaced by w_full.")
      dimnames(n_pp_array)[[2]] <- params@w_full
    }

    if(!is.array(n_pp_array) && !is.matrix(n_pp_array)) n_pp_array <- as.matrix(n_pp_array)

    params <- setResource(params, resource_dynamics = "plankton_forcing")
    # time names of ocean_temp_array might have been modified for thermizer to work, doing the same to n_pp
    dimnames(n_pp_array)[[1]] <- dimnames(ocean_temp_array)[[1]]

    other_params(params)$n_pp_array <- n_pp_array
  }

  # realms
  if(!is.null(vertical_migration_array))
  {
    params <- setVerticality(params, vertical_migration_array, exposure_array)
  } else if(!is.null(exposure_array)){

    vertical_fill <- 1/dim(exposure_array)[1]
    vertical_names <- dimnames(exposure_array)[[1]]
    vertical_dim <- c(dim(exposure_array),
                      length(params@w))

    vertical_migration_array <-
      array(vertical_fill,
            dim = vertical_dim,
            dimnames = list("realm" = vertical_names,
                            "species" = params@species_params$species,
                            "w" = params@w))
    message("exposure_array used to create vertical_migration_array.")
    params <- setVerticality(params, vertical_migration_array, exposure_array)
  } else {
    if(is.null(dim(ocean_temp_array))){
      vertical_fill <- rep(c(rep(c(1,rep(0,length(params@species_params$species))),
                                 length(params@species_params$species)-1),1),length(params@w))
      vertical_names <- params@species_params$species
      vertical_dim <- c(rep(length(params@species_params$species),2), length(params@w))
      message("Assuming one realm per species to create vertical_migration_array.")
    } else if(!isTRUE(all.equal(dimnames(ocean_temp_array)[[2]],params@species_params$species))){
      vertical_fill <- 1/length(dimnames(ocean_temp_array)[[2]])
      vertical_names <- dimnames(ocean_temp_array)[[2]]
      vertical_dim <- c(dim(ocean_temp_array)[2],
                        length(params@species_params$species),
                        length(params@w))
      message("Your realm names don't match your species names, so species are assumed to spend equal time in all realms.")
    } else {
      vertical_fill <- rep(c(rep(c(1,rep(0,length(params@species_params$species))),
                                 length(params@species_params$species)-1),1),length(params@w))
      vertical_names <- params@species_params$species
      vertical_dim <- c(rep(length(params@species_params$species),2), length(params@w))
    }

    vertical_migration_array <-
      array(vertical_fill,
            dim = vertical_dim,
            dimnames = list("realm" = vertical_names,
                            "species" = params@species_params$species,
                            "w" = params@w))
    params <- setVerticality(params, vertical_migration_array)
  }

  ## rate functions
  if(aerobic_effect){
    params <- setRateFunction(params, "Encounter", "therMizerEncounter")
    params <- setRateFunction(params, "PredRate", "therMizerPredRate")
  }

  if(metabolism_effect) params <- setRateFunction(params, "EReproAndGrowth", "therMizerEReproAndGrowth")

  ## time dimension
  other_params(params)$t_idx = - as.numeric(dimnames(other_params(params)$ocean_temp)[[1]][1])

  return(params)
}

#' #' @title Project thermizer object
#' #'
#' #' @description Wrapper function adjusting simulation time and start time
#' #' for the project function
#' #'
#' #' @inheritParams upgradeTherParams
#' #'
#' #' @export
#' #'
#' therProject <- function(params){
#'   sim_times <- c(as.numeric(dimnames(other_params(params)$ocean_temp)[[1]][1]),
#'                  dim(other_params(params)$ocean_temp)[1])
#'
#'   cat(sprintf("The simulation is set to start in %d and will run for %d years.\n",sim_times[1], sim_times[2]))
#'
#'   sim <- project(params, t_start = sim_times[1], t_max = sim_times[2]-1)
#' }
