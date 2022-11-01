#' Vertical integration functions
#'


#' @title Add realms to params object.
#'
#' @description set different species specific realms with different temperatures
#'
#' @inheritParams upgradeTherParams
#'
#' @export
#'

setVerticality <- function(params, vertical_migration_array, exposure_array = NULL){

  species_names <- as.character(params@species_params$species)
  realm_names <- dimnames(vertical_migration_array)$realm
  sizes <- params@w
  ocean_temp_array <- other_params(params)$ocean_temp

  # check if ocean_array has realms
  if(is.null(ocean_temp_array))
    stop("The params object needs to contain an ocean_temp_array.")

  if(is.null(dim(ocean_temp_array))){
    ocean_temp_array <- array(rep(ocean_temp_array, length(dimnames(vertical_migration_array)[[1]])),
                              dim = c(length(ocean_temp_array), length(dimnames(vertical_migration_array)[[1]])),
                              dimnames = list("time" = names(ocean_temp_array),
                                              "realm" = dimnames(vertical_migration_array)[[1]]))
    other_params(params)$ocean_temp <- ocean_temp_array
    message("ocean_temp_array was extended to a matrix with the same names as the vertical_migration_array.")
  }

  if(is.null(exposure_array) & !isTRUE(all.equal(dimnames(ocean_temp_array)[[2]],dimnames(vertical_migration_array)[[1]])))
    stop("The realm names of ocean_temp_array and vertical_migration_array are different.")

  # check if vertical_migration_array is correct
  if(dim(vertical_migration_array)[1] != dim(ocean_temp_array)[2])
    stop("The first dimension of vertical_migration_array must be equal to the number of realms in ocean_temp_array")

  if(dim(vertical_migration_array)[2] != length(species_names))
    stop("The second dimension of vertical_migration_array must be equal to the number of species")

  if(dim(vertical_migration_array)[3] != length(sizes))
    stop("The third dimension of vertical_migration_array must be equal to the number of size classes")

  # And check that all filled size classes sum to 1, no more and no less
  for (iSpecies in 1:length(species_names)) {
    if (!all(apply(vertical_migration_array[ , iSpecies, ],2,sum) == 1)) {
      stop(paste0("Your realm allocations for ", species_names[iSpecies], " don't sum to 1 for all sizes.
		Their time in a given realm is either over- or under-allocated."))
    } # the check above sometimes doesn't work, it may have to do with long decimals and rounding
  }

  other_params(params)$vertical_migration <- vertical_migration_array

  if(is.null(exposure_array)){
    exposure_array <- array(0, dim = (c(length(realm_names), length(species_names))),
                            dimnames = list(realm = realm_names, sp = species_names)) # realm x species

    for (r in seq(1,length(realm_names),1)) {
      for (s in seq(1,length(species_names),1)) {
        if (any(vertical_migration_array[r,s,] > 0)) {
          exposure_array[r,s] = 1
        }
      }
    }
  } else {
    # check if exposure_array is correct
    if(!isTRUE(all.equal(dimnames(vertical_migration_array)[[1]],dimnames(exposure_array)[[1]])))
      stop("The first dimension of exposure_array must have the same names as the first dimension of the vertical_migration_array.")

    if(dim(exposure_array)[2] != length(species_names))
      stop("The second dimension of exposure_array must be equal to the number of species.")

    if(any(exposure_array > 1 & exposure_array <0))
      stop("The values within the exposure array must be between 0 and 1.")
  }
  other_params(params)$exposure <- exposure_array



  return(params)
}
