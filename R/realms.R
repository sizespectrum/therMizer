#' Vertical integration functions
#'


#' @title set verticality
#'
#' @description set different species specifc realms with different temperatures
#'
#' @inheritParams upgradeTherParams
#'
#' @export
#'

setVerticality <- function(params, vertical_migration_array){

  species_names <- as.character(params@species_params$species)
  sizes <- params@w

# check if vertical_migration_array is correct
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

exposure_array <- array(0, dim = (c(length(realm_names), length(species_names))), dimnames = list(realm = realm_names, sp = species_names)) # realm x species

for (r in seq(1,length(realm_names),1)) {
  for (s in seq(1,length(species_names),1)) {
    if (any(vertical_migration_array[r,s,] > 0)) {
      exposure_array[r,s] = 1
    }
  }
}

other_params(params)$exposure <- exposure_array

return(params)
}