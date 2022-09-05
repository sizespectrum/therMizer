#' Vertical integration functions
#'


#' @title set verticality
#'
#' @description set different species specifc realms with different temperatures
#'
#' @param params mizer params object
#' @param real_names A character vector containing the names of the
#' different realms to create.
#'
#' @export
#'

setVerticality <- function(params, realm_names){
# realm_names <- c("upper50m","bottom","DVM_day","DVM_night")
species_names <- as.character(params@species_params$species)
sizes <- params@w

# Create the vertical migration array and fill it
vertical_migration_array <- array(0, dim = (c(length(realm_names), length(species_names), length(sizes))), dimnames = list(realm = realm_names, sp = species_names, w = signif(sizes,3))) # realm x species x size

upp <- which(realm_names == "upper50m") # 0 - 50m average
btm <- which(realm_names == "bottom") # sea floor
DVMd <- which(realm_names == "DVM_day") # 200 - 500m average
DVMn <- which(realm_names == "DVM_night") # 0 - 100m average

# Set all sizes below w_mat for speciesA to "upper50m" and all sizes above w_mat to "bottom
spA <- which(species_names == "speciesA")
vertical_migration_array[upp, spA, sizes < params@species_params$w_mat[spA]] <- 1
vertical_migration_array[btm, spA, sizes >= params@species_params$w_mat[spA]] <- 1

# Have speciesB split its time equally using DVM
spB <- which(species_names == "speciesB")
vertical_migration_array[DVMd, spB, ] <- 0.5
vertical_migration_array[DVMn, spB, ] <- 0.5

# And check that all filled size classes sum to 1, no more and no less
for (s in seq(1,length(species_names),1)) {
  if (!all(apply(vertical_migration_array[ , s, ],2,sum) == 1)) {
    stop(paste("Your realm allocations for ", species_names[s], " don't sum to 1 for all sizes.
		Their time in a given realm is either over- or under-allocated.", sep = ""))
  }
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
