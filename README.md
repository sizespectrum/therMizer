# Introduction to therMizer

The therMizer package is an extension of the mizer package ([Scott et al. 2014](https://doi.org/10.1111/2041-210X.12256)) that allows you to incorporate the effects of temperature on species' metabolic rate and aerobic scope.  These effects can vary by individual body size.

## Installing therMizer
Need to write this section

## Conceptualizing therMizer
Temperature affects species' metabolic rates and aerobic scope.  Yet we lack much species-specific understanding about the relationship between these rates and temperature.  The therMizer package uses approximate relationships based on species' thermal tolerance limits.  

The effects of temperature on metabolism follow a Boltzman factor or Arrhenius relationship.  At the low end of species' thermal tolerance limits, metabolism is the least expensive and more energy can be devoted to growth an reproduction.  At the high end of species thermal tolerance limits metabolism is maximally expensive and less energy is available for growth and reproduction. 

The effects of temperature on aerobic scope are approximated by a generic polynomial rate equation that results in a thermal optimum for each species.  At that temperature, species are able to realize peak aerobic performance and encounter the maximum amount of prey possible.  This falls away at temperatures to either side of the optimum, with reduced predator-prey encounter rates.

The therMizer package attempts to capture these effects based on temperatures to which species are exposed.  This exposure can vary by depth, for example for species that undergo diel vertical migration, as well as by size, for example species that undergo ontogenetic vertical migration.  

Temperature exposure in therMizer is modeled much the same way fishing mortality is modeled in mizer.  `ocean_temp` is analgous to `effort`.  It's assigned to a `realm` the same way effort is assigned to a `gear`.  `exposure` links species and realms much the same way `catchability` links species and gears.  `vertical_migration` brings in body size the same way `selectivity` does.  The following section walks through these parameters in more detail.

## Model parameters and input
This section looks in detail at the parameters and input needed to run therMizer.  There's also some sample code further below for preparing these parameters for use in therMizer.

### Species parameters
In addition to the species parameters required by mizer, you'll also need to supply `temp_min` and `temp_max` parameters for each species.  These represent the lower and upper bounds of species' thermal tolerance limits.  You can find this information in the literature, for example from tagging studies, or in databases such as [rfishbase](https://github.com/ropensci/rfishbase) ([Boettiger et al. 2012](https://doi.org/10.1111/j.1095-8649.2012.03464.x)).  

### Temperature parameters
You'll also need to three additional parameters: `realm`, `vertical_migration`, and `exposure`.

`realm` refers to the vertical realms or depth strata which species inhabit and for which temperatures will be provided.  These could be named something like _epipelagic_, _mesopelagic_, _surface_, _bottom_, whatever you want to call them.  The important thing is that they match the realm names provided in the `ocean_temp` array that you're using.

`vertical_migration` simulates the proportion of time that a given size of a given species spends in a given realm.  It has the dimensions of `realm` $\times$ `sp` $\times$ `w`.  Values can range from 0 to 1, and must sum to 1 across all realms for each species and size.  If values sum to something other than one, it means that either a portion of time is unaccounted for (<1) or that time across realms is over-allocated (>1).  Either way, you'll be modeling something that cannot be replicated in the real ocean.  

`exposure` links `vertical_migration` to `ocean_temp`.  It has the dimensions of `realm` $\times$ `sp`.  The values are 1 for the realms to which a species is exposed and 0 elsewhere.  In theory, you could set all values to 1 and, so long as `vertical_migration` is constructed correctly, get the same results (because when multiplied by `exposure` the result would still be 0 for realms in which species spend no time).  It's up to you whether you'd like to go this route.  However, you do need to ensure that the realm names and order match those used in `vertical_migration` and `ocean_temp`.

### Input
`ocean_temp` is an array that has temperature(s) in degrees Celsius for each `realm`.  It can be a vector, if temperature is constant over time, or an array for dynamic temperatures.  If you're using time-varying temperature, the array will have the dimensions of `time` $\times$ `realm`.  

## Running a scenario
Need to write this section

## Sample code for preparing parameters and input
Here's an example of how to set up the `vertical_migration` array.  We'll assume one species stays in the upper 50 m of the water column until it moves to the bottom at maturity and that all sizes of the other species undergo diel vertical migration (DVM).  This will give us four realms.
```r
realm_names <- c("upper50m","bottom","DVM_day","DVM_night")
species_names <- as.character(params_Temp@species_params$species)
sizes <- params_Temp@w

# Create the vertical migration array and fill it
vertical_migration_array <- array(0, dim = (c(length(realm_names), length(species_names), length(sizes))), dimnames = list(realm = realm_names, sp = species_names, w = signif(sizes,3))) # realm x species x size

upp <- which(realm_names == "upper50m") # 0 - 50m average
btm <- which(realm_names == "bottom") # sea floor
DVMd <- which(realm_names == "DVM_day") # 200 - 500m average
DVMn <- which(realm_names == "DVM_night") # 0 - 100m average

# Set all sizes below w_mat for speciesA to "upper50m" and all sizes above w_mat to "bottom
spA <- which(species_names == "speciesA")
vertical_migration_array[upp, spA, sizes < params_Temp@species_params$w_mat[spA]] <- 1
vertical_migration_array[btm, spA, sizes >= params_Temp@species_params$w_mat[spA]] <- 1

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
```

Using the same scenario, here's an example to set up the `exposure` array.
```r
exposure_array <- array(0, dim = (c(length(realm_names), length(species_names))), dimnames = list(realm = realm_names, sp = species_names)) # realm x species

for (r in seq(1,length(realm_names),1)) {
	for (s in seq(1,length(species_names),1)) {
		if (any(vertical_migration_array[r,s,] > 0)) {
			exposure_array[r,s] = 1
		}
	}
}
```

An example for creating the temperatures for each realm.
```r
# Create temperature array and fill it
times <- 0:500
ocean_temp_array <- array(NA, dim = c(length(times), length(realm_names)), dimnames = list(time = times, realm = realm_names))
temp_inc <- 0
for (i in 1:501) {
  ocean_temp_array[i,] <- c(-4 + temp_inc, -1 + temp_inc, 11 + temp_inc, 14 + temp_inc)
  temp_inc <- temp_inc + 0.01
}
```

## Acknowledgements
Many thanks to Gustav Delius for help with therMizer's code and to Romain Forestier for turning the therMizer code into this handy package.

## References
Boettiger C, Lang DT, Wainwright PC. (2012). rfishbase: exploring, manipulating, and visualizaing FishBase from R. Journal of Fish Biology 81, 2030–2039. https://doi.org/10.1111/j.1095-8649.2012.03464.x

Scott F, Blanchard JL, Andersen KH. (2014) mizer: an R package for
multispecies, trait-based and community size spectrum ecological modelling.
Methods in Ecology and Evolution, 5(10): 1121–1125. doi: 10.1111/2041-210X.12256

For more on what's happening behind the scenes in therMizer, check out the [Temperature-dependent rates in mizer](https://blog.mizer.sizespectrum.org/posts/2022-07-11-thermizer/) blog post.  This post also includes references for the science behind the package.
