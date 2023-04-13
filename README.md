# Introduction to therMizer

The therMizer package is an extension of the mizer package ([Scott et al. 2014](https://doi.org/10.1111/2041-210X.12256)) that allows you to incorporate the effects of temperature on species' metabolic rate and aerobic scope. These effects can vary by individual body size.  The therMizer package also allows you to use a dynamic resource spectrum.

## Installing therMizer

The `remotes` package is needed to install packages hosted on GitHub.

``` r
install.packages("remotes")

remotes::install_github("sizespectrum/therMizer")
```

Finally, load the newly installed package with

``` r
library(therMizer)
```

## Conceptualizing therMizer

Temperature affects species' metabolic rates and aerobic scope. Yet we lack much species-specific understanding about the relationship between these rates and temperature. The therMizer package uses approximate relationships based on species' thermal tolerance limits.

The effects of temperature on metabolism follow a Boltzman factor or Arrhenius relationship. At the low end of species' thermal tolerance limits, metabolism is the least expensive and more energy can be devoted to growth an reproduction. At the high end of species thermal tolerance limits metabolism is maximally expensive and less energy is available for growth and reproduction.

The effects of temperature on aerobic scope are approximated by a generic polynomial rate equation that results in a thermal optimum for each species. At that temperature, species are able to realize peak aerobic performance and encounter the maximum amount of prey possible. This falls away at temperatures to either side of the optimum, with reduced predator-prey encounter rates.

The therMizer package attempts to capture these effects based on temperatures to which species are exposed. This exposure can vary by depth, for example for species that undergo diel vertical migration, as well as by size, for example species that undergo ontogenetic vertical migration.

Temperature exposure in therMizer is modeled much the same way fishing mortality is modeled in mizer. `ocean_temp` is analgous to `effort`. It's assigned to a `realm` the same way effort is assigned to a `gear`. `exposure` links species and realms much the same way `catchability` links species and gears. `vertical_migration` brings in body size the same way `selectivity` does. The following section walks through these parameters in more detail.

The therMizer package also allows users to supply a time-varying background resource. If supplied, this is used in place of the default semi-chemostat resource. This option allows users to simulate changes to the plankton community as a driver of food web change, much the same way time-varying fishing can be used. For example, you could use output from an earth system model to inform a size-structured dynamic resource.

## Model parameters and input

This section looks in detail at the parameters and input needed to run therMizer. There's also some sample code further below for preparing these parameters for use in therMizer.

### Species parameters

In addition to the species parameters required by mizer, you'll also need to supply `temp_min` and `temp_max` parameters for each species. These represent the lower and upper bounds of species' thermal tolerance limits. You can find this information in the literature, for example from tagging studies, or in databases such as [rfishbase](https://github.com/ropensci/rfishbase) ([Boettiger et al. 2012](https://doi.org/10.1111/j.1095-8649.2012.03464.x)).

### Temperature parameters

By default, therMizer provides values for the `realm`, `vertical_migration`, and `exposure` parameters that are appropriate for most use cases. This means that if you don't want to specify custom values for these parameters, you can simply omit them from your code, and therMizer will use the default values instead. This allows you to get started with modelling using just one temperature value at the minimum. Using default values can save you time and effort since you don't need to spend time researching or experimenting with different parameter values. However, if you choose to specify custom values for these parameters, therMizer provides flexibility to do so, allowing you to customize your model to fit your specific research question or study system.

`realm` refers to the vertical realms or depth strata which species inhabit and for which temperatures will be provided. These could be named something like *epipelagic*, *mesopelagic*, *surface*, *bottom*, whatever you want to call them. These realms should be set up as the second dimension (or columns) of `ocean_temp`: one column per realm representing different time series of temperature. By default, if only one vector of temperature is supplied, there will be only one realm.

`vertical_migration` simulates the proportion of time that a given size of a given species spends in a given realm. It has the dimensions of `realm` $\times$ `sp` $\times$ `w`. Values can range from 0 to 1, and must sum to 1 across all realms for each species and size. If values sum to something other than one, it means that either a portion of time is unaccounted for (\<1) or that time across realms is over-allocated (\>1). Either way, you'll be modeling something that cannot be replicated in the real ocean. By default, species are assumed to spend equal time across all realms.

`exposure` links `vertical_migration` to `ocean_temp`. It has the dimensions of `realm` $\times$ `sp`. The values are 1 for the realms to which a species is exposed and 0 elsewhere. In theory, you could set all values to 1 and, so long as `vertical_migration` is constructed correctly, get the same results (because when multiplied by `exposure` the result would still be 0 for realms in which species spend no time). It's up to you whether you'd like to go this route. However, you do need to ensure that the realm names and order match those used in `vertical_migration` and `ocean_temp`. By default, the exposure is set to 1 for all realms and species and therefore has no effects.

### Temperature functions

Temperature affects species within mizer by overwriting mizer's default rate functions and replacing them with custom functions using the new set of parameters. Two new functions `therMizerEncounter` and `therMizerPredRate` affect the encounter and predation rates and one function `therMizerEReproAndGrowth` takes care of the maintenance metabolism. These functions can be disabled by setting the arguments `aerobic_effect` and `metabolism_effect` to `FALSE` for encounter and predation rates and for metabolism, respectively.

These functions can also be overridden by the user using `setRateFunction()`. Example below:

```r

params <- setRateFunction(params,"Encounter","newEncounterFunction")

```

### Input

`ocean_temp` is an array that has temperature(s) in degrees Celsius for each `realm`. It can be a vector, if temperature is constant over time, or an array for dynamic temperatures. If you're using time-varying temperature, the array will have the dimensions of `time` $\times$ `realm`.

`n_pp` is an array that has numerical plankton abundance for each size class. therMizer will convert these abundances to densities for use within mizer. `n_pp` can be a vector, if these abundances are constant over time, or an array for a dynamic resource. If you're using time-varying plankton, the array will have the dimensions of `time` $\times$ `w`.  


## Sample code for preparing parameters and input

Below are examples on how to setup different variables and arrays to use with therMizer. It does not mean you need all of them to run a model. The bare minimum is a mizerParams object, the two temperature range vectors `temp_min` and `temp_max` and finally the temperature array.

Let's create some example species parameters for two fictional species:

```r
species_params = data.frame(species = c("speciesA", "speciesB"), w_inf = c(500, 5000), k_vb = c(0.8, 0.3), w_min = c(0.001, 0.001), w_mat = c(5, 50), beta = c(1000,100), sigma = c(3,3))
species_params$interaction_resource <- c(1,0.5)
params <- newMultispeciesParams(species_params, no_w = 200, kappa = 0.0001) |> 
    steady(tol = 0.001)

# Assign them thermal tolerance limits
species_params(params)$temp_min <- c(-5, 5)
species_params(params)$temp_max <- c(10, 20)
```

Here's an example of how to set up the `vertical_migration` array. We'll assume one species stays in the upper 50 m of the water column until it moves to the bottom at maturity and that all sizes of the other species undergo diel vertical migration (DVM). This will give us four realms.

``` r
realm_names <- c("upper50m","bottom","DVM_day","DVM_night")
species_names <- as.character(params@species_params$species)
sizes <- params@w

# Create the vertical migration array and fill it
vertical_migration_array <- array(0, dim = (c(length(realm_names), 
                                  length(species_names), length(sizes))), 
                                  dimnames = list(realm = realm_names, sp = species_names, 
                                  w = signif(sizes,3))) # realm x species x size

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


```

Using the same scenario, here's an example to set up the `exposure` array.

``` r
exposure_array <- array(0, dim = (c(length(realm_names), length(species_names))), 
                  dimnames = list(realm = realm_names, sp = species_names)) # realm x species

for (r in seq(1,length(realm_names),1)) {
    for (s in seq(1,length(species_names),1)) {
        if (any(vertical_migration_array[r,s,] > 0)) {
            exposure_array[r,s] = 1
        }
    }
}
```

An example for creating the temperatures for each realm.

``` r
# Create temperature array and fill it
times <- 0:500
ocean_temp_array <- array(NA, dim = c(length(times), length(realm_names)), 
                    dimnames = list(time = times, realm = realm_names))
temp_inc <- 0
for (i in 1:501) {
  ocean_temp_array[i,] <- c(-4 + temp_inc, -1 + temp_inc, 11 + temp_inc, 14 + temp_inc)
  temp_inc <- temp_inc + 0.01
}
```

An example for creating a dynamic resource spectra.

``` r
x <- params@w_full
slope <- -1
intercept <- -5

# Create resource array and fill it
n_pp_array <- array(NA, dim = c(length(times), length(x)), 
                    dimnames = list(time = times, w = signif(x,3)))

for (i in 1:501) {
  # Add some noise around the slope and intercept as we fill the array
  n_pp_array[i,] <- (slope * runif(1, min = 0.95, max = 1.05) * log10(x)) + 
                    (intercept * runif(1, min = 0.95, max = 1.05))
}

```

## Running a scenario

The `upgradeTherParams` function combines a standard `mizerParams` object with the therMizer objects described above. The only necessary parameters are `temp_min`, `temp_max` and `ocean_temp_array`.

```r
params <- upgradeTherParams(params, 
                            temp_min = temp_min,
                            temp_max = temp_max,
                            ocean_temp_array = ocean_temp_array,
                            n_pp_array = n_pp_array, 
                            vertical_migration_array = vertical_migration_array,
                            exposure_array = exposure_array, 
                            aerobic_effect = TRUE, 
                            metabolism_effect = TRUE)
                                
```
Note that `ocean_temp_array` defines the specific time frame for the projection in time that needs to occur in the `project()` function. In particular, the `t_start` argument must be within the range of dates covered by `ocean_temp_array`. Similarly, the length of the simulation is determined by `t_max` in years, which should be set within the range of dates of `ocean_temp_array`. If `t_max` exceeds the length of the temperature time series provided by `ocean_temp_array`, therMizer will cycle through the values of `ocean_temp_array` until reaching the end of the simulation.

Finally, it is important to note that mizer is based on a yearly time step. If `ocean_temp_array` does not follow yearly time steps, the `dt` argument in the `project()` function must be adjusted accordingly. The `dt` argument specifies the time step size and should be set to the fraction of a year corresponding to the time step used in `ocean_temp_array`. For example, if `ocean_temp_array` provides monthly temperatures, the `dt` value should be set to 1/12, while daily temperatures would require a `dt` of 1/365 or 1/366 for leap years. By default, `dt` is set to 1/10.

Below is an example to use the `project` function with therMizer:

```r

sim <- project(params, 
               # First date in ocean_temp_array
               t_start = as.numeric(dimnames(other_params(params)$ocean_temp)[[1]][1]),
               # Duration in Years
               t_max = 10,
               # Monthly dates
               dt = 1/12)

```

The `plotThermPerformance` function displays the shape of the thermal performance curves for each species.

```r

plotThermPerformance(params)

```

## Calibrating models with therMizer

Mizer functions used to calibrate models often rely on the `project()` function, which cannot be setup to follow a specific time-frame with `t_start` or `t_max` unless you are coding your own calibration functions. Since therMizer will cycle through the values of `ocean_temp_array` indefinitely, it is possible to set up a unique temperature value (or an array of 1 x realms) as `ocean_temp_array`, **without** including dates, to calibrate a model. That way, therMizer will be date independent. Be mindful that once the end of `ocean_temp_array` is reached, thermizer will pick the first temperature value of `ocean_temp_array` again to continue the simulation.

## Acknowledgements

Many thanks to Gustav Delius for help with therMizer's code and to Romain Forestier for turning the therMizer code into this handy package.

## References

Boettiger C, Lang DT, Wainwright PC. (2012). rfishbase: exploring, manipulating, and visualizaing FishBase from R. Journal of Fish Biology 81, 2030--2039. <https://doi.org/10.1111/j.1095-8649.2012.03464.x>

Scott F, Blanchard JL, Andersen KH. (2014) mizer: an R package for multispecies, trait-based and community size spectrum ecological modelling. Methods in Ecology and Evolution, 5(10): 1121--1125. doi: 10.1111/2041-210X.12256

For more on what's happening behind the scenes in therMizer, check out the [Temperature-dependent rates in mizer](https://blog.mizer.sizespectrum.org/posts/2022-07-11-thermizer/) blog post. This post also includes references for the science behind the package.
