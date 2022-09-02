species_params = data.frame(species = c("speciesA", "speciesB"), w_inf = c(500, 5000), k_vb = c(0.8, 0.3), w_min = c(0.001, 0.001), w_mat = c(5, 50), beta = c(1000,100), sigma = c(3,3))
species_params$interaction_resource <- c(1,0.5)
params <- newMultispeciesParams(species_params, no_w = 200, kappa = 0.0001) |> 
  steady(tol = 0.001)

species_params(params)$temp_min <- c(-5, 5)
species_params(params)$temp_max <- c(10, 20)

# Create temperature array and fill it
times <- 0:500
ocean_temp_array <- array(NA, dim = c(length(times), length(realm_names)), dimnames = list(time = times, realm = realm_names))
temp_inc <- 0
for (i in 1:501) {
  ocean_temp_array[i,] <- c(-5 + temp_inc, -5 + temp_inc, -5 + temp_inc, -5 + temp_inc)
  temp_inc <- temp_inc + 0.1
}

other_params(params)$ocean_temp <- ocean_temp_array

params <- setRateFunction(params, "Encounter", "therMizerEncounter")
params <- setRateFunction(params, "PredRate", "therMizerPredRate")
params <- setRateFunction(params, "EReproAndGrowth", "therMizerEReproAndGrowth")
