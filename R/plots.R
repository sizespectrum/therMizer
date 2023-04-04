# scripts containing all plot functions of the package

# to appease R CMD checks
globalVariables(c("temperature", "scalar","Type"))

#' @title plot thermal performance.
#'
#' @description Take a look at the thermal performance of the species.
#'
#' @inheritParams upgradeTherParams
#'
#' @param return_data Boolean value allowing to return the data frame used for
#' the plot instead of the plot itself. Default is FALSE.
#' @param resolution Numeric value which determines the step-width between each
#' calculation of the species' thermal performance. Default is 0.2.
#'
#' @export

plotThermPerformance <- function(params, resolution = .2, return_data = FALSE){

temp_vec <- seq(min(params@species_params$temp_min), max(params@species_params$temp_max), resolution)
scalar <- NULL
scalarM <- NULL

for(iTemp in temp_vec){

  temp_at_t <- iTemp + 273

  # Encounter rate

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

  scalar <- rbind(scalar,scaled_temp_effect_r)

  # Metabolism

  # Arrhenius equation
  unscaled_temp_effect <- (exp(25.22 - (0.63/((8.62e-5)*(temp_at_t)))))

  # Arrhenius equation scaled to a value between 0 and 1
  temp_effect_metabolism <-
    (unscaled_temp_effect - species_params(params)$metab_min) /
    species_params(params)$metab_range

  # Set temperature effect to 0 if temperatures are outside thermal
  # tolerance limits
  above_max <- (temp_at_t - 273) > species_params(params)$temp_max
  below_min <- (temp_at_t - 273) < species_params(params)$temp_min
  temp_effect_metabolism[above_max | below_min] = 0

  scalarM <- rbind(scalarM, temp_effect_metabolism)

}
  # plot

  rownames(scalar) <- rownames(scalarM) <- temp_vec
  colnames(scalar) <- colnames(scalarM) <- params@species_params$species
  plot_dat <- reshape2::melt(scalar)
  plot_datM <- reshape2::melt(scalarM)
  colnames(plot_dat) <- colnames(plot_datM) <- c("temperature","Species","scalar")
  plot_dat$Type <- "Encounter"
  plot_datM$Type <- "Metabolism"

  plot_dat <- rbind(plot_dat, plot_datM)
  plot_dat$Type <- factor(plot_dat$Type, levels = c("Metabolism","Encounter"))

  # removing scalar = 0 except at the beginning and end of curve (for encounter)
  zero_keep <- NULL
  zero_val <- which(plot_dat$scalar == 0)
  for(i in 1:(length(zero_val)-1)){
    if(zero_val[i+1] != zero_val[i] + 1){
      if(plot_dat$Type[zero_val[i]] == "Encounter") zero_keep <- c(zero_keep, i, i+1)
      else zero_keep <- c(zero_keep, i)
    }
  }
  # zero keep has the position of the row to keep in zero_val
  zero_val <- zero_val[-zero_keep]
  plot_dat <-plot_dat[- zero_val,]

  p <- ggplot(plot_dat)+
    geom_line(aes(x = temperature, y = scalar, color = Type)) +
    scale_x_continuous("Temperature in C")+
    scale_y_continuous("Scalar value") +
    facet_wrap(~Species, scales = "free")

  if(return_data) return(plot_dat) else return(p)
}




