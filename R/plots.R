# scripts containing all plot functions of the package

#' @title plot thermal performance.
#'
#' @description Take a look at the thermal performance of the species.
#'
#' @inheritParams upgradeTherParams
#'
#' @export
#'

plotThermPerformance <- function(params, return_data = FALSE){

  # need temeperature spanning the temp range of all species
  # need run scaled temp with their params
  # plot

temp_vec <- seq(min(params@species_params$temp_min), max(params@species_params$temp_max), .2)
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

  p <- ggplot(plot_dat)+
    geom_line(aes(x = temperature, y = scalar, color = Type)) +
    scale_x_continuous("Temperature in C")+
    scale_y_continuous("Scalar value") +
    facet_wrap(~Species, scales = "free")

  if(return_data) return(plot_dat) else return(p)
}




