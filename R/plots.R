# scripts containing all plot functions of the package

# to appease R CMD checks
globalVariables(c("temperature", "scalar","Type"))

# common theme for the package
therTheme <- function(){
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_line(color = "grey"),
    strip.background = element_blank(),
    legend.key = element_blank()
  )
}

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

plotTherPerformance <- function(params, resolution = .2, return_data = FALSE){

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
    facet_wrap(~Species, scales = "free") +
    therTheme()

  if(return_data) return(plot_dat) else return(p)
}


#' @title plotTherScalar
#'
#' @description Plot the scalar value affecting the encounter rate and metabolsim
#' of each species throughout the provided temperature in ocean_temp_array.
#'
#' @inheritParams plotTherPerformance
#'
#' @param species A character string. Select of specific species to display. It
#' has to correspond to one of the species name in the mizerParams object.
#' Default is NULL.
#' @param species_panel Boolean value. If set to TRUE, the plot will be a panel
#' of the species. Disabled if the argument species is not NULL. Default is TRUE.
#'
#' @export

plotTherScalar <- function(params, species = NULL, species_panel = TRUE, return_data = FALSE){

  params <- validParams(params)
  scalar_en <- NULL
  scalar_met <- NULL
  for(t in as.numeric(dimnames(other_params(params)$ocean_temp)[[1]])){
    # encounter
    scalar_en <- scaled_temp_effect(params,t) %>%
      apply(1,mean) %>%
      rbind(scalar_en,.)

    # Metabolism
    temp_effect_metab_realms <- array(NA, dim = c(dim(other_params(params)$vertical_migration)), dimnames = c(dimnames(other_params(params)$vertical_migration)))

    # Using t+1 to avoid calling ocean_temp[0,] at the first time step
    # Looping through each realm
    nb_realms <- dim(other_params(params)$exposure)[1]
    for (r in seq(1, nb_realms, 1)) {
      index <- which.min(abs(as.numeric(dimnames(other_params(params)$ocean_temp)[[1]]) - t))
      temp_at_t <- other_params(params)$ocean_temp[index,r]
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
    scalar_met <- colSums(temp_effect_metab_realms) %>%
      apply(1,mean) %>%
      rbind(scalar_met,.)
  }
  rownames(scalar_en) <- rownames(scalar_met) <-as.numeric(dimnames(other_params(params)$ocean_temp)[[1]])

  plot_dat2 <- reshape2::melt(scalar_met) %>%
    mutate(Type = "Metabolism")
  plot_dat <- reshape2::melt(scalar_en) %>%
    mutate(Type = "Encounter") %>%
    rbind(plot_dat2) %>%
    rename(Time = Var1, Species = Var2, Scalar = value)

  if(!is.null(species)) plot_dat <- filter(plot_dat, Species == species) %>%
    mutate(Species = as.character(Species))

  p <- ggplot(plot_dat, aes(x = Time, y = Scalar)) +
    scale_y_continuous(limits = c(0,1), name = "Scalar value") +
    scale_x_continuous(name = "Time in years")+
    therTheme()

  if(is.null(species) & species_panel){
    p <- p +
      geom_line(aes(color = Type))+
      facet_wrap(~Species) +
      scale_color_manual(values = c("#619CFF","#F8766D"))
  } else {
    p <- p +
      geom_line(aes(color = Species, linetype = Type)) +
      scale_color_manual(values = params@linecolour[params@species_params$species])
  }
  if(return_data) return(plot_dat) else return(p)
}
