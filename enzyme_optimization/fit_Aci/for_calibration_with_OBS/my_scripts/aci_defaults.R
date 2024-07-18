# Specify some common settings used for A-Ci curve fitting
ELECTRONS_PER_CARBOXYLATION <- 4.5
ELECTRONS_PER_OXYGENATION   <- 5.25

FIX_SOYBEAN_RD <- TRUE

SOYBEAN_FIT_OPTIONS <- if (FIX_SOYBEAN_RD) {
    list(Rd_at_25 = soybean$parameters$Rd,
         alpha_old = 0)
} else {
    list()
}

# Helping function for determining Jmax from J, following equations used in
# BioCro.
get_jmax <- function(
    base_theta,               # dimensionless
    beta_PSII,                # dimensionless
    J,                        # micromol / m^2 / s
    leaf_reflectance,         # dimensionless
    leaf_temperature_celsius, # degrees C
    leaf_transmittance,       # dimensionless
    Qp                        # micromol / m^2 / s
)
{
    # Apply temperature response equations
    dark_adapted_phi_PSII <-
        0.352 + 0.022 * leaf_temperature_celsius -
            3.4 * leaf_temperature_celsius^2 / 10000 # dimensionless

    theta <-
        base_theta + 0.018 * leaf_temperature_celsius -
            3.7e-4 * leaf_temperature_celsius^2 # dimensionless

    # Absorbed light
    Qabs <- Qp * (1.0 - leaf_reflectance - leaf_transmittance)

    # Find useful energy sent to photosystem II
    I2 <- Qabs * dark_adapted_phi_PSII * beta_PSII # micromol / m^2 / s

    # Calculate and return Jmax
    Jmax <- (J * I2 - theta * J^2) / (I2 - J) # micromol / m^2 / s
    
    return(Jmax)
}

#since BioCro's FvCB module has a temperature dependence on the TPU,
#here we use the same temperature function to re-scale the TPU estimated
#by PhotoGEA. Then the FvCB and PhotoGEA will be consistent in TPU
TPU_multiplier <- function(LeafTemperature){
  R     = 8.31446261815324e-3
  TPU_c = 25.5   #BioCro R
  Ha    = 62.99  #BioCro R
  S     = 0.588  #BioCro R
  Hd    = 182.14 #BioCro R
  TPU_rate_sf25 = 306.742 #nomalizer at 25 oC, this is what's used in BioCro

  LeafTemperatureKelvin = LeafTemperature + 273.15  #Leaf temperature in K
  top_term = LeafTemperatureKelvin * exp(TPU_c-Ha/(R*LeafTemperatureKelvin))
  bot_term = 1+ exp((S*LeafTemperatureKelvin - Hd)/(R*LeafTemperatureKelvin))
  TPU_rate_sf = top_term /bot_term / TPU_rate_sf25
  return(TPU_rate_sf)
}
get_TPU_at_25 <- function(TPU,LeafTemperature){
   multiplier  <- TPU_multiplier(LeafTemperature)
   #this is reverse, so we divide the multiplier
   return (as.numeric(TPU)/as.numeric(multiplier))
}
