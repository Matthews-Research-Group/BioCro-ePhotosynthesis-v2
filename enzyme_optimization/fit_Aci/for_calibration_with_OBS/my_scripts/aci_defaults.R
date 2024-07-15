# Specify some common settings used for A-Ci curve fitting
ELECTRONS_PER_CARBOXYLATION <- 4.5
ELECTRONS_PER_OXYGENATION   <- 5.25

FIX_SOYBEAN_RD <- TRUE

SOYBEAN_FIT_OPTIONS <- if (FIX_SOYBEAN_RD) {
    list(Rd_at_25 = soybean$parameters$Rd)
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
