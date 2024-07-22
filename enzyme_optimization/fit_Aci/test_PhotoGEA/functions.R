# Specify some common settings used for A-Ci curve fitting
ELECTRONS_PER_CARBOXYLATION <- 4.5
ELECTRONS_PER_OXYGENATION   <- 5.25
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

Vcmax_multiplier<-function(T_kelvin){
  #Eq. 10 in,
  #Scafaro, A.P., Posch, B.C., Evans, J.R. et al. Rubisco deactivation and chloroplast electron transport rates co-limit photosynthesis above optimal leaf temperature in terrestrial plants. Nat Commun 14, 2820 (2023). https://doi.org/10.1038/s41467-023-38496-4
  Tgrowth = 24
  Ha = 82992-632*Tgrowth
  Tref = 25 + 273.15 #kelvin
  R = 8.314 #gas constant J K-1 mol-1
  deltaS = 668.39 - 1.07 * Tgrowth  #J mol-1
  Hd = 200e3 #J mol-1
  
  term1 = exp(Ha*(T_kelvin - Tref)/(Tref*R*T_kelvin))
  term2 = 1 + exp(Tref*(deltaS - Hd)/(Tref*R))
  term3 = 1 + exp((T_kelvin*deltaS - Hd)/(T_kelvin*R))
  
  multiplier = term1 * term2 / term3
  
  return(multiplier)
}

get_vcmax_jmax <- function(A_Ci_df)
{
  aci_data      = A_Ci_df
  aci_data_exdf = exdf(aci_data)

  aci_data_exdf$units$Tleaf = "degrees C"
  aci_data_exdf$units$Ci = "micromol mol^(-1)"
  aci_data_exdf$units$A = "micromol m^(-2) s^(-1)"
  
  #add a bogus constant Ca
  aci_data_exdf <- set_variable(
    aci_data_exdf,
    'Ca',
    value = 1,
    units = 'micromol mol^(-1)'
  )
  #assume 1 atmos pressure
  aci_data_exdf <- set_variable(
    aci_data_exdf,
    'total_pressure',
    value = 1,
    units = 'bar'
  )
  aci_data_exdf <- set_variable(
    aci_data_exdf,
    'oxygen',
    value = 21,
    units = 'percent'
  )
  #calculate Kc, J_norm, Ko,Vcmax
  aci_data_exdf <- calculate_arrhenius(aci_data_exdf, c3_arrhenius_bernacchi,tleaf_column_name = 'Tleaf')
  #overwrite arrehenius with a new vcmax temperature response
  aci_data_exdf[, 'Vcmax_norm'] <- Vcmax_multiplier(aci_data_exdf[, 'Tleaf']+273.15)

  # Calculate Cc
  aci_data_exdf <- set_variable(
    aci_data_exdf,
    'gmc',
    'mol m^(-2) s^(-1) bar^(-1)',
    'process_soybean_aci',
    Inf
  )
  aci_data_exdf <- apply_gm(aci_data_exdf)
  
  # We can fit just one curve from the data set, although it is rare to do this
  MY_FIT_OPTIONS <-list(Rd_at_25  = soybean$parameters$Rd,
                        Tp        = 10.5,#fixed TPU estimated from LD11
                        alpha_old = 0 ) 

  c3_aci_results <- fit_c3_aci(
    aci_data_exdf,
    a_column_name = 'A',
    Ca_atmospheric = 420,
    fit_options = MY_FIT_OPTIONS,
    atp_use = ELECTRONS_PER_CARBOXYLATION,
    nadph_use = ELECTRONS_PER_OXYGENATION * 2
  )
  
  #here, no average is needed
  c3_aci_averages <- c3_aci_results$parameters
  
  # Compile table of parameter values to use with BioCro
  soybean_ld11_fvcb_parameters <- list(
    cultivar = 'ld11',
    electrons_per_carboxylation = ELECTRONS_PER_CARBOXYLATION,
    electrons_per_oxygenation = ELECTRONS_PER_OXYGENATION,
    Vcmax = c3_aci_averages[, 'Vcmax_at_25'],
    J     = c3_aci_averages[, 'J_at_25'],
    Rd    = c3_aci_averages[, 'Rd_at_25'],
    TPU   = c3_aci_averages[, 'Tp']
  )
  
  # One complication is that PhotoGEA returns a value for J, but not Jmax. In
  # BioCro, values of PhiPSII, Q, and Jmax are used to determine J. Here, we will
  # solve those equations for Jmax.
  leaf_reflectance   = 0.1
  leaf_transmittance = 0.05
  soybean_ld11_fvcb_parameters$Jmax <- get_jmax(
    soybean$parameters$theta,              # dimensionless
    soybean$parameters$beta_PSII,          # dimensionless
    soybean_ld11_fvcb_parameters$J,        # micromol / m^2 / s
    leaf_reflectance,   # dimensionless
    25,                                    # degrees C
    leaf_transmittance, # dimensionless
    mean(aci_data_exdf[, 'PAR'])              # micromol / m^2 / s
  )
  return(c(soybean_ld11_fvcb_parameters$Vcmax,soybean_ld11_fvcb_parameters$Jmax,soybean_ld11_fvcb_parameters$TPU))
}
