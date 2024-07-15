arrhenius_exponential <- function(scaling, activation_energy, temperature_k) {
    ideal_gas_constant <- 8.31446261815324
    exp(scaling - activation_energy / (ideal_gas_constant * temperature_k))
}

Vcmax_new<-function(Vcmax25,T_kelvin){
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
  
  Vcmax = Vcmax25 * term1 * term2 / term3
  
  return(Vcmax)
}

BioCro_FvCB <- function(Qin,Tleaf,Ci,Vcmax,Jmax,Rd,TPU){
  # Make an alternate BioCro comparison; here we use BioCro:FvCB, so we have to
  # manually calculate the temperature response of several parameters using code
  # copied from c3photoC and slightly modified for R
  leaf_reflectance   = 0.1
  leaf_transmittance = 0.05

  ideal_gas_constant <- 8.31446261815324
  
  Tleaf_K <- Tleaf + 273.15
  
  TPU_c = 25.5
  Ha    = 62.99e3
  S     = 0.588e3
  Hd    = 182.14e3
  R     = ideal_gas_constant
  top = Tleaf_K * arrhenius_exponential(TPU_c, Ha, Tleaf_K)
  bot = 1.0 + arrhenius_exponential(S / R, Hd, Tleaf_K)
  TPU_rate_multiplier = (top / bot) / 306.742
  
  Jmax <- Jmax_at_25 * arrhenius_exponential(17.57, 43.54e3, Tleaf_K) # value at 25 C estimated from 2022
  
  theta <- soybean$parameters$theta + 0.018 * Tleaf - 3.7e-4 * Tleaf^2
  
  dark_adapted_phi_PSII <- 0.352 + 0.022 * Tleaf - 3.4 * Tleaf^2 / 1e4
  
  absorbed_ppfd <- Qin * (1 - leaf_transmittance - leaf_reflectance)
  
  I2 <- absorbed_ppfd * dark_adapted_phi_PSII * soybean$parameters$beta_PSII
  
  J <- (Jmax + I2 - sqrt((Jmax + I2)^2 - 4.0 * theta * I2 * Jmax)) / (2.0 * theta)
  
  rc <- module_response_curve(
      'BioCro:FvCB',
      within(soybean$parameters, {
          alpha_TPU = 0
          Oi = O2
          rm(Rd)
      }),
      data.frame(
          Ci = Ci,
          Gstar = arrhenius_exponential(19.02, 37.83e3, Tleaf_K),
          Kc = arrhenius_exponential(38.05, 79.43e3, Tleaf_K),
          Ko = arrhenius_exponential(20.30, 36.38e3, Tleaf_K),
          Rd = Rd_at_25 * arrhenius_exponential(18.72, 46.39e3, Tleaf_K),
          # Vcmax = Vcmax_new(Vcmax_at_25 , Tleaf_K),
          Vcmax = Vcmax_at_25 * arrhenius_exponential(26.35, 65.33e3, Tleaf_K),
          TPU = TPU * TPU_rate_multiplier,
          J = J
      )
  )

  return(rc)
}
