# Load required packages
library(PhotoGEA)
library(lattice)
library(BioCro)
aci_fitting_2021<-function(scaling_factors){
  # Load default A-Ci settings
  source('aci_defaults.R')
  ###
  ### TRANSLATION:
  ### Creating convenient R objects from raw data files
  ###
  
  # Define a vector of paths to the files we wish to load
  file_paths <- c(
      'ACi_AQ_data/licor_raw/2021-08-04 ed ripe1.xlsx',
      'ACi_AQ_data/licor_raw/2021-08-04 ed ripe2.xlsx',
      'ACi_AQ_data/licor_raw/2021-08-04 ed ripe3.xlsx',
      'ACi_AQ_data/licor_raw/2021-08-04 ed ripe4.xlsx'
  )
  
  # Load each file, storing the result in a list
  licor_exdf_list <- lapply(file_paths, function(fpath) {
    read_gasex_file(fpath, 'time')
  })
  
  # Get the names of all columns that are present in all of the Licor files
  columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)
  
  # Extract just these columns
  licor_exdf_list <- lapply(licor_exdf_list, function(x) {
    x[ , columns_to_keep, TRUE]
  })
  
  # Use `rbind` to combine all the data
  licor_data <- do.call(rbind, licor_exdf_list)
  
  # Keep just the soybean data
  licor_data <- licor_data[licor_data[, 'species'] == 'soybean', , TRUE]
  
  # Add columns for year and day of year
  licor_data[, 'year'] <- as.numeric(format(licor_data[, 'time'], '%Y'))
  licor_data[, 'doy'] <- as.numeric(format(licor_data[, 'time'], '%j'))
  
  ###
  ### VALIDATION:
  ### Organizing the data, checking its consistency and quality, cleaning it
  ###
  
  # Create a new identifier column formatted like `instrument - species - plot`
  licor_data[ , 'curve_identifier'] <-
    paste(licor_data[ , 'instrument'], '-', licor_data[ , 'plot'])

  # Make sure the data meets basic requirements
  check_licor_data(licor_data, 'curve_identifier', 16, 'CO2_r_sp')
  
  # Remove points with duplicated `CO2_r_sp` values and order by `Ci`
  licor_data <- organize_response_curve_data(
      licor_data,
      'curve_identifier',
      c(9, 10),
      'Ci'
  )
  
  # Remove curve 'ripe3 - 4' because humidity was not properly controlled; this
  # curve cannot be used.
  licor_data <- remove_points(licor_data, list(curve_identifier = 'ripe3 - 4'))
  
  ###
  ### PROCESSING:
  ### Extracting new pieces of information from the data
  ###
  
  # Calculate total pressure in the Licor chamber
  licor_data <- calculate_total_pressure(licor_data)
  
  # BioCro's C3 photosynthesis module does not consider mesophyl conductance, so
  # we will set it to Inf (which will ensure Cc = Ci).
  licor_data <- set_variable(
      licor_data,
      'gmc',
      'mol m^(-2) s^(-1) bar^(-1)',
      'process_soybean_aci',
      Inf
  )
  
  # Calculate Cc
  licor_data <- apply_gm(licor_data)
  
  # Calculate temperature-dependent values of C3 photosynthetic parameters. Here
  # we use temperature response parameters from Bernacchi because BioCro has these
  # values hard-coded into its C3 photosynthesis code.
  licor_data <- calculate_arrhenius(licor_data, c3_arrhenius_bernacchi)
  
  # Fit the C3 A-Ci curves.
  c3_aci_results <- consolidate(by(
    licor_data,                       # The `exdf` object containing the curves
    licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
    fit_c3_aci,                       # The function to apply to each chunk of `licor_data`
    Ca_atmospheric = 420,
    fit_options = SOYBEAN_FIT_OPTIONS,
    atp_use   = ELECTRONS_PER_CARBOXYLATION,
    nadph_use = ELECTRONS_PER_OXYGENATION * 2
  ))
  
  # Compute the average and standard error of each parameter for each species
  c3_aci_averages <- basic_stats(c3_aci_results$parameters, 'species')
  
  if (!'Rd_at_25_avg' %in% colnames(c3_aci_averages)) {
      # Rd_at_25 was fixed
      c3_aci_averages[, 'Rd_at_25_avg']    <- SOYBEAN_FIT_OPTIONS$Rd
      c3_aci_averages[, 'Rd_at_25_stderr'] <- NA
  }
  
  # Compile table of parameter values to use with BioCro
  soybean_ld11_fvcb_parameters_2021 <- list(
      year = unique(licor_data[, 'year']),
      cultivar = 'ld11',
      electrons_per_carboxylation = ELECTRONS_PER_CARBOXYLATION,
      electrons_per_oxygenation = ELECTRONS_PER_OXYGENATION,
      Vcmax = c3_aci_averages[, 'Vcmax_at_25_avg'],
      J = c3_aci_averages[, 'J_at_25_avg'],
      Rd = c3_aci_averages[, 'Rd_at_25_avg'],
      TPU = c3_aci_averages[, 'Tp_avg']
  )
  
  # One complication is that PhotoGEA returns a value for J, but not Jmax. In
  # BioCro, values of PhiPSII, Q, and Jmax are used to determine J. Here, we will
  # solve those equations for Jmax.
  leaf_reflectance   = 0.1
  leaf_transmittance = 0.05
  soybean_ld11_fvcb_parameters_2021$Jmax <- get_jmax(
      soybean$parameters$theta,              # dimensionless
      soybean$parameters$beta_PSII,          # dimensionless
      soybean_ld11_fvcb_parameters_2021$J,   # micromol / m^2 / s
      leaf_reflectance,   # dimensionless
      25,                                    # degrees C
      leaf_transmittance, # dimensionless
      mean(licor_data[, 'Qin']),             # micromol / m^2 / s
      scaling_factors
  )
  
  return(as.data.frame(soybean_ld11_fvcb_parameters_2021))
} #end function aci_fitting_2021

aci_fitting_2022<-function(scaling_factors){
  # Load default A-Ci settings
  source('aci_defaults.R')
  
  # Make some decisions about what to do
  REMOVE_UNSTABLE_POINTS <- FALSE
  REMOVE_BAD_CURVES <- TRUE
  
  ###
  ### TRANSLATION:
  ### Creating convenient R objects from raw data files
  ###
  
  # Define a vector of paths to the files we wish to load
  file_paths <- c(
      "ACi_AQ_data/licor_raw/2022-07-21 mcgrath1.xlsx",
      "ACi_AQ_data/licor_raw/2022-07-21 mcgrath2.xlsx",
      "ACi_AQ_data/licor_raw/2022-08-03 ed aci ripe2.xlsx",
      "ACi_AQ_data/licor_raw/2022-08-03 ed aci ripe3.xlsx",
      "ACi_AQ_data/licor_raw/2022-08-03 ed aci ripe14.xlsx",
      "ACi_AQ_data/licor_raw/2022-08-03 ed aci ripe15.xlsx",
      "ACi_AQ_data/licor_raw/2022-08-05 ripe2 aci.xlsx",
      "ACi_AQ_data/licor_raw/2022-08-05 ripe3 aci.xlsx",
      "ACi_AQ_data/licor_raw/2022-08-05 ripe14 aci.xlsx",
      "ACi_AQ_data/licor_raw/2022-08-05 ripe15 aci.xlsx"
  )
  
  # Load each file, storing the result in a list
  licor_exdf_list <- lapply(file_paths, function(fpath) {
    read_gasex_file(fpath, 'time')
  })
  
  # Get the names of all columns that are present in all of the Licor files
  columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)
  
  # Extract just these columns
  licor_exdf_list <- lapply(licor_exdf_list, function(x) {
    x[ , columns_to_keep, TRUE]
  })
  
  # Use `rbind` to combine all the data
  licor_data <- do.call(rbind, licor_exdf_list)
  
  # Keep just the soybean data
  licor_data <- licor_data[licor_data[, 'species'] %in% c('ld10', 'ld11'), , TRUE]
  
  # Add columns for year and day of year
  licor_data[, 'year'] <- as.numeric(format(licor_data[, 'time'], '%Y'))
  licor_data[, 'doy'] <- as.numeric(format(licor_data[, 'time'], '%j'))
  
  ###
  ### VALIDATION:
  ### Organizing the data, checking its consistency and quality, cleaning it
  ###
  
  # Create a new identifier column formatted like `species - instrument - replicate`
  licor_data[ , 'curve_identifier'] <-
    paste(licor_data[ , 'species'], licor_data[ , 'instrument'], licor_data[ , 'replicate'], sep = ' - ')
  
  # Check the CO2_r_sp values from one curve to confirm that we have a 16 point
  # curve where points 9 and 10 repeat the initial setpoint value
  print(licor_data[licor_data[, 'curve_identifier'] == licor_data[1, 'curve_identifier'], 'CO2_r_sp'])
  
  # Make sure the data meets basic requirements
  check_licor_data(licor_data, 'curve_identifier', 16, 'CO2_r_sp')
  
  # Remove points with duplicated `CO2_r_sp` values and order by `Ci`
  licor_data <- organize_response_curve_data(
      licor_data,
      'curve_identifier',
      c(9, 10),
      'Ci'
  )
  
  # Some curves have several unstable points, so we may want to exclude those
  # points
  if (REMOVE_UNSTABLE_POINTS) {
      # Only keep points where stability was achieved
      licor_data <- licor_data[licor_data[, 'Stable'] == 2, , TRUE]
  
      # Remove any curves that have fewer than three remaining points
      npts <- by(licor_data, licor_data[, 'curve_identifier'], nrow)
      ids_to_keep <- names(npts[npts > 2])
      licor_data <- licor_data[licor_data[, 'curve_identifier'] %in% ids_to_keep, , TRUE]
  }
  
  # Some curves were not controlled properly, so we may want to exclude them
  if (REMOVE_BAD_CURVES) {
      licor_data <- remove_points(
          licor_data,
          list(
              curve_identifier = 'ld11 - ripe15 - 4' # temperature and humidity controls failed
          )
      )
  }
  
  ###
  ### PROCESSING:
  ### Extracting new pieces of information from the data
  ###
  
  # Calculate total pressure in the Licor chamber
  licor_data <- calculate_total_pressure(licor_data)
  
  # BioCro's C3 photosynthesis module does not consider mesophyl conductance, so
  # we will set it to Inf (which will ensure Cc = Ci).
  licor_data <- set_variable(
      licor_data,
      'gmc',
      'mol m^(-2) s^(-1) bar^(-1)',
      'process_soybean_aci',
      Inf
  )
  
  # Calculate Cc
  licor_data <- apply_gm(licor_data)
  
  # Calculate temperature-dependent values of C3 photosynthetic parameters. Here
  # we use temperature response parameters from Bernacchi because BioCro has these
  # values hard-coded into its C3 photosynthesis code.
  licor_data <- calculate_arrhenius(licor_data, c3_arrhenius_bernacchi)
  
  # Fit the C3 A-Ci curves
  c3_aci_results <- consolidate(by(
    licor_data,                       # The `exdf` object containing the curves
    licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
    fit_c3_aci,                       # The function to apply to each chunk of `licor_data`
    Ca_atmospheric = 420,
    fit_options = SOYBEAN_FIT_OPTIONS,
    atp_use   = ELECTRONS_PER_CARBOXYLATION,
    nadph_use = ELECTRONS_PER_OXYGENATION * 2
  ))
  
  # Compute the average and standard error of each parameter for each species
  c3_aci_averages <- basic_stats(c3_aci_results$parameters, 'species')
  
  if (!'Rd_at_25_avg' %in% colnames(c3_aci_averages)) {
      # Rd_at_25 was fixed to the standard soybean value
      c3_aci_averages[, 'Rd_at_25_avg']    <- SOYBEAN_FIT_OPTIONS$Rd
      c3_aci_averages[, 'Rd_at_25_stderr'] <- NA
  }
  
  soybean_fvcb_parameters <- list(
    year = unique(licor_data[, 'year']),
    cultivar = c3_aci_averages[, 'species'],
    electrons_per_carboxylation = ELECTRONS_PER_CARBOXYLATION,
    electrons_per_oxygenation = ELECTRONS_PER_OXYGENATION,
    Vcmax = c3_aci_averages[, 'Vcmax_at_25_avg'],
    J = c3_aci_averages[, 'J_at_25_avg'],
    Rd = c3_aci_averages[, 'Rd_at_25_avg'],
    TPU = c3_aci_averages[, 'Tp_avg']
  )
  
  # One complication is that PhotoGEA returns a value for J, but not Jmax. In
  # BioCro, values of PhiPSII, Q, and Jmax are used to determine J. Here, we will
  # solve those equations for Jmax.
  leaf_reflectance = 0.1
  leaf_transmittance = 0.05
  soybean_fvcb_parameters$Jmax <- get_jmax(
      soybean$parameters$theta,              # dimensionless
      soybean$parameters$beta_PSII,          # dimensionless
      soybean_fvcb_parameters$J,             # micromol / m^2 / s
      leaf_reflectance,   # dimensionless
      25,                                    # degrees C
      leaf_transmittance, # dimensionless
      mean(licor_data[, 'Qin']),             # micromol / m^2 / s
      scaling_factors
  )
  
  return(as.data.frame(soybean_fvcb_parameters))
} #end function aci_fitting_2022 

fit_AQ<-function(Vcmax, Jmax,TPU,year){
  library(nloptr)
  source("FarquharModel.R")
  #PhiPSII_func<-function(LeafTemperature,sf) {
  #  return(0.352 *sf[1]+ 0.022 * sf[2]* LeafTemperature - sf[3]* 3.4*1e-4*LeafTemperature^2.0)
  #}
  #theta_func<-function(LeafTemperature) {
  #  PhotosynthesisTheta = 0.76
  #  return(PhotosynthesisTheta + 0.01713 * LeafTemperature - 3.75 *LeafTemperature^2.0 / 10000.0)
  #}
  rmse<-function(model, obs){
    return(sqrt(sum((model-obs)^2)))
  }
  R=8.314472E-3 #Gas constant KJ mole^{-1} K^{-1}
  #read in obs An data
  obs_An  = read.csv(paste0('ACi_AQ_data/AQ_DOY217_',year,'.csv'))
  
  Rd25   = 1.28 
  Air_O2 = 210  #
  
  Vcmax25  = Vcmax   #LD11
  Jmax25   = Jmax    #LD11
  Rate_TPu = TPU     #LD11
  
  init_guess = c(1,1)
  best_par = c()
  #here we fit to each replicate separately
  for (curve_id in unique(obs_An$curve_identifier)){
    obs_An_sub = obs_An[obs_An$curve_identifier%in%curve_id,]
    Cis = obs_An_sub$Ci
    Qps = obs_An_sub$Qin
    Tls = obs_An_sub$Tleaf
  #define obj function here  
    obj_func<-function(x){
      An_farquhar = NA*(1:length(Qps))
      for (i in 1:length(Qps)){
        Ci = Cis[i]
        Qp = Qps[i]
        LeafTemperature   = Tls[i] #C
        Rd = Rd25 * exp(18.72 - 46.39 / (R * (LeafTemperature + 273.15)))
        
        output_farquhar  = FarquharModel(LeafTemperature, Ci, Qp, Air_O2,Vcmax25, Jmax25, Rate_TPu,x)
        
        if(length(output_farquhar)>1){
          An_farquhar[i]   = output_farquhar$GA - Rd;
        }
      }
      if(any(is.na(An_farquhar))){
         return(9999)
      }else{
        return(rmse(An_farquhar,obs_An_sub$A))
      }
    }
    # opt_results = hjk(init_guess, obj_func)
    # best_par = rbind(best_par,c(opt_results$par,opt_results$value))
    opt_results = nloptr(x0 = init_guess, eval_f = obj_func, 
                         lb = c(0.5,0.5),
                         ub = c(2.0,2.0),
                         opts = list("algorithm" = "NLOPT_LN_SBPLX","xtol_rel"=1.0e-6))
    best_par = rbind(best_par,c(opt_results$solution,opt_results$objective))
  }
  best_par_mean = colMeans(best_par)
  return(best_par_mean) #par1, par2, rmse
} # end function fit_AQ
