# This script is based on the example in the 'Analyzing C3 A-Ci Curves'
# vignette from the PhotoGEA package, available online at
# https://eloch216.github.io/PhotoGEA/articles/analyzing_c3_aci_curves.html

###
### PRELIMINARIES:
### Loading packages, defining constants, creating helping functions, etc.
###

# Clear the workspace
rm(list=ls())

# Load required packages
library(PhotoGEA)
library(lattice)
library(BioCro)

# Load default A-Ci settings
source(file.path('..', 'all_aci_curves', 'aci_defaults.R'))

# Make some decisions about what to do
MAKE_VALIDATION_PLOTS <- TRUE

###
### TRANSLATION:
### Creating convenient R objects from raw data files
###

# Define a vector of paths to the files we wish to load
file_paths <- c(
    '2021-08-04 ed ripe1.xlsx',
    '2021-08-04 ed ripe2.xlsx',
    '2021-08-04 ed ripe3.xlsx',
    '2021-08-04 ed ripe4.xlsx'
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

# Check the CO2_r_sp values from one curve to confirm that we have a 16 point
# curve where points 9 and 10 repeat the initial setpoint value
print(licor_data[licor_data[, 'curve_identifier'] == 'ripe2 - 1', 'CO2_r_sp'])

# Make sure the data meets basic requirements
check_licor_data(licor_data, 'curve_identifier', 16, 'CO2_r_sp')

# Remove points with duplicated `CO2_r_sp` values and order by `Ci`
licor_data <- organize_response_curve_data(
    licor_data,
    'curve_identifier',
    c(9, 10),
    'Ci'
)

if (MAKE_VALIDATION_PLOTS) {
    # Plot all A-Ci curves in the data set
    dev.new()
    print(xyplot(
      A ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste('Net CO2 assimilation rate [', licor_data$units$A, ']')
    ))

    # Make a plot to check humidity control
    dev.new()
    print(xyplot(
      RHcham + `Humidifier_%` + `Desiccant_%` ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      ylim = c(0, 100),
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']')
    ))

    # Make a plot to check temperature control
    dev.new()
    print(xyplot(
      TleafCnd + Txchg ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      ylim = c(25, 40),
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste0('Temperature (', licor_data$units$TleafCnd, ')')
    ))

    # Make a plot to check CO2 control
    dev.new()
    print(xyplot(
      CO2_s + CO2_r + CO2_r_sp ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste0('CO2 concentration (', licor_data$units$CO2_r, ')')
    ))

    # Make a plot to check stability criteria
    dev.new()
    print(xyplot(
      `A:OK` + `gsw:OK` + Stable ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']')
    ))
}

# Export curves
write.csv(licor_data, 'aci_soybean_ld11_2021.csv', row.names = FALSE)

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

# Plot the C3 A-Ci fits (including limiting rates)
dev.new()
print(xyplot(
  A + Ac + Aj + Ap + A_fit ~ Cc | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto.key = list(space = 'right'),
  grid = TRUE,
  xlab = paste('Chloroplast CO2 concentration [', c3_aci_results$fits$units$Ci, ']'),
  ylab = paste('Net CO2 assimilation rate [', c3_aci_results$fits$units$A, ']'),
  par.settings = list(
    superpose.line = list(col = multi_curve_colors()),
    superpose.symbol = list(col = multi_curve_colors(), pch = 16)
  )
))

# Plot the C3 A-Ci fits
dev.new()
print(xyplot(
  A + A_fit ~ Ci | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste('Intercellular CO2 concentration [', c3_aci_results$fits$units$Ci, ']'),
  ylab = paste('Net CO2 assimilation rate [', c3_aci_results$fits$units$A, ']')
))

# Plot the residuals
dev.new()
print(xyplot(
  A_residuals ~ Ci | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  grid = TRUE,
  xlab = paste('Intercellular CO2 concentration [', c3_aci_results$fits$units$Ci, ']'),
  ylab = paste('Assimilation rate residual (measured - fitted)\n[', c3_aci_results$fits$units$A, ']')
))

###
### SYNTHESIS:
### Using plots and statistics to help draw conclusions from the data
###

# Make bar and box plots of the parameter values
box_bar_plot_param <- list(
    list(
        Y = c3_aci_results$parameters[, 'Vcmax_at_25'],
        X = c3_aci_results$parameters[, 'species'],
        ylim = c(0, 140),
        ylab = paste('Vcmax at 25 degrees C [', c3_aci_results$parameters$units$Vcmax_at_25, ']')
    ),
    list(
        Y = c3_aci_results$parameters[, 'J_at_25'],
        X = c3_aci_results$parameters[, 'species'],
        ylim = c(0, 200),
        ylab = paste('J at 25 degrees C [', c3_aci_results$parameters$units$J_at_25, ']')
    ),
    list(
        Y = c3_aci_results$parameters[, 'Rd_at_25'],
        X = c3_aci_results$parameters[, 'species'],
        ylim = c(0, 1.5),
        ylab = paste('Rd at 25 degrees C [', c3_aci_results$parameters$units$Rd_at_25, ']')
    ),
    list(
        Y = c3_aci_results$parameters[, 'Tp'],
        X = c3_aci_results$parameters[, 'species'],
        ylim = c(0, 16),
        ylab = paste('Tp [', c3_aci_results$parameters$units$Tp, ']')
    )
)

invisible(lapply(box_bar_plot_param, function(x) {
  dev.new()
  print(do.call(bwplot_wrapper, x))

  dev.new()
  print(do.call(barchart_with_errorbars, x))
}))

# Compute the average and standard error of each parameter for each species
c3_aci_averages <- basic_stats(c3_aci_results$parameters, 'species')

if (!'Rd_at_25_avg' %in% colnames(c3_aci_averages)) {
    # Rd_at_25 was fixed
    c3_aci_averages[, 'Rd_at_25_avg']    <- SOYBEAN_FIT_OPTIONS$Rd
    c3_aci_averages[, 'Rd_at_25_stderr'] <- NA
}

# View the averages and errors
columns_to_view <- c(
  'year', 'doy', 'species',
  'Vcmax_at_25_avg', 'Vcmax_at_25_stderr',
  'J_at_25_avg', 'J_at_25_stderr',
  'Rd_at_25_avg', 'Rd_at_25_stderr',
  'Tp_avg', 'Tp_stderr'
)
str(c3_aci_averages[ , columns_to_view, TRUE])

# Wow, the results from these curves are all really consistent with each other!

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
soybean_ld11_fvcb_parameters_2021$Jmax <- get_jmax(
    soybean$parameters$theta,              # dimensionless
    soybean$parameters$beta_PSII,          # dimensionless
    soybean_ld11_fvcb_parameters_2021$J,   # micromol / m^2 / s
    soybean$parameters$leaf_reflectance,   # dimensionless
    25,                                    # degrees C
    soybean$parameters$leaf_transmittance, # dimensionless
    mean(licor_data[, 'Qin'])              # micromol / m^2 / s
)

# Save the parameter values
save(soybean_ld11_fvcb_parameters_2021, file = 'soybean_ld11_fvcb_parameters_2021.RData')
write.csv(as.data.frame(soybean_ld11_fvcb_parameters_2021), file = 'soybean_ld11_fvcb_parameters_2021.csv', row.names = FALSE)
