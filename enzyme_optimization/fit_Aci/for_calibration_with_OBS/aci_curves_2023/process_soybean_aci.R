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
MAKE_INDIVIDUAL_PLOTS <- TRUE
REMOVE_BAD_POINTS <- TRUE
SAVE_TO_PDF <- FALSE

#CURVE_FOR_INDIVIDUAL <- '240 - ripe8 - 1'
#CURVE_FOR_INDIVIDUAL <- '223 - mcgrath1 - 3'
#CURVE_FOR_INDIVIDUAL <- '194 - mcgrath1 - 4'
CURVE_FOR_INDIVIDUAL <- '244 - ripe14 - 2'

# Define a plotting function that automatically names each plot
plot_num <- 1

pdf_print2 <- function(...) {
  pdf_print(..., file = paste0('soybean_', sprintf('%03d', plot_num), '.pdf'), save_to_pdf = SAVE_TO_PDF)
  plot_num <<- plot_num + 1
}

# Define axis limits
a_lim <- c(-10, 70)
dad_lim <- c(0, 45)
j_lim <- c(0, 250)
rd_lim <- c(0, 2.0)
resid_lim <- c(-5, 5)
spad_lim <- c(30, 65)
tp_lim <- c(0, 25)
vcmax_lim <- c(0, 200)

###
### TRANSLATION:
### Creating convenient R objects from raw data files
###

# Define a vector of paths to the files we wish to load
file_paths <- c(
  "time 1 aci/2023-07-13 ed aci mcgrath1.xlsx",
  "time 1 aci/2023-07-13 ed aci mcgrath2.xlsx",
  "time 1 aci/2023-07-13 ed aci ripe2.xlsx",
  "time 1 aci/2023-07-13 ed aci ripe5.xlsx",
  "time 2 aci/2023-08-10 ed mcgrath1 aci.xlsx",
  "time 2 aci/2023-08-10 ed mcgrath2 aci.xlsx",
  "time 2 aci/2023-08-10 ed ripe9 aci.xlsx", # does not include PhiPSII
  "time 2 aci/2023-08-10 ed ripe14 aci.xlsx",
  "time 3 aci/2023-08-28 ed ripe5 aci.xlsx",
  "time 3 aci/2023-08-28 ed ripe6 aci.xlsx",
  "time 3 aci/2023-08-28 ed ripe8 aci.xlsx",
  "time 3 aci/2023-08-28 ed ripe14 aci.xlsx",
  "time 3 aci/2023-09-01 ed ripe5 aci.xlsx",
  "time 3 aci/2023-09-01 ed ripe6 aci.xlsx",
  "time 3 aci/2023-09-01 ed ripe8 aci.xlsx",
  "time 3 aci/2023-09-01 ed ripe14 aci.xlsx"
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

# Only keep the soybean curves (some of the files contain sorghum and maize
# data in addition to soybean)
licor_data <-
  licor_data[licor_data[, 'species'] %in% c('soybean ld10', 'soybean ld11'), , TRUE]

# Add columns for year and day of year
licor_data[, 'year'] <- as.numeric(format(licor_data[, 'time'], '%Y'))
licor_data[, 'doy'] <- as.numeric(format(licor_data[, 'time'], '%j'))

# Fix some of the DOY values; some of the Licors had the incorrect day
licor_data[licor_data[, 'doy'] == 193, 'doy'] <- 194
licor_data[licor_data[, 'doy'] == 222, 'doy'] <- 223
licor_data[licor_data[, 'doy'] == 224, 'doy'] <- 223

# Shorten the species names for plotting purposes
licor_data[licor_data[, 'species'] == 'soybean ld10', 'species'] <- 'ld10'
licor_data[licor_data[, 'species'] == 'soybean ld11', 'species'] <- 'ld11'

# Create new identifier columns
licor_data[ , 'curve_identifier'] <-
  paste(licor_data[, 'doy'], licor_data[, 'instrument'], licor_data[, 'replicate'], sep = ' - ')

licor_data[ , 'species_leaf'] <-
  paste(licor_data[, 'species'], licor_data[, 'leaf'], sep = ' - ')

licor_data[ , 'doy_species_leaf'] <-
  paste(licor_data[, 'doy'], licor_data[, 'species'], licor_data[, 'leaf'], sep = ' - ')

# Add DOY when each leaf was fully developed.
licor_data[licor_data[, 'species_leaf'] == 'ld10 - 2',  'doy_devel'] <- 186
licor_data[licor_data[, 'species_leaf'] == 'ld10 - 3',  'doy_devel'] <- 187
licor_data[licor_data[, 'species_leaf'] == 'ld10 - 4',  'doy_devel'] <- 191
licor_data[licor_data[, 'species_leaf'] == 'ld10 - 8',  'doy_devel'] <- 205
licor_data[licor_data[, 'species_leaf'] == 'ld10 - 10', 'doy_devel'] <- 210
licor_data[licor_data[, 'species_leaf'] == 'ld10 - 11', 'doy_devel'] <- 213 # estimated - LD10 was in V10 on 210 but there is no log entry where LD10 was in V11
licor_data[licor_data[, 'species_leaf'] == 'ld10 - 12', 'doy_devel'] <- 216 # estimated - LD10 was in V10 on 210 and V13 on 220 but there is no log entry where LD10 was in V12
licor_data[licor_data[, 'species_leaf'] == 'ld11 - 3',  'doy_devel'] <- 187
licor_data[licor_data[, 'species_leaf'] == 'ld11 - 4',  'doy_devel'] <- 191
licor_data[licor_data[, 'species_leaf'] == 'ld11 - 8',  'doy_devel'] <- 207
licor_data[licor_data[, 'species_leaf'] == 'ld11 - 11', 'doy_devel'] <- 213 # estimated - LD11 was in V10 on 210 but there is no log entry where LD11 was in V11

# Add developement stage on each day
licor_data[licor_data[, 'doy'] == 194, 'stage'] <- 'V5'
licor_data[licor_data[, 'doy'] == 223, 'stage'] <- 'R4'
licor_data[licor_data[, 'doy'] == 240, 'stage'] <- 'R6'
licor_data[licor_data[, 'doy'] == 244, 'stage'] <- 'R6'

# Add a column for days after development (DAD) for each leaf
licor_data[, 'dad'] <- licor_data[, 'doy'] - licor_data[, 'doy_devel']

# Print curve information
print(unique(licor_data[, c('year', 'doy', 'instrument', 'replicate', 'species', 'plot', 'leaf', 'doy_devel', 'dad')]))

# Factorize some columns to help with plotting later
licor_data <- factorize_id_column(licor_data, 'leaf')

###
### VALIDATION:
### Organizing the data, checking its consistency and quality, cleaning it
###

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

if (MAKE_VALIDATION_PLOTS) {
  # Plot all A-Ci curves in the data set, indicating stability
  pdf_print2(xyplot(
    A ~ Ci | curve_identifier,
    group = Stable,
    data = licor_data$main_data,
    type = 'b',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    ylim = a_lim,
    xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
    ylab = paste('Net CO2 assimilation rate [', licor_data$units$A, ']')
  ))

  if ('PhiPS2' %in% colnames(licor_data)) {
    # Plot all PhiPS2-Ci curves in the data set, indicating stability
    pdf_print2(xyplot(
      PhiPS2 ~ Ci | curve_identifier,
      group = Stable,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste('PhiPS2 [', licor_data$units$PhiPS2, ']')
    ))
  }

  # Plot all A-CO2_r_sp curves in the data set, indicating stability
  pdf_print2(xyplot(
    A ~ CO2_r_sp | curve_identifier,
    group = Stable,
    data = licor_data$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    ylim = a_lim,
    xlab = paste('CO2_r setpoint [', licor_data$units$CO2_r_sp, ']'),
    ylab = paste('Net CO2 assimilation rate [', licor_data$units$A, ']')
  ))

  # Plot all gsw-CO2_r_sp curves in the data set, indicating stability
  pdf_print2(xyplot(
    gsw ~ CO2_r_sp | curve_identifier,
    group = Stable,
    data = licor_data$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlab = paste('CO2_r setpoint [', licor_data$units$CO2_r_sp, ']'),
    ylab = paste('Stomatal conductance to H2O [', licor_data$units$gsw, ']')
  ))

  # Make a plot to check humidity control
  pdf_print2(xyplot(
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
  pdf_print2(xyplot(
    TleafCnd + Txchg ~ Ci | curve_identifier,
    data = licor_data$main_data,
    type = 'b',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
    ylab = paste0('Temperature (', licor_data$units$TleafCnd, ')')
  ))

  # Make a plot to check CO2 control
  pdf_print2(xyplot(
    CO2_s + CO2_r + CO2_r_sp ~ Ci | curve_identifier,
    data = licor_data$main_data,
    type = 'b',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
    ylab = paste0('CO2 concentration (', licor_data$units$CO2_r, ')')
  ))

  # Make a plot to check stability criteria. Here we cannot plot the
  # individual criteria (`A:OK` and `gsw:OK`) because they are not properly
  # included in the file from mcgrath2.
  pdf_print2(xyplot(
    Stable ~ Ci | curve_identifier,
    data = licor_data$main_data,
    type = 'b',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']')
  ))
}

# Export curves
write.csv(licor_data[licor_data[, 'species'] == 'ld10', , TRUE], 'aci_soybean_ld10_2023.csv', row.names = FALSE)
write.csv(licor_data[licor_data[, 'species'] == 'ld11', , TRUE], 'aci_soybean_ld11_2023.csv', row.names = FALSE)

# Some curves have strange points that we might want to exclude
if (REMOVE_BAD_POINTS) {
  licor_data <- remove_points(
    licor_data,
    list(doy = 194, instrument = 'mcgrath1', replicate = 3, CO2_r_sp = c(600, 800, 1000, 1200)),
    list(doy = 194, instrument = 'ripe2',    replicate = 1, CO2_r_sp = c(1200)),
    list(doy = 194, instrument = 'ripe2',    replicate = 3, CO2_r_sp = c()),
    list(doy = 194, instrument = 'ripe5',    replicate = 1, CO2_r_sp = c(600, 800, 1000, 1200, 1500, 1800)),
    list(doy = 223, instrument = 'mcgrath2', replicate = 1, CO2_r_sp = c(1200)),
    list(doy = 223, instrument = 'mcgrath2', replicate = 2, CO2_r_sp = c(1200)),
    list(doy = 223, instrument = 'mcgrath2', replicate = 4, CO2_r_sp = c(600, 800)),
    list(doy = 223, instrument = 'ripe14',   replicate = 1, CO2_r_sp = c(600, 800, 1000, 1200)),
    list(doy = 223, instrument = 'ripe14',   replicate = 2, CO2_r_sp = c(1500)),
    list(doy = 223, instrument = 'ripe14',   replicate = 3, CO2_r_sp = c(600)),
    list(doy = 223, instrument = 'ripe14',   replicate = 4, CO2_r_sp = c(1800)),
    list(doy = 223, instrument = 'ripe9',    replicate = 1, CO2_r_sp = c(800, 1000, 1200, 1500, 1800)),
    list(doy = 223, instrument = 'ripe9',    replicate = 2, CO2_r_sp = c(1800)),
    list(doy = 223, instrument = 'ripe9',    replicate = 3, CO2_r_sp = c(1200, 1500, 1800))
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
  remove_unreliable_param = TRUE
))

# Plot the C3 A-Ci fits (including limiting rates)
pdf_print2(xyplot(
  A + Ac + Aj + Ap + A_fit ~ Ci | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto.key = list(space = 'right', lines = TRUE, points = TRUE),
  grid = TRUE,
  xlab = paste('Intercellular CO2 concentration [', c3_aci_results$fits$units$Ci, ']'),
  ylab = paste('Net CO2 assimilation rate [', c3_aci_results$fits$units$A, ']'),
  par.settings = list(
    superpose.line = list(col = multi_curve_line_colors()),
    superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
  ),
  ylim = a_lim,
  curve_ids = c3_aci_results$fits[, 'curve_identifier'],
  panel = function(...) {
    panel.xyplot(...)
    args <- list(...)
    curve_id <- args$curve_ids[args$subscripts][1]
    fit_param <-
      c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == curve_id, ]
    panel.points(
        fit_param$operating_An_model ~ fit_param$operating_Cc,
        type = 'p',
        col = 'black',
        pch = 1
    )
  }
), width = 12)

# Plot the C3 A-Ci fits
pdf_print2(xyplot(
  A + A_fit ~ Ci | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  ylim = a_lim,
  xlab = paste('Intercellular CO2 concentration [', c3_aci_results$fits$units$Ci, ']'),
  ylab = paste('Net CO2 assimilation rate [', c3_aci_results$fits$units$A, ']')
), width = 12)

# Plot the residuals
pdf_print2(xyplot(
  A_residuals ~ Ci | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  grid = TRUE,
  ylim = resid_lim,
  xlab = paste('Intercellular CO2 concentration [', c3_aci_results$fits$units$Ci, ']'),
  ylab = paste('Assimilation rate residual (measured - fitted)\n[', c3_aci_results$fits$units$A, ']')
), width = 12)

# Plot the residuals
pdf_print2(xyplot(
  A_residuals ~ CO2_r_sp | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  grid = TRUE,
  ylim = resid_lim,
  xlab = paste('CO2_r setpoint [', c3_aci_results$fits$units$CO2_r_sp, ']'),
  ylab = paste('Assimilation rate residual (measured - fitted)\n[', c3_aci_results$fits$units$A, ']')
), width = 12)

###
### SYNTHESIS:
### Using plots and statistics to help draw conclusions from the data
###

print('Max DAD')
print(max(c3_aci_results$parameters[, 'dad'], na.rm = TRUE))
print('Max SPAD')
print(max(c3_aci_results$parameters[, 'spad'], na.rm = TRUE))
print('Max Rd')
print(max(c3_aci_results$parameters[, 'Rd_at_25'], na.rm = TRUE))
print('Max Vcmax')
print(max(c3_aci_results$parameters[, 'Vcmax_at_25'], na.rm = TRUE))
print('Max J')
print(max(c3_aci_results$parameters[, 'J_at_25'], na.rm = TRUE))
print('Max Tp')
print(max(c3_aci_results$parameters[, 'Tp'], na.rm = TRUE))

# Make bar and box plots of the parameter values for each cultivar
box_bar_plot_X <- 'doy_species_leaf'

box_bar_plot_param <- list(
  list(
    Y = c3_aci_results$parameters[, 'Vcmax_at_25'],
    X = c3_aci_results$parameters[, box_bar_plot_X],
    ylim = vcmax_lim,
    ylab = paste('Vcmax at 25 degrees C [', c3_aci_results$parameters$units$Vcmax_at_25, ']')
  ),
  list(
    Y = c3_aci_results$parameters[, 'J_at_25'],
    X = c3_aci_results$parameters[, box_bar_plot_X],
    ylim = j_lim,
    ylab = paste('J at 25 degrees C [', c3_aci_results$parameters$units$J_at_25, ']')
  ),
  list(
    Y = c3_aci_results$parameters[, 'Rd_at_25'],
    X = c3_aci_results$parameters[, box_bar_plot_X],
    ylim = rd_lim,
    ylab = paste('Rd at 25 degrees C [', c3_aci_results$parameters$units$Rd_at_25, ']')
  ),
  list(
    Y = c3_aci_results$parameters[, 'Tp'],
    X = c3_aci_results$parameters[, box_bar_plot_X],
    ylim = tp_lim,
    ylab = paste('Tp [', c3_aci_results$parameters$units$Tp, ']')
  ),
  list(
    Y = c3_aci_results$parameters[, 'spad'],
    X = as.factor(c3_aci_results$parameters[, 'leaf']),
    ylim = c(0, max(spad_lim)),
    ylab = 'SPAD'
  )
)

invisible(lapply(box_bar_plot_param, function(x) {
  pdf_print2(do.call(bwplot_wrapper, x))
  pdf_print2(do.call(barchart_with_errorbars, x))
}))

pdf_print2(xyplot(
  spad ~ leaf | species,
  group = doy,
  data = c3_aci_results$parameters$main_data,
  grid = TRUE,
  auto = TRUE,
  type = 'p',
  pch = 16
), width = 12)

plot_v_spad <- function(exdf_obj) {
  caption <- paste('DOY:', paste(unique(exdf_obj[, 'doy']), collapse = ', '))

  # Plot the fitted parameter values against SPAD for each species
  pdf_print2(xyplot(
    Vcmax_at_25 ~ spad | species,
    group = leaf,
    data = exdf_obj$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = spad_lim,
    ylim = vcmax_lim,
    xlab = 'SPAD',
    ylab = paste('Vcmax at 25 degrees C [', exdf_obj$units$Vcmax_at_25, ']'),
    main = caption
  ))

  pdf_print2(xyplot(
    J_at_25 ~ spad | species,
    group = leaf,
    data = exdf_obj$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = spad_lim,
    ylim = j_lim,
    xlab = 'SPAD',
    ylab = paste('J at 25 degrees C [', exdf_obj$units$J_at_25, ']'),
    main = caption
  ))

  pdf_print2(xyplot(
    Rd_at_25 ~ spad | species,
    group = leaf,
    data = exdf_obj$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = spad_lim,
    ylim = rd_lim,
    xlab = 'SPAD',
    ylab = paste('Rd at 25 degrees C [', exdf_obj$units$Rd_at_25, ']'),
    main = caption
  ))

  pdf_print2(xyplot(
    Tp ~ spad | species,
    group = leaf,
    data = exdf_obj$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = spad_lim,
    ylim = tp_lim,
    xlab = 'SPAD',
    ylab = paste('Tp [', exdf_obj$units$Tp, ']'),
    main = caption
  ))
}

to_plot_v_spad <- list(
  c3_aci_results$parameters[c3_aci_results$parameters[, 'doy'] == 194, , TRUE],
  c3_aci_results$parameters[c3_aci_results$parameters[, 'doy'] == 223, , TRUE],
  c3_aci_results$parameters[c3_aci_results$parameters[, 'doy'] >= 240, , TRUE]
)

lapply(to_plot_v_spad, plot_v_spad)

plot_against <- function(ind, x_lim) {
  # Plot the fitted parameter values against SPAD for each species
  pdf_print2(xyplot(
    Vcmax_at_25 ~ c3_aci_results$parameters[, ind] | species,
    group = doy,
    data = c3_aci_results$parameters$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = x_lim,
    ylim = vcmax_lim,
    xlab = ind,
    ylab = paste('Vcmax at 25 degrees C [', c3_aci_results$parameters$units$Vcmax_at_25, ']')
  ))

  pdf_print2(xyplot(
    J_at_25 ~ c3_aci_results$parameters[, ind] | species,
    group = doy,
    data = c3_aci_results$parameters$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = x_lim,
    ylim = j_lim,
    xlab = ind,
    ylab = paste('J at 25 degrees C [', c3_aci_results$parameters$units$J_at_25, ']')
  ))

  pdf_print2(xyplot(
    Rd_at_25 ~ c3_aci_results$parameters[, ind] | species,
    group = doy,
    data = c3_aci_results$parameters$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = x_lim,
    ylim = rd_lim,
    xlab = ind,
    ylab = paste('Rd at 25 degrees C [', c3_aci_results$parameters$units$Rd_at_25, ']')
  ))

  pdf_print2(xyplot(
    Tp ~ c3_aci_results$parameters[, ind] | species,
    group = doy,
    data = c3_aci_results$parameters$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = x_lim,
    ylim = tp_lim,
    xlab = ind,
    ylab = paste('Tp [', c3_aci_results$parameters$units$Tp, ']')
  ))
}

plot_against('spad', spad_lim)
plot_against('dad',  dad_lim)

if (MAKE_INDIVIDUAL_PLOTS) {
  individual_curve <- c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] == CURVE_FOR_INDIVIDUAL, , TRUE]
  individual_curve_param <- c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == CURVE_FOR_INDIVIDUAL, ]

  # A-Ci
  pdf_print2(xyplot(
    A ~ Ci,
    data = individual_curve$main_data,
    type = 'b',
    pch = 16,
    grid = TRUE,
    xlim = c(-100, 1750),
    ylim = c(-10, 70),
    xlab = paste('Intercellular CO2 concentration [', individual_curve$units$Ci, ']'),
    ylab = paste('Net CO2 assimilation rate [', individual_curve$units$A, ']'),
    panel = function(...) {
        panel.xyplot(...)
        panel.points(
            individual_curve_param$operating_An_model ~ individual_curve_param$operating_Ci,
            type = 'p',
            col = 'red',
            pch = 16
        )
    }
  ))

  # A-Ci with fits
  pdf_print2(xyplot(
    A + Ac + Aj + Ap ~ Ci,
    data = individual_curve$main_data,
    type = 'b',
    pch = 16,
    grid = TRUE,
    xlim = c(-100, 1750),
    ylim = c(-10, 70),
    xlab = paste('Intercellular CO2 concentration [', individual_curve$units$Ci, ']'),
    ylab = paste('Net CO2 assimilation rate [', individual_curve$units$A, ']'),
    par.settings = list(
      superpose.line = list(col = multi_curve_line_colors()),
      superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
    ),
    panel = function(...) {
        panel.xyplot(...)
        panel.points(
            individual_curve_param$operating_An_model ~ individual_curve_param$operating_Ci,
            type = 'p',
            col = 'red',
            pch = 16
        )
    }
  ))

  if ('PhiPS2' %in% colnames(individual_curve)) {
    # Plot all PhiPS2-Ci curves in the data set, indicating stability
    pdf_print2(xyplot(
      PhiPS2 ~ Ci,
      data = individual_curve$main_data,
      type = 'b',
      pch = 16,
      grid = TRUE,
      xlim = c(-100, 1750),
      ylim = c(0, 0.5),
      xlab = paste('Intercellular CO2 concentration [', individual_curve$units$Ci, ']'),
      ylab = 'PhiPS2 [ dimensionless ]'
    ))
  }
}

# Save results
col_to_keep <- c(
  'doy', 'instrument', 'replicate', 'species', 'plot', 'leaf', 'spad',
  'J_at_25', 'Rd_at_25', 'Vcmax_at_25', 'Tp',
  'J_tl_avg', 'Rd_tl_avg', 'Vcmax_tl_avg',
  'operating_Ci', 'operating_An', 'operating_An_model',
  'convergence', 'convergence_msg', 'feval',
  'RSS', 'MSE', 'RMSE', 'RSE'
)

fit_parameters <- c3_aci_results$parameters[, col_to_keep, TRUE]

write.csv(fit_parameters, file = 'soybean_fit_parameters.csv', row.names = FALSE)

pdf_print2(xyplot(
  Tp ~ Vcmax_at_25 | species,
  group = doy,
  data = fit_parameters$main_data,
  type = 'p',
  pch = 16,
  xlim = c(40, 200),
  ylim = c(8, 22),
  auto = TRUE
))

pdf_print2(xyplot(
  operating_An ~ Vcmax_at_25 | species,
  group = doy,
  data = fit_parameters$main_data,
  type = 'p',
  pch = 16,
  xlim = c(40, 200),
  ylim = c(20, 40),
  auto = TRUE,
  grid = TRUE
))

# Plot the fitted parameter values against SPAD for each species
pdf_print2(xyplot(
  Vcmax_at_25 ~ spad | stage,
  group = species,
  data = c3_aci_results$parameters$main_data,
  type = 'p',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlim = spad_lim,
  ylim = vcmax_lim,
  ylab = paste('Vcmax at 25 degrees C [', c3_aci_results$parameters$units$Vcmax_at_25, ']')
))

pdf_print2(xyplot(
  J_at_25 ~ spad | stage,
  group = species,
  data = c3_aci_results$parameters$main_data,
  type = 'p',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlim = spad_lim,
  ylim = j_lim,
  ylab = paste('J at 25 degrees C [', c3_aci_results$parameters$units$J_at_25, ']')
))

pdf_print2(xyplot(
  Rd_at_25 ~ spad | stage,
  group = species,
  data = c3_aci_results$parameters$main_data,
  type = 'p',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlim = spad_lim,
  ylim = rd_lim,
  ylab = paste('Rd at 25 degrees C [', c3_aci_results$parameters$units$Rd_at_25, ']')
))

pdf_print2(xyplot(
  Tp ~ spad | stage,
  group = species,
  data = c3_aci_results$parameters$main_data,
  type = 'p',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlim = spad_lim,
  ylim = tp_lim,
  ylab = paste('Tp [', c3_aci_results$parameters$units$Tp, ']')
))

pdf_print2(xyplot(
  operating_An ~ spad | stage,
  group = species,
  data = c3_aci_results$parameters$main_data,
  type = 'p',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlim = spad_lim,
  ylim = c(0, 45),
  ylab = paste('Operating An [', c3_aci_results$parameters$units$operating_An, ']')
))

pdf_print2(xyplot(
  operating_Ci ~ spad | stage,
  group = species,
  data = c3_aci_results$parameters$main_data,
  type = 'p',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlim = spad_lim,
  ylim = c(220, 340),
  ylab = paste('Operating Ci [', c3_aci_results$parameters$units$operating_Ci, ']')
))
