# This script is based on the example in the 'Analyzing C4 A-Ci Curves'
# vignette from the PhotoGEA package, available online at
# https://eloch216.github.io/PhotoGEA/articles/analyzing_c4_aci_curves.html

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

# Make some decisions about what to do
MAKE_VALIDATION_PLOTS <- TRUE
REMOVE_UNSTABLE_POINTS <- FALSE
REMOVE_BAD_CURVES <- FALSE
SAVE_TO_PDF <- FALSE

# Define a plotting function that automatically names each plot
plot_num <- 1

pdf_print2 <- function(...) {
  pdf_print(..., file = paste0('c4_', sprintf('%03d', plot_num), '.pdf'), save_to_pdf = SAVE_TO_PDF)
  plot_num <<- plot_num + 1
}

# Define axis limits
a_lim <- c(-5, 75)
dad_lim <- c(-5, 45)
rd_lim <- c(0, 3)
resid_lim <- c(-5, 5)
spad_lim <- c(40, 70)
vcmax_lim <- c(0, 40)
vpmax_lim <- c(0, 120)

###
### TRANSLATION:
### Creating convenient R objects from raw data files
###

# Define a vector of paths to the files we wish to load
file_paths <- c(
  "time 1 aci/2023-07-14 ed aci mcgrath1.xlsx",
  "time 1 aci/2023-07-14 ed aci mcgrath2.xlsx",
  "time 1 aci/2023-07-14 ed aci ripe 2.xlsx",
  "time 1 aci/2023-07-14 ed aci ripe5.xlsx",
  "time 2 aci/2023-08-11 ed mcgrath1 aci.xlsx",
  "time 2 aci/2023-08-11 ed mcgrath2 aci.xlsx",
  "time 2 aci/2023-08-11 ed ripe9 aci.xlsx",
  "time 2 aci/2023-08-11 ed ripe14 aci.xlsx",
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

# Only keep the maize and sorghum curves (some of the files contain soybean data
# in addition to maize and sorghum)
licor_data <-
  licor_data[licor_data[, 'species'] %in% c('maize', 'sorghum'), , TRUE]

# Add columns for year and day of year
licor_data[, 'year'] <- as.numeric(format(licor_data[, 'time'], '%Y'))
licor_data[, 'doy'] <- as.numeric(format(licor_data[, 'time'], '%j'))

# Fix some of the DOY values; some of the Licors had the incorrect day
licor_data[licor_data[, 'doy'] == 194, 'doy'] <- 195
licor_data[licor_data[, 'doy'] == 223, 'doy'] <- 224
licor_data[licor_data[, 'doy'] == 225, 'doy'] <- 224

# Create new identifier columns
licor_data[ , 'curve_identifier'] <-
  paste(licor_data[, 'doy'], licor_data[, 'instrument'], licor_data[, 'replicate'], sep = ' - ')

licor_data[ , 'species_leaf'] <-
  paste(licor_data[, 'species'], licor_data[, 'leaf'], sep = ' - ')

licor_data[ , 'doy_species_leaf'] <-
  paste(licor_data[, 'doy'], licor_data[, 'species'], licor_data[, 'leaf'], sep = ' - ')

# Add DOY when each leaf developed a collar. The maize ear leaf is leaf 12, and
# the sorghum flag leaf is leaf 17 (on average).
licor_data[licor_data[, 'species_leaf'] == 'maize - 7',        'doy_devel'] <- 184
licor_data[licor_data[, 'species_leaf'] == 'maize - 10',       'doy_devel'] <- 195
licor_data[licor_data[, 'species_leaf'] == 'maize - ear',      'doy_devel'] <- 202 # leaf 12
licor_data[licor_data[, 'species_leaf'] == 'maize - ear+4',    'doy_devel'] <- 213 # leaf 16; estimated - maize was in V15 on 210 but there is no log entry where maize was in V16
licor_data[licor_data[, 'species_leaf'] == 'sorghum - 9',      'doy_devel'] <- 188
licor_data[licor_data[, 'species_leaf'] == 'sorghum - 12',     'doy_devel'] <- 198
licor_data[licor_data[, 'species_leaf'] == 'sorghum - flag-4', 'doy_devel'] <- 201 # leaf 13; estimated - sorghum was in V12 on 198 but there is no log entry where sorghum was in V13
licor_data[licor_data[, 'species_leaf'] == 'sorghum - flag',   'doy_devel'] <- 210

# Add a column for days after development (DAD) for each leaf
licor_data[, 'dad'] <- licor_data[, 'doy'] - licor_data[, 'doy_devel']

# Print curve information
print(unique(licor_data[, c('year', 'doy', 'instrument', 'replicate', 'species', 'plot', 'leaf', 'doy_devel', 'dad')]))

###
### VALIDATION:
### Organizing the data, checking its consistency and quality, cleaning it
###

# Check the CO2_r_sp values from one curve to confirm that we have a 14 point
# curve where points 8 and 9 repeat the initial setpoint value
print(licor_data[licor_data[, 'curve_identifier'] == licor_data[1, 'curve_identifier'], 'CO2_r_sp'])

# Make sure the data meets basic requirements
check_licor_data(licor_data, 'curve_identifier', 14, 'CO2_r_sp')

# Remove points with duplicated `CO2_r_sp` values and order by `Ci`
licor_data <- organize_response_curve_data(
  licor_data,
  'curve_identifier',
  c(8, 9), # remove the extra points at 400
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
    list(instrument = 'ripe5', replicate = 1) # looks weird
  )
}

# Some points with low Ci produce PCm < 0. These must be removed, or the
# curves with these points cannot be fit.
licor_data <- remove_points(
    licor_data,
    list(seq_num = 7, curve_identifier = c(
        "195 - mcgrath2 - 1",
        "195 - mcgrath2 - 2",
        "195 - mcgrath2 - 3",
        "195 - mcgrath2 - 4",
        "224 - mcgrath2 - 1",
        "224 - mcgrath2 - 2",
        "224 - mcgrath2 - 3",
        "224 - mcgrath2 - 4"
    ))
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
    xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
    ylab = paste('Net CO2 assimilation rate [', licor_data$units$A, ']')
  ))

  # Plot all A-CO2_r_sp curves in the data set, indicating stability
  pdf_print2(xyplot(
    A ~ CO2_r_sp | curve_identifier,
    group = Stable,
    data = licor_data$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
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
    `A:OK` + `gsw:OK` + Stable ~ Ci | curve_identifier,
    data = licor_data$main_data,
    type = 'b',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']')
  ))
}

###
### PROCESSING:
### Extracting new pieces of information from the data
###

# Calculate temperature-dependent values of C4 photosynthetic parameters
licor_data <- calculate_arrhenius(licor_data, c4_arrhenius_von_caemmerer)

# Set mesophyll conductance
licor_data <- set_variable(
  licor_data,
  'gmc',
  'mol m^(-2) s^(-1) bar^(-1)',
  'process_c4_aci',
  Inf
)

# Calculate the total pressure in the Licor chamber
licor_data <- calculate_total_pressure(licor_data)

# Calculate PCm
licor_data <- apply_gm(
  licor_data,
  'C4' # Indicate C4 photosynthesis
)

# Fit the C4 A-Ci curves
c4_aci_results <- consolidate(by(
  licor_data,                       # The `exdf` object containing the curves
  licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
  fit_c4_aci,                       # The function to apply to each chunk of `licor_data`
  Ca_atmospheric = 420,
  alpha = 0,
  gbs = 0,
  Rm_frac = 1
))

# Plot the C4 A-Ci fits
pdf_print2(xyplot(
  A + A_fit ~ Ci | curve_identifier,
  data = c4_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  ylim = a_lim,
  xlab = paste('Intercellular CO2 concentration [', c4_aci_results$fits$units$Ci, ']'),
  ylab = paste('Net CO2 assimilation rate [', c4_aci_results$fits$units$A, ']')
), width = 12)

# Set point colors so that only the first curve has points; the others will all
# have points with alpha = 0 (which is fully transparent)
point_colors <- c(
  multi_curve_colors()[1],
  paste0(multi_curve_colors()[seq(2, length(multi_curve_colors()))], '00')
)

# Set line colors so that the first curve has a line, but the others don't; the
# first line will have alpha = 0 (which is fully transparent)
line_colors <- multi_curve_colors()
line_colors[1] <- paste0(line_colors[1], '00')

# Plot the C4 A-PCm fits (including limiting rates)
pdf_print2(xyplot(
  A + Apc + Ar + A_fit ~ PCm | curve_identifier,
  data = c4_aci_results$fits$main_data,
  type = 'b',
  auto.key = list(space = 'right', lines = TRUE, points = TRUE),
  grid = TRUE,
  ylim = a_lim,
  xlab = paste('Mesophyll CO2 pressure [', c4_aci_results$fits$units$PCm, ']'),
  ylab = paste('Net CO2 assimilation rate [', c4_aci_results$fits$units$A, ']'),
  par.settings = list(
    superpose.line = list(col = line_colors),
    superpose.symbol = list(col = point_colors, pch = 16)
  )
), width = 12)

# Plot the residuals
pdf_print2(xyplot(
  A_residuals ~ Ci | curve_identifier,
  data = c4_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  grid = TRUE,
  ylim = resid_lim,
  xlab = paste('Intercellular CO2 concentration [', c4_aci_results$fits$units$Ci, ']'),
  ylab = paste('Assimilation rate residual (measured - fitted)\n[', c4_aci_results$fits$units$A, ']')
), width = 12)

# Plot the residuals
pdf_print2(xyplot(
  A_residuals ~ CO2_r_sp | curve_identifier,
  data = c4_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  grid = TRUE,
  ylim = resid_lim,
  xlab = paste('CO2_r setpoint [', c4_aci_results$fits$units$CO2_r_sp, ']'),
  ylab = paste('Assimilation rate residual (measured - fitted)\n[', c4_aci_results$fits$units$A, ']')
), width = 12)

###
### SYNTHESIS:
### Using plots and statistics to help draw conclusions from the data
###

print('Max DAD')
print(max(c4_aci_results$parameters[, 'dad']))
print('Max SPAD')
print(max(c4_aci_results$parameters[, 'spad']))
print('Max Rd')
print(max(c4_aci_results$parameters[, 'Rd_at_25']))
print('Max Vcmax')
print(max(c4_aci_results$parameters[, 'Vcmax_at_25']))
print('Max Vpmax')
print(max(c4_aci_results$parameters[, 'Vpmax_at_25']))

# Make bar and box plots of the parameter values for each cultivar
box_bar_plot_X <- 'doy_species_leaf'

box_bar_plot_param <- list(
  list(
    Y = c4_aci_results$parameters[, 'Vcmax_at_25'],
    X = c4_aci_results$parameters[, box_bar_plot_X],
    ylim = vcmax_lim,
    ylab = paste('Vcmax at 25 degrees C [', c4_aci_results$parameters$units$Vcmax_at_25, ']')
  ),
  list(
    Y = c4_aci_results$parameters[, 'Vpmax_at_25'],
    X = c4_aci_results$parameters[, box_bar_plot_X],
    ylim = vpmax_lim,
    ylab = paste('Vpmax at 25 degrees C [', c4_aci_results$parameters$units$Vpmax_at_25, ']')
  ),
  list(
    Y = c4_aci_results$parameters[, 'Rd_at_25'],
    X = c4_aci_results$parameters[, box_bar_plot_X],
    ylim = rd_lim,
    ylab = paste('Rd at 25 degrees C [', c4_aci_results$parameters$units$Rd_at_25, ']')
  ),
  list(
    Y = c4_aci_results$parameters[, 'spad'],
    X = as.factor(c4_aci_results$parameters[, 'leaf']),
    ylim = c(0, max(spad_lim)),
    ylab = 'SPAD'
  )
)

invisible(lapply(box_bar_plot_param, function(x) {
  pdf_print2(do.call(bwplot_wrapper, x))
  pdf_print2(do.call(barchart_with_errorbars, x))
}))

pdf_print2(xyplot(
  spad ~ factor(leaf) | species,
  group = doy,
  data = c4_aci_results$parameters$main_data,
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
    Vpmax_at_25 ~ spad | species,
    group = leaf,
    data = exdf_obj$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = spad_lim,
    ylim = vpmax_lim,
    xlab = 'SPAD',
    ylab = paste('Vpmax at 25 degrees C [', exdf_obj$units$Vpmax_at_25, ']'),
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
}

to_plot_v_spad <- list(
  c4_aci_results$parameters[c4_aci_results$parameters[, 'doy'] == 195, , TRUE],
  c4_aci_results$parameters[c4_aci_results$parameters[, 'doy'] == 224, , TRUE],
  c4_aci_results$parameters[c4_aci_results$parameters[, 'doy'] >= 240, , TRUE]
)

lapply(to_plot_v_spad, plot_v_spad)

# Plot the fitted parameter values against SPAD for each species
plot_against <- function(ind, x_lim) {
  pdf_print2(xyplot(
    Vcmax_at_25 ~ c4_aci_results$parameters[, ind] | species,
    group = doy,
    data = c4_aci_results$parameters$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = x_lim,
    ylim = vcmax_lim,
    xlab = ind,
    ylab = paste('Vcmax at 25 degrees C [', c4_aci_results$parameters$units$Vcmax_at_25, ']')
  ))

  pdf_print2(xyplot(
    Vpmax_at_25 ~ c4_aci_results$parameters[, ind] | species,
    group = doy,
    data = c4_aci_results$parameters$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = x_lim,
    ylim = vpmax_lim,
    xlab = ind,
    ylab = paste('Vpmax at 25 degrees C [', c4_aci_results$parameters$units$Vpmax_at_25, ']')
  ))

  pdf_print2(xyplot(
    Rd_at_25 ~ c4_aci_results$parameters[, ind] | species,
    group = doy,
    data = c4_aci_results$parameters$main_data,
    type = 'p',
    pch = 16,
    auto = TRUE,
    grid = TRUE,
    xlim = x_lim,
    ylim = rd_lim,
    xlab = ind,
    ylab = paste('Rd at 25 degrees C [', c4_aci_results$parameters$units$Rd_at_25, ']')
  ))
}

plot_against('spad', spad_lim)
plot_against('dad',  dad_lim)
