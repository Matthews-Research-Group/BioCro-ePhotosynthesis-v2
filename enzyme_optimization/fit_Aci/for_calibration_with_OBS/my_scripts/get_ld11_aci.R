# Here we just read all the A-Ci Licor files from 2021 and 2022, extract just
# the LD11 curves, and save them to a CSV file

# Clear the workspace
rm(list=ls())

# Load required packages
library(PhotoGEA)
library(BioCro)

# Load tools
source('aci_defaults.R')

source('biocro_FvCB.R') #Vcmax_multiplier 

# Define a vector of paths to the files we wish to load
file_paths <- c(
    file.path('..', 'aci_curves_2021', '2021-08-04 ed ripe1.xlsx'),
    file.path('..', 'aci_curves_2021', '2021-08-04 ed ripe2.xlsx'),
    file.path('..', 'aci_curves_2021', '2021-08-04 ed ripe3.xlsx'),
    file.path('..', 'aci_curves_2021', '2021-08-04 ed ripe4.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-07-21 mcgrath1.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-07-21 mcgrath2.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-08-03 ed aci ripe2.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-08-03 ed aci ripe3.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-08-03 ed aci ripe14.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-08-03 ed aci ripe15.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-08-05 ripe2 aci.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-08-05 ripe3 aci.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-08-05 ripe14 aci.xlsx'),
    file.path('..', 'aci_curves_2022', '2022-08-05 ripe15 aci.xlsx')
)

# Load each file, storing the result in a list
licor_exdf_list <- lapply(file_paths, read_gasex_file)

# Get the names of all columns that are present in all of the Licor files
columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)

# Extract just these columns
licor_exdf_list <- lapply(licor_exdf_list, function(x) {
  x[ , columns_to_keep, TRUE]
})

# Use `rbind` to combine all the data
licor_data <- do.call(rbind, licor_exdf_list)

# Keep just the LD11 data (it was just called 'soybean' in 2021)
licor_data <- licor_data[licor_data[, 'species'] %in% c('soybean', 'ld11'), , TRUE]

# Standardize the species name
licor_data[, 'species'] <- 'ld11'

###
### FITTING
###

# Add a new column for day
licor_data[, 'day'] <- substr(licor_data[, 'date'], 1, 8)

# Add a new column for year
licor_data[, 'year'] <- as.numeric(substr(licor_data[, 'date'], 1, 4))

# Add a curve identifier
licor_data[, 'curve_identifier'] <-
    paste(licor_data[, 'day'], licor_data[, 'instrument'], licor_data[, 'plot'], sep = '-')

# Check and organize the curves
check_response_curve_data(
    licor_data,
    'curve_identifier',
    expected_npts = 16,
    driving_column = 'CO2_r_sp'
)

licor_data <- organize_response_curve_data(
    licor_data,
    'curve_identifier',
    c(9, 10),
    'Ci'
)

# Calculate total pressure in the Licor chamber
licor_data <- calculate_total_pressure(licor_data)

# Specify separate mesophyll conductance values for each species
licor_data <- set_variable(
    licor_data,
    'gmc',
    units = 'mol m^(-2) s^(-1) bar^(-1)',
    value = Inf
)

# Calculate Cc
licor_data <- apply_gm(licor_data)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_data <- calculate_arrhenius(licor_data, c3_arrhenius_bernacchi)
# overwrite arrehenius with a new vcmax temperature response
licor_data[, 'Vcmax_norm'] <- Vcmax_multiplier(licor_data[, 'TleafCnd']+273.15)

# Write to a CSV
write.csv.exdf(licor_data, file = 'ld11_aci.csv')

# Fit the C3 A-Ci curves
c3_aci_results <- consolidate(by(
    licor_data,                       # The `exdf` object containing the curves
    licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
    fit_c3_aci,                       # The function to apply to each chunk of `licor_data`
    Ca_atmospheric = 420,             # Additional argument passed to `fit_c3_aci`
    atp_use = ELECTRONS_PER_CARBOXYLATION,
    nadph_use = 2 * ELECTRONS_PER_OXYGENATION,
    fit_options = SOYBEAN_FIT_OPTIONS
))

# Get Jmax values
c3_aci_results$parameters <- set_variable(
    c3_aci_results$parameters,
    'Jmax_at_25',
    units = c3_aci_results$parameters$units$J_at_25,
    value = sapply(c3_aci_results$parameters[, 'J_at_25'], function(x) {
        get_jmax(
            soybean$parameters$theta,
            soybean$parameters$beta_PSII,
            x,
            soybean$parameters$leaf_reflectance_par,
            25,
            soybean$parameters$leaf_transmittance_par,
            mean(licor_data[, 'Qin'])
        )
    })
)


#get mean Tleaf for each identifier
average_Tleaf = aggregate(licor_data[, 'TleafCnd'],list(licor_data[, 'curve_identifier']),mean)
# Get TPU at 25 C
c3_aci_results$parameters <- set_variable(
    c3_aci_results$parameters,
    'TPU_at_25',
    units = c3_aci_results$parameters$units$Tp,
    value = list(get_TPU_at_25(
             c3_aci_results$parameters[,'Tp'],
            average_Tleaf$x)
            )
)

# Write the parameter values to a CSV file
col_to_write <- c(
    'file_name', 'day', 'year', 'instrument', 'plot',
    'Vcmax_at_25', 'J_at_25', 'Jmax_at_25', 'Rd_at_25', 'TPU_at_25'
)

write.csv.exdf(
    c3_aci_results$parameters[, col_to_write, TRUE],
    file = 'ld11_aci_fit_parameters.csv'
)

# Write average values to a CSV file
param <- c('Vcmax_at_25', 'J_at_25', 'Jmax_at_25', 'Rd_at_25', 'TPU_at_25')

avg_param <- do.call(rbind, lapply(param, function(pn) {
    vals <- c3_aci_results$parameters[, pn]
    vals <- vals[!is.na(vals)]
    data.frame(parameter = pn, mean = mean(vals), stderr = sd(vals) / sqrt(length(vals)))
}))

write.csv(avg_param, row.names = FALSE, file = 'ld11_aci_fit_parameters_avg.csv')


pdf_print(
  plot_c3_aci_fit(
    c3_aci_results,
    'curve_identifier',
    'Ci',
    ylim = c(-10, 50)
  ),
  width = 10,
  save_to_pdf = TRUE,
  file = 'ld11_aci_fits.pdf'
)
