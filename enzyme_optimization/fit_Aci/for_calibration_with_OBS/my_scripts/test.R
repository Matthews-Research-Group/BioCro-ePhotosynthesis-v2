library(PhotoGEA)
library(BioCro)
library(ggplot2)
library(reshape2)
source("biocro_FvCB.R")

parameters <- read.csv.exdf("ld11_aci_fit_parameters.csv")
parameters <- parameters$main_data
curve_identifier_unique = paste(parameters$day,parameters$instrument,parameters$plot,sep='-')

licor_data <- read.csv.exdf("ld11_aci.csv") #this comes from "get_ld11_aci.R"

Tleaf_all  <- licor_data[, 'TleafCnd'] 
Qin_all    <- licor_data[, 'Qin']
Ci_all     <- licor_data[, 'Ci']

date = substr(licor_data[, 'date'], 1, 8)

curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],sep='-')

obs_An = licor_data[,'A']

#check one experiments
i=5
curve_i = curve_identifier_unique[i]
Qin   = Qin_all[curve_identifier == curve_i]
Tleaf = Tleaf_all[curve_identifier == curve_i]
Ci    = Ci_all[curve_identifier == curve_i]
Vcmax_at_25 = parameters$Vcmax_at_25[i]
Jmax_at_25 = parameters$Jmax_at_25[i]
TPU        = parameters$Tp[i]
Rd_at_25   = parameters$Rd_at_25[i]
tmp  = BioCro_FvCB(Qin,Tleaf,Ci,Vcmax_at_25,Jmax_at_25,Rd_at_25,TPU)
tmp  = tmp[order(Ci),]

y1 = tmp$An
y = c3_aci_results$fits$main_data$A_fit
y2 = y[grep("^20210805-ripe4", c3_aci_results$fits$main_data$curve_identifier)]
  
plot(y1,y2,xlab='A by FvCB',ylab='A by fit_c3_aci')