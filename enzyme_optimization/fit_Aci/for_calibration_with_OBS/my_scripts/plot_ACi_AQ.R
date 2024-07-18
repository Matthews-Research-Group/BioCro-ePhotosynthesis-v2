library(PhotoGEA)
library(BioCro)
library(ggplot2)
library(reshape2)
source("biocro_FvCB.R")

plot_type = 2 #1:ACi; 2: AQ
year = 2022

if(plot_type == 1){
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
  
  rc = c()
  for (i in 1:length(curve_identifier_unique)){
    curve_i = curve_identifier_unique[i]
    Qin   = Qin_all[curve_identifier == curve_i]
    Tleaf = Tleaf_all[curve_identifier == curve_i]
    Ci    = Ci_all[curve_identifier == curve_i]
    Vcmax_at_25 = parameters$Vcmax_at_25[i]
    Jmax_at_25 = parameters$Jmax_at_25[i]
    TPU        = parameters$TPU_at_25[i]
    Rd_at_25   = parameters$Rd_at_25[i]
    tmp  = BioCro_FvCB(Qin,Tleaf,Ci,Vcmax_at_25,Jmax_at_25,Rd_at_25,TPU)
    tmp  = tmp[order(Ci),]
    rc   = rbind(rc,tmp)
  }
  
  df = cbind(rc,obs_An,curve_identifier)
  
  df = df[,c("Ci","obs_An","An","curve_identifier")]
  
  df_melt = melt(df,id.vars = c("Ci","curve_identifier"),variable.name = "variable", value.name = "value")
  
  levels(df_melt$variable) = c("Observed","Modelled")
  
  df_sub = df_melt[ grep(paste0("^",year), df_melt$curve_identifier),]
  
  ggplot(df_sub, aes(x = Ci, y = value, color = curve_identifier, linetype = variable,shape = variable)) +
    geom_line() + geom_point() +
    labs(x = "Ci (ppm)", 
         y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
         ) +
    theme_bw() +
    theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
  
}else if(plot_type == 2){

  aci_fit_results <- read.csv("ld11_aci_fit_parameters_avg.csv")
  Vcmax_at_25 <- aci_fit_results$mean[aci_fit_results$parameter=="Vcmax_at_25"]
  Jmax_at_25  <- aci_fit_results$mean[aci_fit_results$parameter=="Jmax_at_25"]
  Rd_at_25    <- 1.28
  TPU          <- aci_fit_results$mean[aci_fit_results$parameter=="TPU_at_25"]
  
  licor_data <- read.csv.exdf("ld11_bb.csv") #this comes from "get_ld11_bb.R"

  Tleaf_all  <- licor_data[, 'TleafCnd'] 
  Qin_all    <- licor_data[, 'Qin']
  Ci_all     <- licor_data[, 'Ci']

  date = substr(licor_data[, 'date'], 1, 8)
  
  curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],licor_data[, 'replicate'],sep='-')
  
  obs_An = licor_data[,'A']
  df = data.frame(Q = Qin_all, Observed=obs_An,Modelled=NA,curve_identifier=curve_identifier)
  for (i in 1:length(Tleaf_all)){
   tmp  = BioCro_FvCB(Qin_all[i],Tleaf_all[i],Ci_all[i],Vcmax_at_25,Jmax_at_25,Rd_at_25,TPU)
   df$Modelled[i] = tmp$An
  }
  df_melt = melt(df,id.vars = c("Q","curve_identifier"),variable.name = "variable", value.name = "value")
  
  df_sub = df_melt[ grep(paste0("^",year), df_melt$curve_identifier),]
  
  ggplot(df_sub, aes(x = Q, y = value, color = curve_identifier, linetype = variable,shape = variable)) +
    geom_line() + geom_point(size=3) +
    labs(x = bquote(Q~(mu*mol ~ m^-2 ~ s^-1)), 
         y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
    ) +
    theme_bw() +
    theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
  
  
}
