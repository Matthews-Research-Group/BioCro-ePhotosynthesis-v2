library(reshape2)
library(ggplot2)
library(BioCro)
library(PhotoGEA)
source("for_calibration_with_OBS/my_scripts/biocro_FvCB.R")
plot_type = 2 #1:ACi; 2: AQ
prefix = c("ACi","AQ")
output_figure_name = paste0("figs/",prefix[plot_type],"_with_OBS_v2_Pi30.pdf")

alpha1_alpha2 = read.csv('ePhotosynthesis_optimal_alpha1_alpha2.csv')
alpha1 = alpha1_alpha2[2]
alpha2 = alpha1_alpha2[3]

aci_fit_results <- read.csv("for_calibration_with_OBS/my_scripts/ld11_aci_fit_parameters_avg.csv")
Vcmax25      <- aci_fit_results$mean[aci_fit_results$parameter=="Vcmax_at_25"]
Jmax25       <- aci_fit_results$mean[aci_fit_results$parameter=="Jmax_at_25"]
Rd25         <- 1.28
TPU25        <- aci_fit_results$mean[aci_fit_results$parameter=="TPU_at_25"]

if(plot_type==1){
  licor_data <- read.csv.exdf("for_calibration_with_OBS/my_scripts/ld11_aci.csv") 
  
  Tleaf_all  <- licor_data[, 'TleafCnd'] 
  Qin_all    <- licor_data[, 'Qin']
  Ci_all     <- licor_data[, 'Ci']
  
  date = substr(licor_data[, 'date'], 1, 8)
  curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],sep='-')
  An_obs = licor_data[,'A']
  
  An_FvCB     = NA* 1:length(Tleaf_all)
  An_ePhoto   = NA* 1:length(Tleaf_all)
   
  for (i in 1:length(Tleaf_all)){ 
    PAR = Qin_all[i]
    Tleaf = Tleaf_all[i]
    Ci  = Ci_all[i]
    #For Farquhar, we need the total Q, as there's a calculation of
    #absorption inside the Farquhar function
    output_farquhar  = BioCro_FvCB(PAR,Tleaf, Ci, Vcmax25, Jmax25, Rd25, TPU25)
    An_FvCB[i]   = output_farquhar$An
  
    #call ephoto c++
    system(paste("./myephoto_single_Pi30.exe",alpha1,alpha2,PAR,Tleaf, Ci))
    #read in ephoto results
    ephoto = read.csv("output.data",header=FALSE)
    Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
    An_ePhoto[i] = ephoto[1] - Rd
  }
  An_ePhoto = as.numeric(An_ePhoto)
  df = data.frame(Ci=Ci_all,An_ePhoto,An_FvCB,An_obs,curve_identifier)

  #plot average line
  unique_identifier = unique(df$curve_identifier)
  xx = 0
  An_ePhoto  = c()
  An_FvCB    = c()
  An_obs     = c()
  for (id in unique_identifier){
    tmp = df[df$curve_identifier == id,1:4]
    xx = xx + tmp
    An_ePhoto = cbind(An_ePhoto,tmp$An_ePhoto)
    An_FvCB   = cbind(An_FvCB  ,tmp$An_FvCB)
    An_obs    = cbind(An_obs   ,tmp$An_obs)
  }
  sd_ePhoto   = apply(An_ePhoto,1,sd)
  sd_FvCB     = apply(An_FvCB,1,sd)
  sd_obs      = apply(An_obs,1,sd)
  df_avg = xx/length(unique_identifier)
  df_melt = melt(df_avg,id.vars = c("Ci"),variable.name = "variable", value.name = "value")
  df_melt = cbind(df_melt,sd = c(sd_ePhoto,sd_FvCB,sd_obs))
  myplot<-ggplot(df_melt, aes(x = Ci, y = value, color = variable)) +
    geom_line() + geom_point(size=3) +
    scale_color_manual(values=c("green","blue", "black"))+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05)) +
    labs(x = bquote(Ci~(ppm)),
         y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
    ) +
    theme_bw() +
    theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
# Save ggplot to PDF
  ggsave(output_figure_name, plot = myplot, width = 6, height = 6)

}else if(plot_type==2){
  licor_data <- read.csv.exdf("for_calibration_with_OBS/my_scripts/ld11_bb.csv") #this comes from "get_ld11_bb.R"

  Tleaf_all  <- licor_data[, 'TleafCnd'] 
  Qin_all    <- licor_data[, 'Qin']
  Ci_all     <- licor_data[, 'Ci']

  date = substr(licor_data[, 'date'], 1, 8)
  
  curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],licor_data[, 'replicate'],sep='-')
  
  An_obs = licor_data[,'A']

  An_FvCB     = NA* 1:length(Tleaf_all)
  An_ePhoto   = NA* 1:length(Tleaf_all)
   
  for (i in 1:length(Tleaf_all)){ 
    PAR = Qin_all[i]
    Tleaf = Tleaf_all[i]
    Ci  = Ci_all[i]
    #For Farquhar, we need the total Q, as there's a calculation of
    #absorption inside the Farquhar function
    output_farquhar  = BioCro_FvCB(PAR,Tleaf, Ci, Vcmax25, Jmax25, Rd25, TPU25)
    An_FvCB[i]   = output_farquhar$An
  
    #call ephoto c++
    system(paste("./myephoto_single_Pi30.exe",alpha1,alpha2,PAR,Tleaf, Ci))
    #read in ephoto results
    ephoto = read.csv("output.data",header=FALSE)
    Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
    An_ePhoto[i] = ephoto[1] - Rd
  }
  An_ePhoto = as.numeric(An_ePhoto)

  df = data.frame(Qin=Qin_all,An_ePhoto,An_FvCB,An_obs,curve_identifier)

  #plot average line
  unique_identifier = unique(df$curve_identifier)
  xx = 0
  An_ePhoto  = c()
  An_FvCB    = c()
  An_obs     = c()
  for (id in unique_identifier){
    tmp = df[df$curve_identifier == id,1:4]
    xx = xx + tmp
    An_ePhoto = cbind(An_ePhoto,tmp$An_ePhoto)
    An_FvCB   = cbind(An_FvCB  ,tmp$An_FvCB)
    An_obs    = cbind(An_obs   ,tmp$An_obs)
  }
  sd_ePhoto   = apply(An_ePhoto,1,sd)
  sd_FvCB     = apply(An_FvCB,1,sd)
  sd_obs      = apply(An_obs,1,sd)
  df_avg = xx/length(unique_identifier)
  df_melt = melt(df_avg,id.vars = c("Qin"),variable.name = "variable", value.name = "value")
  df_melt = cbind(df_melt,sd = c(sd_ePhoto,sd_FvCB,sd_obs))
  myplot<-ggplot(df_melt, aes(x = Qin, y = value, color = variable)) +
    geom_line() + geom_point(size=3) +
    scale_color_manual(values=c("green","blue", "black"))+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05)) +
    labs(x = bquote(Q~(mu*mol ~ m^-2 ~ s^-1)),
         y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
    ) +
    theme_bw() +
    theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
# Save ggplot to PDF
  ggsave(output_figure_name, plot = myplot, width = 6, height = 6)
}
