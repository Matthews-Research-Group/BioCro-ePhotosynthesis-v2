library(reshape2)
library(ggplot2)
library(BioCro)
source("for_calibration_with_OBS/my_scripts/biocro_FvCB.R")
alpha1_alpha2 = read.csv('ePhotosynthesis_optimal_alpha1_alpha2.csv')
alpha1 = alpha1_alpha2[2]
alpha2 = alpha1_alpha2[3]
curve_option = 3 #1: A-Ci;2: A-Q; 3: A-T 
prefix = c("ACi","AQ","AT")
PAR = 400
output_figure_name = paste0("figs/",prefix[curve_option],"_Q",PAR,"_v1.pdf")
#call ephoto c++
system(paste("./myephoto.exe",alpha1,alpha2,PAR,curve_option))
#read in ephoto results
ephoto = read.csv("output.data",header=FALSE)
colnames(ephoto) = c("PAR","Tleaf","Ci","An")

aci_fit_results <- read.csv("for_calibration_with_OBS/my_scripts/ld11_aci_fit_parameters_avg.csv")
Vcmax25      <- aci_fit_results$mean[aci_fit_results$parameter=="Vcmax_at_25"]
Jmax25       <- aci_fit_results$mean[aci_fit_results$parameter=="Jmax_at_25"]
Rd25         <- 1.28
TPU25        <- aci_fit_results$mean[aci_fit_results$parameter=="TPU_at_25"]
#run Farquhar
An_farquhar = NA*(1:dim(ephoto)[1])

for (i in 1:dim(ephoto)[1]){
  #For Farquhar, we need the total Q, as there's a calculation of
  #absorption inside the Farquhar function
  output_farquhar  = BioCro_FvCB(ephoto$PAR[i],ephoto$Tleaf[i], ephoto$Ci[i], Vcmax25, Jmax25, Rd25, TPU25)
  An_farquhar[i]   = output_farquhar$An
  #get ephoto's An
  Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, ephoto$Tleaf[i]+273.15)
  ephoto$An[i] = ephoto$An[i] - Rd
}

df_for_plot = cbind(ephoto,An_FvCB=An_farquhar)
colnames(df_for_plot)[colnames(df_for_plot)=="An"] = "An_ePhoto"
df_for_plot = melt(df_for_plot, id.vars=c("PAR","Tleaf","Ci"))

if(curve_option==1){
  myplot<-ggplot(data=df_for_plot,aes(x=Ci, y=value,colour=variable)) +
       geom_line()  +  geom_point()
}else if(curve_option==2){
  myplot<-ggplot(data=df_for_plot,aes(x=PAR, y=value,colour=variable)) +
       geom_line()  +  geom_point()
}else{
  myplot<-ggplot(data=df_for_plot,aes(x=Tleaf, y=value,colour=variable)) +
       geom_line()  +  geom_point()
}
# Save ggplot to PDF
ggsave(output_figure_name, plot = myplot, width = 6, height = 4)

# df = cbind(Tleaf= ephoto$Tleaf,An = An_farquhar,Aj = Aj_farquhar, Ac = Ac_farquhar, Ap = Ap_farquhar)
# df = as.data.frame(df)
# df = melt(df,id.vars=c("Tleaf"))
# ggplot(data=df,aes(x=Tleaf, y=value,colour=variable)) +
#   geom_line()  +  geom_point()

# k1 = -3E-05
# k2 = 0.0013
# k3=0.0106
# k4=0.8839
# Ru_Act = k1*ephoto$Tleaf^3 + k2*ephoto$Tleaf^2 - k3*ephoto$Tleaf + k4
# plot(ephoto$Tleaf,Ru_Act,type='b')
# y = 2^((ephoto$Tleaf-25)/10)
# plot(ephoto$Tleaf,y,type='b')
# y1 = 1.5^((ephoto$Tleaf-25)/10)
# lines(ephoto$Tleaf,y1,type='b',col='blue')
# y2 = 2.5^((ephoto$Tleaf-25)/10)
# lines(ephoto$Tleaf,y2,type='b',col='red')

# 
# Ci=400
# a1 = 4.5   #electrons_per_carboxylation per C
# a2 = 10.5  #electrons_per_carboxylation per O2
# LeafTemperatureKelvin = ephoto$Tleaf + 273.15
# R=8.31446261815324e-3
# GammaStar = exp(19.02 - 37.83 / (R * LeafTemperatureKelvin))
# beta = (1.0 - GammaStar / Ci) *  Ci /(a1 * Ci + a2 * GammaStar)
# plot(ephoto$Tleaf,beta,type='b')

