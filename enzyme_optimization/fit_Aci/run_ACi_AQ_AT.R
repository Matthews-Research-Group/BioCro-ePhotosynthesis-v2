library(reshape2)
library(ggplot2)
source("for_calibration_with_OBS/FarquharModel.R")
alpha1 = 0.87
alpha2 = 0.98
curve_option = 3 #1: A-Ci;2: A-Q; 3: A-T 
output_figure_name = "figs/AT_Q2000_r4.pdf"
#call ephoto c++
system(paste("./myephoto.exe",alpha1,alpha2,curve_option))
#read in ephoto results
ephoto = read.csv("output.data",header=FALSE)
colnames(ephoto) = c("PAR","Tleaf","Ci","An")

#run Farquhar
An_farquhar = NA*(1:dim(ephoto)[1])
Aj_farquhar = An_farquhar
Ac_farquhar = An_farquhar
Ap_farquhar = An_farquhar
for (i in 1:dim(ephoto)[1]){
  R=8.314472E-3 #Gas constant KJ mole^{-1} K^{-1}
  Rd25   = 1.28 
  Air_O2 = 210  #
  
  Vcmax25  = 123   #LD11
  Jmax25   = 216   #LD11
  Rate_TPu = 14.6  #LD11
  scaling_factors = c(0.72, 0.64)
  #For Farquhar, we need the total Q, as there's a calculation of
  #absorption in the Farquhar function
  output_farquhar  = FarquharModel(ephoto$Tleaf[i], ephoto$Ci[i], ephoto$PAR[i], Air_O2, Vcmax25, Jmax25, Rate_TPu,scaling_factors)
  
  Rd = Rd25 * exp(18.72 - 46.39 / (R * (ephoto$Tleaf[i]+ 273.15)))
  An_farquhar[i]   = output_farquhar$GA - Rd
  Aj_farquhar[i]   = output_farquhar$Aj
  Ac_farquhar[i]   = output_farquhar$Ac
  Ap_farquhar[i]   = output_farquhar$Ap
}

df_for_plot = cbind(ephoto,An_Farquhar=An_farquhar)
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

