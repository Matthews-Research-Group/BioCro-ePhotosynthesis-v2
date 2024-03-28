library(reshape2)
library(ggplot2)
source("FarquharModel.R")
alpha1 = 0.87
alpha2 = 0.98
curve_option = 3 #1: A-Ci;2: A-Q; 3: A-T 
#call ephoto c++
system(paste("./myephoto.exe",alpha1,alpha2,curve_option))
#read in ephoto results
ephoto = read.csv("output.data",header=FALSE)
colnames(ephoto) = c("PAR","Tleaf","Ci","An")

#run Farquhar
An_farquhar = NA*(1:dim(ephoto)[1])
for (i in 1:dim(ephoto)[1]){
  R=8.314472E-3 #Gas constant KJ mole^{-1} K^{-1}
  Rd25   = 1.28 
  Air_O2 = 210  #
  
  Vcmax25  = 123   #LD11
  Jmax25   = 196.5   #LD11
  Rate_TPu = 14.6  #LD11
  scaling_factors = c(1.735, 0.835)
  #For Farquhar, we need the total Q, as there's a calculation of
  #absorption in the Farquhar function
  output_farquhar  = FarquharModel(ephoto$Tleaf[i], ephoto$Ci[i], ephoto$PAR[i], Air_O2, Vcmax25, Jmax25, Rate_TPu,scaling_factors)
  
  Rd = Rd25 * exp(18.72 - 46.39 / (R * (ephoto$Tleaf[i]+ 273.15)))
  An_farquhar[i]   = output_farquhar$GA - Rd;
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
ggsave("figs/AT_Q400_r1.pdf", plot = myplot, width = 6, height = 4)
