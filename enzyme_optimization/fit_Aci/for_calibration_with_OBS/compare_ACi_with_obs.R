compare_ACi_with_obs <- function(year,Vcmax25,Jmax25,Rate_TPu,scaling_factors)
{
  source("FarquharModel.R")
  R=8.314472E-3 #Gas constant KJ mole^{-1} K^{-1}
  #read in obs An data
  obs_An  = read.csv(paste0('ACi_AQ_data/ACi_DOY217_',year,'.csv'))
  unique_obs_curves = unique(obs_An$curve_identifier)
  
  Cis = seq(100,1500,by=100)
  Qps = 2000
  
  Rd25   = 1.28 
  Tls    = 30.4     #C degree, average leaf temperature from ACi observation
  Air_O2 = 210  #
  
  An_farquhar = NA*(1:length(Cis))
  for (i in 1:length(Cis)){
    Ci = Cis[i]
    Qp = Qps
    LeafTemperature   = Tls #C
    Rd = Rd25 * exp(18.72 - 46.39 / (R * (LeafTemperature + 273.15)))
    
    #For Farquhar, we need the total Q, as there's a calculation of
    #absorption in the Farquhar function
    output_farquhar  = FarquharModel(LeafTemperature, Ci, Qp, Air_O2,Vcmax25, Jmax25, Rate_TPu,scaling_factors)
    An_farquhar[i]   = output_farquhar$GA - Rd;
  }
  
  df_for_plot = cbind(A = An_farquhar,Ci = Cis,Tleaf=Tls,Qin=Qps,curve_identifier="Farquhar")
  
  df_for_plot = rbind(obs_An,df_for_plot)
  df_for_plot$A = as.numeric(df_for_plot$A)
  df_for_plot$Ci = as.numeric(df_for_plot$Ci)
  
  p <- ggplot(data=df_for_plot,aes(x=Ci, y=A,colour=curve_identifier)) +
  geom_line()  +  geom_point() + 
  scale_color_manual(values=c("red", rep("blue",length(unique_obs_curves))))+
  theme_bw() +
  theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
  return(p)
}
