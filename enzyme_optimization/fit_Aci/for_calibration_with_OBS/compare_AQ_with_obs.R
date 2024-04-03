compare_AQ_with_obs <- function(year,Vcmax25,Jmax25,Rate_TPu,scaling_factors)
{
  source("FarquharModel.R")
  R=8.314472E-3 #Gas constant KJ mole^{-1} K^{-1}
  #read in obs An data
  obs_An  = read.csv(paste0('ACi_AQ_data/AQ_DOY217_',year,'.csv'))
  unique_obs_curves = unique(obs_An$curve_identifier)
  #to run Farquhar's A-Q, we just need run once with an average of Ci and Q
  Ci = as.numeric(obs_An$Ci)
  Ci = matrix(Ci, nrow=7)   # 7 is the number of measurements for a single A-Q
  Cis = rowMeans(Ci)
  
  Q = as.numeric(obs_An$Qin)
  Q = matrix(Q,nrow=7)
  Qps = rowMeans(Q)
  Tls    = mean(obs_An$Tleaf)     #C degree,average leaf temperature from AQ observation on DOY 217
  
  Rd25   = 1.28 
  Air_O2 = 210  #
  
  An_farquhar = NA*(1:length(Qps))
  for (i in 1:length(Qps)){
    Ci = Cis[i]
    Qp = Qps[i]
    LeafTemperature   = Tls #C
    Rd = Rd25 * exp(18.72 - 46.39 / (R * (LeafTemperature + 273.15)))
    
    #For Farquhar, we need the total Q, as there's a calculation of
    #absorption in the Farquhar function
    output_farquhar  = FarquharModel(LeafTemperature, Ci, Qp, Air_O2,Vcmax25, Jmax25, Rate_TPu,scaling_factors)
    
    if(length(output_farquhar)>1){
       An_farquhar[i]   = output_farquhar$GA - Rd;
    }
  }
  
  df_for_plot = cbind(A = An_farquhar,Qin = Qps, Ci=Cis, Tleaf=Tls, curve_identifier="Farquhar")
  
  df_for_plot     = rbind(obs_An,df_for_plot)
  df_for_plot$A   = as.numeric(df_for_plot$A)
  df_for_plot$Qin = as.numeric(df_for_plot$Qin)

  manual_colors = c(rep("blue",length(unique_obs_curves)),"red")
  if(year==2022){
    manual_colors = c("red",rep("blue",length(unique_obs_curves)))
  }
  
  p<-ggplot(data=df_for_plot,aes(x=Qin, y=A,colour=curve_identifier)) +
  geom_line()  +  geom_point() +
  scale_color_manual(values=manual_colors) +
  theme_bw() +
  theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
  return(p)
}
