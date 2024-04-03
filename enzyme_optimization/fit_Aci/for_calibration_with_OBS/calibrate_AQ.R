library(nloptr)
library(ggplot2)
source("FarquharModel.R")
PhiPSII_func<-function(LeafTemperature,sf) {
  return(0.352 *sf[1]+ 0.022 * sf[2]* LeafTemperature - sf[3]* 3.4*1e-4*LeafTemperature^2.0)
}
theta_func<-function(LeafTemperature) {
  PhotosynthesisTheta = 0.76
  return(PhotosynthesisTheta + 0.01713 * LeafTemperature - 3.75 *LeafTemperature^2.0 / 10000.0)
}
rmse<-function(model, obs){
  return(sqrt(sum((model-obs)^2)))
}
R=8.314472E-3 #Gas constant KJ mole^{-1} K^{-1}
#read in obs An data
obs_An  = read.csv('ACi_AQ_data/AQ_DOY217_2021.csv')

Rd25   = 1.28 
Air_O2 = 210  #

Vcmax25  = 127   #LD11
Jmax25   = 197   #LD11
Rate_TPu = 15  #LD11

init_guess = c(1,1)
best_par = c()
#here we fit to each replicate separately
for (curve_id in unique(obs_An$curve_identifier)){
  obs_An_sub = obs_An[obs_An$curve_identifier%in%curve_id,]
  Cis = obs_An_sub$Ci
  Qps = obs_An_sub$Qin
  Tls = obs_An_sub$Tleaf
#define obj function here  
  obj_func<-function(x){
    An_farquhar = NA*(1:length(Qps))
    for (i in 1:length(Qps)){
      Ci = Cis[i]
      Qp = Qps[i]
      LeafTemperature   = Tls[i] #C
      Rd = Rd25 * exp(18.72 - 46.39 / (R * (LeafTemperature + 273.15)))
      
      output_farquhar  = FarquharModel(LeafTemperature, Ci, Qp, Air_O2,Vcmax25, Jmax25, Rate_TPu,x)
      
      if(length(output_farquhar)>1){
        An_farquhar[i]   = output_farquhar$GA - Rd;
      }
    }
    if(any(is.na(An_farquhar))){
       return(9999)
    }else{
      return(rmse(An_farquhar,obs_An_sub$A))
    }
  }
  # opt_results = hjk(init_guess, obj_func)
  # best_par = rbind(best_par,c(opt_results$par,opt_results$value))
  opt_results = nloptr(x0 = init_guess, eval_f = obj_func, 
                       lb = c(0.5,0.5),
                       ub = c(2.0,2.0),
                       opts = list("algorithm" = "NLOPT_LN_SBPLX","xtol_rel"=1.0e-6))
  best_par = rbind(best_par,c(opt_results$solution,opt_results$objective))
}
best_par_mean = colMeans(best_par)


# 
# df_for_plot = cbind(A = An_farquhar,Qin = obs_mean$Q_mean,curve_identifier="Farquhar",group=5)
# 
# df_for_plot = rbind(obs_An,df_for_plot)
# df_for_plot$A = as.numeric(df_for_plot$A)
# df_for_plot$Qin = as.numeric(df_for_plot$Qin)
# 
# ggplot(data=df_for_plot,aes(x=Qin, y=A,colour=curve_identifier)) +
#   geom_line()  +  geom_point()
leafT = seq(0,40,by=5)
phiPS2_control = PhiPSII_func(leafT,c(1,1,1))
phiPS2_new  = PhiPSII_func(leafT,c(best_par_mean[1:3]))
df = data.frame(leafT,phiPS2_control,phiPS2_new)
df_for_plot = reshape2:: melt(df,id.vars="leafT")
myplot<-ggplot(data=df_for_plot,aes(x=leafT, y=value,colour=variable)) +
  geom_line()  +  geom_point() + ylab("PhiPS2")
myplot
