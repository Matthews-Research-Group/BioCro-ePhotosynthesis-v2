source('aci_and_aq_fitting_functions.R')
year = 2022
iteration_time = 10
#you may check what these scaling factors are applying to in 
#aci_defaults.R
#they are related to scale the PhiPS2 and Theta in get_jmax
# Make sure the the scaling factors usage is consistent in FarquharModel.R and aci_defaults.R !!
# either 2 or 4 parameters to optimize
# length 2: fit to bases only; length 4: fit to Theta base, and all 3 coefs of PhiPS2
# orders: PhiPS2 base, Theta base, linear of PhiPS2, qudratic of PhiPS2
scaling_factors = c(1,1) #initial value; 
lbs = c(0.5,0.5)
ubs = c(2.0,2.0)
fitting_results = c() 
for (i in 1:iteration_time){
  if(year==2021){
    aci = aci_fitting_2021(scaling_factors)
  }else if(year==2022){
    aci = aci_fitting_2022(scaling_factors)
    aci = as.data.frame(aci)
    aci = aci[aci$cultivar=='ld11',]
  }
  aq  = fit_AQ(aci$Vcmax, aci$Jmax,aci$TPU,year,lbs,ubs)
  aq = as.data.frame(matrix(aq,1,length(aq)))
  #the fit_AQ return a vector
  #the last value of the vector is the rmse of AQ fitting
  #the first couple values are the scaling factors that are applied to the bases of PhiPS2 and Theta (order matters)
  #check FarquharModel.R to confirm what these scaling factors do exactly
  
  #save the results for each iteration
  tmp = data.frame(step = i,year = year,Vcmax = aci$Vcmax, Jmax= aci$Jmax,TPU=aci$TPU)
  tmp = cbind(tmp,aq)
  fitting_results = rbind(fitting_results,tmp)

  #update the scaling_factors for next iteration
  scaling_factors = aq[1:length(lbs)]
}

#plot1: parameter evolution
make_plot1<-function(){
  par(mfrow=c(2,2),mar=c(4,4,3,3),mai=c(1,1,0.5,0.5),cex=1.5) #c(bottom, left, top, right)
  linewidth=2
  axislabelspace = 2
  plot(fitting_results$step, fitting_results$Jmax,type='b',
         lwd = linewidth,
         xlab='',ylab='')
  title(xlab = "iterations", line = axislabelspace)           # Add x-axis text
  title(ylab = "Jmax", line = axislabelspace)                 # Add y-axis text
  
  plot(fitting_results$step, fitting_results$V1,type='b',
       lwd = linewidth,
       xlab='',ylab='')
  title(xlab = "iterations", line = axislabelspace)           # Add x-axis text
  title(ylab = "sf (PhiPS2 base)", line = axislabelspace)                 # Add y-axis text
  
  plot(fitting_results$step, fitting_results$V2,type='b',
       lwd = linewidth,
       xlab='',ylab='')
  title(xlab = "iterations", line = axislabelspace)           # Add x-axis text
  title(ylab = "sf (Theta base)", line = axislabelspace)                 # Add y-axis text
  
  plot(fitting_results$step, fitting_results[,ncol(fitting_results)],type='b',
       lwd = linewidth,
       xlab='',ylab='')
  title(xlab = "iterations", line = axislabelspace)           # Add x-axis text
  title(ylab = "RMSE", line = axislabelspace)                 # Add y-axis text
}
make_plot1()

#plot2: compare new and old PhiPSII vs temperature curves
library(ggplot2)
make_plot2<-function(){
  PhiPSII_func<-function(LeafTemperature,sf) {
    return(0.352 *sf[1]+ 0.022 * sf[2]* LeafTemperature - sf[3]* 3.4*1e-4*LeafTemperature^2.0)
  }
  leafT = seq(0,40,by=5)
  phiPS2_control = PhiPSII_func(leafT,c(1,1,1))
  last_result = fitting_results[iteration_time,]
  phiPS2_new  = PhiPSII_func(leafT,c(last_result$V1,1,1))
  if(length(lbs)==4){
    phiPS2_new  = PhiPSII_func(leafT,c(last_result$V1,last_result$V3,last_result$V4))
  }
  df = data.frame(leafT,phiPS2_control,phiPS2_new)
  df_for_plot = reshape2:: melt(df,id.vars="leafT")
  myplot<-ggplot(data=df_for_plot,aes(x=leafT, y=value,colour=variable)) +
    geom_line()  +  geom_point() + ylab("PhiPS2")
  myplot
}
make_plot2()

#plot3: compare ACi and AQ curves between OBS and newly calibrated Farquhar.
last_result = fitting_results[iteration_time,]
scaling_factors = c(last_result[6:(ncol(last_result)-1)]) #exclude the last value, which is RMSE
source('compare_ACi_with_obs.R')
p1=compare_ACi_with_obs(year,last_result$Vcmax,last_result$Jmax,last_result$TPU,scaling_factors)
source('compare_AQ_with_obs.R')
p2=compare_AQ_with_obs(year,last_result$Vcmax,last_result$Jmax,last_result$TPU,scaling_factors)
require(gridExtra)
grid.arrange(p1, p2, ncol=2)

