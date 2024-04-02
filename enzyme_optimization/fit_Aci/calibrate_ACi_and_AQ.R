source('aci_and_aq_fitting_functions.R')
year = 2021
iteration_time = 10
#you may check what these scaling factors are applying to in 
#aci_defaults.R
#they are related to scale the PhiPS2 and Theta in get_jmax
# Make sure the the scaling factors usage is consistent in FarquharModel.R and aci_defaults.R !!
scaling_factors = c(1,1) #initial value
fitting_results = c() 
for (i in 1:iteration_time){
  if(year==2021){
    aci = aci_fitting_2021(scaling_factors)
  }else if(year==2022){
    aci = aci_fitting_2022(scaling_factors)
  }
  aq  = fit_AQ(aci$Vcmax, aci$Jmax,aci$TPU,year)
  #the fit_AQ return a vector
  #the last value of the vector is the rmse of AQ fitting
  #the first couple values are the scaling factors that are applied to the bases of PhiPS2 and Theta (order matters)
  #check FarquharModel.R to confirm what these scaling factors do exactly
  
  #save the results for each iteration
  tmp = data.frame(step = i,year = year,Vcmax = aci$Vcmax, Jmax= aci$Jmax,TPU=aci$TPU,
                   sf1 = aq[1],sf2=aq[2],rmse=aq[3] ) 
  fitting_results = rbind(fitting_results,tmp)

  #update the scaling_factors for next iteration
  scaling_factors = aq[1:2]
}
