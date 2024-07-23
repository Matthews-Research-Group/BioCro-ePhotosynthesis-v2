#This is for fitting the A-Ci curves of ePhotosynthesis
#to match target Vcmax and Jmax
# Clear the workspace
rm(list=ls())
# Load required packages
# library(lattice)
library(nloptr)
library(PhotoGEA)
library(BioCro)
source('for_calibration_with_OBS/my_scripts/aci_defaults.R') #get_jmax function
source('for_calibration_with_OBS/my_scripts/biocro_FvCB.R') #Vcmax_multiplier function
source('aci_functions.R')

obj_func<-function(x){
  # Specify FvCB parameter values
  aci_fit_results <- read.csv("for_calibration_with_OBS/my_scripts/ld11_aci_fit_parameters_avg.csv")
  vcmax_target <- aci_fit_results$mean[aci_fit_results$parameter=="Vcmax_at_25"]
  jmax_target  <- aci_fit_results$mean[aci_fit_results$parameter=="Jmax_at_25"]
  Rd_at_25     <- 1.28
  tpu_target   <- aci_fit_results$mean[aci_fit_results$parameter=="TPU_at_25"]
  PAR = 2000
  system(paste("./myephoto.exe",x[1],x[2],PAR,1))
  ePhoto_result = read.table("output.data",header=FALSE,sep = ",")
  A_Ci_df = as.data.frame(ePhoto_result)
  colnames(A_Ci_df) = c("PAR","Tleaf","Ci","A")
  #ephoto's assimilation is the gross,which needs to substract Rd to get An
  Rd = Rd_at_25 * arrhenius_exponential(18.72, 46.39e3, A_Ci_df$Tleaf+273.15) 
  A_Ci_df$A = A_Ci_df$A - Rd 
  output = get_vcmax_jmax(A_Ci_df)
  vcmax = output[1]
  jmax  = output[2]
  tpu   = output[3]
  
  error = (vcmax - vcmax_target)^2 + (jmax - jmax_target)^2 #+ (tpu - tpu_target)^2
  if(is.na(error)) error=999
  return(error)
}

lbs = c(0.8,0.9) #alpha1, alpha2 
ubs = c(1.0,1.2)
init_guess = rep(1,length(lbs))

opt_results = nloptr(x0 = init_guess, eval_f = obj_func, 
                     lb = lbs,
                     ub = ubs,
                     opts = list("algorithm" = "NLOPT_LN_SBPLX","xtol_rel"=1.0e-6))
alpha1_alpha2 = as.data.frame(opt_results$solution)
rownames(alpha1_alpha2) = c('alpha1','alpha2')
alpha1_alpha2 = t(alpha1_alpha2)
write.csv(alpha1_alpha2,"ePhotosynthesis_optimal_alpha1_alpha2.csv")
