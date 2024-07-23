#This is for fitting the A-Ci curves of ePhotosynthesis
#to match target Vcmax and Jmax
# Clear the workspace
rm(list=ls())
# Load required packages
# library(lattice)
library(PhotoGEA)
library(BioCro)
source('for_calibration_with_OBS/my_scripts/aci_defaults.R') #get_jmax function
source('for_calibration_with_OBS/my_scripts/biocro_FvCB.R') #Vcmax_multiplier function
source('aci_functions.R')

alpha1_or_alpha2 = "alpha2" #fit alpha1 or alpha2
curve_option = 1 #1: A-Ci;2: A-Q; 3: A-T 

if (alpha1_or_alpha2 == "alpha1"){
  alpha1_vec = seq(0.8,0.9,by=0.02)
  alpha2 = 1
}else if (alpha1_or_alpha2 == "alpha2"){
  alpha1 = 0.89
  alpha2_vec = seq(1.06,1.1,by=0.01) 
}else{
  stop('no such option')
}

#note: currently, the PAR, Ci, and Tleaf for myephoto.exe are hard-coded
#double check whether it's high or low light condition

if (alpha1_or_alpha2 == "alpha1"){
  for (alpha1  in alpha1_vec){
    #call ephoto c++
    system(paste("./myephoto.exe",alpha1,alpha2,curve_option))
    ePhoto_result = read.table("output.data",header=FALSE,sep = ",")
    A_Ci_df = as.data.frame(ePhoto_result)
    colnames(A_Ci_df) = c("PAR","Tleaf","Ci","A")
    # A_Ci_df$A1 = BioCro_FvCB(2000,25,A_Ci_df$Ci,134.5,183,1.28,10.5)
    # A_Ci_df$A = A_Ci_df$A1$An
    # A_Ci_df = A_Ci_df[,1:4]
    output = get_vcmax_jmax(A_Ci_df)
    vcmax = output[1]
    jmax  = output[2]
    tpu   = output[3]
    write(c(alpha1,alpha2,vcmax,jmax,tpu),file="VJ_results/VJ_estimate.txt",append=TRUE)
  }
}else if (alpha1_or_alpha2 == "alpha2"){ 
  for (alpha2  in alpha2_vec){
    #call ephoto c++
    system(paste("./myephoto.exe",alpha1,alpha2,curve_option))
    ePhoto_result = read.table("output.data",header=FALSE,sep = ",")
    A_Ci_df = as.data.frame(ePhoto_result)
    colnames(A_Ci_df) = c("PAR","Tleaf","Ci","A")
    output = get_vcmax_jmax(A_Ci_df)
    vcmax = output[1]
    jmax  = output[2]
    tpu   = output[3]
    write(c(alpha1,alpha2,vcmax,jmax,tpu),file="VJ_results/VJ_estimate.txt",append=TRUE)
  }
}

