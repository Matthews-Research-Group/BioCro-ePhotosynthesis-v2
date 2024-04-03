#This is for fitting the A-Ci curves of ePhotosynthesis
#to match target Vcmax and Jmax
# Clear the workspace
rm(list=ls())
# Load required packages
# library(lattice)
library(PhotoGEA)
library(BioCro)
source('for_calibration_with_OBS/aci_defaults.R') #get_jmax function
source('aci_functions.R')

alpha1_or_alpha2 = "alpha2" #fit alpha1 or alpha2
curve_option = 1 #1: A-Ci;2: A-Q; 3: A-T 

if (alpha1_or_alpha2 == "alpha1"){
  alpha1_vec = seq(0.8,0.9,by=0.02)
  alpha2 = 1
}else if (alpha1_or_alpha2 == "alpha2"){
  alpha1 = 0.87 
  alpha2_vec = seq(0.96,1.0,by=0.01) 
}else{
  stop('no such option')
}

Ci = c(100, 150, 200, 250, 300, 400, 500, 600, 800, 1200) 
Tleaf = 25
PAR = 2000

#note: currently, the PAR, Ci, and Tleaf for myephoto.exe are hard-coded
#double check whether it's high or low light condition

if (alpha1_or_alpha2 == "alpha1"){
  for (alpha1  in alpha1_vec){
    #call ephoto c++
    system(paste("./myephoto.exe",alpha1,alpha2,curve_option))
    ePhoto_result = read.table("output.data",header=FALSE,sep = ",")
    A_Ci_df = data.frame(A=ePhoto_result$V4,Ci=Ci,Tleaf=Tleaf,PAR=PAR)
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
    A_Ci_df = data.frame(A=ePhoto_result$V4,Ci=Ci,Tleaf=Tleaf,PAR=PAR)
    output = get_vcmax_jmax(A_Ci_df)
    vcmax = output[1]
    jmax  = output[2]
    tpu   = output[3]
    write(c(alpha1,alpha2,vcmax,jmax,tpu),file="VJ_results/VJ_estimate.txt",append=TRUE)
  }
}

