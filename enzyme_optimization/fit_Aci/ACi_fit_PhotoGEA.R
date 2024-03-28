# Clear the workspace
rm(list=ls())
# Load required packages
# library(lattice)
library(PhotoGEA)
library(BioCro)
source('aci_defaults.R')
source('aci_functions.R')

alpha1_or_alpha2 = "alpha2" #fit alpha1 or alpha2

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

if (alpha1_or_alpha2 == "alpha1"){
  for (alpha1  in alpha1_vec){
    #call ephoto c++
    system(paste("./myephoto.exe",alpha1,alpha2))
    ePhoto_result = read.table("output.data",header=FALSE)
    A_Ci_df = data.frame(A=ePhoto_result$V1,Ci=Ci,Tleaf=Tleaf,PAR=PAR)
    output = get_vcmax_jmax(A_Ci_df)
    vcmax = output[1]
    jmax  = output[2]
    tpu   = output[3]
    write(c(alpha1,alpha2,vcmax,jmax,tpu),file="VJ_results/VJ_estimate.txt",append=TRUE)
  }
}else if (alpha1_or_alpha2 == "alpha2"){ 
  for (alpha2  in alpha2_vec){
    #call ephoto c++
    system(paste("./myephoto.exe",alpha1,alpha2))
    ePhoto_result = read.table("output.data",header=FALSE)
    A_Ci_df = data.frame(A=ePhoto_result$V1,Ci=Ci,Tleaf=Tleaf,PAR=PAR)
    output = get_vcmax_jmax(A_Ci_df)
    vcmax = output[1]
    jmax  = output[2]
    tpu   = output[3]
    write(c(alpha1,alpha2,vcmax,jmax,tpu),file="VJ_results/VJ_estimate.txt",append=TRUE)
  }
}

