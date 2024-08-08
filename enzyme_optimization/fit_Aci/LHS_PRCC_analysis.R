# Load the lhs package
library(lhs)
library(epiR)
rm(list=ls())
number_of_enzymes = 11 #total number of Q10 to be optimized
executable_path <- "./myephoto_single.exe"
# Define the number of samples and the number of variables
n_samples   <- 100  # Number of samples
n_variables <- number_of_enzymes # Number of variables
lower_bound <- 1.5
upper_bound <- 2.5

set.seed(1234)
# Generate the Latin Hypercube Sample
lhs_samples <- randomLHS(n_samples, n_variables)

# Transform samples to the desired range
shifted_samples <- lower_bound + (upper_bound - lower_bound) * lhs_samples

model_run<-function(params){
  source("for_calibration_with_OBS/my_scripts/biocro_FvCB.R")
  source("obj_functions.R")
  curve_type = 1 #1:ACi; 2: AQ
  
  alpha1_alpha2 = read.csv('ePhotosynthesis_optimal_alpha1_alpha2.csv')
  alpha1 = alpha1_alpha2[2]
  alpha2 = alpha1_alpha2[3]
  
  Rd25  = 1.28
  Tleaf = 30
  PAR   = 2000
  Ci    = 600
  
  constants_list = list(curve_type    = curve_type,
                        alpha1        = alpha1,
                        alpha2        = alpha2,
                        Rd25          = Rd25,
                        Tleaf         = Tleaf,
                        PAR           = PAR,
                        Ci            = Ci
  )

  result<-An_obj_function(params, constants_list)
  return(result)

}



model_output <- apply(shifted_samples,1,function(x) model_run(x))

dat = cbind(shifted_samples,model_output)
# Perform PRCC analysis
prcc_results <-  epi.prcc(dat,sided.test = 2, conf.level = 0.95)
importance = prcc_results$est
names(importance) = c("Q10_1",
                      "Q10_2",
                      "Q10_3",
                      "Q10_5",
                      "Q10_6",
                      "Q10_7",
                      "Q10_8",
                      "Q10_9",
                      "Q10_10",
                      "Q10_13",
                      "Q10_23")
bp <- barplot(importance,horiz = FALSE,ylim=c(-1,1),
              ylab = "Partial rank correlation coefficient",
              #main=paste(no_of_samples,"samples"),
              cex.axis = 1.2,cex.names=1.2,cex.lab=1.2)
abline(h=0)
