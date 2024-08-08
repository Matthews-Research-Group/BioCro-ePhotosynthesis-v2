library(BioCro)
library(PhotoGEA)
library(nloptr)
source("for_calibration_with_OBS/my_scripts/biocro_FvCB.R")
source("obj_functions.R")
curve_type = 1 #1:ACi; 2: AQ
prefix = c("ACi","AQ")
outfile_prefix = prefix[curve_type]

all_enzymes = c( 'Q10_1',  'Q10_2', 'Q10_3',
                 'Q10_5',  'Q10_6', 'Q10_7',
                 'Q10_8',  'Q10_9', 'Q10_10',
                 'Q10_13', 'Q10_23')
#the default values of these enzymes
targets_all = c(1.93,rep(2.0,length(all_enzymes)-1))
names(targets_all) = all_enzymes
 
enzymes_to_optimize = c( 'Q10_1',  'Q10_2', 'Q10_3',
                         'Q10_5',  'Q10_6', 'Q10_7',
                         'Q10_8',  'Q10_9', 'Q10_10',
                         'Q10_13', 'Q10_23')
number_of_enzymes = length(enzymes_to_optimize) #total number of Q10 to be optimized
#the default values of these enzymes
targets = targets_all[all_enzymes%in%enzymes_to_optimize]
names(targets) = enzymes_to_optimize

executable_path <- "./myephoto_single.exe"

alpha1_alpha2 = read.csv('ePhotosynthesis_optimal_alpha1_alpha2.csv')
alpha1 = alpha1_alpha2[2]
alpha2 = alpha1_alpha2[3]

aci_fit_results <- read.csv("for_calibration_with_OBS/my_scripts/ld11_aci_fit_parameters_avg.csv")
Vcmax25      <- aci_fit_results$mean[aci_fit_results$parameter=="Vcmax_at_25"]
Jmax25       <- aci_fit_results$mean[aci_fit_results$parameter=="Jmax_at_25"]
Rd25         <- 1.28
TPU25        <- aci_fit_results$mean[aci_fit_results$parameter=="TPU_at_25"]

if(curve_type==1){
  licor_data <- read.csv.exdf("for_calibration_with_OBS/my_scripts/ld11_aci.csv") 
  
  Tleaf_all  <- licor_data[, 'TleafCnd'] 
  Qin_all    <- licor_data[, 'Qin']
  Ci_all     <- licor_data[, 'Ci']
  
  date = substr(licor_data[, 'date'], 1, 8)
  curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],sep='-')
  An_obs = licor_data[,'A']

}else if(curve_type==2){
  licor_data <- read.csv.exdf("for_calibration_with_OBS/my_scripts/ld11_bb.csv") #this comes from "get_ld11_bb.R"

  Tleaf_all  <- licor_data[, 'TleafCnd'] 
  Qin_all    <- licor_data[, 'Qin']
  Ci_all     <- licor_data[, 'Ci']

  date = substr(licor_data[, 'date'], 1, 8)
  
  curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],licor_data[, 'replicate'],sep='-')
  
  An_obs = licor_data[,'A']
}


lbs = rep(1.5,number_of_enzymes) 
ubs = rep(2.5,number_of_enzymes)
init_guess = rep(2.0,length(lbs))

#You can set lambda to 0 to NOT use LASSO
lambda = 0.1

constants_list = list(curve_type    = curve_type,
                 alpha1        = alpha1,
                 alpha2        = alpha2,
                 Vcmax25       = Vcmax25,
                 Jmax25        = Jmax25,
                 Rd25          = Rd25,
                 TPU25         = TPU25,
                 Tleaf_all     = Tleaf_all,
                 Qin_all       = Qin_all,
                 Ci_all        = Ci_all,
                 curve_identifier = curve_identifier,
                 An_obs        = An_obs,
                 targets_all   = targets_all,
                 targets       = targets
                 )
# Define the objective function for nloptr, including constant parameters
objective_function <- function(params) {
  lasso_obj_function(params, constants_list, lambda)
}

opt_results = nloptr(x0 = init_guess, eval_f = objective_function, 
                     lb = lbs,
                     ub = ubs,
                     opts = list("algorithm" = "NLOPT_LN_SBPLX","xtol_rel"=1.0e-4,"print_level"= 1)
                     )

Q10s  = as.data.frame(opt_results$solution)
rownames(Q10s) = names(targets)
Q10s  = t(Q10s)
write.csv(Q10s,paste0("Q10_fitting_results/ePhotosynthesis_optimal_Q10s_",outfile_prefix,"_v5_LASSO.csv"))

if(FALSE){
q10=2
tleaf = seq(5,40,by=5)
y = q10^((tleaf-25)/10)
plot(tleaf,y,type='l',ylab='Q10 multiplier')
q10 = 1.5
y = q10^((tleaf-25)/10)
lines(tleaf,y,col='red')
q10 = 2.5
y = q10^((tleaf-25)/10)
lines(tleaf,y,col='blue')
legend(5,2,legend = c("Q10=2","Q10=1.5","Q10=2.5"),col=c('black','red','blue'),lty = 1)
}
