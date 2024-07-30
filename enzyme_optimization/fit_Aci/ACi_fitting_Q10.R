library(BioCro)
library(PhotoGEA)
source("for_calibration_with_OBS/my_scripts/biocro_FvCB.R")
source("obj_functions.R")
curve_type = 1 #1:ACi; 2: AQ
prefix = c("ACi","AQ")
outfile_prefix = prefix[curve_type]

number_of_enzymes = 11 #total number of Q10 to be optimized

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

  An_FvCB     = NA* 1:length(Tleaf_all)
  An_ePhoto   = NA* 1:length(Tleaf_all)
   
  for (i in 1:length(Tleaf_all)){ 
    PAR = Qin_all[i]
    Tleaf = Tleaf_all[i]
    Ci  = Ci_all[i]
    #For Farquhar, we need the total Q, as there's a calculation of
    #absorption inside the Farquhar function
    output_farquhar  = BioCro_FvCB(PAR,Tleaf, Ci, Vcmax25, Jmax25, Rd25, TPU25)
    An_FvCB[i]   = output_farquhar$An
  
    #call ephoto c++
    args = c(alpha1,alpha2,PAR,Tleaf, Ci)
    ephoto<-system2(executable_path, args = args, stdout = TRUE, stderr = TRUE)
    ephoto =  as.numeric(ephoto)
    Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
    An_ePhoto[i] = ephoto - Rd
  }
  An_ePhoto = as.numeric(An_ePhoto)

  df = data.frame(Qin=Qin_all,An_ePhoto,An_FvCB,An_obs,curve_identifier)

  #plot average line
  unique_identifier = unique(df$curve_identifier)
  xx = 0
  An_ePhoto  = c()
  An_FvCB    = c()
  An_obs     = c()
  for (id in unique_identifier){
    tmp = df[df$curve_identifier == id,1:4]
    xx = xx + tmp
    An_ePhoto = cbind(An_ePhoto,tmp$An_ePhoto)
    An_FvCB   = cbind(An_FvCB  ,tmp$An_FvCB)
    An_obs    = cbind(An_obs   ,tmp$An_obs)
  }
  sd_ePhoto   = apply(An_ePhoto,1,sd)
  sd_FvCB     = apply(An_FvCB,1,sd)
  sd_obs      = apply(An_obs,1,sd)
  df_avg = xx/length(unique_identifier)
  df_melt = melt(df_avg,id.vars = c("Qin"),variable.name = "variable", value.name = "value")
  df_melt = cbind(df_melt,sd = c(sd_ePhoto,sd_FvCB,sd_obs))
}


lbs = rep(1.0,number_of_enzymes) 
ubs = rep(3.0,number_of_enzymes)
init_guess = rep(2.0,length(lbs))

opt_results = nloptr(x0 = init_guess, eval_f = obj_Q10, 
                     lb = lbs,
                     ub = ubs,
                     opts = list("algorithm" = "NLOPT_LN_SBPLX","xtol_rel"=1.0e-4,"print_level"=0),
                     curve_type =  curve_type,
                     alpha1_alpha2 = alpha1_alpha2,
                     Vcmax25 = Vcmax25,
                     Jmax25  = Jmax25,
                     Rd25    = Rd25,
                     TPU25   = TPU25,
                     Tleaf_all = Tleaf_all,Qin_all=Qin_all,
                     Ci_all = Ci_all,curve_identifier=curve_identifier,
                     An_obs = An_obs
                     )
Q10s  = as.data.frame(opt_results$solution)
rownames(Q10s) = c( 'Q10_1',  'Q10_2', 'Q10_3',
                    'Q10_5',  'Q10_6', 'Q10_7',
                    'Q10_8',  'Q10_9', 'Q10_10',
                    'Q10_13', 'Q10_23')
Q10s  = t(Q10s)
write.csv(Q10s,paste0("ePhotosynthesis_optimal_Q10s_",outfile_prefix,"_v1.csv"))

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
