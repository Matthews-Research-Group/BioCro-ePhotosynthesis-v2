# Define the RMSE function
rmse <- function(observed, predicted) {
  # Check if the lengths of the observed and predicted vectors are the same
  if (length(observed) != length(predicted)) {
    stop("The lengths of observed and predicted values must be the same.")
  }
  
  # Calculate the squared differences
  squared_diff <- (observed - predicted)^2
  
  # Calculate the mean of the squared differences
  mean_squared_diff <- mean(squared_diff)
  
  # Calculate the square root of the mean squared differences (RMSE)
  rmse_value <- sqrt(mean_squared_diff)
  
  return(rmse_value)
}

obj_Q10<-function(Q10s,curve_type,alpha1_alpha2,Vcmax25,
                  Jmax25,Rd25,TPU25,
                  Tleaf_all,Qin_all,Ci_all,curve_identifier,An_obs)
{
  if(curve_type==1){
    An_ePhoto   = NA* 1:length(Tleaf_all)
    for (i in 1:length(Tleaf_all)){ 
      PAR = Qin_all[i]
      Tleaf = Tleaf_all[i]
      Ci  = Ci_all[i]
      
      # Arguments for the executable
      args <- c(as.numeric(alpha1), as.numeric(alpha2), PAR, Tleaf,Ci,Q10s)
      
      # Run the executable and capture the output
      ephoto <- system2(executable_path, args = args, stdout = TRUE, stderr = TRUE)
      if(length(ephoto)!=1) {
        print(ephoto)
        ephoto = tail(ephoto,1)
      }
      ephoto = as.numeric(ephoto)
      Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
      An_ePhoto[i] = ephoto - Rd
    }
    An_ePhoto = as.numeric(An_ePhoto)
    df = data.frame(Ci=Ci_all,An_ePhoto,An_obs,curve_identifier)
    #plot average line
    unique_identifier = unique(df$curve_identifier)
    xx = 0
    An_ePhoto  = c()
    An_obs     = c()
    for (id in unique_identifier){
      tmp = df[df$curve_identifier == id,c("Ci","An_ePhoto","An_obs")]
      xx = xx + tmp
      An_ePhoto = cbind(An_ePhoto,tmp$An_ePhoto)
      An_obs    = cbind(An_obs   ,tmp$An_obs)
    }
    df_avg = xx/length(unique_identifier)
    #select only points with Ci larger than some threshold
    df_avg = df_avg[df_avg$Ci>100,] 
    err = rmse(df_avg$An_obs,df_avg$An_ePhoto) 
    return(err)

  }else if(curve_type==2){



  }
}
