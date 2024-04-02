#Farquhar model translated from Matlab
FarquharModel<-function(LeafTemperature, Ci, Radiation_PAR, Air_O2,Vcmax25,Jmax25, Rate_TPu,other_parameters)
{
  R=8.31446261815324e-3
  PhotosynthesisTheta = 0.76
  PhiPS2_base = 0.352
  a1 = 4.5   #electrons_per_carboxylation per C
  a2 = 10.5  #electrons_per_carboxylation per O2
  leaf_reflectance   = 0.1  #BioCro R
  leaf_transmittance = 0.05  #BioCro R
  betaPSII = 0.5 #BioCro R
  TPU_c = 25.5   #BioCro R
  Ha    = 62.99  #BioCro R
  S     = 0.588  #BioCro R
  Hd    = 182.14 #BioCro R
  TPU_rate_sf25 = 306.7509 #nomalizer at 25 oC, this is what's used in BioCro
  
  LeafTemperatureKelvin = LeafTemperature + 273.15  #Leaf temperature in K
  top_term = LeafTemperatureKelvin * exp(TPU_c-Ha/(R*LeafTemperatureKelvin))
  bot_term = 1+ exp((S*LeafTemperatureKelvin - Hd)/(R*LeafTemperatureKelvin))
  TPU_rate_sf = top_term /bot_term / TPU_rate_sf25
  
  GammaStar = exp(19.02 - 37.83 / (R * LeafTemperatureKelvin))
  Ko = exp(20.30 - 36.38 / (R * LeafTemperatureKelvin))
  Kc = exp(38.05 - 79.43 / (R * LeafTemperatureKelvin))	
  Vcmax = Vcmax25 * exp(26.35 - 65.33 / (R * LeafTemperatureKelvin))

#  PhiPS2 = PhiPS2_base*other_parameters[1] + 0.022 * other_parameters[2]*LeafTemperature 
#            - 3.4 * other_parameters[3]*LeafTemperature^2.0 / 10000.0
  PhiPS2 = PhiPS2_base*other_parameters[1] + 0.022 * LeafTemperature 
            - 3.4 * LeafTemperature^2.0 / 10000.0
  
#  PhiPS2 = PhiPS2 * other_parameters[1]
  
  I = Radiation_PAR * PhiPS2 * (1.0 - leaf_reflectance - leaf_transmittance) * betaPSII
  
  ThetaPS2 = PhotosynthesisTheta * other_parameters[2] + 0.01713 * LeafTemperature - 3.75 * LeafTemperature^2.0 / 10000.0
  
#  ThetaPS2 = ThetaPS2 * other_parameters[2]
  Jmax = Jmax25 * exp(17.57 - 43.54 / (R * LeafTemperatureKelvin))
  
  if( (I + Jmax)^2.0 - 4.0 * ThetaPS2 * I * Jmax < 0){
    return(NA)
  }
  
  J = (I + Jmax - sqrt((I + Jmax)^2.0 - 4.0 * ThetaPS2 * I * Jmax)) / (2.0 * ThetaPS2)
  
  LeafAc = (1.0 - GammaStar / Ci) * (Vcmax * Ci) /(Ci + Kc * (1.0 + Air_O2 / Ko)) #Rubisco limited photosynthesis
  LeafAj = (1.0 - GammaStar / Ci) * (J * Ci) /(a1 * Ci + a2 * GammaStar)  #Light limited photosynthesis
  
  if (LeafAj < 0.0){
      LeafAj = 0.0
  }
  
  LeafAp = 3.0 * Rate_TPu * TPU_rate_sf  #TPU limited photosynthesis
  
  if (LeafAp < 0.0){
      LeafAp=0.0
  }
  
  output = list(GA = min(c(LeafAc,LeafAj,LeafAp)), #Minimum of three limitations
                Ac = LeafAc,
                Aj = LeafAj,
                Ap = LeafAp
               )
  return(output)
}
