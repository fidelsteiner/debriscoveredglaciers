################################################################################
# Model Surface Temperature from Air Temperature
# 
# AirSurfTemp.R
#
# ReadMe: Derive Surface Temperatures from a debris-covered surface based on 
# Steiner, J. and Pellicciotti, F. (2016) 'On the variability of air temperature over a debris covered glacier , Nepalese Himalaya', Annals of Glaciology, 57(71), pp. 1-13. doi: 10.3189/2016AoG71A066.
# 
# Output: AirSurfTemp[1]  <- Surface Temperature (timeseries)
#         AirSurfTemp[2]  <- Statistics on model fit (R2, RMSE, MBE)
# Input:  Ta_data         <- Air Temperature (2m measurement height, [degC])
#         TimeString      <- Datetime of timeseries as string ("2012-06-23 12:00:00 +0545")
#         pComb           <- Choice of Parameter set to derive surface from air temperature (1: Steiner&Pellicciotti2016) 
#                   1: based on Steiner&Pellicciotti2016
#         clim            <-  Choice of Climate zone
#                   1: Monsoon dominated Himalaya, uses separation into pre-monsoon, monsoon and post-monsoon as in Steiner&Pellicciotti2016
#
# Necessary dependance: TSStatistics.R
#
# Created:          2018/09/02
# Latest Revision:  2018/09/02
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht | Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 

AirSurfTemp <- function(Ta_data,TimeString,pComb,clim) { 
  
  ###########
  # Model Parameters
  ###########
  if(pComb == 1){
    # all parameters from Steiner&Pellicciotti2016
    a1_preM <- mean(c(0.61,0.56))
    b1_preM <- mean(c(2.55,2.65))
    a2_preM <- mean(c(0.3,0.33))
    b2_preM <- mean(c(4.5,4.67))
    Tc_preM <- mean(c(6.67,7.76))
    Tc_air_preM <- a1_preM*(Tc_preM) + b1_preM
    
    a1_M <- mean(c(0.57,0.52))
    b1_M <- mean(c(3.22,3.67))
    a2_M <- mean(c(0.26,0.28))
    b2_M <- mean(c(5.95,6.5))
    Tc_M <- mean(c(9.99,11.41))
    Tc_air_M <- a1_M*(Tc_M) + b1_M
    
    a1_posM <- mean(c(0.56,0.59))
    b1_posM <- mean(c(2.2,1.84))
    a2_posM <- mean(c(0.33,0.29))
    b2_posM <- mean(c(3.24,4.00))
    Tc_posM <- mean(c(6.30,5.83))
    Tc_air_posM <- a1_posM*(Tc_posM) + b1_posM
    
  }
 
  ###########
  # climate Settings
  ###########
  if(clim == 1){ 
    monPoM <- 9     # month that post-monsoon starts
    dayPoM <- 20    # day that post-monsoon starts
    monPrM <- 6     # month that pre-monsoon ends
    dayPrM <- 15    # day that pre-monsoon ends
  postmonDates <- which(month(TimeString)>monPrM&day(TimeString)>dayPrM|month(TimeString)>(monPrM+1))   # Dates for post-Monsoon
  premonDates <- which(month(TimeString)<(monPrM+1)&day(TimeString)<dayPrM|month(TimeString)<monPrM)    # Dates for pre-Monsoon
  a<-1:length(Ta_data)
  monDates <- a[!(a %in% c(premonDates,postmonDates))]
  
  Dat_preM <- Ta_data[premonDates]
  Dat_M <- Ta_data[monDates]
  Dat_posM <- Ta_data[postmonDates]
  
  Ts_mod_preM <- Dat_preM*0
  Ts_mod_preM[Dat_preM<=Tc_air_preM] <- (Dat_preM[Dat_preM<=Tc_air_preM] - b1_preM) / a1_preM
  Ts_mod_preM[Dat_preM>Tc_air_preM] <- (Dat_preM[Dat_preM>Tc_air_preM] - b2_preM) / a2_preM
  
  Ts_mod_M <- Dat_M*0
  Ts_mod_M[Dat_M<=Tc_air_M] <- (Dat_M[Dat_M<=Tc_air_M] - b1_M) / a1_M
  Ts_mod_M[Dat_M>Tc_air_M] <- (Dat_M[Dat_M>Tc_air_M] - b2_M) / a2_M
  
  Ts_mod_posM <- Dat_posM*0
  Ts_mod_posM[Dat_posM<=Tc_air_M] <- (Dat_posM[Dat_posM<=Tc_air_M] - b1_M) / a1_M
  Ts_mod_posM[Dat_posM>Tc_air_M] <- (Dat_posM[Dat_posM>Tc_air_M] - b2_M) / a2_M
  
  AirSurfTemp2 <- c(Ts_mod_preM,Ts_mod_M,Ts_mod_posM)
  }


  
  # Calculate fit
  Stats <- TSStatistics(Ta_data,AirSurfTemp2) # Calculate basic statistics

  AirSurfTemp <- return(list(v1=AirSurfTemp2,v2=Stats))
}