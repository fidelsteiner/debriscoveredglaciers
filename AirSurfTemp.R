################################################################################
# Model Surface Temperature from Air Temperature
# 
# AirSurfTemp.R
#
# ReadMe: Derive Surface Temperatures from a debris-covered surface based on 
# Steiner, J. and Pellicciotti, F. (2016) 'On the variability of air temperature over a debris covered glacier , Nepalese Himalaya', Annals of Glaciology, 57(71), pp. 1-13. doi: 10.3189/2016AoG71A066.
# 
# Output: Surface Temperature (timeseries)
#         Statistics on model fit
# Input:  Air Temperature (2m measurement height, [degC]) - Ta_data
#         Datetime of timeseries as string ("2012-06-23 12:00:00 +0545") - TimeString
#         Choice of Parameter set to derive surface from air temperature (1: Steiner&Pellicciotti2016) - Param
#         
# Necessary dependance: TSStatistics.R
#
# Created:          2018/09/02
# Latest Revision:  2018/09/02
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht | Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 

AirSurfTemp <- function(Ta_data,TimeString,Param) { 
  
  a1_preM <- Param[1]
  b1_preM <- Param[2]
  a2_preM <- Param[3]
  b3_preM <- Param[4]
  Tc_air_preM <- Param[5]
  a1_M <- Param[6]
  b1_M <- Param[7]
  a2_M <- Param[8]
  b2_M <- Param[9]
  Tc_air_M <- Param[10]
  a1_posM <- Param[11]
  b1_posM <- Param[12]
  a2_posM <- Param[13]
  b2_posM <- Param[14]
  Tc_air_posM <- Param[15]
  
  postmonDates <- which((month(TimeString)>9&day(TimeString)>20)|month(TimeString)>10) # Dates for post-Monsoon
  premonDates <- which((month(TimeString)<7&day(TimeString)<15)|month(TimeString)<6) # Dates for pre-Monsoon
  
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
  
  # Calculate fit
  Stats <- TSStatistics(Ta_data,AirSurfTemp2) # Calculate basic statistics

  AirSurfTemp <- return(list(v1=AirSurfTemp2,v2=Stats))
}