################################################################################
# Produce the diurnal cycle of a variable
# 
# diuCyc.R
#
# ReadMe: produce diurnal cycles of climate variable
#
# Output: diuCyc          <- [mean hourly value, standard deviation, 90th percentile, 10th percentile]
# Input:  TS              <- time series (hourly or higher temporal resolution)
#         time            <- Datetime of timeseries as string ("2012-06-23 12:00:00 +0545")
#
# Necessary dependance: none
#
# Created:          2018/09/02
# Latest Revision:  2018/09/02
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht | Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 


diuCyc <- function(TS,time){
  
  diuCyc_m <- aggregate(TS,list(hour(time)),median,na.rm=T)   # mean values for diurnal cycle
  diuCyc_s <- aggregate(TS,list(hour(time)),sd,na.rm=T)     # standard deviation for diurnal cycle

  diuCyc_90 <- aggregate(TS, list(hour(time)), FUN = quantile, probs  = 0.90,na.rm=T)
  diuCyc_10 <- aggregate(TS, list(hour(time)), FUN = quantile, probs  = 0.10,na.rm=T)
  
  diuCyc <- data.frame(cbind(diuCyc_m,diuCyc_s,diuCyc_90,diuCyc_10))
}