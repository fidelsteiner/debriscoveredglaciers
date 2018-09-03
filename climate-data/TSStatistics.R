################################################################################
# Statistics of Bulk Models
# 
# TSStatsitics.R
#
# ReadMe: calculate Statistics for Model Performance (r2, RMSE, MBE)
#
# Input:  Qmeas - measured time series
#         Qmod - modelled timeseries
#
# Output: TSStatistics - [R2, RMSE, MBE]
# 
# Created:          2018/09/02
# Latest Revision:  2018/09/02
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht | Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 


TSStatistics <- function(Qmeas, Qmod) { 

  TStats <- array(0, c(1, 3));
  TStats[1,1] <- summary(lm(Qmeas ~ Qmod))$r.squared   # R2 of Sensible Heat Flux
  TStats[1,2] <- sqrt(sum((Qmeas - Qmod)^2,na.rm=T) / length(which(is.na((Qmeas - Qmod)^2)==F))); #RMSE
  TStats[1,3] <- sum(Qmeas - Qmod,na.rm=T) / length(which(is.na((Qmeas - Qmod)^2)==F)) # MBE
  
  
  TSStatistics <- TStats
}