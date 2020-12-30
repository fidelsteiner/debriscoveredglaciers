################################################################################
# Statistics of Bulk Models
# 
# LapseStatsitics.R
#
# ReadMe: calculate Statistics for Model Performance
# 
# Created:          2017/11/27
# Latest Revision:  2017/11/27
#
# Jakob F Steiner | PhD candidate | Faculty of Geosciences | Universiteit Utrecht | 
# Heidelberglaan 2, 3584 CS Utrecht | W.C. van Unnik building | Room 124, Zonneveldvleugel | 
# j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org
#

LapseStatsitics <- function(Qmeas, Qmod, match) { 
#browser()
  TStats <- array(0, c(1, 3));
  TStats[1,1] <- summary(lm(Qmeas[match] ~ Qmod[match]))$r.squared   # R2 of Sensible Heat Flux
  TStats[1,2] <- sqrt(sum((Qmeas[match] - Qmod[match])^2,na.rm=T) / length(which(is.na((Qmeas[match] - Qmod[match])^2)==F))); #RMSE
  TStats[1,3] <- sum(Qmeas[match] - Qmod[match],na.rm=T) / length(which(is.na((Qmeas[match] - Qmod[match])^2)==F)) # MBE
  
  
  LapseStatsitics <- TStats
}