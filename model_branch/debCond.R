################################################################################
# Determine debris conductivity from station data
# 
# debCond.R
#
# ReadMe:
# 
#
# Created:          2018/02/05
# Latest Revision:  2017/02/05
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht 
# Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# packages (if not installed yet: install.packages('examplePackage')
library(pacman)
library('lubridate')
library(RColorBrewer)
p_load(rgeos,maptools,raster,foreach,lubridate,compare,colorRamps,data.table,circular,parallel,snowfall,truncnorm,rlecuyer,forecast,rasterVis,R.utils,zoo)

projec<-'+proj=utm +zone=45N +datum=WGS84'
#projec<-'+proj=longlat +datum=WGS84'

path_profiles <- "F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\DebrisData\\Profiles\\"
path <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB'        
path_code <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code'
path_subcode <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes'

path_figs <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Figures'

################# 
# Paths for Extra Codes
file.sources = list.files(path_subcode, 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

# Seasons
monin <- as.numeric(strftime('1970/06/15', format = "%j"))                                         # define beginning of wet season
monout <- as.numeric(strftime('1970/09/19', format = "%j"))                                        # define end of wet season

winin <- as.numeric(strftime('1970/01/01', format = "%j"))                                           # define end of second dry season
winout <- as.numeric(strftime('1970/02/28', format = "%j"))                                         # define beginning of first dry season

ECsite_profile1 <- read.csv(path_profiles&'\\'&'UU_Site\\201810_LirungECSite_TempProfile.csv')
ECsite_profile2 <- read.csv(path_profiles&'\\'&'UU_Site\\201810_LirungECSite2_TempProfile.csv')
ETHSiteprofile <- read.csv(path_profiles&'\\'&'ETH_Site\\CR1000_thermistor1213raw.csv')

ECprofile1 <- list()
Sys.setenv(TZ='Asia/Kathmandu')
ECprofile1$Time_Str <- as.POSIXct(paste(as.character(ECsite_profile1$DATE),as.character(ECsite_profile1$TIME)), format="%m/%d/%Y  %H:%M:%S")
ECprofile1$Time_Num <- as.numeric(ECprofile1$Time_Str)
ECprofile1$surf1 <- TSaggregate(ECsite_profile1$T4,ECprofile1$Time_Str,60,2017,'mean')
ECprofile1$surf2 <- TSaggregate(ECsite_profile1$T3,ECprofile1$Time_Str,60,2017,'mean')
ECprofile1$surf3 <- TSaggregate(ECsite_profile1$T1,ECprofile1$Time_Str,60,2017,'mean')
ECprofile1$surf4 <- TSaggregate(ECsite_profile1$T2,ECprofile1$Time_Str,60,2017,'mean')
ECprofile1$DateTime <- as.POSIXct(ECprofile1$surf4[,1],origin = '1970-01-01')
ECprofile1$Season <- ECprofile1$surf4[,1] * 0
doyMatcher <- as.numeric(strftime(ECprofile1$DateTime, format = "%j"))
ECprofile1$Season[doyMatcher<monin&doyMatcher>=winout] <- 1   # dry season
ECprofile1$Season[doyMatcher>=monin&doyMatcher<monout] <- 2   # wet season
ECprofile1$Season[doyMatcher>=monout] <- 1                    # dry season
ECprofile1$Season[doyMatcher<winout] <- 3                     # cold season

d_ECprofile <- c(0,25,45,65)

cols <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6')
png(file=path_figs&'\\DebrisFreezing.png', res = 160,width=2400,height=600)
par(mfcol=c(1,1), mar=c(4,5,1,2), oma=c(0,0,0,0),cex.axis = 1.5,cex.lab = 1.5)
plot(ECprofile1$DateTime[5200:6300],ECprofile1$surf1[5200:6300,2],type='l',lwd = 1,ylim=c(-5,5),xlab = 'Date', ylab = expression('T'['debris']~~'[°C]'),col=cols[1])
points(ECprofile1$DateTime[5200:6300],ECprofile1$surf2[5200:6300,2],col=cols[2],type='l',lty = 1,lwd = 3)
points(ECprofile1$DateTime[5200:6300],ECprofile1$surf3[5200:6300,2],col=cols[3],type='l',lty = 1,lwd = 3)
points(ECprofile1$DateTime[5200:6300],ECprofile1$surf4[5200:6300,2],col=cols[4],type='l',lty = 2,lwd = 2)
abline(h=seq(-10,20,1),v=seq(ECprofile1$DateTime[5000], ECprofile1$DateTime[6500], by = "day"), col="gray", lty=3)
legend('bottomright',legend = c('surface','0.25 m','0.45 m','0.65 m (ice)'),cex= 1.5,bty = 'n',lty=c(1,1,1,2),lwd = c(1,3,3,2),col=cols)
dev.off()

# Diurnal Cycles
diuProf1 <- list()
diuProf1$surf1_dry <- diuCyc(ECprofile1$surf1[ECprofile1$Season==1,2],as.POSIXct(ECprofile1$surf1[ECprofile1$Season==1,1],origin='1970-01-01'))
diuProf1$surf1_wet <- diuCyc(ECprofile1$surf1[ECprofile1$Season==2,2],as.POSIXct(ECprofile1$surf1[ECprofile1$Season==2,1],origin='1970-01-01'))
diuProf1$surf1_cold <- diuCyc(ECprofile1$surf1[ECprofile1$Season==3,2],as.POSIXct(ECprofile1$surf1[ECprofile1$Season==3,1],origin='1970-01-01'))

diuProf1$surf2_dry <- diuCyc(ECprofile1$surf2[ECprofile1$Season==1,2],as.POSIXct(ECprofile1$surf2[ECprofile1$Season==1,1],origin='1970-01-01'))
diuProf1$surf2_wet <- diuCyc(ECprofile1$surf2[ECprofile1$Season==2,2],as.POSIXct(ECprofile1$surf2[ECprofile1$Season==2,1],origin='1970-01-01'))
diuProf1$surf2_cold <- diuCyc(ECprofile1$surf2[ECprofile1$Season==3,2],as.POSIXct(ECprofile1$surf2[ECprofile1$Season==3,1],origin='1970-01-01'))

diuProf1$surf3_dry <- diuCyc(ECprofile1$surf3[ECprofile1$Season==1,2],as.POSIXct(ECprofile1$surf3[ECprofile1$Season==1,1],origin='1970-01-01'))
diuProf1$surf3_wet <- diuCyc(ECprofile1$surf3[ECprofile1$Season==2,2],as.POSIXct(ECprofile1$surf3[ECprofile1$Season==2,1],origin='1970-01-01'))
diuProf1$surf3_cold <- diuCyc(ECprofile1$surf3[ECprofile1$Season==3,2],as.POSIXct(ECprofile1$surf3[ECprofile1$Season==3,1],origin='1970-01-01'))

diuProf1$surf4_dry <- diuCyc(ECprofile1$surf4[ECprofile1$Season==1,2],as.POSIXct(ECprofile1$surf4[ECprofile1$Season==1,1],origin='1970-01-01'))
diuProf1$surf4_wet <- diuCyc(ECprofile1$surf4[ECprofile1$Season==2,2],as.POSIXct(ECprofile1$surf4[ECprofile1$Season==2,1],origin='1970-01-01'))
diuProf1$surf4_cold <- diuCyc(ECprofile1$surf4[ECprofile1$Season==3,2],as.POSIXct(ECprofile1$surf4[ECprofile1$Season==3,1],origin='1970-01-01'))




dTdt1 <- list()
dT2dz2 <- list()

colRange=colorRampPalette(c("red","yellow", "blue"))
dayhours <- function(TimeVec,Seas,nightday,year){
  if(is.null(year)){
    year <- seq(1,length(TimeVec$Season),1)
  }
  #browser()
  hourofday <- hour(as.POSIXct(TimeVec$surf1[year,][TimeVec$Season[year]==Seas,1],origin='1970-01-01'))
  nighthours_dry <- which(hourofday %in%c(0,1,2,3,4,5,6,7,20,21,22,23))
  dayhours_dry <- which(hourofday %in%c(8,9,10,11,12,13,14,15,16,17,18,19))
  if(nightday==1){
  return(dayhours_dry)
  }
  if(nightday==2){
    return(nighthours_dry)
  }
}

# At Profile 2
dTdt1$p2_dry <- c(NA,diff(ECprofile1$surf2[ECprofile1$Season==1,2]) / 3600)
dTdt1$p2_dry[dTdt1$p2_dry >4*10^-4]<-NA
dTdt1$p2_dry[which(abs(diff(ECprofile1$surf2[ECprofile1$Season==1,2]))<=0.02)] <- NA
dT2dz2$p2_dry <- (ECprofile1$surf1[ECprofile1$Season==1,2] - 2*ECprofile1$surf2[ECprofile1$Season==1,2] + ECprofile1$surf3[ECprofile1$Season==1,2]) / (0.65-0.45)^2
dTdt1$p2_wet <- c(NA,diff(ECprofile1$surf2[ECprofile1$Season==2,2]) / 3600)
dTdt1$p2_wet[which(abs(diff(ECprofile1$surf2[ECprofile1$Season==2,2]))<=0.02)] <- NA
dT2dz2$p2_wet <- (ECprofile1$surf1[ECprofile1$Season==2,2] - 2*ECprofile1$surf2[ECprofile1$Season==2,2] + ECprofile1$surf3[ECprofile1$Season==2,2]) / (0.65-0.45)^2
dTdt1$p2_cold <- c(NA,diff(ECprofile1$surf2[ECprofile1$Season==3,2]) / 3600)
dTdt1$p2_cold[which(abs(diff(ECprofile1$surf2[ECprofile1$Season==3,2]))<=0.02)] <- NA
dT2dz2$p2_cold <- (ECprofile1$surf1[ECprofile1$Season==3,2] - 2*ECprofile1$surf2[ECprofile1$Season==3,2] + ECprofile1$surf3[ECprofile1$Season==3,2]) / (0.65-0.45)^2


#layout(matrix(1:3,ncol=3), width = c(2,2,2),height = c(1,1))
#plot(dT2dz2$p2_dry,dTdt1$p2_dry,col='red')
#points(dT2dz2$p2_dry[dayhours(ECprofile1,1,2,NA)],dTdt1$p2_dry[dayhours(ECprofile1,1,2)])
#abline(lm(dTdt1$p2_dry[dayhours(ECprofile1,1,2)]~dT2dz2$p2_dry[dayhours(ECprofile1,1,2)]))
#plot(dT2dz2$p2_wet ,dTdt1$p2_wet,col='red')
#points(dT2dz2$p2_wet[dayhours(ECprofile1,2,2)],dTdt1$p2_wet[dayhours(ECprofile1,2,2)])
#abline(lm(dTdt1$p2_wet[dayhours(ECprofile1,2,2)]~dT2dz2$p2_wet[dayhours(ECprofile1,2,2)]))
#plot(dT2dz2$p2_cold ,dTdt1$p2_cold,col='red')
#points(dT2dz2$p2_cold[dayhours(ECprofile1,3,2)],dTdt1$p2_cold[dayhours(ECprofile1,3,2)])
#abline(lm(dTdt1$p2_cold[dayhours(ECprofile1,3,2)]~dT2dz2$p2_cold[dayhours(ECprofile1,3,2)]))


# At Profile 3
dTdt1$p3_dry <- c(NA,diff(ECprofile1$surf3[ECprofile1$Season==1,2]) / 3600)
dTdt1$p3_dry[dTdt1$p3_dry>2*10^-4]<-NA
dTdt1$p3_dry[abs(diff(ECprofile1$surf3[ECprofile1$Season==1,2]))<=0.02]<-NA
dT2dz2$p3_dry <- (ECprofile1$surf2[ECprofile1$Season==1,2] - 2*ECprofile1$surf3[ECprofile1$Season==1,2] + ECprofile1$surf4[ECprofile1$Season==1,2]) / (0.2)^2
dTdt1$p3_wet <- c(NA,diff(ECprofile1$surf3[ECprofile1$Season==2,2]) / 3600)
dTdt1$p3_wet[abs(diff(ECprofile1$surf3[ECprofile1$Season==2,2]))<=0.02]<-NA
dT2dz2$p3_wet <- (ECprofile1$surf2[ECprofile1$Season==2,2] - 2*ECprofile1$surf3[ECprofile1$Season==2,2] + ECprofile1$surf4[ECprofile1$Season==2,2]) / (0.2)^2
dTdt1$p3_cold <- c(NA,diff(ECprofile1$surf3[ECprofile1$Season==3,2]) / 3600)
dTdt1$p3_cold[abs(diff(ECprofile1$surf3[ECprofile1$Season==3,2]))<=0.02]<-NA
dT2dz2$p3_cold <- (ECprofile1$surf2[ECprofile1$Season==3,2] - 2*ECprofile1$surf3[ECprofile1$Season==3,2] + ECprofile1$surf4[ECprofile1$Season==3,2]) / (0.2)^2

#plot(dT2dz2$p3_dry[dayhours(ECprofile1,1,1)] ,dTdt1$p3_dry[dayhours(ECprofile1,1,1)],col='red')
#points(dT2dz2$p3_dry[dayhours(ECprofile1,1,2)] ,dTdt1$p3_dry[dayhours(ECprofile1,1,2)])
#plot(dT2dz2$p3_wet[dayhours(ECprofile1,2,1)] ,dTdt1$p3_wet[dayhours(ECprofile1,2,1)],col='red')
#points(dT2dz2$p3_wet[dayhours(ECprofile1,2,2)] ,dTdt1$p3_wet[dayhours(ECprofile1,2,2)])
#plot(dT2dz2$p3_cold[dayhours(ECprofile1,3,1)] ,dTdt1$p3_cold[dayhours(ECprofile1,3,1)],col='red')
#points(dT2dz2$p3_cold[dayhours(ECprofile1,3,2)] ,dTdt1$p3_cold[dayhours(ECprofile1,3,2)])

#Derive diffusivities/conductivities

fit1 <- summary(lm(dTdt1$p2_dry[dayhours(ECprofile1,1,2,NULL)]~dT2dz2$p2_dry[dayhours(ECprofile1,1,2,NULL)]))
fit2 <- summary(lm(dTdt1$p2_wet[dayhours(ECprofile1,2,2,NULL)]~dT2dz2$p2_wet[dayhours(ECprofile1,2,2,NULL)]))
fit3 <- summary(lm(dTdt1$p2_cold[dayhours(ECprofile1,3,2,NULL)]~dT2dz2$p2_cold[dayhours(ECprofile1,3,2,NULL)]))

diff_UU_1 <- c(fit1$coefficients[2,c(1,2)],fit2$coefficients[2,c(1,2)],fit3$coefficients[2,c(1,2)])*10^6

fit1 <- summary(lm(dTdt1$p2_dry~dT2dz2$p2_dry))
fit2 <- summary(lm(dTdt1$p2_wet~dT2dz2$p2_wet))
fit3 <- summary(lm(dTdt1$p2_cold~dT2dz2$p2_cold))

diff_UU_1_whole <- c(fit1$coefficients[2,c(1,2)],fit2$coefficients[2,c(1,2)],fit3$coefficients[2,c(1,2)])*10^6

fit1 <- summary(lm(dTdt1$p3_dry[dayhours(ECprofile1,1,2,NULL)]~dT2dz2$p3_dry[dayhours(ECprofile1,1,2,NULL)]))
fit2 <- summary(lm(dTdt1$p3_wet[dayhours(ECprofile1,2,2,NULL)]~dT2dz2$p3_wet[dayhours(ECprofile1,2,2,NULL)]))
fit3 <- summary(lm(dTdt1$p3_cold[dayhours(ECprofile1,3,2,NULL)]~dT2dz2$p3_cold[dayhours(ECprofile1,3,2,NULL)]))

diff_UU_2 <- c(fit1$coefficients[2,c(1,2)],fit2$coefficients[2,c(1,2)],fit3$coefficients[2,c(1,2)])*10^6

fit1 <- summary(lm(dTdt1$p3_dry~dT2dz2$p3_dry))
fit2 <- summary(lm(dTdt1$p3_wet~dT2dz2$p3_wet))
fit3 <- summary(lm(dTdt1$p3_cold~dT2dz2$p3_cold))

diff_UU_2_whole <- c(fit1$coefficients[2,c(1,2)],fit2$coefficients[2,c(1,2)],fit3$coefficients[2,c(1,2)])*10^6

# Aggregate the data for profile graphs at 00:00, 04:00, 08:00, 12:00. 16:00, 20:00
Tdry_profile <- cbind(diuProf1$surf4_dry[c(1,5,9,13,17,21),2],diuProf1$surf3_dry[c(1,5,9,13,17,21),2],diuProf1$surf2_dry[c(1,5,9,13,17,21),2],diuProf1$surf1_dry[c(1,5,9,13,17,21),2])
Twet_profile <- cbind(diuProf1$surf4_wet[c(1,5,9,13,17,21),2],diuProf1$surf3_wet[c(1,5,9,13,17,21),2],diuProf1$surf2_wet[c(1,5,9,13,17,21),2],diuProf1$surf1_wet[c(1,5,9,13,17,21),2])
Tcold_profile <- cbind(diuProf1$surf4_cold[c(1,5,9,13,17,21),2],diuProf1$surf3_cold[c(1,5,9,13,17,21),2],diuProf1$surf2_cold[c(1,5,9,13,17,21),2],diuProf1$surf1_cold[c(1,5,9,13,17,21),2])

png(file=path_figs&'\\DebrisTempProfile1.png', res = 160,width=2400,height=600)
par(mfcol=c(1,4), mar=c(4,5,1,2), oma=c(0,0,0,0),cex.axis = 1.5,cex.lab = 1.5)
#layout(matrix(c(1,3), nrow = 1, ncol = 3, byrow = TRUE))
plot(diuProf1$surf1_cold[,2],type='l',col='blue',ylim = c(-5,10),xlab = 'Time [h]',ylab = expression('T'['debris']~~'[°C]'),lwd = 2)
points(diuProf1$surf2_cold[,2],type='l',col='blue',lwd = 2,lty = 2)
points(diuProf1$surf3_cold[,2],type='l',col='blue',lwd = 2,lty = 3)
points(diuProf1$surf4_cold[,2],type='l',col='blue',lwd = 2,lty= 4)
grid(NULL,NULL)
plot(diuProf1$surf1_dry[,2],type='l',col='black',ylim = c(0,16),xlab = 'Time [h]',ylab = expression('T'['debris']~~'[°C]'),lwd = 2)
points(diuProf1$surf2_dry[,2],type='l',col='black',lwd = 2,lty = 2)
points(diuProf1$surf3_dry[,2],type='l',col='black',lwd = 2,lty = 3)
points(diuProf1$surf4_dry[,2],type='l',col='black',lwd = 2,lty= 4)
legend('topleft',legend = c('surface','0.25 m','0.45 m','0.65 m (ice)'),cex= 1.5,bty = 'n',lty=c(1,2,3,4),lwd = 2)
grid(NULL,NULL)
plot(diuProf1$surf1_wet[,2],type='l',col='red',ylim = c(0,16),xlab = 'Time [h]',ylab = expression('T'['debris']~~'[°C]'),lwd = 2)
points(diuProf1$surf2_wet[,2],type='l',col='red',lwd = 2,lty = 2)
points(diuProf1$surf3_wet[,2],type='l',col='red',lwd = 2,lty = 3)
points(diuProf1$surf4_wet[,2],type='l',col='red',lwd = 2,lty= 4)
grid(NULL,NULL)
matplot(t(Tdry_profile),rev(d_ECprofile/100),ylim = rev(range(d_ECprofile/100)),type='o',pch = 1:6,xlim=c(-2,17),lty=c(1,2,3,4,5,6),col='black',cex =  1.5,ylab = 'Depth [m]',xlab = expression('T'['debris']~~'[°C]'))
matplot(t(Twet_profile),rev(d_ECprofile/100),type='o',pch = 1:6,xlim=c(-2,17),lty=c(1,2,3,4,5,6),col='red',add=T,cex=2)
matplot(t(Tcold_profile),rev(d_ECprofile/100),type='o',pch = 1:6,xlim=c(-2,17),lty=c(1,2,3,4,5,6),col='blue',add=T,cex=2)
grid(NULL,NULL)
legend('bottomright',legend = c('00:00','04:00','08:00','12:00','16:00','20:00'),cex=  1.5,bty = 'n',lty=c(1,2,3,4,5,6))
dev.off()

lag_transfer_EC_dry <- (24 - which.max(diuProf1$surf1_dry[,2]) + which.max(diuProf1$surf4_dry[,2])) / 0.65
lag_transfer_EC_wet <- (24 - which.max(diuProf1$surf1_wet[,2]) + which.max(diuProf1$surf4_wet[,2])) / 0.65

# At Profile ETH
ETHprofile1 <- list()
Sys.setenv(TZ='Asia/Kathmandu')
ETHprofile1$Time_Str <- as.POSIXct(as.character(ETHSiteprofile$TIMESTAMP), format="%m/%d/%Y  %H:%M")
ETHprofile1$Time_Num <- as.numeric(ETHprofile1$Time_Str)
ETHSiteprofile$T1..C.[ETHSiteprofile$T1..C.>25|ETHSiteprofile$T1..C.<=-5] <- NA
ETHSiteprofile$T1..C.[which(diff(ETHSiteprofile$T1..C.)>5)+1] <- NA
ETHSiteprofile$T2..C.[ETHSiteprofile$T2..C.>25|ETHSiteprofile$T2..C.<=-5] <- NA
ETHSiteprofile$T2..C.[which(diff(ETHSiteprofile$T2..C.)>5)+1] <- NA
ETHSiteprofile$T3..C.[ETHSiteprofile$T3..C.>25|ETHSiteprofile$T3..C.<=-5] <- NA
ETHSiteprofile$T3..C.[which(diff(ETHSiteprofile$T3..C.)>5)+1] <- NA
ETHSiteprofile$T4..C.[ETHSiteprofile$T4..C.>25|ETHSiteprofile$T4..C.<=-5] <- NA
ETHSiteprofile$T4..C.[which(diff(ETHSiteprofile$T4..C.)>5)+1] <- NA
ETHSiteprofile$T5..C.[ETHSiteprofile$T5..C.>25|ETHSiteprofile$T5..C.<=-5] <- NA
ETHSiteprofile$T5..C.[which(diff(ETHSiteprofile$T5..C.)>5)+1] <- NA
ETHSiteprofile$T6..C.[ETHSiteprofile$T6..C.>25|ETHSiteprofile$T6..C.<=-5] <- NA
ETHSiteprofile$T6..C.[which(diff(ETHSiteprofile$T6..C.)>5)+1] <- NA
ETHSiteprofile$T7..C.[ETHSiteprofile$T7..C.>25|ETHSiteprofile$T7..C.<=-5] <- NA
ETHSiteprofile$T7..C.[which(diff(ETHSiteprofile$T7..C.)>5)+1] <- NA
ETHSiteprofile$T8..C.[ETHSiteprofile$T8..C.>25|ETHSiteprofile$T8..C.<=-5] <- NA
ETHprofile1$surf1 <- TSaggregate(ETHSiteprofile$T1..C.,ETHprofile1$Time_Str,60,2017,'mean')
ETHprofile1$surf2 <- TSaggregate(ETHSiteprofile$T2..C.,ETHprofile1$Time_Str,60,2017,'mean')
ETHprofile1$surf3 <- TSaggregate(ETHSiteprofile$T3..C.,ETHprofile1$Time_Str,60,2017,'mean')
ETHprofile1$surf4 <- TSaggregate(ETHSiteprofile$T4..C.,ETHprofile1$Time_Str,60,2017,'mean')
ETHprofile1$surf5 <- TSaggregate(ETHSiteprofile$T5..C.,ETHprofile1$Time_Str,60,2017,'mean')
ETHprofile1$surf6 <- TSaggregate(ETHSiteprofile$T6..C.,ETHprofile1$Time_Str,60,2017,'mean')
ETHprofile1$surf6[,2] <- ETHprofile1$surf6[,2] * NA
ETHprofile1$surf7 <- TSaggregate(ETHSiteprofile$T7..C.,ETHprofile1$Time_Str,60,2017,'mean')
ETHprofile1$surf8 <- TSaggregate(ETHSiteprofile$T8..C.,ETHprofile1$Time_Str,60,2017,'mean')
ETHprofile1$DateTime <- as.POSIXct(ETHprofile1$surf8[,1],origin = '1970-01-01')
ETHprofile1$Season <- ETHprofile1$surf8[,1] * 0
doyMatcher <- as.numeric(strftime(ETHprofile1$DateTime, format = "%j"))
y_2012 <- which(strftime(ETHprofile1$DateTime, format = "%Y")==2012)
y_2013 <- which(strftime(ETHprofile1$DateTime, format = "%Y")==2013)
ETHprofile1$Season[doyMatcher<monin&doyMatcher>=winout] <- 1   # dry season
ETHprofile1$Season[doyMatcher>=monin&doyMatcher<monout] <- 2   # wet season
ETHprofile1$Season[doyMatcher>=monout] <- 1                    # dry season
ETHprofile1$Season[doyMatcher<winout] <- 3                     # cold season

ETHprofile12 <- list()
ETHprofile12$Season <- ETHprofile1$Season[y_2012]
ETHprofile12$surf1 <- ETHprofile1$surf1[y_2012,]
ETHprofile12$surf2 <- ETHprofile1$surf2[y_2012,]
ETHprofile12$surf3 <- ETHprofile1$surf3[y_2012,]
ETHprofile12$surf4 <- ETHprofile1$surf4[y_2012,]
ETHprofile12$surf5 <- ETHprofile1$surf5[y_2012,]
ETHprofile12$surf6 <- ETHprofile1$surf6[y_2012,]
ETHprofile12$surf7 <- ETHprofile1$surf7[y_2012,]
ETHprofile12$surf8 <- ETHprofile1$surf8[y_2012,]
ETHprofile12$DateTime <- ETHprofile1$DateTime[y_2012]

ETHprofile13 <- list()
ETHprofile13$Season <- ETHprofile1$Season[y_2013]
ETHprofile13$surf1 <- ETHprofile1$surf1[y_2013,]
ETHprofile13$surf2 <- ETHprofile1$surf2[y_2013,]
ETHprofile13$surf3 <- ETHprofile1$surf3[y_2013,]
ETHprofile13$surf4 <- ETHprofile1$surf4[y_2013,]
ETHprofile13$surf5 <- ETHprofile1$surf5[y_2013,]
ETHprofile13$surf6 <- ETHprofile1$surf6[y_2013,]
ETHprofile13$surf7 <- ETHprofile1$surf7[y_2013,]
ETHprofile13$surf8 <- ETHprofile1$surf8[y_2013,]
ETHprofile13$DateTime <- ETHprofile1$DateTime[y_2013]

# Depths of sensors in different years
d_2012 = c(70,60,50,40,30,20,10,0) / 100
d_2013 = c(0,2,4,6,10,20,30,45) / 100

# Diurnal Cycles
diuProfETH2012 <- list()
diuProfETH2012$surf1_dry <- diuCyc(ETHprofile12$surf1[ETHprofile12$Season==1,2],as.POSIXct(ETHprofile12$surf1[ETHprofile12$Season==1,1],origin='1970-01-01'))
diuProfETH2012$surf1_wet <- diuCyc(ETHprofile12$surf1[ETHprofile12$Season==2,2],as.POSIXct(ETHprofile12$surf1[ETHprofile12$Season==2,1],origin='1970-01-01'))
diuProfETH2012$surf2_dry <- diuCyc(ETHprofile12$surf2[ETHprofile12$Season==1,2],as.POSIXct(ETHprofile12$surf2[ETHprofile12$Season==1,1],origin='1970-01-01'))
diuProfETH2012$surf2_wet <- diuCyc(ETHprofile12$surf2[ETHprofile12$Season==2,2],as.POSIXct(ETHprofile12$surf2[ETHprofile12$Season==2,1],origin='1970-01-01'))
diuProfETH2012$surf3_dry <- diuCyc(ETHprofile12$surf3[ETHprofile12$Season==1,2],as.POSIXct(ETHprofile12$surf3[ETHprofile12$Season==1,1],origin='1970-01-01'))
diuProfETH2012$surf3_wet <- diuCyc(ETHprofile12$surf3[ETHprofile12$Season==2,2],as.POSIXct(ETHprofile12$surf3[ETHprofile12$Season==2,1],origin='1970-01-01'))
diuProfETH2012$surf4_dry <- diuCyc(ETHprofile12$surf4[ETHprofile12$Season==1,2],as.POSIXct(ETHprofile12$surf4[ETHprofile12$Season==1,1],origin='1970-01-01'))
diuProfETH2012$surf4_wet <- diuCyc(ETHprofile12$surf4[ETHprofile12$Season==2,2],as.POSIXct(ETHprofile12$surf4[ETHprofile12$Season==2,1],origin='1970-01-01'))
diuProfETH2012$surf5_dry <- diuCyc(ETHprofile12$surf5[ETHprofile12$Season==1,2],as.POSIXct(ETHprofile12$surf5[ETHprofile12$Season==1,1],origin='1970-01-01'))
diuProfETH2012$surf5_wet <- diuCyc(ETHprofile12$surf5[ETHprofile12$Season==2,2],as.POSIXct(ETHprofile12$surf5[ETHprofile12$Season==2,1],origin='1970-01-01'))
diuProfETH2012$surf6_dry <- diuCyc(ETHprofile12$surf6[ETHprofile12$Season==1,2],as.POSIXct(ETHprofile12$surf3[ETHprofile12$Season==1,1],origin='1970-01-01'))
diuProfETH2012$surf6_wet <- diuCyc(ETHprofile12$surf6[ETHprofile12$Season==2,2],as.POSIXct(ETHprofile12$surf3[ETHprofile12$Season==2,1],origin='1970-01-01'))
diuProfETH2012$surf7_dry <- diuCyc(ETHprofile12$surf7[ETHprofile12$Season==1,2],as.POSIXct(ETHprofile12$surf5[ETHprofile12$Season==1,1],origin='1970-01-01'))
diuProfETH2012$surf7_wet <- diuCyc(ETHprofile12$surf7[ETHprofile12$Season==2,2],as.POSIXct(ETHprofile12$surf5[ETHprofile12$Season==2,1],origin='1970-01-01'))
diuProfETH2012$surf8_dry <- diuCyc(ETHprofile12$surf8[ETHprofile12$Season==1,2],as.POSIXct(ETHprofile12$surf3[ETHprofile12$Season==1,1],origin='1970-01-01'))
diuProfETH2012$surf8_dry[diuProfETH2012$surf8_dry[,12]<60,2]<-NA
diuProfETH2012$surf8_wet <- diuCyc(ETHprofile12$surf8[ETHprofile12$Season==2,2],as.POSIXct(ETHprofile12$surf3[ETHprofile12$Season==2,1],origin='1970-01-01'))
diuProfETH2012$surf8_wet[diuProfETH2012$surf8_wet[,12]<60,2]<-NA

diuProfETH2013 <- list()
diuProfETH2013$surf1_dry <- diuCyc(ETHprofile1$surf1[y_2013[ETHprofile1$Season[y_2013]==1],2],as.POSIXct(ETHprofile1$surf1[y_2013[ETHprofile1$Season[y_2013]==1],1],origin='1970-01-01'))
diuProfETH2013$surf1_wet <- diuCyc(ETHprofile1$surf1[y_2013[ETHprofile1$Season[y_2013]==2],2],as.POSIXct(ETHprofile1$surf1[y_2013[ETHprofile1$Season[y_2013]==2],1],origin='1970-01-01'))
diuProfETH2013$surf2_dry <- diuCyc(ETHprofile1$surf2[y_2013[ETHprofile1$Season[y_2013]==1],2],as.POSIXct(ETHprofile1$surf2[y_2013[ETHprofile1$Season[y_2013]==1],1],origin='1970-01-01'))
diuProfETH2013$surf2_wet <- diuCyc(ETHprofile1$surf2[y_2013[ETHprofile1$Season[y_2013]==2],2],as.POSIXct(ETHprofile1$surf2[y_2013[ETHprofile1$Season[y_2013]==2],1],origin='1970-01-01'))
diuProfETH2013$surf3_dry <- diuCyc(ETHprofile1$surf3[y_2013[ETHprofile1$Season[y_2013]==1],2],as.POSIXct(ETHprofile1$surf3[y_2013[ETHprofile1$Season[y_2013]==1],1],origin='1970-01-01'))
diuProfETH2013$surf3_wet <- diuCyc(ETHprofile1$surf3[y_2013[ETHprofile1$Season[y_2013]==2],2],as.POSIXct(ETHprofile1$surf3[y_2013[ETHprofile1$Season[y_2013]==2],1],origin='1970-01-01'))
diuProfETH2013$surf4_dry <- diuCyc(ETHprofile1$surf4[y_2013[ETHprofile1$Season[y_2013]==1],2],as.POSIXct(ETHprofile1$surf4[y_2013[ETHprofile1$Season[y_2013]==1],1],origin='1970-01-01'))
diuProfETH2013$surf4_wet <- diuCyc(ETHprofile1$surf4[y_2013[ETHprofile1$Season[y_2013]==2],2],as.POSIXct(ETHprofile1$surf4[y_2013[ETHprofile1$Season[y_2013]==2],1],origin='1970-01-01'))
diuProfETH2013$surf5_dry <- diuCyc(ETHprofile1$surf5[y_2013[ETHprofile1$Season[y_2013]==1],2],as.POSIXct(ETHprofile1$surf5[y_2013[ETHprofile1$Season[y_2013]==1],1],origin='1970-01-01'))
diuProfETH2013$surf5_wet <- diuCyc(ETHprofile1$surf5[y_2013[ETHprofile1$Season[y_2013]==2],2],as.POSIXct(ETHprofile1$surf5[y_2013[ETHprofile1$Season[y_2013]==2],1],origin='1970-01-01'))
diuProfETH2013$surf6_dry <- diuCyc(ETHprofile1$surf6[y_2013[ETHprofile1$Season[y_2013]==1],2],as.POSIXct(ETHprofile1$surf6[y_2013[ETHprofile1$Season[y_2013]==1],1],origin='1970-01-01'))
diuProfETH2013$surf6_wet <- diuCyc(ETHprofile1$surf6[y_2013[ETHprofile1$Season[y_2013]==2],2],as.POSIXct(ETHprofile1$surf6[y_2013[ETHprofile1$Season[y_2013]==2],1],origin='1970-01-01'))
diuProfETH2013$surf7_dry <- diuCyc(ETHprofile1$surf7[y_2013[ETHprofile1$Season[y_2013]==1],2],as.POSIXct(ETHprofile1$surf7[y_2013[ETHprofile1$Season[y_2013]==1],1],origin='1970-01-01'))
diuProfETH2013$surf7_wet <- diuCyc(ETHprofile1$surf7[y_2013[ETHprofile1$Season[y_2013]==2],2],as.POSIXct(ETHprofile1$surf7[y_2013[ETHprofile1$Season[y_2013]==2],1],origin='1970-01-01'))
diuProfETH2013$surf8_dry <- diuCyc(ETHprofile1$surf8[y_2013[ETHprofile1$Season[y_2013]==1],2],as.POSIXct(ETHprofile1$surf8[y_2013[ETHprofile1$Season[y_2013]==1],1],origin='1970-01-01'))
diuProfETH2013$surf8_wet <- diuCyc(ETHprofile1$surf8[y_2013[ETHprofile1$Season[y_2013]==2],2],as.POSIXct(ETHprofile1$surf8[y_2013[ETHprofile1$Season[y_2013]==2],1],origin='1970-01-01'))

diuProfETH2013$surf1_dry[diuProfETH2013$surf1_dry[,12]<30,2]<-NA
diuProfETH2013$surf2_dry[diuProfETH2013$surf2_dry[,12]<30,2]<-NA
diuProfETH2013$surf3_dry[diuProfETH2013$surf3_dry[,12]<30,2]<-NA
diuProfETH2013$surf4_dry[diuProfETH2013$surf4_dry[,12]<30,2]<-NA
diuProfETH2013$surf5_dry[diuProfETH2013$surf5_dry[,12]<30,2]<-NA
diuProfETH2013$surf6_dry[diuProfETH2013$surf6_dry[,12]<30,2]<-NA
diuProfETH2013$surf7_dry[diuProfETH2013$surf7_dry[,12]<30,2]<-NA
diuProfETH2013$surf8_dry[diuProfETH2013$surf8_dry[,12]<30,2]<-NA
diuProfETH2013$surf1_wet[diuProfETH2013$surf1_wet[,12]<30,2]<-NA
diuProfETH2013$surf2_wet[diuProfETH2013$surf2_wet[,12]<30,2]<-NA
diuProfETH2013$surf3_wet[diuProfETH2013$surf3_wet[,12]<30,2]<-NA
diuProfETH2013$surf4_wet[diuProfETH2013$surf4_wet[,12]<30,2]<-NA
diuProfETH2013$surf5_wet[diuProfETH2013$surf5_wet[,12]<30,2]<-NA
diuProfETH2013$surf6_wet[diuProfETH2013$surf6_wet[,12]<30,2]<-NA
diuProfETH2013$surf7_wet[diuProfETH2013$surf7_wet[,12]<30,2]<-NA
diuProfETH2013$surf8_wet[diuProfETH2013$surf8_wet[,12]<30,2]<-NA

lag_transfer_2012_dry <- (24-which.max(diuProfETH2012$surf7_dry[,2]) + which.max(diuProfETH2012$surf2_dry[,2])) / 0.5
lag_transfer_2012_wet <- (24 - which.max(diuProfETH2012$surf7_wet[,2]) + which.max(diuProfETH2012$surf2_wet[,2])) / 0.5

lag_transfer_2013_dry <- (which.max(diuProfETH2013$surf8_dry[,2]) - which.max(diuProfETH2013$surf1_dry[,2])) / 0.45
lag_transfer_2013_wet <- (which.max(diuProfETH2013$surf8_wet[,2]) - which.max(diuProfETH2013$surf1_wet[,2])) / 0.45

lag_carenzo <- c(0,3,5,7,10)
lag_carenzo_depth <- c(0.1,0.2,0.3,0.4,0.5)

# Aggregate the data for profile graphs at 00:00, 04:00, 08:00, 12:00. 16:00, 20:00
Tdry_profile_2012 <- cbind(diuProfETH2012$surf8_dry[c(1,5,9,13,17,21),2],diuProfETH2012$surf7_dry[c(1,5,9,13,17,21),2],diuProfETH2012$surf6_dry[c(1,5,9,13,17,21),2],diuProfETH2012$surf5_dry[c(1,5,9,13,17,21),2],diuProfETH2012$surf4_dry[c(1,5,9,13,17,21),2],diuProfETH2012$surf3_dry[c(1,5,9,13,17,21),2],diuProfETH2012$surf2_dry[c(1,5,9,13,17,21),2],diuProfETH2012$surf1_dry[c(1,5,9,13,17,21),2])
Twet_profile_2012 <- cbind(diuProfETH2012$surf8_wet[c(1,5,9,13,17,21),2],diuProfETH2012$surf7_wet[c(1,5,9,13,17,21),2],diuProfETH2012$surf6_wet[c(1,5,9,13,17,21),2],diuProfETH2012$surf5_wet[c(1,5,9,13,17,21),2],diuProfETH2012$surf4_wet[c(1,5,9,13,17,21),2],diuProfETH2012$surf3_wet[c(1,5,9,13,17,21),2],diuProfETH2012$surf2_wet[c(1,5,9,13,17,21),2],diuProfETH2012$surf1_wet[c(1,5,9,13,17,21),2])

Tdry_profile_2013 <- cbind(diuProfETH2013$surf1_dry[c(1,5,9,13,17,21),2],diuProfETH2013$surf2_dry[c(1,5,9,13,17,21),2],diuProfETH2013$surf3_dry[c(1,5,9,13,17,21),2],diuProfETH2013$surf4_dry[c(1,5,9,13,17,21),2],diuProfETH2013$surf5_dry[c(1,5,9,13,17,21),2],diuProfETH2013$surf6_dry[c(1,5,9,13,17,21),2],diuProfETH2013$surf7_dry[c(1,5,9,13,17,21),2],diuProfETH2013$surf8_dry[c(1,5,9,13,17,21),2])
Twet_profile_2013 <- cbind(diuProfETH2013$surf1_wet[c(1,5,9,13,17,21),2],diuProfETH2013$surf2_wet[c(1,5,9,13,17,21),2],diuProfETH2013$surf3_wet[c(1,5,9,13,17,21),2],diuProfETH2013$surf4_wet[c(1,5,9,13,17,21),2],diuProfETH2013$surf5_wet[c(1,5,9,13,17,21),2],diuProfETH2013$surf6_wet[c(1,5,9,13,17,21),2],diuProfETH2013$surf7_wet[c(1,5,9,13,17,21),2],diuProfETH2013$surf8_wet[c(1,5,9,13,17,21),2])

png(file=path_figs&'\\DebrisTempProfile2012.png', res = 160,width=2400,height=600)
par(mfcol=c(1,4), mar=c(4,5,1,2), oma=c(0,0,0,0),cex.axis = 1.5,cex.lab = 1.5)
#layout(matrix(c(1,3), nrow = 1, ncol = 3, byrow = TRUE))
plot(diuProfETH2012$surf8_dry[,2],type='l',col='black',ylim = c(0,18),xlab = 'Time [h]',ylab = expression('T'['debris']~~'[°C]'),lwd = 2)
points(diuProfETH2012$surf7_dry[,2],type='l',col='black',lwd = 2,lty = 2)
#points(diuProfETH2012$surf6_dry[,2],type='l',col='black',lwd = 2,lty = 3)
points(diuProfETH2012$surf5_dry[,2],type='l',col='black',lwd = 2,lty= 3)
points(diuProfETH2012$surf4_dry[,2],type='l',col='black',lwd = 2,lty= 4)
points(diuProfETH2012$surf3_dry[,2],type='l',col='black',lwd = 2,lty= 5)
points(diuProfETH2012$surf2_dry[,2],type='l',col='black',lwd = 2,lty= 6)
points(diuProfETH2012$surf1_dry[,2],type='l',col='black',lwd = 2,lty= 7)
legend('topleft',legend = c('surface','0.1 m','0.3 m','0.4 m','0.5 m','0.6 m','0.7 m (ice)'),cex= 1.5,bty = 'n',lty=c(1,2,3,4,5,6,7),lwd = 2)
grid(NULL,NULL)
plot(diuProfETH2012$surf8_wet[,2],type='l',col='red',ylim = c(0,18),xlab = 'Time [h]',ylab = expression('T'['debris']~~'[°C]'),lwd = 2)
points(diuProfETH2012$surf7_wet[,2],type='l',col='red',lwd = 2,lty = 2)
#points(diuProfETH2012$surf6_wet[,2],type='l',col='black',lwd = 2,lty = 3)
points(diuProfETH2012$surf5_wet[,2],type='l',col='red',lwd = 2,lty= 3)
points(diuProfETH2012$surf4_wet[,2],type='l',col='red',lwd = 2,lty= 4)
points(diuProfETH2012$surf3_wet[,2],type='l',col='red',lwd = 2,lty= 5)
points(diuProfETH2012$surf2_wet[,2],type='l',col='red',lwd = 2,lty= 6)
points(diuProfETH2012$surf1_wet[,2],type='l',col='red',lwd = 2,lty= 7)
grid(NULL,NULL)
matplot(t(Tdry_profile_2012),rev(d_2012), ylim = rev(range(d_2012)),type='o',pch = 1:6,xlim=c(-2,17),lty=c(1,2,3,4,5,6),col='black',cex =  1.5,ylab = 'Depth [m]',xlab = expression('T'['debris']~~'[°C]'))
matplot(t(Twet_profile_2012),rev(d_2012),type='o',pch = 1:6,xlim=c(-2,17),lty=c(1,2,3,4,5,6),col='red',add=T,cex=2)
legend('bottomright',legend = c('00:00','04:00','08:00','12:00','16:00','20:00'),cex=  1.5,bty = 'n',lty=c(1,2,3,4,5,6))
grid(NULL,NULL)
dev.off()

png(file=path_figs&'\\DebrisTempProfile2013.png', res = 160,width=2400,height=600)
par(mfcol=c(1,4), mar=c(4,5,1,2), oma=c(0,0,0,0),cex.axis = 1.5,cex.lab = 1.5)
#layout(matrix(c(1,3), nrow = 1, ncol = 3, byrow = TRUE))
plot(diuProfETH2013$surf1_dry[,2],type='l',col='black',ylim = c(0,20),xlab = 'Time [h]',ylab = expression('T'['debris']~~'[°C]'),lwd = 2)
points(diuProfETH2013$surf2_dry[,2],type='l',col='black',lwd = 2,lty = 2)
points(diuProfETH2013$surf3_dry[,2],type='l',col='black',lwd = 2,lty= 3)
points(diuProfETH2013$surf4_dry[,2],type='l',col='black',lwd = 2,lty= 4)
points(diuProfETH2013$surf5_dry[,2],type='l',col='black',lwd = 2,lty= 5)
points(diuProfETH2013$surf7_dry[,2],type='l',col='black',lwd = 2,lty= 6)
points(diuProfETH2013$surf8_dry[,2],type='l',col='black',lwd = 2,lty= 7)
legend('topleft',legend = c('surface','0.02 m','0.04 m','0.06 m','0.1 m','0.3 m','0.45 m'),cex= 1.5,bty = 'n',lty=c(1,2,3,4,5,6,7),lwd = 2)
grid(NULL,NULL)
plot(diuProfETH2013$surf1_wet[,2],type='l',col='red',ylim = c(0,20),xlab = 'Time [h]',ylab = expression('T'['debris']~~'[°C]'),lwd = 2)
points(diuProfETH2013$surf2_wet[,2],type='l',col='red',lwd = 2,lty = 2)
points(diuProfETH2013$surf3_wet[,2],type='l',col='red',lwd = 2,lty= 3)
points(diuProfETH2013$surf4_wet[,2],type='l',col='red',lwd = 2,lty= 4)
points(diuProfETH2013$surf5_wet[,2],type='l',col='red',lwd = 2,lty= 5)
points(diuProfETH2013$surf7_wet[,2],type='l',col='red',lwd = 2,lty= 6)
points(diuProfETH2013$surf8_wet[,2],type='l',col='red',lwd = 2,lty= 7)
grid(NULL,NULL)

matplot(t(Tdry_profile_2013),d_2013, ylim = rev(range(d_2013)),type='o',pch = 1:6,xlim=c(-2,17),lty=c(1,2,3,4,5,6),col='black',cex =  1.5,ylab = 'Depth [m]',xlab = expression('T'['debris']~~'[°C]'))
matplot(t(Twet_profile_2013),d_2013,type='o',pch = 1:6,xlim=c(-2,17),lty=c(1,2,3,4,5,6),col='red',add=T,cex=2)
legend('bottomright',legend = c('00:00','04:00','08:00','12:00','16:00','20:00'),cex=  1.5,bty = 'n',lty=c(1,2,3,4,5,6))
grid(NULL,NULL)
dev.off()

dTdt1_ETH <- list()
dT2dz2_ETH <- list()

# get diffusivity for different profiles parts
k_cal <- function(ListTrunk,profilecenter,profileup,profiledown,year,Season,dz,dT2){
  dT1_ETH <- c(NA,diff(profilecenter[ListTrunk$Season==Season,2]) / 3600)
  #browser()
  dT1_ETH[dT1_ETH>5*10^-5|dT1_ETH<=-5*10^-5]<-NA
  #dT1_ETH[abs(diff(profilecenter[year[ListTrunk$Season[year]==Season],2]))<=0.02]<-NA
  dT2_ETH <- (profileup[ListTrunk$Season==Season,2] - 2*profilecenter[ListTrunk$Season==Season,2] + profiledown[ListTrunk$Season==Season,2]) / (dz)^2
  if(dT2 == 1){
    return(dT1_ETH)
  }
  if(dT2 == 2){
    return(dT2_ETH)
  }
}

dTdt1_ETH$p1_dry <- k_cal(ETHprofile12,ETHprofile12$surf2,ETHprofile12$surf3,ETHprofile12$surf1,y_2012,1,d_2012[3]-d_2012[1],1)
dT2dz2_ETH$p1_dry <- k_cal(ETHprofile12,ETHprofile12$surf2,ETHprofile12$surf3,ETHprofile12$surf1,y_2012,1,d_2012[3]-d_2012[1],2)
dTdt1_ETH$p1_wet <- k_cal(ETHprofile12,ETHprofile12$surf2,ETHprofile12$surf3,ETHprofile12$surf1,y_2012,2,d_2012[3]-d_2012[1],1)
dT2dz2_ETH$p1_wet <- k_cal(ETHprofile12,ETHprofile12$surf2,ETHprofile12$surf3,ETHprofile12$surf1,y_2012,2,d_2012[3]-d_2012[1],2)
dTdt1_ETH$p2_dry <- k_cal(ETHprofile12,ETHprofile12$surf3,ETHprofile12$surf4,ETHprofile12$surf2,y_2012,1,d_2012[4]-d_2012[2],1)
dT2dz2_ETH$p2_dry <- k_cal(ETHprofile12,ETHprofile12$surf3,ETHprofile12$surf4,ETHprofile12$surf2,y_2012,1,d_2012[4]-d_2012[2],2)
dTdt1_ETH$p2_wet <- k_cal(ETHprofile12,ETHprofile12$surf3,ETHprofile12$surf4,ETHprofile12$surf2,y_2012,2,d_2012[4]-d_2012[2],1)
dT2dz2_ETH$p2_wet <- k_cal(ETHprofile12,ETHprofile12$surf3,ETHprofile12$surf4,ETHprofile12$surf2,y_2012,2,d_2012[4]-d_2012[2],2)
dTdt1_ETH$p3_dry <- k_cal(ETHprofile12,ETHprofile12$surf4,ETHprofile12$surf5,ETHprofile12$surf3,y_2012,1,d_2012[5]-d_2012[3],1)
dT2dz2_ETH$p3_dry <- k_cal(ETHprofile12,ETHprofile12$surf4,ETHprofile12$surf5,ETHprofile12$surf3,y_2012,1,d_2012[5]-d_2012[3],2)
dTdt1_ETH$p3_wet <- k_cal(ETHprofile12,ETHprofile12$surf4,ETHprofile12$surf5,ETHprofile12$surf3,y_2012,2,d_2012[5]-d_2012[3],1)
dT2dz2_ETH$p3_wet <- k_cal(ETHprofile12,ETHprofile12$surf4,ETHprofile12$surf5,ETHprofile12$surf3,y_2012,2,d_2012[5]-d_2012[3],2)
dTdt1_ETH$p4_dry <- k_cal(ETHprofile12,ETHprofile12$surf5,ETHprofile12$surf7,ETHprofile12$surf4,y_2012,1,d_2012[7]-d_2012[4],1)
dT2dz2_ETH$p4_dry <- k_cal(ETHprofile12,ETHprofile12$surf5,ETHprofile12$surf7,ETHprofile12$surf4,y_2012,1,d_2012[7]-d_2012[4],2)
dTdt1_ETH$p4_wet <- k_cal(ETHprofile12,ETHprofile12$surf5,ETHprofile12$surf7,ETHprofile12$surf4,y_2012,2,d_2012[7]-d_2012[4],1)
dT2dz2_ETH$p4_wet <- k_cal(ETHprofile12,ETHprofile12$surf5,ETHprofile12$surf7,ETHprofile12$surf4,y_2012,2,d_2012[7]-d_2012[4],2)
dTdt1_ETH$p5_dry <- k_cal(ETHprofile12,ETHprofile12$surf6,ETHprofile12$surf7,ETHprofile12$surf5,y_2012,1,d_2012[7]-d_2012[5],1)
dT2dz2_ETH$p5_dry <- k_cal(ETHprofile12,ETHprofile12$surf6,ETHprofile12$surf7,ETHprofile12$surf5,y_2012,1,d_2012[7]-d_2012[5],2)
dTdt1_ETH$p5_wet <- k_cal(ETHprofile12,ETHprofile12$surf6,ETHprofile12$surf7,ETHprofile12$surf5,y_2012,2,d_2012[7]-d_2012[5],1)
dT2dz2_ETH$p5_wet <- k_cal(ETHprofile12,ETHprofile12$surf6,ETHprofile12$surf7,ETHprofile12$surf5,y_2012,2,d_2012[7]-d_2012[5],2)
dTdt1_ETH$p6_dry <- k_cal(ETHprofile12,ETHprofile12$surf7,ETHprofile12$surf8,ETHprofile12$surf5,y_2012,1,d_2012[8]-d_2012[5],1)
dT2dz2_ETH$p6_dry <- k_cal(ETHprofile12,ETHprofile12$surf7,ETHprofile12$surf8,ETHprofile12$surf5,y_2012,1,d_2012[8]-d_2012[5],2)
dTdt1_ETH$p6_wet <- k_cal(ETHprofile12,ETHprofile12$surf7,ETHprofile12$surf8,ETHprofile12$surf5,y_2012,2,d_2012[8]-d_2012[5],1)
dT2dz2_ETH$p6_wet <- k_cal(ETHprofile12,ETHprofile12$surf7,ETHprofile12$surf8,ETHprofile12$surf5,y_2012,2,d_2012[8]-d_2012[5],2)

par(mfcol=c(1,6), mar=c(4,5,1,2), oma=c(0,0,0,0),cex.axis = 1.5,cex.lab = 1.5)
plot(dT2dz2_ETH$p1_dry,dTdt1_ETH$p1_dry)
points(dT2dz2_ETH$p1_wet,dTdt1_ETH$p1_wet,col='red')
plot(dT2dz2_ETH$p2_dry,dTdt1_ETH$p2_dry)
points(dT2dz2_ETH$p2_wet,dTdt1_ETH$p2_wet,col='red')
plot(dT2dz2_ETH$p3_dry,dTdt1_ETH$p3_dry,xlim=c(-20,20))
points(dT2dz2_ETH$p3_wet,dTdt1_ETH$p3_wet,col='red')
abline(lm(dTdt1_ETH$p3_dry~dT2dz2_ETH$p3_dry))
abline(lm(dTdt1_ETH$p3_wet~dT2dz2_ETH$p3_wet),col='red')
plot(dT2dz2_ETH$p4_dry,dTdt1_ETH$p4_dry)
points(dT2dz2_ETH$p4_wet,dTdt1_ETH$p4_wet,col='red')
plot(dT2dz2_ETH$p5_dry,dTdt1_ETH$p5_dry)
points(dT2dz2_ETH$p5_wet,dTdt1_ETH$p5_wet,col='red')
plot(dT2dz2_ETH$p6_dry,dTdt1_ETH$p6_dry)
points(dT2dz2_ETH$p6_wet,dTdt1_ETH$p6_wet,col='red')
abline(lm(dTdt1_ETH$p6_dry~dT2dz2_ETH$p6_dry))
abline(lm(dTdt1_ETH$p6_wet~dT2dz2_ETH$p6_wet),col='red')

fit1 = summary(lm(dTdt1_ETH$p3_wet[dayhours(ETHprofile12,2,2,NULL)]~dT2dz2_ETH$p3_wet[dayhours(ETHprofile12,2,2,NULL)]))
fit3 = summary(lm(dTdt1_ETH$p4_wet[dayhours(ETHprofile12,2,2,NULL)]~dT2dz2_ETH$p4_wet[dayhours(ETHprofile12,2,2,NULL)]))
fit4 = summary(lm(dTdt1_ETH$p6_wet[dayhours(ETHprofile12,2,2,NULL)]~dT2dz2_ETH$p6_wet[dayhours(ETHprofile12,2,2,NULL)]))
fit2 = summary(lm(dTdt1_ETH$p2_dry[dayhours(ETHprofile12,1,2,NULL)]~dT2dz2_ETH$p2_dry[dayhours(ETHprofile12,1,2,NULL)]))

diff_ETH12wet <- c(fit1$coefficients[2,c(1,2)],fit3$coefficients[2,c(1,2)],fit4$coefficients[2,c(1,2)])*10^6

fit1 = summary(lm(dTdt1_ETH$p3_wet~dT2dz2_ETH$p3_wet))
fit3 = summary(lm(dTdt1_ETH$p4_wet~dT2dz2_ETH$p4_wet))
fit4 = summary(lm(dTdt1_ETH$p6_wet~dT2dz2_ETH$p6_wet))

fit2 = summary(lm(dTdt1_ETH$p2_dry~dT2dz2_ETH$p2_dry))

diff_ETH12wet_whole <- c(fit1$coefficients[2,c(1,2)],fit3$coefficients[2,c(1,2)],fit4$coefficients[2,c(1,2)])*10^6

diff_ETH12 <- c(fit1$coefficients[2,c(1,2)],fit2$coefficients[2,c(1,2)])*10^6

dTdt1_ETH13 <- list()
dT2dz2_ETH13 <- list()

dTdt1_ETH13$p1_dry <- k_cal(ETHprofile13,ETHprofile13$surf2,ETHprofile13$surf1,ETHprofile13$surf3,y_2013,1,d_2013[3]-d_2013[1],1)
dT2dz2_ETH13$p1_dry <- k_cal(ETHprofile13,ETHprofile13$surf2,ETHprofile13$surf1,ETHprofile13$surf3,y_2013,1,d_2013[3]-d_2013[1],2)

dTdt1_ETH13$p1_wet <- k_cal(ETHprofile13,ETHprofile13$surf2,ETHprofile13$surf1,ETHprofile13$surf3,y_2013,2,d_2013[3]-d_2013[1],1)
dT2dz2_ETH13$p1_wet <- k_cal(ETHprofile13,ETHprofile13$surf2,ETHprofile13$surf1,ETHprofile13$surf3,y_2013,2,d_2013[3]-d_2013[1],2)
dTdt1_ETH13$p2_dry <- k_cal(ETHprofile13,ETHprofile13$surf3,ETHprofile13$surf2,ETHprofile13$surf4,y_2013,1,d_2013[4]-d_2013[2],1)
dT2dz2_ETH13$p2_dry <- k_cal(ETHprofile13,ETHprofile13$surf3,ETHprofile13$surf2,ETHprofile13$surf4,y_2013,1,d_2013[4]-d_2013[2],2)
dTdt1_ETH13$p2_wet <- k_cal(ETHprofile13,ETHprofile13$surf3,ETHprofile13$surf2,ETHprofile13$surf4,y_2013,2,d_2013[4]-d_2013[2],1)
dT2dz2_ETH13$p2_wet <- k_cal(ETHprofile13,ETHprofile13$surf3,ETHprofile13$surf2,ETHprofile13$surf4,y_2013,2,d_2013[4]-d_2013[2],2)
dTdt1_ETH13$p3_dry <- k_cal(ETHprofile13,ETHprofile13$surf4,ETHprofile13$surf3,ETHprofile13$surf5,y_2013,1,d_2013[5]-d_2013[3],1)
dT2dz2_ETH13$p3_dry <- k_cal(ETHprofile13,ETHprofile13$surf4,ETHprofile13$surf3,ETHprofile13$surf5,y_2013,1,d_2013[5]-d_2013[3],2)
dTdt1_ETH13$p3_wet <- k_cal(ETHprofile13,ETHprofile13$surf4,ETHprofile13$surf3,ETHprofile13$surf5,y_2013,2,d_2013[5]-d_2013[3],1)
dT2dz2_ETH13$p3_wet <- k_cal(ETHprofile13,ETHprofile13$surf4,ETHprofile13$surf3,ETHprofile13$surf5,y_2013,2,d_2013[5]-d_2013[3],2)
dTdt1_ETH13$p4_dry <- k_cal(ETHprofile13,ETHprofile13$surf5,ETHprofile13$surf4,ETHprofile13$surf7,y_2013,1,d_2013[7]-d_2013[4],1)
dT2dz2_ETH13$p4_dry <- k_cal(ETHprofile13,ETHprofile13$surf5,ETHprofile13$surf4,ETHprofile13$surf7,y_2013,1,d_2013[7]-d_2013[4],2)
dTdt1_ETH13$p4_wet <- k_cal(ETHprofile13,ETHprofile13$surf5,ETHprofile13$surf4,ETHprofile13$surf7,y_2013,2,d_2013[7]-d_2013[4],1)
dT2dz2_ETH13$p4_wet <- k_cal(ETHprofile13,ETHprofile13$surf5,ETHprofile13$surf4,ETHprofile13$surf7,y_2013,2,d_2013[7]-d_2013[4],2)
dTdt1_ETH13$p5_dry <- k_cal(ETHprofile13,ETHprofile13$surf6,ETHprofile13$surf5,ETHprofile13$surf7,y_2013,1,d_2013[7]-d_2013[5],1)
dT2dz2_ETH13$p5_dry <- k_cal(ETHprofile13,ETHprofile13$surf6,ETHprofile13$surf5,ETHprofile13$surf7,y_2013,1,d_2013[7]-d_2013[5],2)
dTdt1_ETH13$p5_wet <- k_cal(ETHprofile13,ETHprofile13$surf6,ETHprofile13$surf5,ETHprofile13$surf7,y_2013,2,d_2013[7]-d_2013[5],1)
dT2dz2_ETH13$p5_wet <- k_cal(ETHprofile13,ETHprofile13$surf6,ETHprofile13$surf5,ETHprofile13$surf7,y_2013,2,d_2013[7]-d_2013[5],2)
dTdt1_ETH13$p6_dry <- k_cal(ETHprofile13,ETHprofile13$surf7,ETHprofile13$surf5,ETHprofile13$surf8,y_2013,1,d_2013[5]-d_2013[8],1)
dT2dz2_ETH13$p6_dry <- k_cal(ETHprofile13,ETHprofile13$surf7,ETHprofile13$surf5,ETHprofile13$surf8,y_2013,1,d_2013[5]-d_2013[8],2)
dTdt1_ETH13$p6_wet <- k_cal(ETHprofile13,ETHprofile13$surf7,ETHprofile13$surf5,ETHprofile13$surf8,y_2013,2,d_2013[5]-d_2013[8],1)
dT2dz2_ETH13$p6_wet <- k_cal(ETHprofile13,ETHprofile13$surf7,ETHprofile13$surf5,ETHprofile13$surf8,y_2013,2,d_2013[5]-d_2013[8],2)
dT2dz2_ETH13$p6_dry[dT2dz2_ETH13$p6_dry>10] <- NA
dT2dz2_ETH13$p6_wet[dT2dz2_ETH13$p6_wet>10] <- NA
dTdt1_ETH13$p6_wet[dTdt1_ETH13$p6_wet>-1.5*10^(-5)] <- NA
dTdt1_ETH13$p6_dry[dTdt1_ETH13$p6_dry>-1.5*10^(-5)] <- NA

par(mfcol=c(1,6), mar=c(4,5,1,2), oma=c(0,0,0,0),cex.axis = 1.5,cex.lab = 1.5)
plot(dT2dz2_ETH13$p1_dry,dTdt1_ETH13$p1_dry)
points(dT2dz2_ETH13$p1_wet,dTdt1_ETH13$p1_wet,col='red')
plot(dT2dz2_ETH13$p2_dry,dTdt1_ETH13$p2_dry)
points(dT2dz2_ETH13$p2_wet,dTdt1_ETH13$p2_wet,col='red')
plot(dT2dz2_ETH13$p3_dry,dTdt1_ETH13$p3_dry)
points(dT2dz2_ETH13$p3_wet,dTdt1_ETH13$p3_wet,col='red')
plot(dT2dz2_ETH13$p4_dry,dTdt1_ETH13$p4_dry)
points(dT2dz2_ETH13$p4_wet,dTdt1_ETH13$p4_wet,col='red')
plot(dT2dz2_ETH13$p5_dry,dTdt1_ETH13$p5_dry)
points(dT2dz2_ETH13$p5_wet,dTdt1_ETH13$p5_wet,col='red')
plot(dT2dz2_ETH13$p6_dry,dTdt1_ETH13$p6_dry)
points(dT2dz2_ETH13$p6_wet,dTdt1_ETH13$p6_wet,col='red')
abline(lm(dTdt1_ETH13$p6_wet~dT2dz2_ETH13$p6_wet),col='red')
abline(lm(dTdt1_ETH13$p6_dry~dT2dz2_ETH13$p6_dry))

fit1 = summary(lm(dTdt1_ETH13$p6_wet[dayhours(ETHprofile13,2,2,NULL)]~dT2dz2_ETH13$p6_wet[dayhours(ETHprofile13,2,2,NULL)]))
fit2 = summary(lm(dTdt1_ETH13$p6_dry[dayhours(ETHprofile13,1,2,NULL)]~dT2dz2_ETH13$p6_dry[dayhours(ETHprofile13,1,2,NULL)]))

diff_ETH13 <- c(fit1$coefficients[2,c(1,2)],fit2$coefficients[2,c(1,2)])*10^6

#############
# Read additional data for calculation of conductivity
#############

DebrisDataelse <- read.csv(path&'\\'&'RawData\\DebrisData\\AggregatedDebrisDataRanges.csv')

meanBulk_rho <- mean(DebrisDataelse$bulk.density..kg.m3.,na.rm=T)
sdBulk_rho <- sd(DebrisDataelse$bulk.density..kg.m3.,na.rm=T)

meanmoist <- mean(DebrisDataelse$volumetric.soil.moisture..m3.m3.,na.rm=T)
sdmoist <- sd(DebrisDataelse$volumetric.soil.moisture..m3.m3.,na.rm=T)

meanpor <- mean(DebrisDataelse$calculated.porosity....,na.rm=T)
sdpor <- sd(DebrisDataelse$calculated.porosity....,na.rm=T)

png(file=path_figs&'\\SoilMoisture.png', res = 160,width=700,height=600)
par(mfcol=c(1,1), mar=c(4,5,1,2), oma=c(0,0,0,0),cex.axis = 1,cex.lab = 1)
plot(DebrisDataelse$volumetric.soil.moisture..m3.m3.,DebrisDataelse$measured.Depth/DebrisDataelse$Total.Depth,ylim = rev(range(DebrisDataelse$measured.Depth/DebrisDataelse$Total.Depth,na.rm=T)),xlab = expression('soil moisture ['~m^{3}~''~ m^{-3} ~']'),ylab = expression('relative depth [-]'),pch=25, bg="blue",col='red')
grid(NULL,NULL)
dev.off()

# convert from diffusivity to conductivity
phi_deb <- meanpor
phiDeb_lit <- c(0.2,0.3,0.43,0.33,0.3) # Nicholson2006 (2x), Popovnin2002, Nicholson2012, Conway and Rasmussen2000
c_deb <- 948
rho_d <- meanBulk_rho       # 19 samples 

rho_rock <- 2650# kg m-3
c_rock <- c(750,948,804,900,890,948) # J kg-1 K-1

c_w <- 4181.3 
c_air <- 1003.5
rho_w <<- 999.7
rho_air <<- 0.819

volum_heat_cap_void <- rho_w * c_w * meanmoist/0.2 + rho_air * c_air * (1-meanmoist/0.2)

volum_heat_cap_wet <- rho_w * c_w
volum_heat_cap_dry <- rho_air * c_air

juenDiff <- mean(c(0.393,0.364,0.304))
k_diff_lit_mean <- c(0.3,0.38,juenDiff,0.95,0.689,0.6,0.9,0.82,0.89,0.82,0.6,0.59,0.42,0.61) / 10^6 # m2 s-1

k_diff_lit_sd <- c(0.05,0.02,0.01,0.1,0.1,0.01,0.01,0.01,0.01,0.03,0.05) / 10^6

k_diff_lit <- mean(c(0.3,0.38,0.95,0.689,0.6,0.9))

k_diff_lir_wet <- c(diff_UU_1_whole[3],diff_UU_2_whole[3])
k_diff_lir_dry <- c(diff_UU_1_whole[1],diff_UU_2_whole[1],0.61)

cond_lit_dry <- c(0.585,0.637,1.29)
cond_lit_wet <- c(1.669,1.776)
cond_lit <- c(0.94,0.96,1.33,1.62,1.523,1.29,1.28,0.85) # Reid2010
cond_tot_distrib <- k_diff_lit_mean * (rho_rock * 890 * (1-phi_deb) + volum_heat_cap_void * phi_deb)
cond_tot_phi0.3 <- k_diff_lit_mean * (rho_rock * 890 * (1-0.3) + volum_heat_cap_void * 0.3)

cond_tot_Lirung_dry1 <- mean(k_diff_lir_dry[1])/10^6 * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)
cond_tot_Lirung_dry2 <- mean(k_diff_lir_dry[2])/10^6 * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)
cond_tot_Lirung_dry3 <- mean(k_diff_lir_dry[3])/10^6 * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)
cond_tot_Lirung_wet1 <- mean(k_diff_lir_wet[1])/10^6 * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)
cond_tot_Lirung_wet2 <- mean(k_diff_lir_wet[2])/10^6 * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)


cond_tot_wet <- k_diff_lit_mean * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_wet * phi_deb)
cond_tot_dry <- k_diff_lit_mean * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_dry * phi_deb)
cond_tot <- k_diff_lit_mean * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)


library(scales)
png(file=path_figs&'\\Conductivity.png', res = 160,width=700,height=600)
par(mfcol=c(1,1), mar=c(4,5,1,2), oma=c(0,0,0,0),cex.axis = 1,cex.lab = 1)
boxplot(cbind(cond_tot_distrib,cond_tot_wet,cond_tot_dry),names=c('mix','wet','dry'),ylab = expression('thermal conductivity [J '~m^{-1}~''~ K^{-1} ~']'))
  myjitter<-jitter(rep(2, length(cond_lit_wet)), amount=0.2)
  points(myjitter, cond_lit_wet, pch=20, col=alpha('blue',0.4) ,cex=2,xlab='',xaxt='n')
  myjitter<-jitter(rep(3, length(cond_lit_dry)), amount=0.2)
  points(myjitter, cond_lit_dry, pch=20, col=alpha('red',0.4) ,cex=2,xlab='',xaxt='n')
  myjitter<-jitter(rep(1, length(cond_lit)), amount=0.2)
  points(myjitter, cond_lit, pch=20, col=alpha('grey',0.4) ,cex=2,xlab='',xaxt='n')
  
  myjitter<-jitter(rep(1, length(cond_tot_Lirung_dry1)), amount=0.2)
  points(myjitter, cond_tot_Lirung_dry1 , pch=10, col='red' ,cex=2,xlab='',xaxt='n')
  myjitter<-jitter(rep(1, length(cond_tot_Lirung_dry2)), amount=0.2)
  points(myjitter, cond_tot_Lirung_dry2 , pch=10, col='red' ,cex=2,xlab='',xaxt='n')
  myjitter<-jitter(rep(1, length(cond_tot_Lirung_dry3)), amount=0.2)
  points(myjitter, cond_tot_Lirung_dry3 , pch=10, col='red' ,cex=2,xlab='',xaxt='n')
  myjitter<-jitter(rep(1, length(cond_tot_Lirung_wet1)), amount=0.2)
  points(myjitter, cond_tot_Lirung_wet1 , pch=10, col='blue' ,cex=2,xlab='',xaxt='n')
  myjitter<-jitter(rep(1, length(cond_tot_Lirung_wet2)), amount=0.2)
  points(myjitter, cond_tot_Lirung_wet2 , pch=10, col='blue' ,cex=2,xlab='',xaxt='n')
  
  dev.off()