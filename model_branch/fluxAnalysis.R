################################################################################
# fluxAnalysis - overview over all fluxes on Lirung stations
# 
# fluxAnalysis.R
#
# ReadMe:
#
# Created:          2019/03/11
# Latest Revision:  2019/03/11
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
path_data <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData'
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\TSaggregate.R")

climData2012 <- 'AWSLirung_onglacier_2012_Data.csv'
climData2013 <- 'AWSLirung_onglacier_2013_Data.csv'
climData2014 <- 'AWSLirung_onglacier_2014_Data.csv'
climData2016 <- 'Lirung_low_frequency_data.csv'
# Open AWS Data Sheet
climInput2012 <- read.csv(path_data&'\\'&climData2012,header = T)
Sys.setenv(TZ='Asia/Kathmandu')
climInput2012$Time_Str <- as.POSIXct(paste(as.character(climInput2012$Date),as.character(climInput2012$Time)), format="%m/%d/%Y  %H:%M:%S")
climInput2012$Time_Num <- as.numeric(climInput2012$Time_Str)

climInput2013 <- read.csv(path_data&'\\'&climData2013,header = T)
Sys.setenv(TZ='Asia/Kathmandu')
climInput2013$Time_Str <- as.POSIXct(paste(as.character(climInput2013$DATE),as.character(climInput2013$TIME)), format="%m/%d/%Y  %H:%M:%S")
climInput2013$Time_Num <- as.numeric(climInput2013$Time_Str)

climInput2014 <- read.csv(path_data&'\\'&climData2014,header = T)
Sys.setenv(TZ='Asia/Kathmandu')
climInput2014$Time_Str <- as.POSIXct(paste(as.character(climInput2014$Date),as.character(climInput2014$Time)), format="%m/%d/%Y  %H:%M:%S")
climInput2014$Time_Num <- as.numeric(climInput2014$Time_Str)

hourly2014 <- list()
climInput2014$SWin[climInput2014$SWin<0] <- 0
hourly2014$SWin <- TSaggregate(climInput2014$SWin,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$SWout <- TSaggregate(climInput2014$SWout,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$LWin <- TSaggregate(climInput2014$LWin,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$LWout <- TSaggregate(climInput2014$LWout,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$Ta <- TSaggregate(climInput2014$Ta,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$RH <- TSaggregate(climInput2014$rH,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$ws <- TSaggregate(climInput2014$ws,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$wd <- TSaggregate(climInput2014$wd,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$Ts <- TSaggregate(climInput2014$Ts,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2014$Date <- as.POSIXct(TSaggregate(climInput2014$SWin,climInput2014$Time_Str,60,'2014','mean')[,1],origin='1970-01-01')

Tab2014 <- cbind(data.frame(hourly2014$Date,hourly2014$Date,hourly2014$Ta,hourly2014$Ts,hourly2014$ws,hourly2014$wd,hourly2014$SWin,hourly2014$SWout,hourly2014$LWin,hourly2014$LWout, stringsAsFactors = FALSE))
colnames(Tab2014) <-c('DATE','TIME','TAIR','TS','WSPD','WDIR','KINC','KOUT','LINC','LOUT')
write.csv(Tab2014,path_data&'\\AWSLirung_onglacier_2014_Data_hourly.csv')

climInput2016 <- read.csv(path_data&'\\'&climData2016,header = T)
Sys.setenv(TZ='Asia/Kathmandu')
climInput2016$Time_Str <- as.POSIXct(paste(as.character(climInput2016$datetime),as.character(climInput2016$datetime.1)), format="%Y-%m-%d  %H:%M:%S")
climInput2016$Time_Num <- as.numeric(climInput2016$Time_Str)

hourly2016 <- list()
climInput2016$SWinc_Avg[climInput2016$SWinc_Avg<0] <- 0
hourly2016$SWin <- TSaggregate(climInput2016$SWinc_Avg,climInput2016$Time_Str,60,'2016','mean')[,2]
hourly2016$SWout <- TSaggregate(climInput2016$SWout_Avg,climInput2016$Time_Str,60,'2016','mean')[,2]
hourly2016$LWin <- TSaggregate(climInput2016$LWinc_co_Avg,climInput2016$Time_Str,60,'2016','mean')[,2]
hourly2016$LWout <- TSaggregate(climInput2016$LWout_co_Avg,climInput2016$Time_Str,60,'2016','mean')[,2]
hourly2016$Ta <- TSaggregate(climInput2016$T_air_Avg,climInput2016$Time_Str,60,'2016','mean')[,2]
hourly2016$RH <- TSaggregate(climInput2016$RH_air_Avg,climInput2016$Time_Str,60,'2016','mean')[,2]
#hourly2016$ws <- TSaggregate(climInput2016$,climInput2014$Time_Str,60,'2014','mean')[,2]
#hourly2016$wd <- TSaggregate(climInput2014$wd,climInput2014$Time_Str,60,'2014','mean')[,2]
#hourly2016$Ts <- TSaggregate(climInput2014$Ts,climInput2014$Time_Str,60,'2014','mean')[,2]
hourly2016$Date <- as.POSIXct(TSaggregate(climInput2016$SWinc_Avg,climInput2016$Time_Str,60,'2016','mean')[,1],origin='1970-01-01')

Tab2016 <- cbind(data.frame(hourly2016$Date,hourly2016$Date,hourly2016$Ta,hourly2016$SWin,hourly2016$SWout,hourly2016$LWin,hourly2016$LWout, stringsAsFactors = FALSE))
colnames(Tab2016) <-c('DATE','TIME','TAIR','KINC','KOUT','LINC','LOUT')
write.csv(Tab2016,path_data&'\\AWSLirung_onglacier_2016_Data_hourly.csv')


####################################
# Albedo in all years
####################################

albedo2012 <- climInput2012$SWout/climInput2012$SWin
albedo2012[albedo2012>1] <- NA
albedo2012[albedo2012<0] <- NA
albedo2012_d <- TSaggregate(albedo2012,climInput2012$Time_Str,3600,'2012','median')[,2]

albedo2013 <- climInput2013$KUPW/climInput2013$KINC
albedo2013[albedo2013>1] <- NA
albedo2013[albedo2013<0] <- NA
albedo2013_d <- TSaggregate(albedo2013,climInput2013$Time_Str,3600,'2013','median')[,2]

albedo2014 <- climInput2014$SWout/climInput2014$SWin
albedo2014[albedo2014>1] <- NA
albedo2014[albedo2014<0] <- NA
albedo2014_d <- TSaggregate(albedo2014,climInput2014$Time_Str,3600,'2014','median')[,2]

albedo2016 <- climInput2016$SWout/climInput2016$SWin
albedo2016[albedo2016>1] <- NA
albedo2016[albedo2016<0] <- NA
albedo2016_d <- TSaggregate(albedo2016,climInput2016$Time_Str,3600,'2016','median')[,2]