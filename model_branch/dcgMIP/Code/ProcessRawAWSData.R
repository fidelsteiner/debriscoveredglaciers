

# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# packages (if not installed yet: install.packages('examplePackage')
library(rgdal)
library(rgeos)
library(maptools)
library(raster)
library(ggplot2)
library(foreach)
library(lubridate)
library(compare)
library(colorRamps)
library(data.table)
library(circular)
library(parallel)
library(snowfall)
library(truncnorm)
library(raster)
library(rlecuyer)
library(forecast)
library(rasterVis)
library(R.utils)
library(zoo)

path_subcode <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\'
# Paths for Extra Codes
file.sources = list.files(path_subcode, 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

LirungRaw <- read.csv('F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\Lirung2014_DCWGdata_5min.csv',header = T)

Sys.setenv(TZ='Asia/Kathmandu')
LirungRaw$Time_Str <- as.POSIXct(paste(LirungRaw$Date,LirungRaw$Time), format="%m/%d/%Y  %H:%M")

LirungData <- list()
LirungData$Ta <- TSaggregate(LirungRaw$Ta,LirungRaw$Time_Str,60,2014,'mean')
LirungData$RH <- TSaggregate(LirungRaw$rH,LirungRaw$Time_Str,60,2014,'mean')
LirungData$ws <- TSaggregate(LirungRaw$ws,LirungRaw$Time_Str,60,2014,'mean')
LirungData$wd <- TSaggregate(LirungRaw$wd,LirungRaw$Time_Str,60,2014,'mean')
LirungData$sd <- TSaggregate(LirungRaw$height,LirungRaw$Time_Str,60,2014,'mean')
LirungData$SWIN <- TSaggregate(LirungRaw$SWin,LirungRaw$Time_Str,60,2014,'mean')
LirungData$SWOUT <- TSaggregate(LirungRaw$SWout,LirungRaw$Time_Str,60,2014,'mean')
LirungData$LWIN <- TSaggregate(LirungRaw$LWin,LirungRaw$Time_Str,60,2014,'mean')
LirungData$LWOUT <- TSaggregate(LirungRaw$LWout,LirungRaw$Time_Str,60,2014,'mean')
LirungData$TimeStr <- as.POSIXlt(LirungData$LWOUT[,1],origin = '1970-01-01')

Lirung_hourlyData <- cbind(data.frame(LirungData$TimeStr,LirungData$Ta[,2],LirungData$RH[,2],LirungData$ws[,2],LirungData$wd[,2],LirungData$sd[,2],LirungData$SWIN[,2],LirungData$SWOUT[,2],LirungData$LWIN[,2],LirungData$LWOUT[,2]))
colnames(Lirung_hourlyData) <-c('Time','Tair','RH','WS','WD','SnowDepth','SWin','SWout','LWin','LWout')
write.table(Lirung_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\Lirung2014_DCWGdata_60min.csv.csv',sep = ",",row.names = FALSE)



ARORaw <- read.csv('F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\1_ARO_2010\\ARO_2010\\ARO_2010_60.csv',header = T)

Sys.setenv(TZ='UTC')
ARORaw$TimeStr <- as.POSIXct(paste(paste(ARORaw$YEAR,'/',ARORaw$MONTH,'/',ARORaw$DAY,sep=""),ARORaw$HOUR),format ='%Y/%m/%d %H')
ARO_hourlyData <- cbind(data.frame(ARORaw$TimeStr,ARORaw$T,ARORaw$RH,ARORaw$FF,ARORaw$DIR,ARORaw$SNOW,ARORaw$SWIN,ARORaw$SWOUT,ARORaw$LWIN,ARORaw$LWOUT))

colnames(ARO_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(ARO_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\ARO_hourly_FINAL.csv',sep = ",",row.names = FALSE)

CNURaw <- read.csv('F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\2_CN-5470_20152017\\CN-5470_20152017\\CN-5470_20152017_hourly.csv',header = T)

Sys.setenv(TZ='UTC')
CNURaw$TimeStr <- as.POSIXct(paste(paste(CNURaw$YEAR,'/',CNURaw$MONTH,'/',CNURaw$DAY,sep=""),CNURaw$HOUR),format ='%Y/%m/%d %H')
CNU_hourlyData <- cbind(data.frame(CNURaw$TimeStr,CNURaw$T,CNURaw$RH,CNURaw$FF,CNURaw$DIR,CNURaw$SNOW,CNURaw$SWIN,CNURaw$SWOUT,CNURaw$LWIN,CNURaw$LWOUT))

colnames(CNU_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(CNU_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\CNU_hourly_FINAL.csv',sep = ",",row.names = FALSE)

DJARaw <- read.csv('F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\3_DJA_20072008\\DJA_20072008\\DJA_20072008_60.csv',header = T)

Sys.setenv(TZ='UTC')
DJARaw$TimeStr <- as.POSIXct(paste(paste(DJARaw$YEAR,'/',DJARaw$MONTH,'/',DJARaw$DAY,sep=""),DJARaw$HOUR),format ='%Y/%m/%d %H')
DJA_hourlyData <- cbind(data.frame(DJARaw$TimeStr,DJARaw$T,DJARaw$RH,DJARaw$FF,DJARaw$DIR,DJARaw$SNOW,DJARaw$SWIN,DJARaw$SWOUT,DJARaw$LWIN,DJARaw$LWOUT))

colnames(DJA_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(DJA_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\DJA_hourly_FINAL.csv',sep = ",",row.names = FALSE)

LIRRaw <- read.csv("F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\4_LIR_2014\\LIR_2014_60.csv",header = T)

Sys.setenv(TZ='UTC')
LIRRaw$TimeStr <- as.POSIXct(paste(paste(LIRRaw$YEAR,'/',LIRRaw$MONTH,'/',LIRRaw$DAY,sep=""),LIRRaw$HOUR),format ='%Y/%m/%d %H')
LIR_hourlyData <- cbind(data.frame(LIRRaw$TimeStr,LIRRaw$T,LIRRaw$RH,LIRRaw$FF,LIRRaw$DIR,LIRRaw$SNOW,LIRRaw$SWIN,LIRRaw$SWOUT,LIRRaw$LWIN,LIRRaw$LWOUT))

colnames(LIR_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(LIR_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\LIR_hourly_FINAL.csv',sep = ",",row.names = FALSE)

MIARaw <- read.csv("F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\5_MIA_2005\\MIA_2005_60.csv",header = T)

Sys.setenv(TZ='UTC')
MIARaw$TimeStr <- as.POSIXct(paste(paste(MIARaw$YEAR,'/',MIARaw$MONTH,'/',MIARaw$DAY,sep=""),MIARaw$HOUR),format ='%Y/%m/%d %H')
MIA_hourlyData <- cbind(data.frame(MIARaw$TimeStr,MIARaw$T,MIARaw$RH,MIARaw$FF,MIARaw$DIR,MIARaw$SNOW,MIARaw$SWIN,MIARaw$SWOUT,MIARaw$LWIN,MIARaw$LWOUT))

colnames(MIA_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(MIA_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\MIA_hourly_FINAL.csv',sep = ",",row.names = FALSE)


PIRRaw <- read.csv("F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\6_PIR_20142015\\PIR_20142015\\PIR_20142015_60.csv",header = T)

Sys.setenv(TZ='UTC')
PIRRaw$TimeStr <- as.POSIXct(paste(paste(PIRRaw$YEAR,'/',PIRRaw$MONTH,'/',PIRRaw$DAY,sep=""),PIRRaw$HOUR),format ='%Y/%m/%d %H')
PIR_hourlyData <- cbind(data.frame(PIRRaw$TimeStr,PIRRaw$T,PIRRaw$RH,PIRRaw$FF,PIRRaw$DIR,PIRRaw$SNOW,PIRRaw$SWIN,PIRRaw$SWOUT,PIRRaw$LWIN,PIRRaw$LWOUT))

colnames(PIR_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(PIR_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\PIR_hourly_FINAL.csv',sep = ",",row.names = FALSE)

SDFRaw <- read.csv("F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\7_SDF_2016\\SDF_2016\\SDF_2016_60.csv",header = T)

Sys.setenv(TZ='UTC')
SDFRaw$TimeStr <- as.POSIXct(paste(paste(SDFRaw$YEAR,'/',SDFRaw$MONTH,'/',SDFRaw$DAY,sep=""),SDFRaw$HOUR),format ='%Y/%m/%d %H')
SDF_hourlyData <- cbind(data.frame(SDFRaw$TimeStr,SDFRaw$T,SDFRaw$RH,SDFRaw$FF,SDFRaw$DIR,SDFRaw$SNOW,SDFRaw$SWIN,SDFRaw$SWOUT,SDFRaw$LWIN,SDFRaw$LWOUT))

colnames(SDF_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(SDF_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\SDF_hourly_FINAL.csv',sep = ",",row.names = FALSE)

TAPRaw <- read.csv("F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\8_TAP_20132015\\TAP_20132015\\TAP_20142015_60.csv",header = T)

Sys.setenv(TZ='UTC')
TAPRaw$TimeStr <- as.POSIXct(paste(paste(TAPRaw$YEAR,'/',TAPRaw$MONTH,'/',TAPRaw$DAY,sep=""),TAPRaw$HOUR),format ='%Y/%m/%d %H')
TAP_hourlyData <- cbind(data.frame(TAPRaw$TimeStr,TAPRaw$T,TAPRaw$RH,TAPRaw$FF,TAPRaw$DIR,TAPRaw$SNOW,TAPRaw$SWIN,TAPRaw$SWOUT,TAPRaw$LWIN,TAPRaw$LWOUT))

colnames(TAP_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(TAP_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\TAP_hourly_FINAL.csv',sep = ",",row.names = FALSE)

TASRaw <- read.csv("F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\9_TAS_20172018\\TAS_20172018\\TAS_20172018_60.csv",header = T)

Sys.setenv(TZ='UTC')
TASRaw$TimeStr <- as.POSIXct(paste(paste(TASRaw$YEAR,'/',TASRaw$MONTH,'/',TASRaw$DAY,sep=""),TASRaw$HOUR),format ='%Y/%m/%d %H')
TAS_hourlyData <- cbind(data.frame(TASRaw$TimeStr,TASRaw$T,TASRaw$RH,TASRaw$FF,TASRaw$DIR,TASRaw$SNOW,TASRaw$SWIN,TASRaw$SWOUT,TASRaw$LWIN,TASRaw$LWOUT))

colnames(TAS_hourlyData) <-c('Time','Tair','rH','ws','wd','SnowDepth','SWin','SWout','LWin','LWout')
write.table(TAS_hourlyData,'F:\\PhD\\Research\\EB_DCG\\DCWG\\Data\\AWSPrep\\TAS_hourly_FINAL.csv',sep = ",",row.names = FALSE)
