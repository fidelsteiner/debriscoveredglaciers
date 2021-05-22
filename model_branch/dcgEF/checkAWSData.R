################################################################################
# D1EB - point-scale debris Energy Balance Model
# 
# d1Deb.R
#
# ReadMe:
# Based on Reid and Brock, 2010
# Code adapted for global debris EB point-scale model
#
# % NOTE: To run this successfully you will need the files: 
#% 'isleapyear.m', 'debristemp.m', 'Lup.m', 'G.m', 'H.m', 'LE.m', 'P.m', 'loadEBdata.m' and 'diucyc.m'
#% They must all be saved in the same directory as this routine
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

##################
# packages (if not installed yet: install.packages('examplePackage')
##################
library(rgdal)
library(rgeos)
library(maptools)
library(raster)
library(foreach)
library(lubridate)
library(compare)
library(colorRamps)
library(data.table)
library(circular)
library(parallel)
library(snowfall)
library(truncnorm)
library(rlecuyer)
library(forecast)

# Create and specify the path of all files for the Model on your station. The folder needs the follwoing subfolders:
# 'Code','Figures','Data','Temp', 'Output'

path <- 'F:\\PHD\\Research\\Collaborations\\EMiles\\CaseStudies_EB'
path_code <- 'F:\\PHD\\Research\\Collaborations\\EMiles\\CaseStudies_EB\\Code'
path_figs <- 'F:\\PHD\\Research\\Collaborations\\EMiles\\CaseStudies_EB\\Figures'
path_data <- 'F:\\PHD\\Research\\Collaborations\\EMiles\\CaseStudies_EB\\Data'
path_temp <- 'F:\\PHD\\Research\\Collaborations\\EMiles\\CaseStudies_EB\\Temp'
path_output <- 'F:\\PHD\\Research\\Collaborations\\EMiles\\CaseStudies_EB\\Output'
path_subcode <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\'
################# 
# Paths for Extra Codes
#library(R.utils)
#sourceDirectory(path_code)
# Paths for Extra Codes
file.sources = list.files(path_subcode, 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\EBModCore.R")
##################
# File Paths/File Names and basic Settings
##################

climData <- 'AWS_ChangriNup_debris.csv'       # Climate Data File

# Open AWS Data Sheet
climInput <- read.csv(path_data&'\\'&climData,header = T)

Sys.setenv(TZ='UTC')
climInput$Time_Str <- as.POSIXct(paste(as.character(climInput$Date),as.character(climInput$Time)), format="%m/%d/%Y  %H:%M:%S")
climInput$Time_Num <- as.numeric(climInput$Time_Str)

diuSWin <- diuCyc(climInput$SWin[4500:5500],climInput$Time_Str[4500:5500])
diuSWout <- diuCyc(climInput$SWout,climInput$Time_Str)
diuLWin <- diuCyc(climInput$LWin,climInput$Time_Str)
diuLWout <- diuCyc(climInput$LWout,climInput$Time_Str)
diuTA <- diuCyc(climInput$Ta,climInput$Time_Str)
diurH <- diuCyc(climInput$rH,climInput$Time_Str)
diuws <- diuCyc(climInput$ws,climInput$Time_Str)

climInput$SWin[which(is.na(climInput$SWin))] <- diuSWin[hour(climInput$Time_Str)[which(is.na(climInput$SWin))]+1,2]
climInput$SWout[which(is.na(climInput$SWout))] <- diuSWout[hour(climInput$Time_Str)[which(is.na(climInput$SWout))]+1,2]
climInput$LWin[which(is.na(climInput$LWin))] <- diuLWin[hour(climInput$Time_Str)[which(is.na(climInput$LWin))]+1,2]
climInput$LWout[which(is.na(climInput$LWout))] <- diuLWout[hour(climInput$Time_Str)[which(is.na(climInput$LWout))]+1,2]
climInput$Ta[which(is.na(climInput$Ta))] <- diuTA[hour(climInput$Time_Str)[which(is.na(climInput$Ta))]+1,2]
climInput$rH[which(is.na(climInput$rH))] <- diurH[hour(climInput$Time_Str)[which(is.na(climInput$rH))]+1,2]
climInput$ws[which(is.na(climInput$ws))] <- diuws[hour(climInput$Time_Str)[which(is.na(climInput$ws))]+1,2]



write.csv(climInput,path_data&'\\NEW'&climData,header = T)
