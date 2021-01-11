################################################################################
# D2EB - distributed debris Energy Balance Model
# 
# dDeb.R
#
# ReadMe:
# Based on Reid and Brock, 2010
# Only edit within the blocks (a) PARAMETERS, (b) VARIABLES and (c) VISUALIZE
# (a)   - Paths: Make sure all paths are correctly pointed to the directory you work in
#       - year (change to specific year, numeric)
#       - StakeName (change to specific Stake name, see Folder StakeData)
#       - AWS_loc (change to 'on' if you use on-glacier data and 'off' if you use off-glacier data)
#       - StationName (change according to the Meteo data file you use, see folder MeteoData)
#
# % NOTE: To run this successfully you will need the files: 
#% 'isleapyear.m', 'debristemp.m', 'Lup.m', 'G.m', 'H.m', 'LE.m', 'P.m', 'loadEBdata.m' and 'diucyc.m'
#% They must all be saved in the same directory as this routine
# 
# Created:          2018/02/05
# Latest Revision:  2017/02/05
#
# Jakob F Steiner | PhD candidate | Faculty of Geosciences | Universiteit Utrecht | 
# Heidelberglaan 2, 3584 CS Utrecht | W.C. van Unnik building | Room 124, Zonneveldvleugel | 
# j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org
#
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

##################
# File Paths/File Names and basic Settings
##################

debThickness <- 'debris_thickness_mean_f1_f3_1000r.tif'
cliff_outlines <- 'Cliffs_May13.shp'
pond_outlines <- 'Ponds_May13.shp'
ModRes <- 30    # Spatial Resolution of Model [m]
#spDom <-'AWSLirung_onglacier_2012_Data.csv'
climData <- 'Kyangjin_ICIMOD.csv'

location<-'Lirung'

# set start and end date of calculation of EB

EBstart <- '2013/05/18 00:00:00'
EBend <- '2013/10/22 00:00:00'

seasonsDyn <- 'on';                     # choose variable seasons (1 dry, 1 wet, 1 dry)

monin <<- '06/15'       # define beginning of wet season
monout <<- '09/19'      # define end of wet season

winin <<- '1/1'        # define end of second dry season
winout <<- '02/28'      # define beginning of first dry season

albedoDyn <<- 'const';                   # choose whether a constant debris albedo value ('const') or a time series is used ('var')

StationName <- 'Lirung_'&toString(year);          # on-glacier Station
StationName_off <- 'Kyanjing_'&toString(year);    # off-glacier Station

# Paths PC
path<-'F:\\PHD\\Research\\EB_DCG\\DistributedEB'
path_code<-'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code'
path_figs<-'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Figures'
path_data<-'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData'


################# 
# Paths for Extra Codes
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\TS_Check.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\diuCyc.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\diu_cycle.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\TSaggregate.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\BulkFlux_Sensible.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\BulkFlux_Latent.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\calc_footprint_FFP_climatology.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\Bulk_Litt.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\Bulk_Reid.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\Bulk_Nicholson.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\EBModCore.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\debristemp.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\Q_H.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\Q_LE.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\Q_G.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\Lup.R")
#source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\DEMderiv.R")

################
projec<-'+proj=utm +zone=45N +datum=WGS84'
#projec<-'+proj=longlat +datum=WGS84'

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
library('raster')

# ========================
#  % PARAMETERS/VARIABLES (may change across sites; EDIT ACCORDING TO FIELDSITE)
# ========================

# CONSTANTS (not likely to change across different sites)

grav <<- 9.81;               # Gravitational acceleration (m s^-2)
k_vk <<- 0.41;                # Von Karman's constant
sigma_sb <<- 5.67e-8;            # Stefan-Boltzmann constant (W m^-2 K^-4)
Rgas <<- 8.31447;             # Gas constant (J mol^-1 K^-1)
Mair <<- 0.0289644;           # Molar mass of dry air (kg mol^-1)
p0 <<- 101325;                # Standard sea level pressure (Pa)
T0 <<- 288.15;                # Standard sea level temperature (K)
L_v <<- 2476000;              # Latent heat of vaporization of water (J kg^-1)
L_f <<- 334000;               # Latent heat of fusion of water (J kg^-1)
rho_w <<- 999.7;              # Density of water (kg m^-3)
rho_i <<- 915;                 # Density of ice (kg m^-3)
c_w <<- 4181.3;               # Specific heat capacity of water (J kg^-1 K^-1)
c_ad <<- 1005;                # Specific heat capacity of air at constant humidity (J kg^-1 K^-1)
T_f <<- 273.15;                    # Freezing point of water (C)

# VARIABLES (to be edited; EDIT ACCORDING TO FIELDSITE)
timestep <<- 3600;           # Model timestep (s)
altitude <<- 3980;           # Altitude of measurement site (m)
z_a <<- 2.31;                # % Height of air temp / wind / humidity measurements (m)
z_0 <<- 0.03;               # % Debris aerodynamic roughness length (m)
Lapse <<- 0.0065;            # % Temperature lapse rate (K m^-1)
k_d_wet <<- 1.55;             # Debris thermal conductivity wet season (W m^-1 K^-1)
k_d_wet_range <- c(1.40,1.70);# possible range for k_d_wet (based on Nicholson & Benn 2012)
k_d_dry <<- 1.04;                 # Debris thermal conductivity dry season (W m^-1 K^-1)
k_d_dry_range <- c(0.94,1.14);    # possible range for k_d (based on Nicholson & Benn 2012)
rho_d <<- 1496;               # Debris density (kg m^-3)
rho_d_range <- c(1000,1900);  # range for rho_d; 1496 (Reid & Brock 2010), 2700 (Nicholson & Benn 2012)
c_d <<- 948;                  # Debris specific heat capacity (J kg^-1 K^-1)
c_d_range <- c(850,1045);     # range for c_d
epsilon_d <<- 0.94;           # Debris emissivity
epsilon_d_range <- c(0.92,0.97); # range for epsilon_d
albedo_const <- 0.15;        # Debris albedo (will be calculated if dynamic version is chosen)

# Calculate air pressure in Pa based on altitude
p_a <- p0*((1-(Lapse*altitude/T0))^(grav*Mair/(Rgas*Lapse)));


# ========================
#  Load DATA for each Grid Cell
# ========================

# Read Model Domain / Number of Cells

# Debris Thickness (FOR THE TIME BEING DEFINES THE MODEL DOMAIN)
Hdeb <- raster(path_data&'\\DebrisData\\'&debThickness)

# resample to standard model resolution
Hdeb_ModRes <- aggregate(Hdeb, fact = floor(ModRes/res(Hdeb)) , median)
Hdeb_dataframe <- as.data.frame(Hdeb_ModRes)

# Number of Model cells
DomDim <- dim(Hdeb_ModRes)
CellNo <- DomDim[1] * DomDim[2]

# Read in Surface Features with different melt properties
ogrInfo(path_data&'\\Outlines_Features\\'&cliff_outlines)
cliff_mask<-readOGR(dsn=path_data&'\\Outlines_Features\\'&cliff_outlines)
projection(cliff_mask)<-projec
r_cliffmask <- rasterize(cliff_mask, Hdeb_ModRes)

ogrInfo(path_data&'\\Outlines_Features\\'&pond_outlines)
pond_mask<-readOGR(dsn=path_data&'\\Outlines_Features\\'&pond_outlines)
projection(pond_mask)<-projec
r_pondmask <- rasterize(pond_mask, Hdeb_ModRes)

Hdeb_masks <- Hdeb_ModRes
Hdeb_masks[!is.na(r_pondmask)] <- 10
Hdeb_masks[!is.na(r_cliffmask)] <- 20

png(file=path_figs&'\\Thickness_Shift.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
plot(Hdeb,ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list(text=expression('Debris Thickness [m]'), side=2, font=1, line=0.3, cex=1))
plot(pond_mask,add=T,col = 'blue')
plot(cliff_mask, add = T)
grid(nx=NULL, ny=NULL)
dev.off()

Hdeb2 <- aggregate(Hdeb,fact = 2,median)
# Read DEM from same year as thickness map
DEM_201605<-raster(path_data&'/DEMs/'&'20160430_lirung_dem_20cm.tif')
DEM_201605_dom <- aggregate(crop(DEM_201605,extent(Hdeb2)),fact = 15,median)
DEM_201605_dom <- mask(DEM_201605_dom,Hdeb2)

plot(terrain(DEM_201605_dom, opt=c('slope'), unit='degrees'),Hdeb2)

# Read ASTER GDEM for surrounding area
DEM_ASTER<-raster(path_data&'/DEMs/'&'LangtangASTERGDEM.tif')

#####FROM HERE #####

# Open AWS Data Sheet
climInput <- read.csv(path_data&'\\'&climData,header = T)

Sys.setenv(TZ='Asia/Kathmandu')
climInput$Time_Str <- as.POSIXct(paste(as.character(climInput$DATE),as.character(climInput$TIME)), format="%m/%d/%Y  %H:%M:%S")
climInput$Time_Num <- as.numeric(climInput$Time_Str)

# find IDs to limit the time series of calculation
modstart <- as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S")
modend <- as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S")
IDstart <- which(!is.na(match(climInput$Time_Str,modstart)))
IDend <- which(!is.na(match(climInput$Time_Str,modend)))


######### 
# TEST SCENARIO DEBUGGING in SINGLE RUN MODE
# #######
# Calculate d2EB at each individual grid cell with available time series
EBCoreInput <- list()

# Gap Filling
GAP_IDs <- which(is.na(climInput$KINC))

gapFill <- function(TS_withHoles){    # Gap Filling Function
  for (k in 25:length(TS_withHoles)){
    if(is.na(TS_withHoles[k])){TS_withHoles[k] <- TS_withHoles[k-24]}
  }
  for (k in 1:24){
    if(is.na(TS_withHoles[k])){TS_withHoles[k] <- TS_withHoles[k+24]}
  }
  return(TS_withHoles)
}

EBCoreInput$timeline_str <- climInput$Time_Str[IDstart:IDend]
timeline_num <- climInput$Time_Num[IDstart:IDend]
EBCoreInput$SWin <- gapFill(climInput$KINC[IDstart:IDend])            # Downwelling shortwave radiation (W m^-2)
EBCoreInput$SWout <- gapFill(climInput$KUPW[IDstart:IDend])             # Upwelling shortwave radiation (W m^-2)
EBCoreInput$LWin <- gapFill(climInput$LINC[IDstart:IDend])              # Downwelling longwave radiation (W m^-2)
EBCoreInput$LWout <- gapFill(climInput$LUPW[IDstart:IDend])             # Upwelling longwave radiation (W m^-2)
EBCoreInput$T_s_data <- (EBCoreInput$LWout/5.67/10^(-8)/1)^(1/4);	    # Surface temperature (K)
EBCoreInput$T_a <- gapFill(climInput$TAIR[IDstart:IDend] + 273.15);      	    # Air temperature (K)
EBCoreInput$u <- gapFill(climInput$WSPD[IDstart:IDend])                     # Wind speed (m s^-1)
EBCoreInput$wd <- gapFill(climInput$WDIR[IDstart:IDend])
EBCoreInput$RH_a <- gapFill(climInput$RH[IDstart:IDend])          # Air relative humidity (%)
EBCoreInput$p_a <- gapFill(climInput$PRES[IDstart:IDend])          # Air pressure (kPa)  

# Albedo Model (Constant vs. Variable)
switch(albedoDyn,
       const = {                       # constant Debris albedo
         albedo_d = albedo_const + EBCoreInput$SWin*0 ;},            
       var = {                         # variable debris albedo
         albedo_d = EBCoreInput$SWout / EBCoreInput$SWin
         # remove non-sensical values
         albedo_d[albedo_d<0 | albedo_d>1 | is.infinite(albedo_d) | albedo_d==0] <- NA;
         # remove albedo at night time (few values, not reasonable, scatter effects at low solar angle)
         albedo_d[hour(EBCoreInput$timeline_str) < 9 | hour(EBCoreInput$timeline_str) > 18] = NaN;})
EBCoreInput$albedo_d <- albedo_d


# Debris Layers
d <- Hdeb_dataframe[91,1]

range <- 0.5;                   # Range for central difference calculation of the derivative of surface flux when calculating surface temperature T_s.
N <<- floor(d/0.10);             # Number of layers used for debris temperature calculation (chosen to make each layer 1cm deep)
if(N<10){N <- 10}               # Set minimum number of layers = 10
#N <- 4
h <<- d/N;                       # Size of each layer in calculation


# View Factors for Radiation

DEM_dataframe <- cbind(seq(1,length(DEM_res[]),1),DEM_res[],coordinates(DEM_res))
#DEM_dataframe_new <- DEM_dataframe[which(!is.na(DEM_dataframe[,2])),]
old <- Sys.time()
#testFunc <- function(i){
vf <- viewFactors_core(DEM_dataframe,DEM_r,res(DEM_res)[1],res(DEM_res)[1],NULL,NULL)
EBModRes <-  EBModCore(EBCoreInput);


png(file=path_figs&'\\ModelFit _ PointScale.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
#par(mfrow=c(3,1))
axDate <- round(range(EBCoreInput$timeline_str,na.rm=T), "days")
layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
plot(EBCoreInput$timeline_str,EBModRes[,1]*1000,
     xlab ='Date', ylab = expression('melt'['mod']~~'[mm]'),
     ylim=range(c(0,1)),
     type='l',col='black', xaxt = 'n',frame.plot=TRUE,xlim = c(as.numeric(axDate[1]), as.numeric(axDate[2])))
axis.POSIXct(1, at = seq(axDate[1], axDate[2], by = "month"), format = "%m/%Y")
abline(h=seq(0,1,0.2),v=seq(axDate[1], axDate[2], by = "month"), col="gray", lty=3)
plot(EBCoreInput$timeline_str,EBModRes[,2],
     xlab ='Date', ylab = expression('melt'['cum']~~'[m]'),
     ylim=range(c(0,2)),
     type='l',col='black', xaxt = 'n',frame.plot=TRUE,xlim = c(as.numeric(axDate[1]), as.numeric(axDate[2])))
axis.POSIXct(1, at = seq(axDate[1], axDate[2], by = "month"), format = "%m/%Y")
abline(h=seq(0,2,0.2),v=seq(axDate[1], axDate[2], by = "month"), col="gray", lty=3)
plot(EBCoreInput$timeline_str,EBCoreInput$T_a-273.15,
     xlab ='Date', ylab = expression('T'['air']~~'[°C]'),
     ylim=range(c(-5,20)),
     type='l',col='black', xaxt = 'n',frame.plot=TRUE,xlim = c(as.numeric(axDate[1]), as.numeric(axDate[2])))
axis.POSIXct(1, at = seq(axDate[1], axDate[2], by = "month"), format = "%m/%Y")
abline(h=seq(-5,20,5),v=seq(axDate[1], axDate[2], by = "month"), col="gray", lty=3)
plot(EBCoreInput$timeline_str,EBModRes[,4],
     xlab ='Date', ylab = expression('Debris Conductivity [W '~ m^{-1}~K^{-1} ~']'),
     ylim=range(c(0,2)),
     type='l',col='black', xaxt = 'n',frame.plot=TRUE,xlim = c(as.numeric(axDate[1]), as.numeric(axDate[2])))
axis.POSIXct(1, at = seq(axDate[1], axDate[2], by = "month"), format = "%m/%Y")
abline(h=seq(0,2,0.2),v=seq(axDate[1], axDate[2], by = "month"), col="gray", lty=3)
plot(EBCoreInput$timeline_str,EBCoreInput$SWin - EBCoreInput$SWout + EBCoreInput$LWin - EBCoreInput$LWout,
     xlab ='Date', ylab = expression('Net Flux [W '~ m^{-2} ~']'),
     ylim=range(c(-1000,1000)),
     type='l',col='black', xaxt = 'n',frame.plot=TRUE,xlim = c(as.numeric(axDate[1]), as.numeric(axDate[2])))
points(EBCoreInput$timeline_str,EBModRes[,10] + EBModRes[,11], type = 'l', col = 'red')
axis.POSIXct(1, at = seq(axDate[1], axDate[2], by = "month"), format = "%m/%Y")
abline(h=seq(-1000,1000,100),v=seq(axDate[1], axDate[2], by = "month"), col="gray", lty=3)
plot(EBModRes[,3] - 273.15,EBCoreInput$T_s_data - 273.15,xlim=c(-10,40),ylim=c(-10,40),ylab=expression('T'['surf-meas']~~'[°C]'),xlab=expression('T'['surf-mod']~~'[°C]'))
grid(nx=NULL, ny=NULL)
abline(0,1)
dev.off()

###### Run function parallel multicore #####

# define simple test function
testFunc <- function(i){
  EBCoreInput <- list()
  
  # Gap Filling
  GAP_IDs <- which(is.na(EBCoreInput$SWin))
  
  gapFill <- function(TS_withHoles){    # Gap Filling Function
    for (k in 25:length(TS_withHoles)){
      if(is.na(TS_withHoles[k])){TS_withHoles[k] <- TS_withHoles[k-24]}
    }
    for (k in 1:24){
      if(is.na(TS_withHoles[k])){TS_withHoles[k] <- TS_withHoles[k+24]}
    }
    return(TS_withHoles)
  }
  
  EBCoreInput$timeline_str <- climInput$Time_Str[IDstart:IDend]
  timeline_num <- climInput$Time_Num[IDstart:IDend]
  EBCoreInput$SWin <- gapFill(climInput$KINC[IDstart:IDend])            # Downwelling shortwave radiation (W m^-2)
  EBCoreInput$SWout <- gapFill(climInput$KUPW[IDstart:IDend])             # Upwelling shortwave radiation (W m^-2)
  EBCoreInput$LWin <- gapFill(climInput$LINC[IDstart:IDend])              # Downwelling longwave radiation (W m^-2)
  EBCoreInput$LWout <- gapFill(climInput$LUPW[IDstart:IDend])             # Upwelling longwave radiation (W m^-2)
  EBCoreInput$T_s_data <- gapFill(climInput$TSOIL[IDstart:IDend] + 273.15);	    # Surface temperature (K)
  EBCoreInput$T_a <- gapFill(climInput$TAIR[IDstart:IDend] + 273.15);      	    # Air temperature (K)
  EBCoreInput$u <- gapFill(climInput$WSPD[IDstart:IDend])                     # Wind speed (m s^-1)
  EBCoreInput$wd <- gapFill(climInput$WDIR[IDstart:IDend])
  EBCoreInput$RH_a <- gapFill(climInput$RH[IDstart:IDend])          # Air relative humidity (%)
  EBCoreInput$p_a <- gapFill(climInput$PRES[IDstart:IDend])          # Air pressure (kPa)  
  
  # Albedo Model (Constant vs. Variable)
  switch(albedoDyn,
         const = {                       # constant Debris albedo
           albedo_d = albedo_const + EBCoreInput$SWin*0 ;},            
         var = {                         # variable debris albedo
           albedo_d = EBCoreInput$SWout / EBCoreInput$SWin
           # remove non-sensical values
           albedo_d[albedo_d<0 | albedo_d>1 | is.infinite(albedo_d) | albedo_d==0] <- NA;
           # remove albedo at night time (few values, not reasonable, scatter effects at low solar angle)
           albedo_d[hour(EBCoreInput$timeline_str) < 9 | hour(EBCoreInput$timeline_str) > 18] = NaN;})
  EBCoreInput$albedo_d <- albedo_d
  
  # ========
  # MODEL CORE
  # ========
  #browser()
  # Debris Layers
  d <- i
  
  range <- 0.5;                   # Range for central difference calculation of the derivative of surface flux when calculating surface temperature T_s.
  N <<- floor(d/0.01);             # Number of layers used for debris temperature calculation (chosen to make each layer 1cm deep)
  if(N<10){N <- 10}               # Set minimum number of layers = 10
  h <<- d/N;                       # Size of each layer in calculation
  
  EBModRes <-  EBModCore(EBCoreInput);
  return(EBModRes[,2])

}

# number of workers (leave one to keep system more responsive)
nworkers <- detectCores()-3

# run entire function parallelized
sfInit(parallel=T,cpu=nworkers, type="SOCK")                                                # initialize cluster
loadlib <- lapply(c('raster','truncnorm'),                                                  # load required packages in cluster
                  function(x) sfLibrary(x, character.only=T))
sfExportAll()                                                                               # transfer all variables to clusters
sfClusterSetupRNG()                                                                         # set up random number generator

tout.multi <- system.time(
  out <- sfLapply(Hdeb_dataframe[!is.na(Hdeb_dataframe[,1]),1], testFunc)
)
sfStop() 

save(out, file=path_data&'//temp_ModOutput')
load(path_data&'//temp_ModOutput')
###########
#Test Model Plots
allPix <- lapply(out, max);         # Total melt after complete model period
maxMelt <- Hdeb_dataframe[,1]*0
maxMelt[!is.na(Hdeb_dataframe[,1])] <- allPix

Hdeb_ModRes2 <- Hdeb_ModRes
Hdeb_ModRes2[1:length(Hdeb_dataframe[,1])] <- unlist(maxMelt)   # make raster with melt values

# Load DEMs for Validation Data
DEM_201305<-raster(path_data&'/DEMs/'&'201305_V4_DEM_20cm.tif')
DEM_201305_dom <- aggregate(crop(DEM_201305,extent(Hdeb_ModRes)),fact = floor(res(Hdeb_ModRes)/res(DEM_201305)),median)
DEM_201305_dom <- mask(DEM_201305_dom,Hdeb_ModRes)

DEM_201310<-raster(path_data&'/DEMs/'&'201310_V4_DEM_20cm.tif')
DEM_201310_dom <- aggregate(crop(DEM_201310,extent(Hdeb_ModRes)),fact = (res(Hdeb_ModRes)/res(DEM_201310)),median)
DEM_201310_dom <- mask(DEM_201310_dom,Hdeb_ModRes)

DEM_201405<-raster(path_data&'/DEMs/'&'201405_V6_Lirung_DEM_20cm_BilinearSnap.tif')
DEM_201405_dom <- aggregate(crop(DEM_201405,extent(Hdeb_ModRes)),fact = (res(Hdeb_ModRes)/res(DEM_201405)),median)
DEM_201405_dom <- resample(DEM_201405_dom,Hdeb_ModRes,'bilinear')
DEM_201405_dom <- mask(DEM_201405_dom,Hdeb_ModRes)

DEMdates <- c('05/18/2013','10/22/2013')
DEMini <- DEM_201305_dom
DEMend <- DEM_201310_dom
dh <- DEMend - DEMini
dt <- (as.numeric(as.POSIXct(DEMdates[2], format="%m/%d/%Y")) - as.numeric(as.POSIXct(DEMdates[1], format="%m/%d/%Y")))/3600/24
#DEM_201410<-raster(path_data&'/DEMs/'&'General_DEM_model_recalculated_after_merging_chunks.tif')
#DEM_201410_dom <- aggregate(crop(DEM_201410,extent(Hdeb_ModRes)),fact = floor(res(Hdeb_ModRes)/res(DEM_201410)),median)
#DEM_201410_dom <- mask(DEM_201410_dom,Hdeb_ModRes)

#b <- layerize(DEM_201410)
#fact <- round(dim(r_hr)[1:2] / dim(r_lr)[1:2])
#a <- aggregate(b, fact)
#x <- resample(a, r_lr)

d_oestrem <- c(0,5,10,15,20,25,30)
m_oestrem <- c(4.5,3,2,1,0.8,0.6,0.5)
approx_ostrem <- cbind(d_oestrem,m_oestrem)
png(file=path_figs&'\\ModelPerformance.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
#par(mfrow=c(3,1))
layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
plot(DEMini,ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list(text=expression('Elevation [m]'), side=2, font=1, line=0.3, cex=1))
grid(nx=NULL, ny=NULL)
plot(Hdeb_ModRes,ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list(text=expression('debris thickness [m]'), side=2, font=1, line=0.3, cex=1))
grid(nx=NULL, ny=NULL)
plot(Hdeb_ModRes2/length(EBCoreInput$SWin)*24*1000,zlim=c(0,20),ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list(text=expression('melt mod. [mm/day]'), side=2, font=1, line=0.3, cex=1.1))
grid(nx=NULL, ny=NULL)
plot(Hdeb_ModRes*100,Hdeb_ModRes2/length(EBCoreInput$SWin)*24*1000,xlim=c(0,100),ylim=c(0,15),ylab='melt [mm/day]',xlab='deb thickness [cm]')
grid(nx=NULL, ny=NULL)
points(approx_ostrem[,1],approx_ostrem[,2],col='red')
legend('topleft',c('mod','Oestrem'),pch=c(1,1),col=c('black','red'))
plot(dh/dt*-1000,zlim=c(0,20),ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list(text=expression('dh/dt [mm/day]'), side=2, font=1, line=0.3, cex=1.1))
grid(nx=NULL, ny=NULL)
plot(dh/dt*-1000,Hdeb_ModRes2/length(EBCoreInput$SWin)*24*1000,xlim=c(0,20),ylim=c(0,20),ylab='model [mm/day]',xlab='observed [mm/day]')
grid(nx=NULL, ny=NULL)
abline(0,1)
dev.off()



