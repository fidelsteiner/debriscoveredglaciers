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

climData <- 'AWSLirung_onglacier_2014_Data.csv'       # Climate Data File
paramData <- 'Params_5000.csv'              # Parameter Data File

##### only for 
glacCharacFile <- 'GlacierCharac.csv'
abbrev <- 'LIR14'
glacCharac <- read.csv(path_data&'\\'&glacCharacFile,header = T)
idCharac <- which(glacCharac$GLA==abbrev)

location<-abbrev                        # Location Name for outputs
StationName <- abbrev;                  # on-glacier Station
yearTS <- 2014                             # Year of measurements (if more than one give first year, will be disabled)
lat <- glacCharac$LAT[idCharac]
lon <- glacCharac$LON[idCharac]

# Define whether station is located on a clean ice ('clean'), debris-covered ('debris') glacier or off-glacier ('off')
station_loc <- 'debris'

# specify height of T_a and ws sensor (assumed height is 2 m)
sh_Ta <- glacCharac$H_T[idCharac];
sh_ws <- glacCharac$H_W[idCharac];

# set start and end date of calculation of EB
EBstart <- glacCharac$SIMSTART[idCharac]
EBend <- glacCharac$SIMEND[idCharac]

# for the case of wet/dry seasons turn the following on and specifiy end/start dates
seasonsDyn <- 'off';                        # choose variable seasons (1 dry, 1 wet, 1 dry)
monin <<- '06/15'                           # define beginning of wet season
monout <<- '09/19'                          # define end of wet season

winin <<- '1/1'                             # define end of second dry season
winout <<- '02/28'                          # define beginning of first dry season

albedoDyn <<- 'const';                      # choose whether a constant albedo value ('const') or a time series is used ('var')

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
altitude <<- glacCharac$ELE[idCharac];           # Altitude of measurement site (m)
z_a <<- 2.00;                # % Height of air temp / wind / humidity measurements (m)
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

# Variable parameter space

paramInput <- read.csv(path_data&'\\'&paramData,header = T)

# ========================
#  Load DATA 
# ========================

# Open AWS Data Sheet
climInput <- read.csv(path_data&'\\'&climData,header = T)

Sys.setenv(TZ='UTC')
climInput$Time_Str <- as.POSIXct(paste(as.character(climInput$Date),as.character(climInput$Time)), format="%m/%d/%Y  %H:%M:%S")
climInput$Time_Num <- as.numeric(climInput$Time_Str)

# find IDs to limit the time series of calculation
modstart <- as.POSIXct(EBstart, format="%m/%d/%Y  %H:%M")
modend <- as.POSIXct(EBend, format="%m/%d/%Y  %H:%M")
IDstart <- which(!is.na(match(climInput$Time_Str,modstart)))
IDend <- which(!is.na(match(climInput$Time_Str,modend)))


######### 
# Read all data into Structure
# #######

EBCoreInput <- list()

# scale wind speeds to different sensor height
climInput$ws_corr <- climInput$ws*(log(sh_ws)/log(2))

climInput$SWin[climInput$SWin>1361] <- 1361
EBCoreInput$SWin <- TSaggregate(climInput$SWin,climInput$Time_Str,60,yearTS,'mean')[,2]            # Downwelling shortwave radiation (W m^-2)

EBCoreInput$SWinhour(climInput$Time_Str)

EBCoreInput$SWin[EBCoreInput$SWin<10] <- 0
EBCoreInput$SWout <- TSaggregate(climInput$SWout,climInput$Time_Str,60,yearTS,'mean')[,2]             # Upwelling shortwave radiation (W m^-2)
EBCoreInput$SWout[EBCoreInput$SWout>EBCoreInput$SWin] <- EBCoreInput$SWin[EBCoreInput$SWout>EBCoreInput$SWin]

EBCoreInput$LWin <-TSaggregate(climInput$LWin,climInput$Time_Str,60,yearTS,'mean')[,2]              # Downwelling longwave radiation (W m^-2)
EBCoreInput$LWout <- TSaggregate(climInput$LWout,climInput$Time_Str,60,yearTS,'mean')[,2]             # Upwelling longwave radiation (W m^-2)
EBCoreInput$T_s_data <- (EBCoreInput$LWout/5.67/10^(-8)/1)^(1/4);	    # Surface temperature (K)
EBCoreInput$T_a <- TSaggregate(climInput$Ta,climInput$Time_Str,60,yearTS,'mean')[,2] + 273.15;      	    # Air temperature (K)
EBCoreInput$u <- TSaggregate(climInput$ws_corr,climInput$Time_Str,60,yearTS,'mean')[,2]                     # Wind speed (m s^-1)
EBCoreInput$wd <- TSaggregate(climInput$wd,climInput$Time_Str,60,yearTS,'mean')[,2]
EBCoreInput$RH_a <- TSaggregate(climInput$rH,climInput$Time_Str,60,yearTS,'mean')[,2]          # Air relative humidity (%)
EBCoreInput$RH_a[EBCoreInput$RH_a>100] <- 100
EBCoreInput$p_a <- p_a +  EBCoreInput$RH_a* 0         # Air pressure (kPa)  
if(!is.null(climInput$pa)){
EBCoreInput$p_a <- TSaggregate(climInput$pa,climInput$Time_Str,60,yearTS,'mean')[,2] * 100
}
EBCoreInput$mod <- station_loc;                             # type of model (debris/clean ice/pond/cliff)
EBCoreInput$lat <- lat;
EBCoreInput$lon <- lon;

EBCoreInput$T_a[1] <- EBCoreInput$T_a[2]
EBCoreInput$LWout[1] <- EBCoreInput$LWout[2]
EBCoreInput$LWin[1] <- EBCoreInput$LWin[2]
EBCoreInput$RH_a[1] <- EBCoreInput$RH_a[2]
EBCoreInput$SWout[1] <- EBCoreInput$SWout[2]
EBCoreInput$u[1] <- EBCoreInput$u[2]
EBCoreInput$timeline_str <- as.POSIXct(TSaggregate(climInput$SWin,climInput$Time_Str,60,yearTS,'mean')[,1],origin = "1970-01-01")

EBCoreInput$mod <- 'debris'
EBModType <<- 'no_TS'
CDrange <<- 0.5;                   # Range for central difference calculation of the derivative of surface flux when calculating surface temperature T_s.


######### 
# Single core Model Run for Debugging
########
for(i in 1:10){
z_0 <<- paramInput$z_0_d..m.[i];               # % Debris aerodynamic roughness length (m)
k_d_wet <<- paramInput$k_d..W.m...1..K...1..[i];             # Debris thermal conductivity wet season (W m^-1 K^-1)
k_d_dry <<- k_d_wet;                 # Debris thermal conductivity dry season (W m^-1 K^-1)
rho_d <<- paramInput$rho_r..kg.m...3..[i];               # Debris density (kg m^-3)
c_d <<- paramInput$c_d..J.kg...1..K...1..[i];                  # Debris specific heat capacity (J kg^-1 K^-1)
epsilon_d <<- paramInput$eps_d....[i];           # Debris emissivity
EBCoreInput$albedo_d <- paramInput$alpha_d....[i];        # Debris albedo (will be calculated if dynamic version is chosen)

EBCoreInput$T_s_data <- (EBCoreInput$LWout/5.67/10^(-8)/epsilon_d)^(1/4)

# Debris Layers
d <- paramInput$d..m.[i]
range <- 0.5;                   # Range for central difference calculation of the derivative of surface flux when calculating surface temperature T_s.
N <<- floor(d/0.01);             # Number of layers used for debris temperature calculation (chosen to make each layer 1cm deep)
if(N<5){N <<- 5}               # Set minimum number of layers = 10
#N <- 10
h <<- d/N;                       # Size of each layer in calculation

EBModRes <-  EBModCore(EBCoreInput);
}

plot(EBModRes[1:944,3]-273.15,type='l',ylim=c(-10,40),ylab ='T_surface')
points(EBCoreInput$T_s_data[1:944]-273.15,type='l',col='red')
legend('topleft',legend=c('mod','meas'),lty=1,lwd=1,col=c('black','red'),bty='n')


plot(diuCyc(EBModRes[1:944,3]-273.15,EBCoreInput$timeline_str[1:944])[,2],type='l',ylim=c(-5,30),ylab ='T_surface')
points(diuCyc(EBCoreInput$T_s_data[1:944]-273.15,EBCoreInput$timeline_str[1:944])[,2],type='l',col='red')
legend('topleft',legend=c('mod','meas'),lty=1,lwd=1,col=c('black','red'),bty='n')

# Flux at each layer:

lowerT_s <- dim(EBModRes)[2]-2
Q <- EBModRes[,12:lowerT_s] * 0
Q <- Q[,-1]

for(pk in 1:dim(EBModRes)[1]){
for(k in 13:(lowerT_s)){
  Q[pk,k-12] <- (EBModRes[pk,(k-1)]-EBModRes[pk,k])/h*k_d_wet}
}

Qp <- Q*0
for(pk in 2:dim(EBModRes)[1]){
  for(k in 2:(lowerT_s-13)){
Qp[pk,k] <- Q[pk-1,k]+(Q[pk-1,k-1]-Q[pk-1,k])-(Q[pk-1,k]-Q[pk-1,k+1]) - Q[pk,k]

  }
}

png(file=path_figs&'\\tempFit.png', res = 160,width=1800,height=1800)
par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
plot(EBCoreInput$T_s_data-273.15,EBModRes[,3]-273.15,xlim = c(-10,35),ylim = c(-10,35),xlab='measured T_surface [C]',ylab='modelled T_surface [C]')
abline(0,1)
grid(NULL,NULL)
dev.off()

dT_s <- (EBCoreInput$T_s_data-273.15)-(EBModRes[,3]-273.15)

png(file=path_figs&'\\Q_turb.png', res = 160,width=1800,height=1800)
par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
plot(EBModRes[,'Q_H_v']+EBModRes[,'Q_LE_v'],ylab='turbulent fluxes [W m-2]',type='l')

grid(NULL,NULL)
dev.off()

png(file=path_figs&'\\dT_ice.png', res = 160,width=1800,height=1800)
par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
plot(EBModRes[,42]-EBModRes[,41],ylab='dT at ice interface',type='l')

grid(NULL,NULL)
dev.off()


png(file=path_figs&'\\dT_ice.png', res = 160,width=1800,height=1800)
par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
plot(EBModRes[2020,12:(dim(EBModRes)[2]-2)]-273.15,seq(d,0,-h),ylab='debris depth [m]',xlab='debris T [C]',type='l')
for(i in 1:240){
points(EBModRes[2020+i,12:(dim(EBModRes)[2]-2)]-273.15,seq(d,0,-h),add=T,type='l')
}
grid(NULL,NULL)
dev.off()


######### 
# Multicore Model Runs
########

MultiCore_EB <- function(i){

  z_0 <<- paramInput$z_0_d..m.[i];               # % Debris aerodynamic roughness length (m)
  k_d_wet <<- paramInput$k_d..W.m...1..K...1..[i];             # Debris thermal conductivity wet season (W m^-1 K^-1)
  k_d_dry <<- k_d_wet;                 # Debris thermal conductivity dry season (W m^-1 K^-1)
  rho_d <<- paramInput$rho_r..kg.m...3..[i];               # Debris density (kg m^-3)
  c_d <<- paramInput$c_d..J.kg...1..K...1..[i];                  # Debris specific heat capacity (J kg^-1 K^-1)
  epsilon_d <<- paramInput$eps_d....[i];           # Debris emissivity
  EBCoreInput$albedo_d <- paramInput$alpha_d....[i];        # Debris albedo (will be calculated if dynamic version is chosen)
  
  # Debris Layers
  d <- paramInput$d..m.[i]
  range <- 0.5;                   # Range for central difference calculation of the derivative of surface flux when calculating surface temperature T_s.
  N <<- floor(d/0.01);             # Number of layers used for debris temperature calculation (chosen to make each layer 1cm deep)
  if(N<5){N <<- 5}               # Set minimum number of layers = 10
  #N <- 4
  h <<- d/N;                       # Size of each layer in calculation
  
  EBModRes <-  EBModCore(EBCoreInput);
    print(paste("Actual Progress:", floor(i/10*100),"%"))

  write((EBModRes[,3]),path_temp&'\\T_s_run_'&i&'.txt',ncolumns=1)    # temporary save file for T_surf
  write((EBModRes[,1]),path_temp&'\\melt_run_'&i&'.txt',ncolumns=1)   # temporary save file for melt
  #write((EBModRes[,dim(EBModRes)[2]]),path_temp&'\\refreeze_run_'&i&'.txt',ncolumns=1)   # temporary save file for refreeze
  #write((EBModRes[,7]),path_temp&'\\Snet_run_'&i&'.txt',ncolumns=1)   # temporary save file for net solar
  #write((EBModRes[,8]),path_temp&'\\Lnet_run_'&i&'.txt',ncolumns=1)   # temporary save file for net longwave
  write((EBModRes[,10]+EBModRes[,11]),path_temp&'\\Turb_run_'&i&'.txt',ncolumns=1)   # temporary save file for turb
  #write((EBModRes[,10]),path_temp&'\\QH_run_'&i&'.txt',ncolumns=1)   # temporary save file for sensible
  #write((EBModRes[,11]),path_temp&'\\QLE_run_'&i&'.txt',ncolumns=1)   # temporary save file for latent
}

# number of workers (leave one to keep system more responsive)
nworkers <- detectCores()-1

# run entire function parallelized
sfInit(parallel=T,cpu=nworkers, type="SOCK")                                                # initialize cluster
loadlib <- lapply(c('raster','truncnorm'),                                                  # load required packages in cluster
                  function(x) sfLibrary(x, character.only=T))
sfExportAll()                                                                               # transfer all variables to clusters
sfClusterSetupRNG()                                                                         # set up random number generator

tout.multi <- system.time(
  out <- sfLapply(paramInput$RunNo....[1:5000], MultiCore_EB)
)
sfStop() 

EB_Ts <- matrix(, nrow = length(EBCoreInput$timeline_str), ncol = 5001,byrow=T)
EB_melt <- matrix(, nrow = length(EBCoreInput$timeline_str), ncol = 5001,byrow=T) 
#EB_Snet <- matrix(, nrow = length(EBCoreInput$timeline_str), ncol = 1001,byrow=T) 
#EB_Lnet <- matrix(, nrow = length(EBCoreInput$timeline_str), ncol = 1001,byrow=T) 
EB_Turb <- matrix(, nrow = length(EBCoreInput$timeline_str), ncol = 5001,byrow=T) 
#EB_refreeze <- matrix(, nrow = length(EBCoreInput$timeline_str), ncol = 1001,byrow=T) 
for(i in 2:5001){
  EB_Ts[,i] <- read.table(path_temp&'\\T_s_run_'&i-1&'.txt')[,1]
  EB_melt[,i] <- read.table(path_temp&'\\melt_run_'&i-1&'.txt')[,1]
  #EB_refreeze[,i] <- read.table(path_temp&'\\refreeze_run_'&i-1&'.txt')[,1]
  #EB_Snet[,i] <- read.table(path_temp&'\\Snet_run_'&i-1&'.txt')[,1]
  #EB_Lnet[,i] <- read.table(path_temp&'\\Lnet_run_'&i-1&'.txt')[,1]
  EB_Turb[,i] <- read.table(path_temp&'\\Turb_run_'&i-1&'.txt')[,1]
}


  EB_Ts[,1] <- as.character(EBCoreInput$timeline_str)
  EB_melt[,1] <- as.character(EBCoreInput$timeline_str)
  #EB_refreeze[,1] <- as.character(EBCoreInput$timeline_str)
  #EB_Snet[,1] <- as.character(EBCoreInput$timeline_str)
  #EB_Lnet[,1] <- as.character(EBCoreInput$timeline_str)
  EB_Turb[,1] <- as.character(EBCoreInput$timeline_str)
  
  colnames(EB_Ts) <- c('Date',1:5000)
  colnames(EB_melt)<- c('Date',1:5000)
  #colnames(EB_refreeze) <- c('Date',1:5000)
  #colnames(EB_Snet) <- c('Date',1:5000)
  #colnames(EB_Lnet) <- c('Date',1:5000)
  colnames(EB_Turb) <- c('Date',1:5000)
  
###### Plotting and Outputs #####

write.csv(EB_Ts, file = path_output&'\\Tsurf.csv', row.names=FALSE)
write.csv(EB_melt, file = path_output&'\\Melt.csv', row.names=FALSE)
#write.csv(EB_Snet, file = path_output&'\\Snet.csv', row.names=FALSE)
#write.csv(EB_Lnet, file = path_output&'\\Lnet.csv', row.names=FALSE)
write.csv(EB_Turb, file = path_output&'\\Turb.csv', row.names=FALSE)
#write.csv(EB_refreeze, file = path_output&'\\refreeze.csv', row.names=FALSE)



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
     xlab ='Date', ylab = expression('T'['air']~~'[�C]'),
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
plot(EBModRes[,3] - 273.15,EBCoreInput$T_s_data - 273.15,xlim=c(-10,40),ylim=c(-10,40),ylab=expression('T'['surf-meas']~~'[�C]'),xlab=expression('T'['surf-mod']~~'[�C]'))
grid(nx=NULL, ny=NULL)
abline(0,1)
dev.off()