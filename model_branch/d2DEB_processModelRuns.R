################################################################################
# D2EB - Visualize Model runs from the d2EB
# 
# d2DEB_processModelRuns.R
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
#install.packages('pacman')
library(pacman)
p_load(tidyverse,rgdal,rgeos,maptools,raster,foreach,lubridate,compare,colorRamps,data.table,circular,parallel,snowfall,truncnorm,rlecuyer,forecast,rasterVis,R.utils,zoo,rlist)

projec<-'+proj=utm +zone=45N +datum=WGS84'
#projec<-'+proj=longlat +datum=WGS84'

##################
# File Paths/File Names and basic Settings
##################

# Paths PC1
path <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB'        
path_code <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code'
path_subcode <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\'
path_figs <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Figures'
path_data <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData'
path_dhdt <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt'
path_shade <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\localShading'

glac_outline <- 'Lirung_2015.shp'                   # Outline of glacier tongue to be investigated
dhdt_map <- '2013Lirung_dhdt_total.tif'                  # dh/dt map for calibration/validation
dhdt_uncertainty <- sqrt(0.25^2+0.25^2)
cl_thres <- 25                                            # threshold from which pixel is classified as cliff (%, relative area cover)
po_thres <- 25                                            # threshold from which pixel is classified as pond (%, relative area cover)
ModRes <- 10                                              # Spatial Resolution of Model [m]
DEM_domain <- '20130518_lirung_dem_20cm_quakeshift.tif'              # DEM for model domain (for lapsing)

############
# optional additional data
############
pond_outlines <- 'Ponds_May13.shp'                        # Outlines for pond features (write NULL if not applicable)
cliff_outlines <- 'Cliffs_May13.shp'                      # Outlines for cliff features (write NULL if not applicable)
pond_raster <- '2013pondsRasterized.tif'                  # Rasterized pond outlines (use featureRasterize_d2EB.R)
cliff_raster <- '2013cliffsRasterized.tif'                # Rasterized cliff outlines (use featureRasterize_d2EB.R)
dhdt_map_emer <- '2013Lirung_dhdtEmergence.tif'           # map of emergence velocity, if not available set to 'NA'
deb_thick_map <- 'modelledThickness_OestremModel.tif'     # debris thickness map if available; if not available set to 'NA'
pixID <- 'MeasuredThickness_1m.tif'                # Location where thickness measurements exist and initial runs can be executed

debrisview_domain <- 'debrisView_2013.tif'               # debris view raster for topographic shading
cliffFileName_out <- '2013maycliffsRasterized.tif'
pondFileName_out <- '2013maypondsRasterized.tif'

# parameter ranges + debris thickness range

paramSize <- 50

cond_debris_file <- path&'\\MCRanges\\debris_conductivity_n'&paramSize&'.csv'             # debris conductivity range for Lirung Glacier
por_debris_file <- path&'\\MCRanges\\por_lognormal_n'&paramSize&'.csv'                    # debris porosity range for Lirung Glacier
den_debris_file <- path&'\\MCRanges\\den_lognormal_n'&paramSize&'.csv'                    # debris porosity range for Lirung Glacier
z0_debris_file <- path&'\\MCRanges\\z0_lognormal_n'&paramSize&'.csv'                      # surface roughness range for Lirung Glacier
alpha_debris_file <- path&'\\MCRanges\\alpha_n'&paramSize&'.csv'             # debris albedo range for Lirung Glacier
thickness_debris_file <- path&'\\MCRanges\\DebrisThickness_lognormal_n'&paramSize&'.csv'  # debris thickness range for Lirung Glacier

# uncertainty range of climate variables
tair_file <- path&'\\MCRanges\\TAir_n'&paramSize&'.csv'                                   # Air Temperature
rh_file <- path&'\\MCRanges\\RH_n'&paramSize&'.csv'                                   # Air Temperature
ws_file <- path&'\\MCRanges\\ws_n'&paramSize&'.csv'                                   # Wind speed
LW_file <- path&'\\MCRanges\\LW_n'&paramSize&'.csv'                                   # incoming longwave radiation

climData <- 'AWSLirung_onglacier_2013_Data.csv'           # Climate data file
AWSelev <- 4058                                           # Elevation of climate data station

location <- 'Lirung'                                      # name of model location for file outputs
lat <- 28.23259709                                        # Latitude (decimal degrees)
lon <- 85.5621322                                         # Longitude (decimal degrees)

# Define whether station is located on a clean ice ('clean'), debris-covered ('debris') glacier or off-glacier ('off')
station_loc <- 'debris'

# Specify height of T_a and ws sensor (assumed height is 2 m). Used to scale to standard height.
sh_Ta <- 2;
sh_ws <- 2;

# set start and end date of calculation of EB
EBstart <- '2013/05/18 00:00:00'
EBend <- '2013/10/22 00:00:00'

seasonsDyn <- 'on';                                       # choose variable seasons (1 dry, 1 wet, 1 dry)

monin <<- '06/15'                                         # define beginning of wet season
monout <<- '09/19'                                        # define end of wet season

winin <<- '1/1'                                           # define end of second dry season
winout <<- '02/28'                                        # define beginning of first dry season

albedoDyn <<- 'const';                                    # choose whether a constant debris albedo value ('const') or a time series is used ('var')

#StationName <- 'Lirung_'&toString(year);                 # on-glacier Station
#StationName_off <- 'Kyanjing_'&toString(year);           # off-glacier Station



################# 
# Paths for Extra Codes
file.sources = list.files(path_subcode, 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\EBModCore.R")
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\mOm.R")
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\TSaggregate.R")
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
altitude <<- 4058;           # Altitude of measurement site (m)
z_a <<- 2.00;                # % Height of air temp / wind / humidity measurements (m)
z_0 <<- 0.03;               # % Debris aerodynamic roughness length (m)
Lapse <<- 0.0065;            # % Temperature lapse rate (K m^-1)
k_d_wet <<- 2.03;             # Debris thermal conductivity wet season (W m^-1 K^-1)
k_d_wet_range <- c(1.40,1.70);# possible range for k_d_wet (based on Nicholson & Benn 2012)
k_d_dry <<- 0.85;                 # Debris thermal conductivity dry season (W m^-1 K^-1)
k_d_dry_range <- c(0.94,1.14);    # possible range for k_d (based on Nicholson & Benn 2012)
rho_d <<- 1588;               # Debris density (kg m^-3)
rho_d_range <- c(1000,1900);  # range for rho_d; 1496 (Reid & Brock 2010), 2700 (Nicholson & Benn 2012)
c_d <<- 948;                  # Debris specific heat capacity (J kg^-1 K^-1)
c_d_range <- c(850,1045);     # range for c_d
epsilon_d <<- 1;           # Debris emissivity
epsilon_d_range <- c(0.92,0.97); # range for epsilon_d
albedo_const <- 0.13;        # Debris albedo (will be calculated if dynamic version is chosen)

# Calculate air pressure in Pa based on altitude
p_a <- p0*((1-(Lapse*altitude/T0))^(grav*Mair/(Rgas*Lapse)));


# ========================
#  % DISTRIBUTED DATA (load distributed data including thickness, dhdt map, ponds/cliffs, etc)
# ========================

# Read Model Domain / Number of Cells

# Glacier Outline
ogrInfo(path_data&'\\Outlines\\'&glac_outline)
glac_mask<-readOGR(dsn=path_data&'\\Outlines\\'&glac_outline)
projection(glac_mask)<-projec

raster_domain <- raster()
extent(raster_domain) <- extent(glac_mask)
res(raster_domain) <- 1
raster_domain[] <- 0
projection(raster_domain)<-projec
raster_domain <- mask(raster_domain,glac_mask)

# Validation/Calibration domain (in this case available dh/dt map, once with and without emergence)
dhdt <- raster(path_data&'\\dhdt\\'&dhdt_map)
projection(dhdt)<-projec
dhdt <- crop(dhdt,raster_domain)
dhdt[dhdt<0]<-NA
#remove outliers from raster (mean annual melt rates exceeding 3cm/day are removed)
cutoff <- (lubridate::yday(EBend)-lubridate::yday(EBstart)) * 0.03
dhdt[dhdt>cutoff] <- NA
dhdt_res <- raster::resample(dhdt, raster_domain,method='bilinear')
dhdt_res <- mask(dhdt_res,glac_mask)      # Available calib/valid data


# Debris Thickness 
debThickVal <- NA
if(!is.na(deb_thick_map)){
  debThickVal <- raster(path_data&'\\DebrisThickness\\'&deb_thick_map)
  projection(debThickVal) <- projec
  debThickVal <- crop(debThickVal,dhdt_res)
  debThickVal <- raster::resample(debThickVal,dhdt_res,method='bilinear')
}

debThickMeas <- NA
if(!is.na(pixID)){
  debThickMeas <- raster(path_data&'\\DebrisThickness\\'&pixID)
  projection(debThickMeas) <- projec
  debThickMeas <- crop(debThickMeas,dhdt_res)
  debThickMeas <- raster::resample(debThickMeas,dhdt_res,method='bilinear')
}


# Read in Surface Features with different melt properties
raster_domain_masks <- raster_domain
cliffStack <- raster(path_data&'\\Outlines_Features\\'&cliffFileName_out)
cliffStack$typeID <- raster(path_data&'\\Outlines_Features\\'&cliffFileName_out,band=1)
cliffStack$indID <- raster(path_data&'\\Outlines_Features\\'&cliffFileName_out,band=2)
pondStack <- raster(path_data&'\\Outlines_Features\\'&pondFileName_out)
pondStack$typeID <- raster(path_data&'\\Outlines_Features\\'&pondFileName_out,band=1)
pondStack$indID <- raster(path_data&'\\Outlines_Features\\'&pondFileName_out,band=2)

if(!is.null(cliffStack)){
  raster_domain_masks[cliffStack$indID>0] <- NA
  visCliff <- as.factor(cliffStack$typeID)
  rat <- levels(visCliff)[[1]]
  rat[["pond"]] <- c("-pond","+pond")
  levels(visCliff) <- rat
}

if(!is.null(pondStack)){
  raster_domain_masks[pondStack$indID>0] <- NA
  visPond <- as.factor(pondStack$typeID)
  rat <- levels(visPond)[[1]]
  rat[["pond"]] <- c("-cliff","+cliff")
  levels(visPond) <- rat
}

# Read DEM from same year as thickness map
DEM_model<-raster(path_data&'/DEMs/'&DEM_domain)
DEM_model <- raster::resample(DEM_model, raster_domain,method='bilinear')
DEM_model <- mask(DEM_model,raster_domain)
DEM_model <- mask(DEM_model,dhdt_res)

# resample to standard model resolution
raster_domain <- aggregate(raster_domain_masks, fact = floor(ModRes/res(raster_domain_masks)) , median)

# if provided read in pixels where model should be performed
if(!is.na(pixID)){
  ModPix <- raster(path_data&'\\DebrisThickness\\'&pixID)
  ModPix <- raster::resample(ModPix, raster_domain,method='bilinear')
}

dhdt_res[which(is.na(raster_domain[]))] <- NA
##### ####### #######
dhdt_res <- aggregate(dhdt_res, fact = floor(ModRes/res(dhdt_res)) , median)
DEM_model <- aggregate(DEM_model, fact = floor(ModRes/res(DEM_model)) , median)
debThick <- aggregate(debThickVal, fact = floor(ModRes/res(debThickVal)) , median)

# Dataframe of DEM for later computations
DEM_dataframe <- cbind(seq(1,length(DEM_model[]),1),DEM_model[],coordinates(DEM_model))
#DEM_dataframe_new <- DEM_dataframe[which(!is.na(DEM_dataframe[,2])),]

DomDim <- dim(dhdt_res)
CellNo <- DomDim[1] * DomDim[2]

# View Factors for Radiation (debris view factor)
vf <- raster(path_shade&'\\'&debrisview_domain)
vf_deb <- aggregate(vf, fact = floor(ModRes/res(vf)) , median)
vf_deb <- raster::resample(vf_deb, raster_domain,method='bilinear')
vf_deb <- mask(vf_deb,raster_domain)
vf_sky <- 1 - vf_deb

if(station_loc=='debris'){
  # get raster ID of AWS on DEM raster
  AWSloc <- data.frame(lon=lon, lat=lat)
  coordinates(AWSloc) <- c("lon", "lat")
  proj4string(AWSloc) <-  CRS("+init=epsg:4326")
  coord <- spTransform(AWSloc,projec)
  ID_AWS <- raster::extract(DEM_model,coord,cellnumbers=T)[1]
}

# relative skyview factor vs the location of the AWS
vf_sky_actual <- vf_sky / vf_sky[ID_AWS]

writeRaster(raster_domain_masks, filename = path_figs&'\\DEBRISDOMAIN'&ModRes&'m.tiff', 'GTiff',overwrite=TRUE)
#projection(visPond) <- projection(raster_domain_masks)
#projection(visCliff) <- projection(raster_domain_masks)
writeRaster(merge(raster_domain_masks,visPond, visCliff), filename = path_figs&'\\ALLDOMAIN'&ModRes&'m.tiff', 'GTiff',overwrite=TRUE)

#png(file=path_figs&'\\ModelDomain_res'&ModRes&'m.png', res = 160,width=1800,height=1800)
#par(mar=c(7,7,4,4),cex.lab=1.5,cex.axis=1.5)
#layout(matrix(c(1,2), nrow = 2, ncol = 1,byrow=TRUE))
#par(xpd=FALSE)
#plot(merge(raster_domain_masks,visPond, visCliff),col=c('grey',terrain.colors(4)),legend=F,xlab = 'Easting [m]', ylab ='Northing [m]')
#plot(cliff_mask,add=T, border='red',xaxt='n', yaxt='n')
#plot(pond_mask,add=T, border='blue',xaxt='n', yaxt='n')
#grid(NULL,NULL)
#par(xpd=TRUE)
#legend('bottom',legend=c('debris',"cliff - pond","cliff + pond","pond - cliff","pond + cliff"),fill=c('grey',terrain.colors(4)),horiz = TRUE,inset = -0.14)
#plot(DEM_model,col=terrain.colors(200),legend=F,xlab = 'Easting [m]', ylab ='Northing [m]')
#dev.off()

# ========================
#  % CLIMATE DATA
# ========================

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


# ========================
# TEST SCENARIO DEBUGGING in SINGLE RUN MODE
# ========================
# Calculate d2EB at each individual grid cell with available time series
EBCoreInput <- list()

# Gap Filling
gapFill <- function(TS_withHoles){    # Gap Filling Function
  for (k in 25:length(TS_withHoles)){
    if(is.na(TS_withHoles[k])){TS_withHoles[k] <- TS_withHoles[k-24]}
  }
  for (k in 1:24){
    if(is.na(TS_withHoles[k])){TS_withHoles[k] <- TS_withHoles[k+24]}
  }
  return(TS_withHoles)
}

# Lapse climate variables to model pixels

LR_TAIR_off <- read.csv(path_data&'\\'&'LapseoffGlacier.csv')[IDstart:IDend,3]
LR_TAIR_on <- read.csv(path_data&'\\'&'LapseonGlacier.csv')[IDstart:IDend,3]
LR_RH <- read.csv(path_data&'\\'&'LR_Langtang_RH.csv')[IDstart:IDend,3]

# Initial Debris Thickness
d_raster <- debThick

# ========================
# MULTI CORE RUN
# ========================

# Prepare parameter space for Lirung case
cond <- read.csv(cond_debris_file, header = T)
k_wet_sample <- cond$conductivity_tot
k_dry_sample <- cond$conductivity_tot
rho_d_sample <- read.csv(den_debris_file, header = T)$x
z0_d_sample <- read.csv(z0_debris_file, header = T)$x
alpha_d_sample <- read.csv(alpha_debris_file, header = T)$x
d_sample <- read.csv(thickness_debris_file, header = T)$x

# Prepare climate variable uncertainty
tair_range <- read.csv(tair_file, header = T)$x
rh_range <- read.csv(rh_file, header = T)$x
ws_range <- read.csv(ws_file, header = T)$x
LW_range <- read.csv(LW_file, header = T)$x

# Gap Filling
gapFill <- function(TS_withHoles){    # Gap Filling Function
  for (k in 25:length(TS_withHoles)){
    if(is.na(TS_withHoles[k])){TS_withHoles[k] <- TS_withHoles[k-24]}
  }
  for (k in 1:24){
    if(is.na(TS_withHoles[k])){TS_withHoles[k] <- TS_withHoles[k+24]}
  }
  return(TS_withHoles)
}

# load climate data from station
EBCoreInput <- list()
EBCoreInput$timeline_str <- climInput$Time_Str[IDstart:IDend]
timeline_num <- climInput$Time_Num[IDstart:IDend]
EBCoreInput$SWout <- gapFill(climInput$KUPW[IDstart:IDend])             # Upwelling shortwave radiation (W m^-2)
EBCoreInput$SWin <- gapFill(climInput$KINC[IDstart:IDend])            # Downwelling shortwave radiation (W m^-2)
EBCoreInput$SWin[EBCoreInput$SWin<0] <- 0
SW_orig <- EBCoreInput$SWin
EBCoreInput$LWin <- gapFill(climInput$LINC[IDstart:IDend]) # Downwelling longwave radiation (W m^-2)
LWinorig <- EBCoreInput$LWin
EBCoreInput$T_a <- gapFill(climInput$TAIR[IDstart:IDend] + 273.15);      	    # Air temperature (K)
Ta_orig <- EBCoreInput$T_a
EBCoreInput$u <- gapFill(climInput$WSPD[IDstart:IDend])                     # Wind speed (m s^-1)
uorig <- EBCoreInput$u
EBCoreInput$wd <- gapFill(climInput$WDIR[IDstart:IDend])
EBCoreInput$RH_a <- gapFill(climInput$RH[IDstart:IDend])          # Air relative humidity (%)
RHa_orig <- EBCoreInput$RH_a

# Check possible missing data
if(is.null(climInput$LUPW[IDstart:IDend])){
  EBCoreInput$LWout <- gapFill((climInput$TS[IDstart:IDend] + 273.15)^4*5.67*10^(-8)*epsilon_d)      # use outgoing LW for surface temperature
}else{
  EBCoreInput$LWout <- gapFill(climInput$LUPW[IDstart:IDend])             # Upwelling longwave radiation (W m^-2)
}
if(is.null(climInput$TS[IDstart:IDend])){
  EBCoreInput$T_s_data <- (EBCoreInput$LWout/5.67/10^(-8)/epsilon_d)^(1/4)      # use outgoing LW for surface temperature
}else{
  EBCoreInput$T_s_data <- gapFill(climInput$TS[IDstart:IDend] + 273.15);	    # Surface temperature (K)
}

# Albedo Model (Constant vs. Variable)
switch(albedoDyn,
       const = {                       # constant Debris albedo
         albedo_d = albedo_const + EBCoreInput$SWin*0},            
       var = {                         # variable debris albedo
         albedo_d = EBCoreInput$SWout / EBCoreInput$SWin
         # remove non-sensical values
         albedo_d[albedo_d<0 | albedo_d>1 | is.infinite(albedo_d) | albedo_d==0] <- NA
         # remove albedo at night time (few values, not reasonable, scatter effects at low solar angle)
         albedo_d[hour(EBCoreInput$timeline_str) < 9 | hour(EBCoreInput$timeline_str) > 18] = NaN})
EBCoreInput$albedo_d <- albedo_d

EBCoreInput$mod <- 'debris'
EBModType <<- 'no_TS'

nMC <- paramSize
#Modcount <- 1
dhdt_param_mod <- vector()
dhdt_clim_mod <- vector()
meltCD <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
tairvec <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
Snet <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
Lnet <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
QH <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
QLE <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
refreeze <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)

CDrange <<- 0.5;                   # Range for central difference calculation of the derivative of surface flux when calculating surface temperature T_s.
library(lubridate)

# Read Output files and create Figures

ModelRunSeries <- DEM_dataframe[which(!is.na(dhdt_res[])&!is.na(ModPix[])),1];
bestselect = 25
mat_opt_param <- list()
sort_depth <- sort(debThick[which(!is.na(dhdt_res[])&!is.na(ModPix[]))],index.return=T); sort_depth <- sort_depth$ix


#####
# Validation of the point scale model at the AWS site
#####
resultsfolder <- path&'\\ModelOutput\\results_pointscaleAWS'
modThick <- dhdt_res * NA # thickness as determined from the Monte Carlo
modThick_sd <- dhdt_res * NA # thickness standard deviation as determined from the Monte Carlo
modCond <- dhdt_res * NA # conductivity as determined from the Monte Carlo
modCond_sd <- dhdt_res * NA # conductivity standard deviation as determined from the Monte Carlo
depth_matrix <- matrix(NA,nrow=bestselect, ncol=1)
cond_matrix <- matrix(NA,nrow=bestselect, ncol=1)
rho_matrix <- matrix(NA,nrow=bestselect, ncol=1)

modmelt_matrix <- matrix(NA,nrow=bestselect, ncol=1)

#read melt curves
meltruns <- read.csv(file = resultsfolder&'\\melt'&ID_AWS&'.csv') # melt from all different runs
tsurf <- read.csv(file = resultsfolder&'\\tempsurf'&ID_AWS&'.csv') # tsurf from all different runs

qsens <- read.csv(file = resultsfolder&'\\QH'&ID_AWS&'.csv') # sensible heat from all different runs
lesens <- read.csv(file = resultsfolder&'\\QLE'&ID_AWS&'.csv') # latent heat from all different runs
lwsens <- read.csv(file = resultsfolder&'\\Lnet'&ID_AWS&'.csv') # longwave radiation from all different runs
swsens <- EBCoreInput$SWin - EBCoreInput$SWout

cumMelt <- cumsum(meltruns[,1:nMC+1])
maxCum <- max(cumMelt)
totMelt <- unlist(cumMelt[dim(meltruns)[1],1:nMC]) # maximum melt at end of model period
measMassLoss <- dhdt_res[ID_AWS] # measured mass Loss for pixel

mOmout <- vector()
for(Ts_cor in 1:nMC){
mOmout[Ts_cor] <- mOm(EBCoreInput$T_s_data,tsurf[,Ts_cor+1])
}

bestFitTemp <- sort(mOmout,index.return = TRUE)$ix[1:bestselect]
bestfit <- which(rank(abs(totMelt - measMassLoss), ties.method='min') <=  bestselect) # model runs that match closest
leastfit <- which(rank(abs(totMelt - measMassLoss), ties.method='max') >  nMC-bestselect)

diuObs <- diuCyc(EBCoreInput$T_s_data,EBCoreInput$timeline_str)
diuCalc <- matrix(NA,nrow=24,ncol=nMC)
for(diuFit in 1:nMC){
diuCalc[,diuFit] <- diuCyc(tsurf[,diuFit+1],EBCoreInput$timeline_str)[,2]
}

diumOm <- vector()
for(Ts_cor in 1:nMC){
  diumOm[Ts_cor] <- mOm(diuObs[,2],diuCalc[,Ts_cor])
}

rmseT <- vector()
mbeT <- vector()
for(RMSERun in 1: bestselect){
rmseT[RMSERun] <- sqrt(1/length(tsurf[,bestFitTemp[RMSERun]+1])*sum((tsurf[,bestFitTemp[RMSERun]+1]-EBCoreInput$T_s_data)^2,na.rm=T)); # RMSE
mbeT[RMSERun] <- 1/length(tsurf[,bestFitTemp[RMSERun]+1])*sum(tsurf[,bestFitTemp[RMSERun]+1]-EBCoreInput$T_s_data,na.rm=T);
}

# Plot example for surface temperature fit of point scale model
library(matrixStats)
png(file=path_figs&'\\ModelResultsAWSLoc.png', res = 300,width=3600,height=1800)

par(mar=c(1,2,1,1),mai = c(0.3, 0.7, 0.1, 0.2),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1,1,2,2,3,3,4,4), nrow = 2, ncol = 4, byrow = T))

datAx1 <- as.POSIXct(EBCoreInput$timeline_str[springTime],origin = "1970-01-01")
datAx2 <- as.POSIXct(EBCoreInput$timeline_str[monsoonTime],origin = "1970-01-01")
springTime <- c(150:500)
monsoonTime <- c(1800:2100)

plot(EBCoreInput$timeline_str[springTime],EBCoreInput$T_s_data[springTime] - 273.15,type='l',lwd=2,ylim = c(0,41),xlab ='', ylab = expression('T'['surf']~~'[°C]'),xaxt='n')
meanTS_T <- rowMeans(tsurf[springTime,bestFitTemp+1])-273.15
#sdTS <- rowSds(as.matrix(tsurf[springTime,bestFitTemp+1]))
#points(EBCoreInput$timeline_str[springTime],meanTS,type='l',col='red')
#polygon(c(EBCoreInput$timeline_str[springTime], rev(EBCoreInput$timeline_str[springTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(1, 0, 0,0.3),border="red")
meanTS <- rowMeans(tsurf[springTime,bestfit+1])-273.15
sdTS <- rowSds(as.matrix(tsurf[springTime,bestfit+1]))
points(EBCoreInput$timeline_str[springTime],meanTS_T,type='l',col='red',lwd=2)
points(EBCoreInput$timeline_str[springTime],meanTS,type='l',col='blue')
polygon(c(EBCoreInput$timeline_str[springTime], rev(EBCoreInput$timeline_str[springTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 0, 1,0.3),border="blue")


abline(h = seq(0,40,10),v = seq(datAx1[1], datAx1[length(datAx1)], by = "day"), col="gray", lty=3)

plot(EBCoreInput$timeline_str[monsoonTime],EBCoreInput$T_s_data[monsoonTime] - 273.15,type='l',lwd=2,ylim = c(0,41),yaxt='n',xaxt='n',ylab = '',xlab ='')
meanTS_T <- rowMeans(tsurf[monsoonTime,bestFitTemp+1])-273.15
#sdTS <- rowSds(as.matrix(tsurf[monsoonTime,bestFitTemp+1]))
#points(EBCoreInput$timeline_str[monsoonTime],meanTS,type='l',col='red')
#polygon(c(EBCoreInput$timeline_str[monsoonTime], rev(EBCoreInput$timeline_str[monsoonTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(1, 0, 0,0.3),border="red" )
meanTS <- rowMeans(tsurf[monsoonTime,bestfit+1])-273.15
sdTS <- rowSds(as.matrix(tsurf[monsoonTime,bestfit+1]))
points(EBCoreInput$timeline_str[monsoonTime],meanTS_T,type='l',col='red',lwd=2)
points(EBCoreInput$timeline_str[monsoonTime],meanTS,type='l',col='blue')
polygon(c(EBCoreInput$timeline_str[monsoonTime], rev(EBCoreInput$timeline_str[monsoonTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 0, 1,0.3),border="blue")


abline(h = seq(0,40,10),v = seq(datAx2[1], datAx2[length(datAx2)], by = "day"), col="gray", lty=3)
legend("topright",c("measured","mod UAV","mod T"),lty=c(1,1,1), col=c("black","blue","red"),lwd=c(2,1,2),bty = "n",cex = 1.2)


plot(EBCoreInput$timeline_str[springTime],swsens[springTime],type='l',xlab ='',ylab = expression('Flux [W '~ m^{-2}~']'),ylim=c(-350,1000),xaxt='n')
meanTS <- rowMeans(lwsens[springTime,bestfit+1])
sdTS <- rowSds(as.matrix(lwsens[springTime,bestfit+1]))
points(EBCoreInput$timeline_str[springTime],meanTS,type='l',col='red')
polygon(c(EBCoreInput$timeline_str[springTime], rev(EBCoreInput$timeline_str[springTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(1, 0, 0,0.3),border="red" )
meanTS <- rowMeans(qsens[springTime,bestfit+1])
sdTS <- rowSds(as.matrix(qsens[springTime,bestfit+1]))
points(EBCoreInput$timeline_str[springTime],meanTS,type='l',col='green')
polygon(c(EBCoreInput$timeline_str[springTime], rev(EBCoreInput$timeline_str[springTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 1, 0,0.3),border="green" )
meanTS <- rowMeans(lesens[springTime,bestfit+1])
sdTS <- rowSds(as.matrix(lesens[springTime,bestfit+1]))
points(EBCoreInput$timeline_str[springTime],meanTS,type='l',col='blue')
polygon(c(EBCoreInput$timeline_str[springTime], rev(EBCoreInput$timeline_str[springTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 0, 1,0.3),border="blue" )
par(new = TRUE)
plot(EBCoreInput$timeline_str[springTime],rowMeans(meltruns[springTime,bestfit+1])*1000,type="l", lty=3,col='red',lwd=1,xaxt='n',yaxt='n',ylab='')
axis(side=4, at = seq(0.2,1.3,0.2),col='red')
mtext("melt [mm]", side = 4, line = 3,col='red')
axis.POSIXct(1, at = seq(datAx1[1], datAx1[length(datAx1)], by = "day"), format = "%b-%d")
abline(h = seq(-300,1000,100),v = seq(datAx1[1], datAx1[length(datAx1)], by = "day"), col="gray", lty=3)

plot(EBCoreInput$timeline_str[monsoonTime],swsens[monsoonTime],type='l',xlab ='',ylab = '',ylim=c(-350,1000),xaxt='n')
meanTS <- rowMeans(lwsens[monsoonTime,bestfit+1])
sdTS <- rowSds(as.matrix(lwsens[monsoonTime,bestfit+1]))
points(EBCoreInput$timeline_str[monsoonTime],meanTS,type='l',col='red')
polygon(c(EBCoreInput$timeline_str[monsoonTime], rev(EBCoreInput$timeline_str[monsoonTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(1, 0, 0,0.3),border="red" )

meanTS <- rowMeans(qsens[monsoonTime,bestfit+1])
sdTS <- rowSds(as.matrix(qsens[monsoonTime,bestfit+1]))
points(EBCoreInput$timeline_str[monsoonTime],meanTS,type='l',col='green')
polygon(c(EBCoreInput$timeline_str[monsoonTime], rev(EBCoreInput$timeline_str[monsoonTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 1, 0,0.3),border="green" )

meanTS <- rowMeans(lesens[monsoonTime,bestfit+1])
sdTS <- rowSds(as.matrix(lesens[monsoonTime,bestfit+1]))
points(EBCoreInput$timeline_str[monsoonTime],meanTS,type='l',col='blue')
polygon(c(EBCoreInput$timeline_str[monsoonTime], rev(EBCoreInput$timeline_str[monsoonTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 0, 1,0.3),border="blue" )
par(new = TRUE)
plot(EBCoreInput$timeline_str[monsoonTime],rowMeans(meltruns[monsoonTime,bestfit+1])*1000,type="l", lty=3,col='red',lwd=1,xaxt='n',yaxt='n',ylab='')
axis.POSIXct(1, at = seq(datAx2[1], datAx2[length(datAx2)], by = "day"), format = "%b-%d")
abline(h = seq(-300,1000,100),v = seq(datAx2[1], datAx2[length(datAx2)], by = "day"), col="gray", lty=3)
legend("topright",c("SW","LW","H","LE","melt"),lty=c(1,1,1,1,3), col=c("black","red","green","blue","red"),lwd=c(1,1,1,1,1),bty = "n",cex = 1.2)
dev.off()

# Plot example for mass loss fit

lirungUDG <- read.csv('F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\massloss2013Lirung.csv',header = T)
UDGtimeStr <- as.POSIXct(lirungUDG$timVec,origin='1970-01-01')
UDGheight <- max(lirungUDG$X.1,na.rm=T) - lirungUDG$X.1
matchDataStart <- which(UDGtimeStr==EBCoreInput$timeline_str[1])
UDGheight <- UDGheight - UDGheight[matchDataStart]

png(file=path_figs&'\\AWSLoc_massLoss.png', res = 300,width=3600,height=1800)
par(mar=c(6,3,1,1),mai = c(0.5, 1, 0.1, 0.1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1), nrow = 1, ncol = 1, byrow = T))
datAx1 <- as.POSIXct(EBCoreInput$timeline_str,origin = "1970-01-01")
meanTS <- rowMeans(cumMelt[,bestFitTemp])
sdTS <- rowSds(as.matrix(cumMelt[,bestFitTemp]))
plot(EBCoreInput$timeline_str,meanTS,type='l',lwd=2,ylim = c(0,2.3),xlab ='', ylab = expression('melt [m]'),xaxt='n',col='red')
polygon(c(EBCoreInput$timeline_str, rev(EBCoreInput$timeline_str)), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(1, 0, 0,0.3),border="red")
axis.POSIXct(1, at = seq(datAx1[1], datAx1[length(datAx1)], by = "month"), format = "%b-%d")
meanTS <- rowMeans(cumMelt[,bestfit])
sdTS <- rowSds(as.matrix(cumMelt[,bestfit]))
points(EBCoreInput$timeline_str,meanTS,type='l',lwd=2,ylim = c(0,2.3),xlab ='', ylab = 'n',col='blue')
polygon(c(EBCoreInput$timeline_str, rev(EBCoreInput$timeline_str)), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 0, 1,0.3),border="blue")
abline(h = seq(0,3,0.2),v = seq(datAx1[1], datAx1[length(datAx1)], by = "month"), col="gray", lty=3)
points(datAx1[length(datAx1)],measMassLoss,pch=19)
arrows(datAx1[length(datAx1)], measMassLoss-dhdt_uncertainty/2, datAx1[length(datAx1)], measMassLoss+dhdt_uncertainty/2, length=0.05, angle=90, code=3,lwd=2)
points(UDGtimeStr[matchDataStart :length(UDGtimeStr)],UDGheight[matchDataStart:length(UDGtimeStr)],type='l',lwd=2,ylim = c(0,2.3),xlab ='', ylab = 'n',col='black')
dev.off()

png(file=path_figs&'\\AWSLoc_massLoss_inset.png', res = 300,width=3600,height=1800)
par(mar=c(6,3,1,1),mai = c(0.5, 1, 0.1, 0.1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1), nrow = 1, ncol = 1, byrow = T))
datAx1 <- as.POSIXct(EBCoreInput$timeline_str,origin = "1970-01-01")
meanTS <- rowMeans(cumMelt[,bestFitTemp])
sdTS <- rowSds(as.matrix(cumMelt[,bestFitTemp]))
plot(EBCoreInput$timeline_str,meanTS,type='l',lwd=2,xlim=c(EBCoreInput$timeline_str[1],EBCoreInput$timeline_str[700]),ylim = c(0,0.5),xlab ='', ylab = expression('melt [m]'),xaxt='n',col='red')
polygon(c(EBCoreInput$timeline_str, rev(EBCoreInput$timeline_str)), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(1, 0, 0,0.3),border="red")
axis.POSIXct(1, at = seq(datAx1[1], datAx1[length(datAx1)], by = "week"), format = "%b-%d")
meanTS <- rowMeans(cumMelt[,bestfit])
sdTS <- rowSds(as.matrix(cumMelt[,bestfit]))
points(EBCoreInput$timeline_str,meanTS,type='l',lwd=2,ylim = c(0,2.3),xlab ='', ylab = 'n',col='blue')
polygon(c(EBCoreInput$timeline_str, rev(EBCoreInput$timeline_str)), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 0, 1,0.3),border="blue")
abline(h = seq(0,2,0.1),v = seq(datAx1[1], datAx1[length(datAx1)], by = "week"), col="gray", lty=3)
points(datAx1[length(datAx1)],measMassLoss,pch=19)
arrows(datAx1[length(datAx1)], measMassLoss-dhdt_uncertainty/2, datAx1[length(datAx1)], measMassLoss+dhdt_uncertainty/2, length=0.05, angle=90, code=3,lwd=2)
points(UDGtimeStr[matchDataStart :length(UDGtimeStr)],UDGheight[matchDataStart:length(UDGtimeStr)],type='l',lwd=2,ylim = c(0,2.3),xlab ='', ylab = 'n',col='black')
dev.off()

png(file=path_figs&'\\AWSLoc_conductivity.png', res = 300,width=900,height=1800)
par(mar=c(2,1,1,1),mai = c(0.1, 1, 0.1, 0.1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1), nrow = 1, ncol = 1, byrow = T))
boxplot(cbind(k_wet_sample[bestFitTemp],k_wet_sample[bestfit]),outline=F,border = c('red','blue'),col='white',horizontal=F,xaxt='n',ylab = expression('conductivity [W '~ m^{-1}~ K^{-1}~']'))
myjitter<-jitter(rep(1, length(k_wet_sample[bestFitTemp])), amount=0.2)
points(myjitter, k_wet_sample[bestFitTemp], pch=20, col='red',cex=1.5)
myjitter<-jitter(rep(2, length(k_wet_sample[bestfit])), amount=0.2)
points(myjitter, k_wet_sample[bestfit], pch=20, col='blue',cex=1.5)
dev.off()

png(file=path_figs&'\\AWSLoc_density.png', res = 300,width=900,height=1800)
par(mar=c(2,1,1,1),mai = c(0.1, 1, 0.1, 0.1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1), nrow = 1, ncol = 1, byrow = T))
boxplot(cbind(rho_d_sample[bestFitTemp],rho_d_sample[bestfit]),outline=F,border = c('red','blue'),col='white',horizontal=F,xaxt='n',ylab = expression('debris density [kg '~ m^{-3}~']'))
myjitter<-jitter(rep(1, length(rho_d_sample[bestFitTemp])), amount=0.2)
points(myjitter, rho_d_sample[bestFitTemp], pch=20, col='red',cex=1.5)
myjitter<-jitter(rep(2, length(rho_d_sample[bestfit])), amount=0.2)
points(myjitter, rho_d_sample[bestfit], pch=20, col='blue',cex=1.5)
dev.off()

png(file=path_figs&'\\AWSLoc_z0.png', res = 300,width=900,height=1800)
par(mar=c(2,1,1,1),mai = c(0.1, 1, 0.1, 0.1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1), nrow = 1, ncol = 1, byrow = T))
boxplot(cbind(z0_d_sample[bestFitTemp],z0_d_sample[bestfit]),outline=F,col='white',border = c('red','blue'),horizontal=F,xaxt='n',ylab = expression('z'['0']~'[m]'~''))
myjitter<-jitter(rep(1, length(z0_d_sample[bestFitTemp])), amount=0.2)
points(myjitter, z0_d_sample[bestFitTemp], pch=20, col='red',cex=1.5)
myjitter<-jitter(rep(2, length(z0_d_sample[bestfit])), amount=0.2)
points(myjitter, z0_d_sample[bestfit], pch=20, col='blue',cex=1.5)
dev.off()


meanTS <- rowMeans(meltruns[springTime,bestFitTemp+1]*1000)
sdTS <- rowSds(as.matrix(meltruns[springTime,bestFitTemp+1]*1000))

plot(EBCoreInput$timeline_str[springTime],meanTS,type='l',xlab ='',ylab = expression('melt [mm]'),ylim=c(0,1.5),xaxt='n',col='red')
polygon(c(EBCoreInput$timeline_str[springTime], rev(EBCoreInput$timeline_str[springTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(1, 0, 0,0.3),border="red")

meanTS <- rowMeans(meltruns[springTime,bestfit+1]*1000)
sdTS <- rowSds(as.matrix(meltruns[springTime,bestfit+1]*1000))

points(EBCoreInput$timeline_str[springTime],meanTS,type='l',xlab ='',ylab = '',ylim=c(0,1.5),xaxt='n',col='blue')
polygon(c(EBCoreInput$timeline_str[springTime], rev(EBCoreInput$timeline_str[springTime])), c(meanTS+sdTS ,rev(meanTS-sdTS)), col = rgb(0, 0, 1,0.3),border="blue")



modThick[ID_AWS] <- median(d_sample[bestfit])
modThick_sd[ID_AWS] <- sd(d_sample[bestfit])

modCond[ID_AWS] <- median(k_wet_sample[bestfit])
modCond_sd[ID_AWS] <- sd(k_wet_sample[bestfit])

modmelt_matrix <- totMelt[bestfit]

tsurf_best <- tsurf[,bestfit]

depth_matrix <- d_sample[bestfit]
cond_matrix <- k_wet_sample[bestfit]
rho_matrix <- rho_d_sample[bestfit]
z0_matrix <- z0_d_sample[bestfit]
alphad_matrix <- alpha_d_sample[bestfit]

tair_matrix <- tair_range[bestfit]
rh_matrix <- rh_range[bestfit]
ws_matrix <- ws_range[bestfit]
LW_matrix <- LW_range[bestfit]



#########
# sensitivity to debris and climate variables at the same time
#########
resultsfolder <- path&'\\ModelOutput\\results_allsens_1000_fixedthickness'
modThick <- dhdt_res * NA # thickness as determined from the Monte Carlo
modThick_sd <- dhdt_res * NA # thickness standard deviation as determined from the Monte Carlo
modCond <- dhdt_res * NA # conductivity as determined from the Monte Carlo
modCond_sd <- dhdt_res * NA # conductivity standard deviation as determined from the Monte Carlo
z0_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))
alphad_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))
tair_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))
rh_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))
ws_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))
LW_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))
depth_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))
cond_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))
rho_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))

modmelt_matrix <- matrix(NA,nrow=bestselect, ncol=length(ModelRunSeries))

for(i in 1:length(ModelRunSeries)){
  #read melt curves
  meltruns <- read.csv(file = resultsfolder&'\\melt'&ModelRunSeries[i]&'.csv') # melt from all different runs
  cumMelt <- cumsum(meltruns[,1:nMC+1])
  maxCum <- max(cumMelt)
  totMelt <- unlist(cumMelt[dim(meltruns)[1],1:nMC]) # maximum melt at end of model period
  measMassLoss <- dhdt_res[ModelRunSeries[i]] # measured mass Loss for pixel
  
  bestfit <- which(rank(abs(totMelt - measMassLoss), ties.method='min') <=  bestselect) # model runs that match closest
  leastfit <- which(rank(abs(totMelt - measMassLoss), ties.method='max') >  nMC-bestselect)
  
  modThick[ModelRunSeries[i]] <- median(d_sample[bestfit])
  modThick_sd[ModelRunSeries[i]] <- sd(d_sample[bestfit])
  
  modCond[ModelRunSeries[i]] <- median(k_wet_sample[bestfit])
  modCond_sd[ModelRunSeries[i]] <- sd(k_wet_sample[bestfit])
  
  mat_opt_param[[i]] <- cbind(matrix(1,5,1)*measMassLoss,matrix(1,5,1)*ModPix[ModelRunSeries[i]],d_sample[bestfit],k_dry_sample[bestfit],k_wet_sample[bestfit],rho_d_sample[bestfit])
  
  modmelt_matrix[,i] <- totMelt[bestfit]
  
  depth_matrix[,i] <- d_sample[bestfit]
  cond_matrix[,i] <- k_wet_sample[bestfit]
  rho_matrix[,i] <- rho_d_sample[bestfit]
  
  z0_matrix[,i] <- z0_d_sample[bestfit]
  alphad_matrix[,i] <- alpha_d_sample[bestfit]
  
  tair_matrix[,i] <- tair_range[bestfit]
  rh_matrix[,i] <- rh_range[bestfit]
  ws_matrix[,i] <- ws_range[bestfit]
  LW_matrix[,i] <- LW_range[bestfit]
  
}

# find locations where model does not match observations
modSuccess <- vector()
for(i in 1:length(ModelRunSeries)){
  if(min(modmelt_matrix[,sort_depth][,i])<dhdt_res[which(!is.na(dhdt_res[])&!is.na(ModPix[]))][sort_depth][i]&&max(modmelt_matrix[,sort_depth][,i])>dhdt_res[which(!is.na(dhdt_res[])&!is.na(ModPix[]))][sort_depth][i]){modSuccess[i]<-1}
}

thickrangelow <- which(debThick[which(!is.na(dhdt_res[])&!is.na(ModPix[]))][sort_depth]<=0.5)
thickrangelow <- thickrangelow[!thickrangelow %in% which(is.na(modSuccess))]
thickrangemid <- which(debThick[which(!is.na(dhdt_res[])&!is.na(ModPix[]))][sort_depth]<=1&debThick[which(!is.na(dhdt_res[])&!is.na(ModPix[]))][sort_depth]>0.5)
thickrangemid <- thickrangemid[!thickrangemid %in% which(is.na(modSuccess))]
thickrangehigh <- which(debThick[which(!is.na(dhdt_res[])&!is.na(ModPix[]))][sort_depth]>1)
thickrangehigh <- thickrangehigh[!thickrangehigh %in% which(is.na(modSuccess))]

# get conductivities for different thicknesses

cond_thin <- median(cond_matrix[,sort_depth][,thickrangelow])
cond_thick <- median(cond_matrix[,sort_depth][,thickrangehigh])
cond_mid <- median(cond_matrix[,sort_depth][,thickrangemid])

png(file=path_figs&'\\ModelSensitivity_targetPixels.png', res = 300,width=3600,height=1800)
par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
layout(matrix(c(1,1,2,3,4,5,6,7,8,9), nrow =1, ncol = 10, byrow = T))

mycol <- rgb(0, 0, 255, max = 255, alpha = 70, names = "blue50")

boxplot(modmelt_matrix[,sort_depth],horizontal=T,yaxt='n',xlab = 'mass loss [m]',ylim=c(0,2.5))
points(dhdt_res[which(!is.na(dhdt_res[])&!is.na(ModPix[]))][sort_depth],seq(1,64,1),col='blue',ylim=c(0,2.5))
par(new = TRUE)
points(debThick[which(!is.na(dhdt_res[])&!is.na(ModPix[]))][sort_depth],seq(1,64,1),col='red',type = "o",ylim=c(0.2,1.8))
axis(side=3, at = seq(0.2,1.8,0.2),yaxt='n',col='red')
mtext("debris thickness [m]", side = 3, line = 3,col='red')
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))
grid(NULL,NULL)

boxplot(cond_matrix[,sort_depth],horizontal=T,xlab =expression('k [W '~ m^{-1}~ K^{-1}~']'), axes=T,yaxt='n',xaxt='n',ylim=c(0.6,2))
axis(1, at = seq(0.6, 2, by = 1.4), las=2)
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))

boxplot(rho_matrix[,sort_depth],horizontal=T,xlab = expression(paste(rho,' [kg '~ m^{-3}~']')), axes=T,xaxt='n',yaxt='n',ylim=c(1200,2000))
axis(1, at = seq(1200, 2000, by = 800), las=2)
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))

boxplot(z0_matrix[,sort_depth],horizontal=T,xlab = expression('z'['0']~~'[m]'), axes=T,xaxt='n',yaxt='n',ylim=c(0.005,0.35),log = "x")
axis(1, at = seq(0.005, 0.35, by = 0.345), las=2)
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))

boxplot(alphad_matrix[,sort_depth],horizontal=T,xlab = 'albedo [-]',axes=T,xaxt='n',yaxt='n',ylim=c(0.06,0.2))
axis(1, at = seq(0.06, 0.2, by = 0.14), las=2)
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))

boxplot(tair_matrix[,sort_depth],horizontal=T,xlab = expression('dT'['air']~~'[°C]'),axes=T,xaxt='n',yaxt='n',ylim=c(-2,2))
axis(1, at = seq(-2, 2, by = 4), las=2)
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))

boxplot(rh_matrix[,sort_depth],horizontal=T, xlab = expression('dRH'['air']~~'[%]'),axes=T,xaxt='n',yaxt='n',ylim=c(-20,20))
axis(1, at = seq(-20, 20, by = 40), las=2)
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))

boxplot(ws_matrix[,sort_depth],horizontal=T, xlab = expression('dws [m '~ s^{-1}~']'),axes=T,xaxt='n',yaxt='n',ylim=c(0.4,2.5))
axis(1, at = seq(0.4, 2.5, by = 2.1), las=2)
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))

boxplot(LW_matrix[,sort_depth],horizontal=T, xlab = expression('dLW [W '~ m^{-2}~']'),axes=T,xaxt='n',yaxt='n',ylim=c(-50,50))
axis(1, at = seq(-50, 50, by = 100), las=2)
abline(h=c(which(is.na(modSuccess))), col=c("grey"), lty=c(1), lwd=c(2))
abline(h=43.5, col=c("red"), lty=c(1), lwd=c(2))
abline(h=17.5, col=c("red"), lty=c(1), lwd=c(2))

dev.off()

# average parameter values
rhmean <- apply(rh_matrix, 2, mean)
tairmean <- apply(tair_matrix, 2, mean)
wsmean <- apply(ws_matrix, 2, mean)
LWmean <- apply(LW_matrix, 2, mean)
albedomean <- apply(alphad_matrix, 2, median)
z0mean <- apply(z0_matrix, 2, median)
rhomean <- apply(rho_matrix, 2, median)
condmean <- apply(cond_matrix, 2, median)

####################################
# All Results from complete glacier surface
####################################
resultsfolder <- path&'\\ModelOutput\\results_distributedscale'
ModelRunSeries <- DEM_dataframe[which(!is.na(dhdt_res[])),1];
modCond <- dhdt_res * NA # conductivity as determined from the Monte Carlo
modmassloss <- dhdt_res * NA # modelled mass loss
modmasslossSD <- dhdt_res * NA # modelled mass loss standard deviation
modmasslossfixcond <- dhdt_res * NA
modmasslossmean <- dhdt_res * NA
modmasslossmedian <- dhdt_res * NA
cond_matrix <- matrix(NA,nrow=paramSize, ncol=1)
nMC <- paramSize

modmelt_matrix <- matrix(NA,nrow=50, ncol=1)
tothourlymelt <- matrix(NA,nrow = length(EBCoreInput$T_a),ncol=length(ModelRunSeries))
actualTa <- matrix(NA,nrow = length(EBCoreInput$T_a),ncol=length(ModelRunSeries))
actualSW <- matrix(NA,nrow = length(EBCoreInput$SWin),ncol=length(ModelRunSeries))
vf_sky_actual[ModelRunSeries][which(is.na(vf_sky_actual[ModelRunSeries]))] <-1

#read melt curves
for(i in 1:length(ModelRunSeries)){
  #read melt curves
  if(file.exists(resultsfolder&'\\melt'&ModelRunSeries[i]&'.csv') == TRUE){
    
  meltruns <- read.csv(file = resultsfolder&'\\melt'&ModelRunSeries[i]&'.csv') # melt from all different runs
  cumMelt <- cumsum(meltruns[,1:nMC+1])
  maxCum <- max(cumMelt)
  totMelt <- unlist(cumMelt[dim(meltruns)[1],1:nMC]) # maximum melt at end of model period
  measMassLoss <- dhdt_res[ModelRunSeries[i]] # measured mass Loss for pixel
  
  bestfit <- which(rank(abs(totMelt - measMassLoss), ties.method='min') <=  1) # model runs that match closest

  
  modCond[ModelRunSeries[i]] <- k_wet_sample[bestfit]
  
  modmassloss[ModelRunSeries[i]] <- totMelt[bestfit]
  modmasslossmean[ModelRunSeries[i]] <- mean(totMelt,na.rm=T)
  modmasslossmedian[ModelRunSeries[i]] <- median(totMelt,na.rm=T)
  modmasslossSD[ModelRunSeries[i]] <- sd(totMelt)
  modmasslossfixcond[ModelRunSeries[i]] <- totMelt[which.min(abs(k_wet_sample-median(k_wet_sample)))]

  # Data for the TI model
  tothourlymelt[,i] <- meltruns[,bestfit+1]
  actualTa[,i] <- LR_TAIR_on * (DEM_model[ModelRunSeries[i]] - AWSelev) + EBCoreInput$T_a
  actualSW[,i] <- SW_orig * vf_sky_actual[ModelRunSeries[i]]
  }
}


# Draw ROI
r1NaM <- is.na(as.matrix(dhdt_res))
colNotNA <- which(colSums(r1NaM) != nrow(dhdt_res))
rowNotNA <- which(rowSums(r1NaM) != ncol(dhdt_res))
r3Extent <- extent(dhdt_res, rowNotNA[1], rowNotNA[length(rowNotNA)],
                   colNotNA[1], colNotNA[length(colNotNA)])

# Cliff/Pond Pixels for pre and post model period
cliff10m <-  raster::resample(cliffStack$X2013maycliffsRasterized,DEM_model)

cliffoct <- raster::resample(raster(path_data&'\\Outlines_Features\\2013octcliffsRasterized.tif'),DEM_model)
pondoct <- raster::resample(raster(path_data&'\\Outlines_Features\\2013octpondsRasterized.tif'),DEM_model)

pond10m <-  raster::resample(pondStack$X2013maypondsRasterized,DEM_model)

dhdt_res2 <- dhdt_res
dhdt_res2[!is.na(cliff10m)] <- NA
dhdt_res2[!is.na(pond10m)] <- NA
#dhdt_res2[!is.na(cliffoct)] <- NA
#dhdt_res2[!is.na(pondoct)] <- NA
modmassloss[!is.na(cliff10m)] <- NA
modmassloss[!is.na(pond10m)] <- NA
#modmassloss[!is.na(cliffoct)] <- NA
#modmassloss[!is.na(pondoct)] <- NA
modmasslossmean[!is.na(cliff10m)] <- NA
modmasslossmean[!is.na(pond10m)] <- NA
#modmasslossmean[!is.na(cliffoct)] <- NA
#modmasslossmean[!is.na(pondoct)] <- NA
modmasslossmedian[!is.na(cliff10m)] <- NA
modmasslossmedian[!is.na(pond10m)] <- NA
#modmasslossmedian[!is.na(cliffoct)] <- NA
#modmasslossmedian[!is.na(pondoct)] <- NA
modmasslossfixcond[!is.na(cliff10m)] <- NA
modmasslossfixcond[!is.na(pond10m)] <- NA
#modmasslossfixcond[!is.na(cliffoct)] <- NA
#modmasslossfixcond[!is.na(pondoct)] <- NA
modCond[!is.na(cliff10m)] <- NA
modCond[!is.na(pond10m)] <- NA
#modCond[!is.na(cliffoct)] <- NA
#modCond[!is.na(pondoct)] <- NA

par(mar=c(1,5,2,2),cex.lab=1.5,cex.axis=1.5)
#par(mfrow=c(1,4))
layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = T))
plot(crop(dhdt_res2,r3Extent),zlim=c(0,4))
plot(crop(modmassloss,r3Extent),zlim=c(0,4))
plot(crop(modmasslossmedian,r3Extent),zlim=c(0,4))

# convert to [m/d] with measured time period
#(lubridate::yday(EBend)-lubridate::yday(EBstart))

writeRaster(crop(modmassloss,r3Extent),path&'\\ModelOutput\\spatialResults\\modelledMassLoss2013.tif',format='GTiff',overwrite=T)
writeRaster(crop(dhdt_res2,r3Extent),path&'\\ModelOutput\\spatialResults\\measuredMassLoss2013.tif',format='GTiff',overwrite=T)
writeRaster(crop(modmasslossfixcond,r3Extent),path&'\\ModelOutput\\spatialResults\\modelledMassLoss_fixedkd2013.tif',format='GTiff',overwrite=T)

writeRaster(crop(modmasslossSD/modmassloss*100,r3Extent),path&'\\ModelOutput\\spatialResults\\errorpercent2013.tif',format='GTiff',overwrite=T)
writeRaster(crop(modCond,r3Extent),path&'\\ModelOutput\\spatialResults\\modCond2013.tif',format='GTiff',overwrite=T)
writeRaster(crop(debThick,r3Extent),path&'\\ModelOutput\\spatialResults\\thickness2013.tif',format='GTiff',overwrite=T)


# Bin mass loss by elevation band

coords_extent <- coordinates(DEM_model)

binned_matrix_model <- as.data.frame(cbind(coords_extent[which(!is.na(modmassloss)[]),2],DEM_model[!is.na(modmassloss)],modmassloss[!is.na(modmassloss)]))

binnedmodel <- binned_matrix_model %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 2)

binned_matrix_model_fix <- as.data.frame(cbind(coords_extent[which(!is.na(modmassloss)[]),2],DEM_model[!is.na(modmassloss)],modmasslossfixcond[!is.na(modmasslossfixcond)]))

binnedmodelfix <- binned_matrix_model_fix %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 1)

binned_matrix_meas <- as.data.frame(cbind(coords_extent[which(!is.na(modmassloss)[]),2],DEM_model[!is.na(dhdt_res2)],dhdt_res[!is.na(dhdt_res2)]))

binnedmeas <- binned_matrix_meas %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 3)

binnedresults <- rbind(binnedmodelfix,binnedmodel,binnedmeas)


# RMSE
dDays <- (lubridate::yday(EBend)-lubridate::yday(EBstart))
rmsemod <- sqrt(sum((binnedmeas[,3] - binnedmodel[,3])^2,na.rm=T) / length(which(is.na((binnedmeas[,3] - binnedmodel[,3])^2)==F))) / dDays * 100
rmsemodfix <- sqrt(sum((binnedmeas[,3] - binnedmodelfix[,3])^2,na.rm=T) / length(which(is.na((binnedmeas[,3] - binnedmodelfix[,3])^2)==F)))/ dDays * 100

mbemod <- sum(binnedmeas[,3] - binnedmodel[,3],na.rm=T) / length(which(is.na((binnedmeas[,3] - binnedmodel[,3])^2)==F))/ dDays * 100 # MBE
mbemodfix <- sum(binnedmeas[,3] - binnedmodelfix[,3],na.rm=T) / length(which(is.na((binnedmeas[,3] - binnedmodelfix[,3])^2)==F))/ dDays * 100 # MBE



png(file=path_figs&'\\Distributed_Boxplots.png', res = 300,width=900,height=1800)
binnedresults %>%
ggplot(aes(y = bin_dist, x = V3,fill=factor(type),colour = factor(type))) +
  geom_boxplot(outlier.colour = NULL) + 
  theme_bw() + 
  scale_y_discrete(breaks=seq(31233,31247,1),labels=seq(100,1500,100),position = "right") +
  labs(x = expression(paste("mass loss [m]")), y = "distance to terminus [m]") +
  theme(legend.position="none")
dev.off()

thickQuantiles <- quantile(debThick[],probs = seq(0,1,0.10),na.rm=T)

percmelt <- matrix(NA,nrow=length(tothourlymelt[,1]),ncol=length(seq(0,1,0.10)))
percTa <- matrix(NA,nrow=length(tothourlymelt[,1]),ncol=length(seq(0,1,0.10)))
percSW <- matrix(NA,nrow=length(tothourlymelt[,1]),ncol=length(seq(0,1,0.10)))
for(k in 1:length(tothourlymelt[,1])){
  meltClean <- tothourlymelt[k,]
  meltClean[which(tothourlymelt[3769,]==0)] <- NA
  percmelt[k,] <- quantile(meltClean,probs = seq(0,1,0.10),na.rm=T)
  percTa[k,] <- quantile(actualTa[k,],probs = seq(0,1,0.10),na.rm=T)
  percSW[k,] <- quantile(actualSW[k,],probs = seq(0,1,0.10),na.rm=T)
  
}

tothourlymelt[,which(tothourlymelt[3769,]==0)] <- NA
sort(abs(debThick[ModelRunSeries]-thickQuantiles[2]),descending = T)

thinDeb <- which(!is.na(match(abs(debThick[ModelRunSeries]-thickQuantiles[2]),sort(abs(debThick[ModelRunSeries]-thickQuantiles[2]))[1:10])))
midDeb <- which(!is.na(match(abs(debThick[ModelRunSeries]-thickQuantiles[6]),sort(abs(debThick[ModelRunSeries]-thickQuantiles[6]))[1:10])))
higDeb <- which(!is.na(match(abs(debThick[ModelRunSeries]-thickQuantiles[10]),sort(abs(debThick[ModelRunSeries]-thickQuantiles[10]))[1:10])))

permelt_low <- rowMeans(tothourlymelt[,higDeb])
permelt_mid <- rowMeans(tothourlymelt[,midDeb])
permelt_hig <- rowMeans(tothourlymelt[,thinDeb])

timXAxis <- as.numeric(EBCoreInput$timeline_str)
daterange=c(as.POSIXct(timXAxis[1]-17*3600*24,origin='1970-01-01'), as.POSIXct(timXAxis[length(timXAxis)]-17*3600*24,origin='1970-01-01'))
library(RColorBrewer)

png(file=path_figs&'\\climateforcingmelt.png', res = 160,width=1800,height=900)
par(mar=c(3,7,2,5),cex.lab=1.2,cex.axis=1.2)
layout(matrix(c(1,1), nrow = 1, ncol = 1, byrow = FALSE))
plot(timXAxis,EBCoreInput$T_a-273.15,col=brewer.pal(3, 'Reds')[3],type='l',lwd=1,ylim=c(-15,15),axes = F, xlab = "", ylab = "")
axis(side=4,at = seq(0,20,5),labels=seq(0,20,5), col=brewer.pal(3, 'Reds')[3],col.axis=brewer.pal(3, 'Reds')[3])
axis.POSIXct(1, at=seq(daterange[1], daterange[2], by="month"),labels=format(seq(daterange[1], daterange[2], by="month"),"%b"), format="%b")

mtext(text=expression('                                   T'['air']~~'[°C]'), side=4,line=3, col=brewer.pal(3, 'Reds')[3],cex=1.2)
abline(h=seq(0,20,5),v=seq(daterange[1], daterange[2], by="month"),col='grey')
par(new = TRUE)
plot(timXAxis,EBCoreInput$SWin,xaxt='n',col=brewer.pal(3, 'Greys')[2],type='l',lwd=1,ylim=c(0,2400),axes = F, xlab = "", ylab = "")
abline(h=seq(0,1000,200),v=seq(daterange[1], daterange[2], by="month"),col='grey')
axis(side=4,at = seq(0,1000,200),labels=seq(0,1000,200), col=brewer.pal(3, 'Greys')[3],col.axis=brewer.pal(3, 'Greys')[3])
mtext(text=expression('SW/LW'['in']~~'[W '~ m^{-2}~ ']                                      '), side=4,line=3, col=brewer.pal(3, 'Greys')[3],cex=1.2)

par(new = TRUE)
plot(timXAxis,EBCoreInput$LWin,xaxt='n',col=brewer.pal(3, 'Greys')[3],type='l',lwd=1,ylim=c(0,2400),axes = F, xlab = "", ylab = "")


par(new = TRUE)
plot(timXAxis,permelt_hig*100*24,xaxt='n',type='l',ylim=c(0,1.8),col=brewer.pal(3, 'Blues')[2],lwd=1.5,yaxt='n',ylab='')
axis(side=2,at = seq(0,1.8,0.1),labels=seq(0,1.8,0.1), col=brewer.pal(3, 'Blues')[3],col.axis=brewer.pal(3, 'Blues')[3])
mtext(text=expression('melt [cm '~ d^{-1}~ ']'), side=2,line=3, col=brewer.pal(3, 'Blues')[3],cex=1.2)
points(timXAxis,permelt_mid*100*24,type='l',col=brewer.pal(3, 'Blues')[3],lwd=1.5)
points(timXAxis,permelt_low*100*24,type='l',col=brewer.pal(3, 'Blues')[2],lwd=0.5)
dev.off()


focalTime <- seq(500,1000,1)
png(file=path_figs&'\\climateforcingmelt_focus.png', res = 160,width=1300,height=900)
par(mar=c(3,7,2,5),cex.lab=1.2,cex.axis=1.2)
layout(matrix(c(1,1), nrow = 1, ncol = 1, byrow = FALSE))
plot(timXAxis[focalTime],EBCoreInput$T_a[focalTime]-273.15,col=brewer.pal(3, 'Reds')[3],type='l',lwd=1,ylim=c(-15,15),axes = F, xlab = "", ylab = "")
axis(side=4,at = seq(0,20,5),labels=seq(0,20,5), col=brewer.pal(3, 'Reds')[3],col.axis=brewer.pal(3, 'Reds')[3])
axis.POSIXct(1, at=seq(daterange[1], daterange[2], by="week"),labels=format(seq(daterange[1], daterange[2], by="week"),"%d%b"), format="%d%b")

mtext(text=expression('                                   T'['air']~~'[°C]'), side=4,line=3, col=brewer.pal(3, 'Reds')[3],cex=1.2)
abline(h=seq(0,20,5),v=seq(daterange[1], daterange[2], by="week"),col='grey')

par(new = TRUE)
plot(timXAxis[focalTime],EBCoreInput$SWin[focalTime],xaxt='n',col=brewer.pal(3, 'Greys')[2],type='l',lwd=1,ylim=c(0,2400),axes = F, xlab = "", ylab = "")
axis(side=4,at = seq(0,1000,200),labels=seq(0,1000,200), col=brewer.pal(3, 'Greys')[3],col.axis=brewer.pal(3, 'Greys')[3])
mtext(text=expression('SW/LW'['in']~~'[W '~ m^{-2}~ ']                                      '), side=4,line=3, col=brewer.pal(3, 'Greys')[3],cex=1.2)
abline(h=seq(0,1000,200),v=seq(daterange[1], daterange[2], by="week"),col='grey')

par(new = TRUE)
plot(timXAxis[focalTime],EBCoreInput$LWin[focalTime],xaxt='n',col=brewer.pal(3, 'Greys')[3],type='l',lwd=1,ylim=c(0,2400),axes = F, xlab = "", ylab = "")


par(new = TRUE)
plot(timXAxis[focalTime],permelt_hig[focalTime]*100*24,xaxt='n',type='l',ylim=c(0,1.8),col=brewer.pal(3, 'Blues')[2],lwd=1.5,yaxt='n',ylab='')
axis(side=2,at = seq(0,1.8,0.1),labels=seq(0,1.8,0.1), col=brewer.pal(3, 'Blues')[3],col.axis=brewer.pal(3, 'Blues')[3])
mtext(text=expression('melt [cm '~ d^{-1}~ ']'), side=2,line=3, col=brewer.pal(3, 'Blues')[3],cex=1.2)

points(timXAxis[focalTime],permelt_mid[focalTime]*100*24,type='l',col=brewer.pal(3, 'Blues')[3],lwd=1.5)
points(timXAxis[focalTime],permelt_low[focalTime]*100*24,type='l',col=brewer.pal(3, 'Blues')[2],lwd=0.5)
#legend('bottomright',c('melt ~ 45 cm', 'melt ~ 115 cm','melt ~ 163 cm','air temperature', 'shortwave radiation'),lty=1,col=c(brewer.pal(3, 'Blues')[3],brewer.pal(3, 'Blues')[2],brewer.pal(3, 'Blues')[3],brewer.pal(3, 'Reds')[3],brewer.pal(3, 'Greys')[3]),bty='n')
dev.off()

diumeltlow <- diuCyc(percmelt[,2]*100*24,EBCoreInput$timeline_str)
diumeltmid <- diuCyc(percmelt[,6]*100*24,EBCoreInput$timeline_str)
diumelthig <- diuCyc(percmelt[,10]*100*24,EBCoreInput$timeline_str)

diuSW <- diuCyc(percSW[,6],EBCoreInput$timeline_str)
diuTa <- diuCyc(percTa[,6]-273.15,EBCoreInput$timeline_str)

png(file=path_figs&'\\climateforcingmelt_diurnal.png', res = 160,width=1800,height=900)
par(mar=c(3,7,2,5),cex.lab=1.2,cex.axis=1.2)
layout(matrix(c(1,1), nrow = 1, ncol = 1, byrow = FALSE))
plot(seq(1,24,1),diuTa[,2],col=brewer.pal(3, 'Reds')[3],type='l',lwd=1,ylim=c(-15,15),axes = F, xlab = "", ylab = "")
axis(side=4,at = seq(0,20,5),labels=seq(0,20,5), col=brewer.pal(3, 'Reds')[3],col.axis=brewer.pal(3, 'Reds')[3])
axis(1, at=seq(6, 24, 6),labels=seq(6,24,6))

mtext(text=expression('                                   T'['air']~~'[°C]'), side=4,line=3, col=brewer.pal(3, 'Reds')[3],cex=1.2)
par(new = TRUE)
plot(seq(1,24,1),diuSW[,2],xaxt='n',col=brewer.pal(3, 'Greys')[2],type='l',lwd=1,ylim=c(0,2200),axes = F, xlab = "", ylab = "")
axis(side=4,at = seq(0,1000,200),labels=seq(0,1000,200), col=brewer.pal(3, 'Greys')[3],col.axis=brewer.pal(3, 'Greys')[3])
mtext(text=expression('SW'['in']~~'[W '~ m^{-2}~ ']                                      '), side=4,line=3, col=brewer.pal(3, 'Greys')[2],cex=1.2)
par(new = TRUE)
plot(seq(1,24,1),diumeltlow[,2],xaxt='n',type='l',ylim=c(0,1.8),col=brewer.pal(3, 'Blues')[2],lwd=2,yaxt='n',ylab='')
axis(side=2,at = seq(0,1.8,0.1),labels=seq(0,1.8,0.1), col=brewer.pal(3, 'Blues')[3],col.axis=brewer.pal(3, 'Blues')[3])
mtext(text=expression('melt [cm '~ d^{-1}~ ']'), side=2,line=3, col=brewer.pal(3, 'Blues')[3],cex=1.2)
points(seq(1,24,1),diumeltmid[,2],type='l',col=brewer.pal(3, 'Blues')[3],lwd=2)
points(seq(1,24,1),diumelthig[,2],type='l',col=brewer.pal(3, 'Blues')[2],lwd=2)
dev.off()

#########

# TI model setup
timelag <- 2
TAir_all <- vector()
SW_all <- vector()
melt_all <- tothourlymelt
depth_all <- debThick[]

for(rF in 1:length(ModelRunSeries)){
  TI_TAirSeries <- c(seq(1,timelag,1)*0+visDataArray[[rF]]$V5[1:timelag]-273.15,visDataArray[[rF]]$V5-273.15)[1:length(visDataArray[[rF]]$V5)] # shifted temperature time series
  TAir_all <- c(TAir_all,TI_TAirSeries)
  ETI_SWinSeries <- c(seq(1,timelag,1)*0+visDataArray[[rF]]$V6[1:timelag],visDataArray[[rF]]$V6)[1:length(visDataArray[[rF]]$V6)]
  SW_all <- c(SW_all,ETI_SWinSeries)
  
  melt_all <- c(melt_all,visDataArray[[rF]]$V4)
  actualDepth <- mean(d_sample[visDataArray[[rF]]$V13],na.rm=T)
  depth_all <- c(depth_all,actualDepth + visDataArray[[rF]]$V4*0)
}

# melt runs to exclude because best run is zero melt over the whole period.
which(tothourlymelt[3769,]==0)

# get lag between T_air/SW and melt by finding peaks for each day
lagMatrix <- matrix(NA,nrow=dim(tothourlymelt)[2],ncol=3)
for(lagF in 1:dim(tothourlymelt)[2]){
if(!is.na(mean(tothourlymelt[,lagF],na.rm=T))){
  xyccf <- ccf(tothourlymelt[,lagF],actualTa[,lagF],lag.max=50)

  lag_Ta <- which(xyccf$acf[xyccf$lag>0]==max(xyccf$acf[xyccf$lag>0]))

  xyccf <- ccf(tothourlymelt[,lagF],actualSW[,lagF],lag.max=50)
  
  lag_SW <- which(xyccf$acf[xyccf$lag>0]==max(xyccf$acf[xyccf$lag>0]))
  
lagMatrix[lagF,1] <- debThick[ModelRunSeries[lagF]]
lagMatrix[lagF,2] <- lag_Ta
lagMatrix[lagF,3] <- lag_SW
}
  else{lagMatrix[lagF,] <- c(NA,NA,NA)}
}

lagMatrix[which(tothourlymelt[3769,]==0),] <- NA
(a * debThick^b) * T_air [t - c*debThick]

f <- function(x) sum((melt_all - x[1] * depth_all^x[2] * TAir_all*depth_all)^2,na.rm=T) / sum((melt_all - mean(x[1] * depth_all^x[2] * TAir_all*depth_all,na.rm=T))^2,na.rm=T)
DTI_p <- optim(c(0.001,0.1),f)

f <- function(x) sum((melt_all - (x[1] * debThMod^x[2] * TI_TAirSeries*debThMod + x[3] * (1 - EBCoreInput$albedo_d) * ETI_SWinSeries))^2,na.rm=T) / sum((meltMod - mean(x[1] * debThMod^x[2] * TI_TAirSeries*debThMod + x[3] * (1 - EBCoreInput$albedo_d) * ETI_SWinSeries))^2,na.rm=T)
DETI_p <- optim(c(0,0,0),f)


DTI_melt <- DTI_p$par[1] * depth_all^DTI_p$par[2] * TAir_all*depth_all
DTI_melt[DTI_melt<=0] <- 0

DETI_melt <- DETI_p$par[1] * depth_all^DETI_p$par[2] * TI_TAirSeries * depth_all + DETI_p$par[3] * (1 - EBCoreInput$albedo_d[1]) * SW_all
DETI_melt[DETI_melt<=0] <- 0

plot(cumsum(melt_all))
points(cumsum(DTI_melt),col='red')
points(cumsum(DETI_melt),col='green')

plot(DTI_melt[2000:2500]*1000,type='l',ylim=c(0,4),ylab = 'melt [mm h-1]')
points(melt_all[2000:2500]*1000,type='l',col='red')
#points(DETI_melt[2000:2500]*1000,type='l',col='green')

TI_daily <- TSaggregate(DTI_melt,rep(EBCoreInput$timeline_str,234),3600,2013,'sum')
EB_daily <- TSaggregate(melt_all,rep(EBCoreInput$timeline_str,234),3600,2013,'sum')

TI_diu <- diuCyc(TI_melt*1000,EBCoreInput$timeline_str)
EB_diu <- diuCyc(out[[2]]$melt*1000,EBCoreInput$timeline_str)


plot(TI_daily[,2],EB_daily[,2], xlab ='TI model daily [mm d-1]',ylab ='EB model daily [mm d-1]')
abline(0,1)

plot(TI_diu[,2],type='l',ylim=c(0,2),ylab = 'melt [mm h-1]')
points(EB_diu[,2],type='l',col='red')

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

png(file=path_figs&'\\ReconstructedThickness.png', res = 160,width=1800,height=600)
par(mar=c(4,5,2,2),cex.lab=1.5,cex.axis=1.5)
#par(mfrow=c(3,1))
layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = FALSE))
plot(d_raster,zlim=c(0,2.5),main='measured thickness [m]')
plot(glac_mask,add=T)
grid(NULL,NULL)
plot(d_raster_new,zlim=c(0,2.5),legend=FALSE,main='modelled thickness [m]')
plot(glac_mask,add=T)
grid(NULL,NULL)
boxplot(d_raster_new[]-d_raster[],ylim=c(-0.5,0.5),ylab =c('error [m]'),main = 'offset')
myjitter<-jitter(rep(1, length(which(!is.na(d_raster_new[]-d_raster[])))), amount=0.2)
points(myjitter, (d_raster_new[]-d_raster[])[which(!is.na(d_raster_new[]-d_raster[]))], pch=1, col=alpha('black',0.8) )
grid(NULL,NULL)
dev.off()

dThickordered <- which(!is.na(debThick[]))[1:length(out)]
png(file=path_figs&'\\ModelPerformance_indivPixel.png', res = 160,width=3600,height=1800)
ModelRunDays <- as.numeric(strftime(EBend, format = "%j"))-as.numeric(strftime(EBstart, format = "%j"))
par(mar=c(4,5,2,2),cex.lab=1.2,cex.axis=1.2)
#par(mfrow=c(3,1))
layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = FALSE))
resStat <- boxplot(out2[order(debThick[dThickordered])],ylim=c(0,4),xlab = 'pixels @ '&ModRes&'m resolution',ylab = 'mass loss [m]')
kl <- 1
for(i in order(debThick[dThickordered])){
  myjitter<-jitter(rep(kl, length(out2[order(dThickordered)][[i]])), amount=0.2)
  points(myjitter, out2[[i]], pch=20, col=alpha('black',0.2) )
  kl <- kl + 1
}
points(dhdt_res[dThickordered][order(debThick[dThickordered])],col='red',pch=2)
grid(NULL,NULL)
legend('bottomleft',c('mod','obs'),pch=c(NA,2),fill = c(NA,NA),border=c(gray.colors(1),NA),col=c('black','red'),bty='n')
#bplot <- barplot(debThick[dThickordered][order(debThick[dThickordered])],axisnames=T,space = 0, axes = FALSE,col="gray",xlab = 'pixels @ '&ModRes&'m resolution',ylab = 'thickness [m]')
bplot <- barplot(debOpt_mu,axisnames=T,space = 0, axes = FALSE,col="gray",xlab = 'pixels @ '&ModRes&'m resolution',ylab = 'thickness [m]',ylim=c(0,2))
error.bar(bplot,debOpt_mu, debOpt_sd)
at_tick <- seq_len(length(out) + 1)
axis(side = 1, at = at_tick - 1, labels = FALSE)
axis(side=1, at = seq_along(seq(1,length(out),1)) - 0.5, tick = FALSE,labels = seq(1,length(out),1))
axis(side = 2, at = seq(0,5,0.2), labels = seq(0,5,0.2))
box()
grid(NULL,NULL)
plot(debThick[dThickordered][order(debThick[dThickordered])],(resStat$stats[3,]-dhdt_res[dThickordered][order(debThick[dThickordered])]) / ModelRunDays * 100,ylab = 'error [cm/day]',xlab = 'thickness [m]')
points(debThick[dThickordered][order(debThick[dThickordered])],(optdhdt[order(debThick[dThickordered])]-dhdt_res[dThickordered][order(debThick[dThickordered])]) / ModelRunDays * 100,ylab = 'error [cm/day]',xlab = 'thickness [m]',col='red')

grid(NULL,NULL)
plot(debThick[dThickordered][order(debThick[dThickordered])],dhdt_res[dThickordered][order(debThick[dThickordered])] / ModelRunDays * 100,ylab = 'mass loss measured [cm/day]',xlab = 'thickness [m]')
grid(NULL,NULL)
dev.off()

save(out, file=path_data&'//temp_ModOutput')
load(path_data&'//temp_ModOutput')
###########
#Test Model Plots
allPix <- lapply(out, max);         # Total melt after complete model period
maxMelt <- DEM_dataframe[!is.na(DEM_dataframe[,2]),1]*0
maxMelt <- allPix

Model_Results_Raster <- DEM_model * 0
Model_Results_Raster[!is.na(Model_Results_Raster[])] <- unlist(maxMelt)   # make raster with melt values

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
