################################################################################
# D2EB - distributed debris Energy Balance Model
# 
# d2Deb.R
#
# ReadMe:
# Based on Reid and Brock, 2010
# Only edit within the blocks (a) PARAMETERS, (b) VARIABLES and (c) VISUALIZE
# (a)   - Paths:  (*) Make sure all paths are correctly pointed to the directory you work in. 
#                 (*) This code (d2Deb.R) has to be located in the path specified as 'path_code'
#                 All necessary subcodes need to be located in 'path_subcode'
#                 (*) The Path 'path_data' needs a subfolder called 'Outline_Features'
# (b)   - extra Files:(*) compulsory extra codes in 'path_subcode' are TSaggregate.R', 'TS_Check.R', 'Q_P.R', 'Q_LE.R', 'Q_H.R', 'Q_G.R', 'Lup.R', 'EBModCore.R', 'debristemp.R' and 'diuCyc.R'
#                     (*) If available you can provide cliff and pond outlines for the glacier surface.They need to be polygons saved as .shp files. 
#                         Each of the mapped polygons ideally has an ID, to define its properties (see Steiner et al. 2019). 
#                         For cliffs, ID=1 means no pond is associated to it and ID=2 signifies an adjacent pond
#                         For ponds, ID=7 signifies no cliff and ID=8 signifies a cliff.
#       - StakeName (change to specific Stake name, see Folder StakeData)
#       - AWS_loc (change to 'on' if you use on-glacier data and 'off' if you use off-glacier data)
#       - StationName (change according to the Meteo data file you use, see folder MeteoData)
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
p_load(rgdal,rgeos,maptools,raster,foreach,lubridate,compare,colorRamps,data.table,circular,parallel,snowfall,truncnorm,rlecuyer,forecast,rasterVis,R.utils,zoo,rlist)

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

paramSize <- 1000

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
AWSelev <- 4075.84                                           # Elevation of climate data station

location <- 'Lirung'                                      # name of model location for file outputs
lat <- 28.23259709                                        # Latitude (decimal degrees)
lon <- 85.5621322                                         # Longitude (decimal degrees)

# Define whether station is located on a clean ice ('clean'), debris-covered ('debris') glacier or off-glacier ('off')
station_loc <- 'debris'

# Specify height of T_a and ws sensor (assumed height is 2 m). Used to scale to standard height.
sh_Ta <- 2.00;
sh_ws <- 2.00;

# set start and end date of calculation of EB
EBstart <- '2013/05/18 00:00:00'
EBend <- '2013/10/22 00:00:00'

seasonsDyn <- 'on';                                       # choose variable seasons (1 dry, 1 wet, 1 dry)

monin <<- '06/15'                                         # define beginning of wet season
monout <<- '09/19'                                        # define end of wet season

winin <<- '1/1'                                           # define end of second dry season
winout <<- '02/28'                                        # define beginning of first dry season

albedoDyn <<- 'var';                                    # choose whether a constant debris albedo value ('const') or a time series is used ('var')

#StationName <- 'Lirung_'&toString(year);                 # on-glacier Station
#StationName_off <- 'Kyanjing_'&toString(year);           # off-glacier Station



################# 
# Paths for Extra Codes
file.sources = list.files(path_subcode, 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\EBModCore.R")
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
altitude <<- 4075;           # Altitude of measurement site (m)
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
epsilon_d <<- 0.94;           # Debris emissivity
epsilon_d <<- 1;
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
tsurfvec <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
Snet <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
Lnet <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
QH <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
QLE <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
refreeze <- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)

t002<- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
t004<- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
t006<- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
t01<- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
t03<- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)
t045<- matrix(nrow=length(EBCoreInput$timeline_str),ncol=nMC)





CDrange <<- 0.5;                   # Range for central difference calculation of the derivative of surface flux when calculating surface temperature T_s.
library(lubridate)


dfix <- 1       # 1 if measured debris thickness is used, 0 if thickness is determined from MC
sensrun <- 3    # 1 for sensitivity to debris properties, 2 climate properties, 3 all properties
#sensrun <- 4    # 4 Model run with final parameter results for complete domain

# define simple test function
CoreFunc <- function(i){

  # pixel specific climate data
  EBCoreInput$SWin <- SW_orig * vf_sky_actual[i]
  EBCoreInput$T_a <- Ta_orig
  EBCoreInput$RH_a <- RHa_orig
  if(is.null(climInput$PRES[IDstart:IDend])){
    p_0 <- 101325   # air pressure at level [Pa]
    T_0 <- 288.15   # Temperature at Sea Level
    g <- 9.81       # gravitational acceleration [m s-2]
    M_a <- 0.02896  # Molar Mass of Air [kg mol-1]
    R <- 8.31       # Universal Gas Constant
    L <- 0.0065     # Temperature Lapse Rate
    EBCoreInput$p_a <- p_0* ((1-(L*DEM_model[i]/T_0))^(g*M_a/(R*L))) + EBCoreInput$T_a * 0     # lapse Air Pressure 
  }  else{
    EBCoreInput$p_a <- gapFill(climInput$PRES[IDstart:IDend])          # Air pressure (kPa)  
  }
  
  # Lapse Station data to specific grid location
  if(station_loc == 'debris'){
    DEM_model[i] <- AWSelev
  EBCoreInput$T_a <- LR_TAIR_on * (DEM_model[i] - AWSelev) + EBCoreInput$T_a
  EBCoreInput$T_s_data <- AirSurfTemp(EBCoreInput$T_a-273.15,EBCoreInput$timeline_str,1,1)$v1
  }
  if(station_loc == 'off'){
  EBCoreInput$T_a <- LR_TAIR_off * (min(DEM_model[],na.rm=TRUE) - AWSelev) + LR_TAIR_on * (DEM_model[i] - min(DEM_model[],na.rm=TRUE)) + EBCoreInput$T_a  
  EBCoreInput$RH_a <- EBCoreInput$RH_a * LR_RH
  EBCoreInput$T_s_data <- AirSurfTemp(EBCoreInput$T_a-273.15,EBCoreInput$timeline_str,1,1)$v1
  EBCoreInput$LWin <- LWMod_HMA(EBCoreInput$T_a,EBCoreInput$RH_a,EBCoreInput$SWin,EBCoreInput$timeline_str)
  }
  
  # ========
  # MODEL CORE
  # ========

  for(mc in 1:nMC){
    if(sensrun==1){
    k_d_dry <<- k_dry_sample[mc]
    k_d_wet <<- k_wet_sample[mc]
    rho_d <<- rho_d_sample[mc]
    if(dfix==1){d <- debThick[i]} else if(dfix==0){d <- d_sample[mc]}
    
    }else if(sensrun==2){
    d <- debThick[i]
    z_0 <<- z0_d_sample[mc]
    EBCoreInput$albedo_d <<- alpha_d_sample[mc]
    EBCoreInput$T_a <- Ta_orig + tair_range[mc]
    EBCoreInput$u <- uorig * ws_range[mc]
    EBCoreInput$RH_a <- RHa_orig + rh_range[mc]
    EBCoreInput$RH_a[EBCoreInput$RH_a>100]<-100
    EBCoreInput$LWin <- LWinorig + LW_range[mc]
    
    }else if(sensrun==3){
      k_d_dry <<- k_dry_sample[mc]
      k_d_wet <<- k_wet_sample[mc]
      rho_d <<- rho_d_sample[mc]
      if(dfix==1){d <- debThick[i]} else if(dfix==0){d <- d_sample[mc]}
      z_0 <<- z0_d_sample[mc]
      EBCoreInput$albedo_d <- alpha_d_sample[mc]
      EBCoreInput$T_a <- Ta_orig + tair_range[mc]
      EBCoreInput$u <- uorig * ws_range[mc]
      EBCoreInput$RH_a <- RHa_orig + rh_range[mc]
      EBCoreInput$RH_a[EBCoreInput$RH_a>100]<-100
      EBCoreInput$LWin <- LWinorig + LW_range[mc]
      }

    #EBCoreInput$SWout <- EBCoreInput$SWin * albedo_d
    # Debris Layers
   

    N <<- floor(d/0.01);             # Number of layers used for debris temperature calculation (chosen to make each layer 1cm deep)
    if(N<10){N <<- 10}               # Set minimum number of layers = 5
    h <<- d/N;                       # Size of each layer in calculation

    EBModRes <-  EBModCore(EBCoreInput)
    #dhdt_param_mod[mc] <- max(cumsum(EBModRes[,dim(EBModRes)[2]-1]))
    meltCD[,mc] <- EBModRes[,dim(EBModRes)[2]-1]
    tairvec[,mc] <- EBCoreInput$T_a
    Snet[,mc] <- EBModRes[,7]
    Lnet[,mc] <- EBModRes[,8]
    QH[,mc] <- EBModRes[,10]
    QLE[,mc] <- EBModRes[,11]
    refreeze[,mc] <- EBModRes[,dim(EBModRes)[2]]
    tsurfvec[,mc] <-  EBModRes[,3]
    t002[,mc] <- EBModRes[,14]
    t004[,mc] <- EBModRes[,16]
    t006[,mc] <- EBModRes[,18]
    t01[,mc] <- EBModRes[,22]
    t03[,mc] <- EBModRes[,41]
    t045[,mc] <- EBModRes[,49]

  }

  #write(i, file = path_code&'\\Temp\\Run'&i&'.txt',append=T)
  # Find the Model run most closely matching the measurement (best 25 versions)
  #MCcutoff <- 100
  #ID_fit <- which(rank(abs(dhdt_param_mod - dhdt_res[i]), ties.method='min') <=  MCcutoff)
  
  outPut <- list()
  #outPut$melt <- meltCD[,ID_fit[1: MCcutoff]]
  #outPut$temp <- tairvec[,ID_fit[1: MCcutoff]]
  #outPut$Snet <- Snet[,ID_fit[1: MCcutoff]]
  #outPut$Lnet <- Lnet[,ID_fit[1: MCcutoff]]
  #outPut$QH <- QH[,ID_fit[1: MCcutoff]]
  #outPut$QLE <- QLE[,ID_fit[1: MCcutoff]]
  #outPut$Sin <- EBCoreInput$SWin
  #outPut$refreeze <- refreeze[,ID_fit[1: MCcutoff]]
  #outPut$dhdt <- dhdt_param_mod
  #outPut$ID <- ID_fit + dhdt_param_mod*0
  write.csv(meltCD[,], file = path_code&'\\Temp\\melt'&i&'.csv')
  write.csv(tairvec[,], file = path_code&'\\Temp\\temp'&i&'.csv')
  write.csv(Snet[,], file = path_code&'\\Temp\\Snet'&i&'.csv')
  write.csv(Lnet[,], file = path_code&'\\Temp\\Lnet'&i&'.csv')
  write.csv(QH[,], file = path_code&'\\Temp\\QH'&i&'.csv')
  write.csv(QLE[,], file = path_code&'\\Temp\\QLE'&i&'.csv')
  write.csv(refreeze[,], file = path_code&'\\Temp\\refreeze'&i&'.csv')
  write.csv(tsurfvec[,], file = path_code&'\\Temp\\tempsurf'&i&'.csv')
  write.csv(t002[,], file = path_code&'\\Temp\\temp002'&i&'.csv')
  write.csv(t004[,], file = path_code&'\\Temp\\temp004'&i&'.csv')
  write.csv(t006[,], file = path_code&'\\Temp\\temp006'&i&'.csv')
  write.csv(t01[,], file = path_code&'\\Temp\\temp01'&i&'.csv')
  write.csv(t03[,], file = path_code&'\\Temp\\temp03'&i&'.csv')
  write.csv(t045[,], file = path_code&'\\Temp\\temp045'&i&'.csv')
}

  
# number of workers (leave one to keep system more responsive)
nworkers <- detectCores()-4

# run entire function parallelized
sfInit(parallel=T,cpu=nworkers, type="SOCK")                                                # initialize cluster
loadlib <- lapply(c('raster','truncnorm'),                                                  # load required packages in cluster
                  function(x) sfLibrary(x, character.only=T))
sfExportAll()                                                                               # transfer all variables to clusters
sfClusterSetupRNG()                                                                         # set up random number generator

# identify locations of dhdt measurements that are NA or below the accuracy threshold

#ModelRunSeries <- DEM_dataframe[which(!is.na(dhdt_res[])&!is.na(ModPix[])),1];
ModelRunSeries <- c(ID_AWS)
tout.multi <- system.time(
  out <- sfLapply(ModelRunSeries, CoreFunc)
)
sfStop() 