################################################################################
# D2EB - mass loss with TI model
# 
# d2DETI.R
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
dhdt_map_emer <- '2013Lirung_dhdtEmergence.tif'           # map of emergence velocity, if not available set to 'NA'
deb_thick_map <- 'modelledThickness_OestremModel.tif'     # debris thickness map if available; if not available set to 'NA'
pixID <- 'MeasuredThickness_1m.tif'                # Location where thickness measurements exist and initial runs can be executed

debrisview_domain <- 'debrisView_2013.tif'               # debris view raster for topographic shading
cliffFileName_out <- '2013octcliffsRasterized.tif'
pondFileName_out <- '2013octpondsRasterized.tif'

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
xStart <- as.numeric(as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S"))
xEnd <- as.numeric(as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S"))


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
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\mOmdTI.R")
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\TSStatistics.R")
# ========================
#  % PARAMETERS/VARIABLES (may change across sites; EDIT ACCORDING TO FIELDSITE)
# ========================

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
cutoff <- (xEnd - xStart)/3600/24 * 0.03
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

DomDim <- dim(dhdt_res)
CellNo <- DomDim[1] * DomDim[2]

if(station_loc=='debris'){
  # get raster ID of AWS on DEM raster
  AWSloc <- data.frame(lon=lon, lat=lat)
  coordinates(AWSloc) <- c("lon", "lat")
  proj4string(AWSloc) <-  CRS("+init=epsg:4326")
  coord <- spTransform(AWSloc,projec)
  ID_AWS <- raster::extract(DEM_model,coord,cellnumbers=T)[1]
}


# off glacier T data

Kyanjing_TData <- read.csv("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\KyanjingTData_formassloss.csv",header = T)       # Lapse Rate from Heynen et al. 2016

xStart <- as.numeric(as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S"))
xEnd <- as.numeric(as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S"))

idIN <- which(Kyanjing_TData[,2]==xStart)
idOUT <- which(Kyanjing_TData[,2]==xEnd)

# TI model parameters
TF1 <- 0.029
TF2 <- -0.919

minElevLir <- 4000    # Elevation of Lirung Glacier Snout
ElevKya <- 3862

ID_pix <- which(!is.na(dhdt_res[]))
meltCum <- dhdt_res*NA
for(i in 1:length(ID_pix)){
  
  shiftTim <- seq(idIN - round(17.7 *debThick[ID_pix[i]]), idOUT - round(17.7 *debThick[ID_pix[i]]),1)
  
  lapsedT <- Kyanjing_TData[shiftTim,3] +
    (minElevLir- ElevKya) * Kyanjing_TData[shiftTim,5] +
    (DEM_model[ID_pix[i]] - minElevLir) * Kyanjing_TData[shiftTim,6]
    
  meltCum[ID_pix[i]] <- cumsum(TF1 * debThick[ID_pix[i]]^TF2 * lapsedT)[length(lapsedT)] / 0.915 / 1000
  
}

coords_extent <- coordinates(DEM_model)

binned_matrix_model <- as.data.frame(cbind(coords_extent[which(!is.na(dhdt_res)[]),2],DEM_model[!is.na(dhdt_res)],meltCum[!is.na(meltCum)]))

binnedmodel <- binned_matrix_model %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 2)

binned_matrix_meas <- as.data.frame(cbind(coords_extent[which(!is.na(dhdt_res)[]),2],DEM_model[!is.na(dhdt_res)],dhdt_res[!is.na(dhdt_res)]))

binnedmeas <- binned_matrix_meas %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 3)

binnedresults <- rbind(binnedmodel,binnedmeas)


binnedresults %>%
  ggplot(aes(y = bin_dist, x = V3,fill=factor(type),colour = factor(type))) +
  geom_boxplot(outlier.colour = NULL) + 
  theme_bw() + 
  scale_y_discrete(breaks=seq(31233,31247,1),labels=seq(100,1500,100),position = "right") +
  labs(x = expression(paste("mass loss [m]")), y = "distance to terminus [m]") +
  theme(legend.position="none")

StatsTI_2013 <- TSStatistics(dhdt_res[],meltCum[])

dhdt_2013 <- dhdt_res
modmelt_2013 <- meltCum
days_2013 <- (xEnd - xStart)/3600/24


########################
# 2015 to 2016 data

glac_outline <- 'Lirung_2015.shp'                   # Outline of glacier tongue to be investigated
dhdt_map <- '2015-2016Lirung_dhdt_total.tif'                  # dh/dt map for calibration/validation
dhdt_uncertainty <- sqrt(0.25^2+0.25^2)
cl_thres <- 25                                            # threshold from which pixel is classified as cliff (%, relative area cover)
po_thres <- 25                                            # threshold from which pixel is classified as pond (%, relative area cover)
ModRes <- 10                                              # Spatial Resolution of Model [m]
DEM_domain <- '20151018_lirung_dem_20cm.tif'              # DEM for model domain (for lapsing)

############
# optional additional data
############
deb_thick_map <- 'modelledThickness_OestremModel.tif'     # debris thickness map if available; if not available set to 'NA'
pixID <- 'MeasuredThickness_1m.tif'                # Location where thickness measurements exist and initial runs can be executed

cliffFileName_out <- '2015octcliffsRasterized.tif'
pondFileName_out <- '2015octpondsRasterized.tif'


location <- 'Lirung'                                      # name of model location for file outputs
lat <- 28.23259709                                        # Latitude (decimal degrees)
lon <- 85.5621322                                         # Longitude (decimal degrees)

# Define whether station is located on a clean ice ('clean'), debris-covered ('debris') glacier or off-glacier ('off')
station_loc <- 'debris'

# set start and end date of calculation of EB
EBstart <- '2015/10/18 00:00:00'
EBend <- '2016/04/30 00:00:00'
xStart <- as.numeric(as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S"))
xEnd <- as.numeric(as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S"))

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
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\mOmdTI.R")
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\TSStatistics.R")
# ========================
#  % PARAMETERS/VARIABLES (may change across sites; EDIT ACCORDING TO FIELDSITE)
# ========================

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
cutoff <- (xEnd - xStart)/3600/24 * 0.03
dhdt[dhdt>cutoff] <- NA
dhdt_res <- raster::resample(dhdt, raster_domain,method='bilinear')
dhdt_res <- mask(dhdt_res,glac_mask)      # Available calib/valid data



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

DomDim <- dim(dhdt_res)
CellNo <- DomDim[1] * DomDim[2]

if(station_loc=='debris'){
  # get raster ID of AWS on DEM raster
  AWSloc <- data.frame(lon=lon, lat=lat)
  coordinates(AWSloc) <- c("lon", "lat")
  proj4string(AWSloc) <-  CRS("+init=epsg:4326")
  coord <- spTransform(AWSloc,projec)
  ID_AWS <- raster::extract(DEM_model,coord,cellnumbers=T)[1]
}


# off glacier T data

Kyanjing_TData <- read.csv("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\KyanjingTData_formassloss.csv",header = T)       # Lapse Rate from Heynen et al. 2016

xStart <- as.numeric(as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S"))
xEnd <- as.numeric(as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S"))

idIN <- which(Kyanjing_TData[,2]==xStart)
idOUT <- which(Kyanjing_TData[,2]==xEnd)

# TI model parameters
TF1 <- 0.029
TF2 <- -0.919

minElevLir <- 4000    # Elevation of Lirung Glacier Snout
ElevKya <- 3862

ID_pix <- which(!is.na(debThick[]))
meltCum <- dhdt_res*NA
for(i in 1:length(ID_pix)){
  
  shiftTim <- seq(idIN - round(17.7 *debThick[ID_pix[i]]), idOUT - round(17.7 *debThick[ID_pix[i]]),1)
  
  lapsedT <- Kyanjing_TData[shiftTim,3] +
    (minElevLir- ElevKya) * Kyanjing_TData[shiftTim,5] +
    (DEM_model[ID_pix[i]] - minElevLir) * Kyanjing_TData[shiftTim,6]
  
  meltVec <- TF1 * debThick[ID_pix[i]]^TF2 * lapsedT
  meltVec[lapsedT<=0] <- 0 
  meltCum[ID_pix[i]] <- cumsum(meltVec)[length(lapsedT)] / 0.915 / 1000

}

coords_extent <- coordinates(DEM_model)

binned_matrix_model <- as.data.frame(cbind(coords_extent[which(!is.na(dhdt_res)[]),2],DEM_model[!is.na(dhdt_res)],meltCum[!is.na(meltCum)]))

binnedmodel <- binned_matrix_model %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 2)

binned_matrix_meas <- as.data.frame(cbind(coords_extent[which(!is.na(dhdt_res)[]),2],DEM_model[!is.na(dhdt_res)],dhdt_res[!is.na(dhdt_res)]))

binnedmeas <- binned_matrix_meas %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 3)

binnedresults <- rbind(binnedmodel,binnedmeas)


binnedresults %>%
  ggplot(aes(y = bin_dist, x = V3,fill=factor(type),colour = factor(type))) +
  geom_boxplot(outlier.colour = NULL) + 
  theme_bw() + 
  scale_y_discrete(breaks=seq(31233,31247,1),labels=seq(100,1500,100),position = "right") +
  labs(x = expression(paste("mass loss [m]")), y = "distance to terminus [m]") +
  theme(legend.position="none")

StatsTI <- TSStatistics(dhdt_res[],meltCum[])

# Draw ROI
r1NaM <- is.na(as.matrix(dhdt_res))
colNotNA <- which(colSums(r1NaM) != nrow(dhdt_res))
rowNotNA <- which(rowSums(r1NaM) != ncol(dhdt_res))
r3Extent <- extent(dhdt_res, rowNotNA[1], rowNotNA[length(rowNotNA)],
                   colNotNA[1], colNotNA[length(colNotNA)])

StatsTI_201516 <- TSStatistics(dhdt_res[],meltCum[])

dhdt_201516 <- dhdt_res
modmelt_201516 <- meltCum
days_201516 <- (xEnd - xStart)/3600/24



########################
# 2016 data

glac_outline <- 'Lirung_2015.shp'                   # Outline of glacier tongue to be investigated
dhdt_map <- '2016Lirung_dhdt_total.tif'                  # dh/dt map for calibration/validation
dhdt_uncertainty <- sqrt(0.25^2+0.25^2)
cl_thres <- 25                                            # threshold from which pixel is classified as cliff (%, relative area cover)
po_thres <- 25                                            # threshold from which pixel is classified as pond (%, relative area cover)
ModRes <- 10                                              # Spatial Resolution of Model [m]
DEM_domain <- '20160430_lirung_dem_20cm.tif'              # DEM for model domain (for lapsing)

############
# optional additional data
############
deb_thick_map <- 'modelledThickness_OestremModel.tif'     # debris thickness map if available; if not available set to 'NA'
pixID <- 'MeasuredThickness_1m.tif'                # Location where thickness measurements exist and initial runs can be executed

cliffFileName_out <- '2016aprcliffsRasterized.tif'
pondFileName_out <- '2016aprpondsRasterized.tif'


location <- 'Lirung'                                      # name of model location for file outputs
lat <- 28.23259709                                        # Latitude (decimal degrees)
lon <- 85.5621322                                         # Longitude (decimal degrees)

# Define whether station is located on a clean ice ('clean'), debris-covered ('debris') glacier or off-glacier ('off')
station_loc <- 'debris'

# set start and end date of calculation of EB
EBstart <- '2016/04/30 00:00:00'
EBend <- '2016/10/06 00:00:00'
xStart <- as.numeric(as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S"))
xEnd <- as.numeric(as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S"))

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
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\mOmdTI.R")
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\TSStatistics.R")
# ========================
#  % PARAMETERS/VARIABLES (may change across sites; EDIT ACCORDING TO FIELDSITE)
# ========================

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
cutoff <- (xEnd - xStart)/3600/24 * 0.03
dhdt[dhdt>cutoff] <- NA
dhdt_res <- raster::resample(dhdt, raster_domain,method='bilinear')
dhdt_res <- mask(dhdt_res,glac_mask)      # Available calib/valid data



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

DomDim <- dim(dhdt_res)
CellNo <- DomDim[1] * DomDim[2]

if(station_loc=='debris'){
  # get raster ID of AWS on DEM raster
  AWSloc <- data.frame(lon=lon, lat=lat)
  coordinates(AWSloc) <- c("lon", "lat")
  proj4string(AWSloc) <-  CRS("+init=epsg:4326")
  coord <- spTransform(AWSloc,projec)
  ID_AWS <- raster::extract(DEM_model,coord,cellnumbers=T)[1]
}


# off glacier T data

Kyanjing_TData <- read.csv("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\KyanjingTData_formassloss.csv",header = T)       # Lapse Rate from Heynen et al. 2016

xStart <- as.numeric(as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S"))
xEnd <- as.numeric(as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S"))

idIN <- which(Kyanjing_TData[,2]==xStart)
idOUT <- which(Kyanjing_TData[,2]==xEnd)

# TI model parameters
TF1 <- 0.029
TF2 <- -0.919

minElevLir <- 4000    # Elevation of Lirung Glacier Snout
ElevKya <- 3862

ID_pix <- which(!is.na(debThick[]))
meltCum <- dhdt_res*NA
for(i in 1:length(ID_pix)){
  
  shiftTim <- seq(idIN - round(17.7 *debThick[ID_pix[i]]), idOUT - round(17.7 *debThick[ID_pix[i]]),1)
  
  lapsedT <- Kyanjing_TData[shiftTim,3] +
    (minElevLir- ElevKya) * Kyanjing_TData[shiftTim,5] +
    (DEM_model[ID_pix[i]] - minElevLir) * Kyanjing_TData[shiftTim,6]
  
  meltVec <- TF1 * debThick[ID_pix[i]]^TF2 * lapsedT
  meltVec[lapsedT<=0] <- 0 
  meltCum[ID_pix[i]] <- cumsum(meltVec)[length(lapsedT)] / 0.915 / 1000
  
}

coords_extent <- coordinates(DEM_model)

binned_matrix_model <- as.data.frame(cbind(coords_extent[which(!is.na(dhdt_res)[]),2],DEM_model[!is.na(dhdt_res)],meltCum[!is.na(meltCum)]))

binnedmodel <- binned_matrix_model %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 2)

binned_matrix_meas <- as.data.frame(cbind(coords_extent[which(!is.na(dhdt_res)[]),2],DEM_model[!is.na(dhdt_res)],dhdt_res[!is.na(dhdt_res)]))

binnedmeas <- binned_matrix_meas %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 3)

binnedresults <- rbind(binnedmodel,binnedmeas)


binnedresults %>%
  ggplot(aes(y = bin_dist, x = V3,fill=factor(type),colour = factor(type))) +
  geom_boxplot(outlier.colour = NULL) + 
  theme_bw() + 
  scale_y_discrete(breaks=seq(31233,31247,1),labels=seq(100,1500,100),position = "right") +
  labs(x = expression(paste("mass loss [m]")), y = "distance to terminus [m]") +
  theme(legend.position="none")

StatsTI <- TSStatistics(dhdt_res[],meltCum[])

StatsTI_2016 <- TSStatistics(dhdt_res[],meltCum[])

dhdt_2016 <- dhdt_res
modmelt_2016 <- meltCum
days_2016 <- (xEnd - xStart)/3600/24




########################
# 2017 - 2018 data

glac_outline <- 'Lirung_2015.shp'                   # Outline of glacier tongue to be investigated
dhdt_map <- '2017-2018Lirung_dhdt_total.tif'                  # dh/dt map for calibration/validation
dhdt_uncertainty <- sqrt(0.25^2+0.25^2)
cl_thres <- 25                                            # threshold from which pixel is classified as cliff (%, relative area cover)
po_thres <- 25                                            # threshold from which pixel is classified as pond (%, relative area cover)
ModRes <- 10                                              # Spatial Resolution of Model [m]
DEM_domain <- '20171019_lirung_dem_20cm.tif'              # DEM for model domain (for lapsing)

############
# optional additional data
############
deb_thick_map <- 'modelledThickness_OestremModel.tif'     # debris thickness map if available; if not available set to 'NA'
pixID <- 'MeasuredThickness_1m.tif'                # Location where thickness measurements exist and initial runs can be executed

cliffFileName_out <- '2017octcliffsRasterized.tif'
pondFileName_out <- '2017octpondsRasterized.tif'


location <- 'Lirung'                                      # name of model location for file outputs
lat <- 28.23259709                                        # Latitude (decimal degrees)
lon <- 85.5621322                                         # Longitude (decimal degrees)

# Define whether station is located on a clean ice ('clean'), debris-covered ('debris') glacier or off-glacier ('off')
station_loc <- 'debris'

# set start and end date of calculation of EB
EBstart <- '2017/10/19 00:00:00'
EBend <- '2018/04/28 00:00:00'
xStart <- as.numeric(as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S"))
xEnd <- as.numeric(as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S"))

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
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\mOmdTI.R")
source("F:\\PhD\\BasicCode\\timeseriesAnalysis\\TSStatistics.R")
# ========================
#  % PARAMETERS/VARIABLES (may change across sites; EDIT ACCORDING TO FIELDSITE)
# ========================

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
cutoff <- (xEnd - xStart)/3600/24 * 0.03
dhdt[dhdt>cutoff] <- NA
dhdt_res <- raster::resample(dhdt, raster_domain,method='bilinear')
dhdt_res <- mask(dhdt_res,glac_mask)      # Available calib/valid data



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

DomDim <- dim(dhdt_res)
CellNo <- DomDim[1] * DomDim[2]

if(station_loc=='debris'){
  # get raster ID of AWS on DEM raster
  AWSloc <- data.frame(lon=lon, lat=lat)
  coordinates(AWSloc) <- c("lon", "lat")
  proj4string(AWSloc) <-  CRS("+init=epsg:4326")
  coord <- spTransform(AWSloc,projec)
  ID_AWS <- raster::extract(DEM_model,coord,cellnumbers=T)[1]
}


# off glacier T data

Kyanjing_TData <- read.csv("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\KyanjingTData_formassloss.csv",header = T)       # Lapse Rate from Heynen et al. 2016

xStart <- as.numeric(as.POSIXct(EBstart, format="%Y/%m/%d  %H:%M:%S"))
xEnd <- as.numeric(as.POSIXct(EBend, format="%Y/%m/%d  %H:%M:%S"))

idIN <- which(Kyanjing_TData[,2]==xStart)
idOUT <- which(Kyanjing_TData[,2]==xEnd)
filler <- frollmean(Kyanjing_TData[idIN:idOUT,3],n=12,na.rm=T)
Kyanjing_TData[seq(idIN-72,idOUT+72,1),3][which(is.na(Kyanjing_TData[seq(idIN-72,idOUT+72,1),3]))] <- mean(Kyanjing_TData[idIN:idOUT,3],na.rm=T)

# TI model parameters
TF1 <- 0.029
TF2 <- -0.919

minElevLir <- 4000    # Elevation of Lirung Glacier Snout
ElevKya <- 3862

ID_pix <- which(!is.na(debThick[]))
meltCum <- dhdt_res*NA
for(i in 1:length(ID_pix)){
  
  shiftTim <- seq(idIN - round(17.7 *debThick[ID_pix[i]]), idOUT - round(17.7 *debThick[ID_pix[i]]),1)
  
  lapsedT <- Kyanjing_TData[shiftTim,3] +
    (minElevLir- ElevKya) * Kyanjing_TData[shiftTim,5] +
    (DEM_model[ID_pix[i]] - minElevLir) * Kyanjing_TData[shiftTim,6]
  
  meltVec <- TF1 * debThick[ID_pix[i]]^TF2 * lapsedT
  meltVec[lapsedT<=0] <- 0 
  meltCum[ID_pix[i]] <- cumsum(meltVec)[length(lapsedT)] / 0.916 / 1000
  
}

coords_extent <- coordinates(DEM_model)

binned_matrix_model <- as.data.frame(cbind(coords_extent[which(!is.na(dhdt_res)[]),2],DEM_model[!is.na(dhdt_res)],meltCum[!is.na(meltCum)]))

binnedmodel <- binned_matrix_model %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 2)

binned_matrix_meas <- as.data.frame(cbind(coords_extent[which(!is.na(dhdt_res)[]),2],DEM_model[!is.na(dhdt_res)],dhdt_res[!is.na(dhdt_res)]))

binnedmeas <- binned_matrix_meas %>% 
  mutate(bin_dist = factor(V1%/%100)) %>%
  mutate(type = 3)

binnedresults <- rbind(binnedmodel,binnedmeas)


binnedresults %>%
  ggplot(aes(y = bin_dist, x = V3,fill=factor(type),colour = factor(type))) +
  geom_boxplot(outlier.colour = NULL) + 
  theme_bw() + 
  scale_y_discrete(breaks=seq(31233,31247,1),labels=seq(100,1500,100),position = "right") +
  labs(x = expression(paste("mass loss [m]")), y = "distance to terminus [m]") +
  theme(legend.position="none")

StatsTI <- TSStatistics(dhdt_res[],meltCum[])

StatsTI_201718 <- TSStatistics(dhdt_res[],meltCum[])

dhdt_201718 <- dhdt_res
modmelt_201718 <- meltCum
days_201718 <- (xEnd - xStart)/3600/24


png(file=path_figs&'\\Distributed_Boxplots_TIModel_summer.png', res = 300,width=900,height=1800)
par(mar=c(3,5,1,2),mai = c(0.6, 0.1, 0.1, 0.1),cex.lab = 1.5, cex.axis = 1.5)
#par(pty="s")
layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = T))

plot.new()
boxplot(cbind(dhdt_2013[],modmelt_2013[]),ylim=c(0,5),xaxt='n',ylab = 'mass loss [m]',col=c('grey','white'))
mtext(text=expression('mass loss [m]'), side=2,line=3, col='black',cex=1)
boxplot(cbind(dhdt_2016[],modmelt_2016[]),ylim=c(0,5),xaxt='n',yaxt='n',col=c('grey','white'))

dev.off()

png(file=path_figs&'\\Distributed_Boxplots_TIModel_winter.png', res = 300,width=900,height=1800)
par(mar=c(3,5,1,2),mai = c(0.6, 0.1, 0.1, 0.1),cex.lab = 1.5, cex.axis = 1.5)
#par(pty="s")
layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = T))

plot.new()
boxplot(cbind(dhdt_201516[],modmelt_201516[]),ylim=c(0,2),xaxt='n',col=c('grey','white'))
boxplot(cbind(dhdt_201718[],modmelt_201718[]),ylim=c(0,2),yaxt='n',xaxt='n',col=c('grey','white'))
dev.off()

png(file=path_figs&'\\Distributed_Boxplots_TIModel_summer2013_violin.png', res = 300,width=900,height=1800)

par(mar=c(3,5,1,2),mai = c(0.6, 0.1, 0.1, 0.1),cex.lab = 2, cex.axis = 2)

factVio <- as.factor(c(rep(1, length(dhdt_2013)),rep(2,length(modmelt_2013[]))))
factCol <- as.factor(c(rep('UAV', length(dhdt_2013)),rep('dTI',length(modmelt_2013[]))))
voPl <- data.frame(factVio,c(dhdt_2013[]/days_2013*365,modmelt_2013[]/days_2013*365),factCol)
voPl <- voPl[-which(is.na(voPl[,2])),]
names(voPl)[1] <- 'type'
names(voPl)[2] <- 'ml'
p <- ggplot(voPl,aes(x=type,y=ml,fill=factCol)) + 
  geom_violin(trim=F) + geom_boxplot(width=0.1) + scale_fill_grey() +theme_classic() + ylim(0,11.5)
p + labs( y = expression('melt [m '~ yr^{-1}~ ']')) + scale_fill_grey(start=0.6,end=0.9) + theme(text = element_text(size=18)) +  theme(legend.position = "none") + theme(axis.title.x = element_blank(),axis.ticks = element_blank(),axis.text.x=element_blank())

dev.off()

png(file=path_figs&'\\Distributed_Boxplots_TIModel_summer2016_violin.png', res = 300,width=900,height=1800)

par(mar=c(3,5,1,2),mai = c(0.6, 0.1, 0.1, 0.1),cex.lab = 2, cex.axis = 2)

factVio <- as.factor(c(rep(1, length(dhdt_2016)),rep(2,length(modmelt_2016[]))))
factCol <- as.factor(c(rep('UAV', length(dhdt_2016)),rep('dTI',length(modmelt_2016[]))))
voPl <- data.frame(factVio,c(dhdt_2016[]/days_2016*365,modmelt_2016[]/days_2016*365),factCol)
voPl <- voPl[-which(is.na(voPl[,2])),]
names(voPl)[1] <- 'type'
names(voPl)[2] <- 'ml'
p <- ggplot(voPl,aes(x=type,y=ml,fill=factCol)) + 
  geom_violin(trim=F) + geom_boxplot(width=0.1) + scale_fill_grey() +theme_classic() + ylim(0,11.5)
p + labs( y = "") + scale_fill_grey(start=0.6,end=0.9) + theme(text = element_text(size=18)) +  theme(legend.position = "none") + theme(axis.title.x = element_blank(),axis.ticks = element_blank(),axis.text.x=element_blank())

dev.off()



png(file=path_figs&'\\Distributed_Boxplots_TIModel_winter201516_violin.png', res = 300,width=900,height=1800)

par(mar=c(3,5,1,2),mai = c(0.6, 0.1, 0.1, 0.1),cex.lab = 2, cex.axis = 2)

factVio <- as.factor(c(rep(1, length(dhdt_201516)),rep(2,length(modmelt_201516[]))))
factCol <- as.factor(c(rep('UAV', length(dhdt_201516)),rep('dTI',length(modmelt_201516[]))))
voPl <- data.frame(factVio,c(dhdt_201516[]/days_201516*365,modmelt_201516[]/days_201516*365),factCol)
voPl <- voPl[-which(is.na(voPl[,2])),]
names(voPl)[1] <- 'type'
names(voPl)[2] <- 'ml'
p <- ggplot(voPl,aes(x=type,y=ml,fill=factCol)) + 
  geom_violin(trim=F) + geom_boxplot(width=0.1) + scale_fill_grey() +theme_classic() + ylim(0,11.5)
p + labs( y = "") + scale_fill_grey(start=0.6,end=0.9) + theme(text = element_text(size=18)) +  theme(legend.position = "none") + theme(axis.title.x = element_blank(),axis.ticks = element_blank(),axis.text.x=element_blank())

dev.off()

png(file=path_figs&'\\Distributed_Boxplots_TIModel_winter201718_violin.png', res = 300,width=900,height=1800)

par(mar=c(3,5,1,2),mai = c(0.6, 0.1, 0.1, 0.1),cex.lab = 2, cex.axis = 2)

factVio <- as.factor(c(rep(1, length(dhdt_201718)),rep(2,length(modmelt_201718[]))))
factCol <- as.factor(c(rep('UAV', length(dhdt_201718)),rep('dTI',length(modmelt_201718[]))))
voPl <- data.frame(factVio,c(dhdt_201718[]/days_201718*365,modmelt_201718[]/days_201718*365),factCol)
voPl <- voPl[-which(is.na(voPl[,2])),]
names(voPl)[1] <- 'type'
names(voPl)[2] <- 'ml'
p <- ggplot(voPl,aes(x=type,y=ml,fill=factCol)) + 
  geom_violin(trim=F) + geom_boxplot(width=0.1) + scale_fill_grey() +theme_classic() + ylim(0,11.5)
p + labs( y = "") + scale_fill_grey(start=0.6,end=0.9) + theme(text = element_text(size=18)) +  theme(legend.position = "none") + theme(axis.title.x = element_blank(),axis.ticks = element_blank(),axis.text.x=element_blank())

dev.off()

