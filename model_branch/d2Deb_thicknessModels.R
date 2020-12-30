################################################################################
# debHModels
# 
# d2Deb_thicknessModels.R
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
#install.packages('pacman')
library(pacman)
library(gstat)
p_load(rgdal,rgeos,maptools,raster,foreach,lubridate,compare,colorRamps,data.table,circular,parallel,snowfall,truncnorm,rlecuyer,forecast,rasterVis,R.utils,zoo,rlist,ggplot2)

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

# Paths PC2
#path <- 'C:\\Users\\Stein023\\DistributedEB'        
#path_code <- 'C:\\Users\\Stein023\\DistributedEB\\Code'
#path_subcode <- 'C:\\Users\\Stein023\\DistributedEB\\Code\\subCodes\\'
#path_figs <- 'C:\\Users\\Stein023\\DistributedEB\\Figures'
#path_data <- 'C:\\Users\\Stein023\\DistributedEB\\RawData'
#path_dhdt <- 'C:\\Users\\Stein023\\DistributedEB\\RawData\\dhdt'
#path_shade <- 'C:\\Users\\Stein023\\DistributedEB\\RawData\\localShading'

debThickness <- 'debris_thickness_mean_f1_f3_1000r.tif'   # Debris thickness raster
glac_outline <- 'Lirung_2015.shp'                   # Outline of glacier tongue to be investigated
dhdt_map <- '2013Lirung_dhdt_total.tif'                  # dh/dt map for calibration/validation
dhdt_uncertainty <- sqrt(0.25^2+0.2^25)
cl_thres <- 25                                            # threshold from which pixel is classified as cliff (%, relative area cover)
po_thres <- 25                                            # threshold from which pixel is classified as pond (%, relative area cover)
ModRes <- 10                                              # Spatial Resolution of Model [m]
DEM_domain <- '20160430_lirung_dem_20cm.tif'              # DEM for model domain (for lapsing)

############
# optional additional data
############
pond_outlines <- 'Ponds_May13.shp'                        # Outlines for pond features (write NULL if not applicable)
cliff_outlines <- 'Cliffs_May13.shp'                      # Outlines for cliff features (write NULL if not applicable)
pond_raster <- '2013maypondsRasterized.tif'                  # Rasterized pond outlines (use featureRasterize_d2EB.R)
cliff_raster <- '2013maycliffsRasterized.tif'                # Rasterized cliff outlines (use featureRasterize_d2EB.R)
dhdt_map_emer <- '2013Lirung_dhdtEmergence.tif'           # map of emergence velocity, if not available set to 'NA'




deb_thick_map <- 'GPRThickness_2015.tif'                  # debris thickness map if available; if not available set to 'NA'
pit_thick_map <- 'PitThickness_2016.tif'                  # debris thickness map if available; if not available set to 'NA'
debrisview_domain <- 'debrisView_2013.tif'                # debris view raster for topographic shading
cliffFileName_out <- '2013maycliffsRasterized.tif'
pondFileName_out <- '2013maypondsRasterized.tif'

# parameter ranges + debris thickness range
cond_debris_file <- path&'\\MCRanges\\debris_conductivity_n100.csv'             # debris conductivity range for Lirung Glacier
por_debris_file <- path&'\\MCRanges\\por_lognormal_n100.csv'                    # debris porosity range for Lirung Glacier
den_debris_file <- path&'\\MCRanges\\den_lognormal_n100.csv'                    # debris porosity range for Lirung Glacier
z0_debris_file <- path&'\\MCRanges\\z0_lognormal_n100.csv'                      # surface roughness range for Lirung Glacier
alpha_debris_file <- path&'\\MCRanges\\alpha_n100.csv'             # debris albedo range for Lirung Glacier
thickness_debris_file <- path&'\\MCRanges\\DebrisThickness_lognormal_n100.csv'  # debris thickness range for Lirung Glacier

# uncertainty range of climate variables
tair_file <- path&'\\MCRanges\\TAir_n100.csv'                                   # Air Temperature
rh_file <- path&'\\MCRanges\\RH_n100.csv'                                   # Air Temperature
ws_file <- path&'\\MCRanges\\ws_n100.csv'                                   # Wind speed
LW_file <- path&'\\MCRanges\\LW_n100.csv'                                   # incoming longwave radiation

#spDom <-'AWSLirung_onglacier_2012_Data.csv'
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

dhdt_res <- raster::resample(dhdt, raster_domain,method='bilinear')
dhdt_res <- mask(dhdt_res,glac_mask)      # Available calib/valid data

# Debris Thickness (GPR and Pits)
debThickVal <- NA
if(!is.na(deb_thick_map)){
  debThickVal <- raster(path_data&'\\DebrisThickness\\'&deb_thick_map)
  projection(debThickVal) <- projec
  debThickVal <- crop(debThickVal,dhdt_res)
  debThickVal <- raster::resample(debThickVal,dhdt_res,method='bilinear')
}

debThickVal_pits <- NA
if(!is.na(pit_thick_map)){
  debThickVal_pits <- raster(path_data&'\\DebrisThickness\\'&pit_thick_map)
  projection(debThickVal_pits) <- projec
  debThickVal_pits <- crop(debThickVal_pits,dhdt_res)
  debThickVal_pits <- raster::resample(debThickVal_pits,dhdt_res,method='bilinear')
}

# Read in Surface Features with different melt properties (need to be preprocessed with featureRasterize_d2EB.R)

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
DEM_model <- raster(path_data&'/DEMs/'&DEM_domain)
DEM_model <- raster::resample(DEM_model, raster_domain,method='bilinear')
DEM_model <- mask(DEM_model,raster_domain)
DEM_model <- mask(DEM_model,dhdt_res)

# resample to standard model resolution
raster_domain <- aggregate(raster_domain_masks, fact = floor(ModRes/res(raster_domain_masks)) , median)
dhdt_res[which(is.na(raster_domain_masks[]))] <- NA
#remove outliers from raster (mean annual melt rates exceeding 3cm/day are removed); also mass loss below 0.2m is removed to avoid outliers
cutoff <- (lubridate::yday(EBend)-lubridate::yday(EBstart)) * 0.03
dhdt_res[dhdt_res[]>cutoff] <- NA
dhdt_res <- aggregate(dhdt_res, fact = floor(ModRes/res(dhdt_res)) , median)
DEM_model <- aggregate(DEM_model, fact = floor(ModRes/res(DEM_model)) , median)
debThick <- aggregate(debThickVal, fact = floor(ModRes/res(debThickVal)) , median)
debThick_pits <- aggregate(debThickVal_pits, fact = floor(ModRes/res(debThickVal_pits)) , median)
debAll <-merge(debThick,debThick_pits)
writeRaster(debAll,path_data&'\\DebrisThickness\\measuredThickness_10mRaster.tif','GTiff',overwrite=T)



# Dataframe of DEM for later computations
DEM_dataframe <- cbind(seq(1,length(DEM_model[]),1),DEM_model[],coordinates(DEM_model))
#DEM_dataframe_new <- DEM_dataframe[which(!is.na(DEM_dataframe[,2])),]

DomDim <- dim(dhdt_res)
CellNo <- DomDim[1] * DomDim[2]



resultFiles <- list.files(path = path_code&'\\Temp\\', pattern = '.txt',full.names = T)

mod_dhdt_res <- dhdt_res * NA
dhdt_res2 <- dhdt_res
mod2_dhdt_res <- dhdt_res * NA
d_thick <- dhdt_res * NA
d_thick2 <- dhdt_res * NA
ID_pos <- which(dhdt_res[]>=0)
ID_neg <- which(dhdt_res[]<0)

d_sample <- read.csv(thickness_debris_file, header = T)$x
for(i in 1: length(resultFiles)){
  sp1 <- strsplit(resultFiles[i],'\\',fixed=TRUE)
  sp1 <- sp1[[1]][length(sp1[[1]])]
  sp2 <- strsplit(sp1,'cycle',fixed=TRUE)[[1]][2]
  sp3 <- strsplit(sp2,'ID',fixed=TRUE)[[1]]
  sp4 <- strsplit(sp3,'.',fixed=TRUE)[[2]][1]
  finID <- as.numeric(sp3[1]) + as.numeric(sp4[1]) - 1
  
  resData <- read.table(resultFiles[finID],skip = 1)
  
  IDmin <- which.min(abs(dhdt_res[DEM_dataframe[which(!is.na(dhdt_res[])),1][finID]] - resData$V12[resData$V13[1:25]]))
  
  mod_dhdt_res[DEM_dataframe[which(!is.na(dhdt_res[])),1][finID]] <- max(cumsum(resData$V4))
  mod2_dhdt_res[DEM_dataframe[which(!is.na(dhdt_res[])),1][finID]] <- resData$V12[resData$V13[1:25]][IDmin]
  d_thick[DEM_dataframe[which(!is.na(dhdt_res[])),1][finID]] <- d_sample[resData$V13[1:25]][IDmin]
  d_thick2[DEM_dataframe[which(!is.na(dhdt_res[])),1][finID]] <- mean(d_sample[resData$V13[1:25]])
  
  dhdt_res[DEM_dataframe[which(!is.na(dhdt_res[])),1][finID]]
  
  #DEM_dataframe[which(!is.na(dhdt_res[])),1]
  #dhdt_res[which(!is.na(dhdt_res[])),1]
}

mod_dhdt_res[ID_neg] <- NA
mod2_dhdt_res[ID_neg] <- NA
d_thick[ID_neg] <- NA
d_thick2[ID_neg] <- NA
dhdt_res2[ID_neg] <- NA

# Model Debris Thickness based on Oestrem Curve
# Load different model parameters

OCparams <- read.table('F:\\PhD\\Research\\EB_DCG\\OestremLiterature\\OestremParameters.csv',header=T,sep=',')

GPRDatID <- which(!is.na(debThick[]))
PitDatID <- which(!is.na(debThick_pits[]))





# Make dhdt map for all years, wherever data is available

dhdtList <- c('2013-2014Lirung_dhdtEmergence.tif','2014-2015Lirung_dhdtEmergence.tif','2015-2016Lirung_dhdtEmergence.tif','2016Lirung_dhdtEmergence.tif',
              '2016-2017Lirung_dhdtEmergence.tif','2017Lirung_dhdtEmergence.tif','2017-2018Lirung_dhdtEmergence.tif')

outName <- c('2013-2014Lirung_dhdt_nofeatures.tif','2014-2015Lirung_dhdt_nofeatures.tif','2015-2016Lirung_dhdt_nofeatures.tif','2016Lirung_dhdt_nofeatures.tif',
             '2016-2017Lirung_dhdt_nofeatures.tif','2017Lirung_dhdt_nofeatures.tif','2017-2018Lirung_dhdt_nofeatures.tif')

cliffList <- c('2013octcliffsRasterized.tif','2014maycliffsRasterized.tif','2015octcliffsRasterized.tif','2015octcliffsRasterized.tif','2016aprcliffsRasterized.tif','2016octcliffsRasterized.tif','2017aprcliffsRasterized.tif','2017octcliffsRasterized.tif','2018aprcliffsRasterized.tif')
pondList <- c('2013octpondsRasterized.tif','2014maypondsRasterized.tif','2015octpondsRasterized.tif','2015octpondsRasterized.tif','2016aprpondsRasterized.tif','2016octpondsRasterized.tif','2017aprpondsRasterized.tif','2017octpondsRasterized.tif','2018aprpondsRasterized.tif')

dhdt_nofeatures <- function(i,feat){
  
  # Validation/Calibration domain (in this case available dh/dt map, once with and without emergence)
  raster_domain <- raster()
  extent(raster_domain) <- extent(glac_mask)
  res(raster_domain) <- 1
  raster_domain[] <- 0
  projection(raster_domain)<-projec
  raster_domain <- mask(raster_domain,glac_mask)
  
  dhdt <- raster(path_data&'\\dhdt\\'&dhdtList[i])
  projection(dhdt)<-projec
  dhdt <- crop(dhdt,raster_domain)
  
  dhdt_res <- raster::resample(dhdt, raster_domain,method='bilinear')
  dhdt_res <- mask(dhdt_res,glac_mask)      # Available calib/valid data
  dhdt_res_withfeatures <- dhdt_res
  
  # Read in Surface Features with different melt properties (need to be preprocessed with featureRasterize_d2EB.R)
  cliffFileName_out <- cliffList[i]
  pondFileName_out <- pondList[i]
  
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
   # rat <- levels(visCliff)[[1]]
  #  rat[["pond"]] <- c("-pond","+pond")
   # levels(visCliff) <- rat
  }
  
  if(!is.null(pondStack)){
    raster_domain_masks[pondStack$indID>0] <- NA
    visPond <- as.factor(pondStack$typeID)
    #rat <- levels(visPond)[[1]]
    #rat[["pond"]] <- c("-cliff","+cliff")
    #levels(visPond) <- rat
  }
  
  # resample to standard model resolution
  raster_domain <- aggregate(raster_domain_masks, fact = floor(ModRes/res(raster_domain_masks)) , median)
  dhdt_res[which(is.na(raster_domain_masks[]))] <- NA
  #remove outliers from raster (mean annual melt rates exceeding 3cm/day are removed); also mass loss below 0.2m is removed to avoid outliers
  cutoff <- (lubridate::yday(EBend)-lubridate::yday(EBstart)) * 0.03
  dhdt_res[dhdt_res[]>cutoff] <- NA
  dhdt_res <- aggregate(dhdt_res, fact = floor(ModRes/res(dhdt_res)) , median)
  dhdt_res_withfeatures <- aggregate(dhdt_res_withfeatures, fact = floor(ModRes/res(dhdt_res_withfeatures)) , median) 
  writeRaster(dhdt_res,path_data&'\\dhdt\\'&outName[i],format='GTiff',overwrite=TRUE )
  if(feat == 1){return(dhdt_res_withfeatures)}
  if(feat == 2){return(dhdt_res)}
}

dhdt1314 <- dhdt_nofeatures(1,2)
dhdt1415 <- dhdt_nofeatures(2,2)
dhdt1516 <- dhdt_nofeatures(3,2)
dhdt16 <- dhdt_nofeatures(4,2)
dhdt16_features <- dhdt_nofeatures(4,1)
dhdt1617 <- dhdt_nofeatures(5,2)
dhdt17 <- dhdt_nofeatures(6,2)
dhdt1718 <- dhdt_nofeatures(7,2)

options(stringsAsFactors = FALSE)
flightDates <- read.csv('F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt\\UAVflightDates.csv',header = T)
library(fields) 
library(gridExtra)
library(automap)

debThickMod <- function(dhdtMap,dhdtMap_features,item,flightDates,OCparams,crit_melt,modelType,crit_melt_std){
  
  if(item==9){
    deltaDays <- 166
  }
  else{deltaDays <- as.numeric((as.Date(as.character(flightDates$EndDate[item]), format="%m/%d/%Y") - as.Date(as.character(flightDates$StartDate[item]), format="%m/%d/%Y")))
  }
  #browser()
  ID_neg <- which(dhdtMap[]<0)
  dhdtMap[ID_neg] <- NA
  dhdtMap[dhdtMap[]>quantile(dhdtMap[], c(.95),na.rm=T) ] <- NA
  dhdtMap[dhdtMap[]<quantile(dhdtMap[], c(.05),na.rm=T) ] <- NA
  dhdtWE <- dhdtMap / deltaDays
  dhdt_std <- dhdt_uncertainty / deltaDays
  dhdtWE[dhdtWE[]<dhdt_std]<-NA # kick out all data points that may be below the accuracy of the product
  
  thickInterpolate <- function(debModel,aggrDim){
    aggrMod <- aggregate(debModel,fact = aggrDim)
    xy <- data.frame(xyFromCell(aggrMod , 1:ncell(aggrMod )))
    v <- getValues(aggrMod )
    xy <- xy[which(v>0),]
    v <- v[which(v>0)]
    
    grd <- expand.grid(x=xy[,1],y=xy[,2])
    test <- cbind(v,xy)
    coordinates(test) = ~x+y
    
    m = SpatialPixelsDataFrame(points = grd[c("x", "y")], data = grd)
    
    idw <- idw(formula = v ~ 1, locations = test, newdata = m) # apply idw model for the data
    debMod_all_interp <- raster(idw, layer=1, values=TRUE)
    debMod_all_interp <- raster::resample(debMod_all_interp,debModel)
    debMod_all_interp[which(is.na(dhdtMap_features[]))] <- NA
    debMod_all_interp[which(!is.na(debModel[]))] <- debModel[which(!is.na(debModel[]))]
    
    return(debMod_all_interp)
  }
  
  dhdtWEInterp <- thickInterpolate(dhdtWE,2)

  debMod_all <- ((dhdtWEInterp / crit_melt) / OCparams[modelType,'a'])^(1/OCparams[modelType,'b'])
  debMod_all_low <- (((dhdtWEInterp)/ (crit_melt - crit_melt_std)) / (OCparams["allCurves_noMod",'a'] - OCparams["allCurvesStd_noMod",'a']))^(1/(OCparams["allCurves_noMod",'b'] + OCparams["allCurvesStd_noMod",'b']))
  debMod_all_high <- (((dhdtWEInterp) / (crit_melt + crit_melt_std)) / (OCparams["allCurves_noMod",'a'] + OCparams["allCurvesStd_noMod",'a']))^(1/(OCparams["allCurves_noMod",'b'] - OCparams["allCurvesStd_noMod",'b']))

  #debMod_all[which((dhdtWE[] / crit_melt)<0.05)] <- NA
  debMod_all[debMod_all[]>2] <- NA
  debMod_all_low[debMod_all_low[]>2] <- NA
  debMod_all_high[debMod_all_high[]>2] <- NA
  # Interpolate the missing values
  
  debMod_all_interp <- thickInterpolate(debMod_all,2)
  debMod_all_interp_low <- thickInterpolate(debMod_all_low,2)
  debMod_all_interp_high <- thickInterpolate(debMod_all_high,2)

  return(stack(debMod_all_interp,debMod_all,dhdtWEInterp,debMod_all_interp_low,debMod_all_interp_high,dhdtWEInterp))
  
}

crit_melt <- 0.04
crit_melt_std <- 0.01
debMod16_noMod <- debThickMod(dhdt16,dhdt16_features,5,flightDates,OCparams,crit_melt,"allCurves_noMod",crit_melt_std)

# Uncertainty Tests

sigma_mb <- sqrt(0.25^2+0.25^2)
sigma_cm <- 0.01
sigma_a <- 0.002
sigma_b <- 0.006

cm_range <-rnorm(1000,cm,sigma_cm)
a_range <- rnorm(1000,a,sigma_a)
b_range <- rnorm(1000,b,sigma_b)
mb_range <- rnorm(1000,0,sigma_mb)
a <- 0.13
b <- -0.52
cm <- 0.04
mb_in <- debMod16_noMod[[3]][c(GPRDatID,PitDatID)]

delta_dmed <- vector()
for(k in 1:length(mb_in)){
delta_d <- vector()
for(mc_dmodel in 1:1000){
  a<-a_range[mc_dmodel]
  b<-b_range[mc_dmodel]
  mb<-mb_in[k] + mb_range[mc_dmodel]
  cm<-cm_range[mc_dmodel]
delta_d[mc_dmodel] <- sqrt(
  (((a^(-1/b) * cm^(-1/b) * (mb)^(1/b-1)))/b)^2*sigma_mb^2 +
    (((a^(-1/b) * cm^(-1/b - 1)* (mb)^(1/b-1))/b))^2*sigma_cm^2 +
    (-((a^(-1/b-1)*cm(-1/b) * (mb)^(1/b))/b))^2*sigma_a^2 +
    (-(((mb)/a/cm)^(1/b)*log((mb)/a/cm))/b^2)^2*sigma_b^2)

}
delta_dmed[k] <- median(delta_d,na.rm=T)
}



debMod16_noMod_min <- debMod16_noMod[[4]]
debMod16_noMod_min[which(debMod16_noMod_min[]>=debMod16_noMod[[1]][])] <- debMod16_noMod[[1]][which(debMod16_noMod_min[]>=debMod16_noMod[[1]][])]
debMod16_noMod_max <- debMod16_noMod[[5]]
debMod16_noMod_max[which(debMod16_noMod_max[]<=debMod16_noMod[[1]][])] <- debMod16_noMod[[1]][which(debMod16_noMod_max[]<=debMod16_noMod[[1]][])]

df <- data.frame(X = debMod16_noMod[[1]][GPRDatID], errX_p = delta_dmed[1:length(GPRDatID)],errX_m = delta_dmed[1:length(GPRDatID)], Y = debThick[GPRDatID], errY = 0.04/2)
#df <- df[-which(is.na(df[,1])),]
df2 <- data.frame(X = debMod16_noMod[[1]][PitDatID], errX_p = delta_dmed[(length(GPRDatID)+1):length(delta_dmed)],errX_m = delta_dmed[(length(GPRDatID)+1):length(delta_dmed)], Y = debThick_pits[PitDatID ], errY = 0.05/2)
#df2 <- df2[-which(is.na(df2[,1])),]


RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

RMSE(debMod16_noMod[[1]][c(GPRDatID,PitDatID)],debAll[c(GPRDatID,PitDatID)])

mbe_dh <- 1/length(debMod16_noMod[[1]][c(GPRDatID,PitDatID)])*sum(debMod16_noMod[[1]][c(GPRDatID,PitDatID)]-debAll[c(GPRDatID,PitDatID)])
mae_dh <- sum(abs(debMod16_noMod[[1]][c(GPRDatID,PitDatID)]-debAll[c(GPRDatID,PitDatID)]))/length(debMod16_noMod[[1]][c(GPRDatID,PitDatID)])

png(file=path_figs&'\\GeneralModel_2016_Lirung.png', res = 300,width=1800,height=1800)
ggplot(data = df, aes(x = X, y = Y)) + geom_point(data = df, aes(x = X, y = Y)) + #main graph
 geom_errorbar(aes(ymin = Y-errY, ymax = Y+errY)) + 
 geom_errorbarh(aes(xmin = X-errX_m, xmax = X+errX_p)) +
 geom_point(data = df2,aes(x = X, y = Y),color='red') +
  geom_errorbar(data = df2,aes(ymin = df2$Y-df2$errY, ymax = df2$Y+df2$errY),color='red') + 
 geom_errorbarh(data = df2,aes(xmin = df2$X - df2$errX_m, xmax = df2$X + df2$errX_p),color='red') +
geom_abline() +
 coord_cartesian(xlim=c(0, 2),ylim=c(0, 2), expand = c(0, 0)) +
   labs(x = "modelled thickness [m]", y = "measured thickness [m]") + 
  theme_bw() +
  theme(text = element_text(size=20))+
   theme(plot.margin=unit(c(1,1,0,0),"cm"))

dev.off()

writeRaster(debMod16_noMod[[1]],path_data&'\\DebrisThickness\\modelledThickness_OestremModel.tif','GTiff',overwrite=T)
