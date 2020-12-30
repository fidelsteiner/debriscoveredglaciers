################################################################################
# View Factors for DEMs used in energy balance models
# 
# vEB.R
#
# ReadMe:
# Based on the original codes developed for the cliff models (Steiner et al. 2015, Buri et al. 2016)
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

DEM_file <- 'dem1m_final.tif'
temp_file <- 'F3_2016_K.tif'
temp_file <- 'tsurf_f1_e_bias.tif'
emis_file <- 'emis_201610.tif'
ModRes <- 1.5    # Spatial Resolution of Model [m]

location<-'Lirung'
projec<-'+proj=utm +zone=45N +datum=WGS84'

# Paths PC
path<-'F:\\PHD\\Research\\EB_DCG\\DistributedEB'
path_code<-'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code'
path_figs<-'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Figures'
path_data<-'F:\\PhD\\Research\\Collaborations\\PleunBonekamp\\LES\\SW'

# packages (if not installed yet: install.packages('examplePackage')
library(rgdal)
library(rgeos)
library(maptools)
library(raster)
library(ggplot2)
library(foreach)
library(GISTools)
library(prettymapr)
library(circular)
library(parallel)
library(snowfall)
library(truncnorm)
library(zoo)

source("F:\\PHD\\Research\\Collaborations\\PleunBonekamp\\LES\\viewFactors_core.R")

##################
# Load Input Data
##################

# read DEM
DEM_r<-raster(path_data&'/'&DEM_file)
projection(DEM_r)<-projec
#DEMclip_r<-crop(DEM_r,extent(glacier_p))

projec<-'+proj=utm +zone=45N +datum=WGS84'
GOutline <- ogrInfo(path_data&'/final_area_30m.shp')
GOutline<-readOGR(path_data&'/final_area_30m.shp')
projection(GOutline)<-projec

# Read Footprint
FP1130 <- raster('F://PhD//Research//EB_DCG//EC_DataProcessing//'&"multilayer_FP.grd",band=691)
FP1100 <- raster('F://PhD//Research//EB_DCG//EC_DataProcessing//'&"multilayer_FP.grd",band=690)
FP1030 <- raster('F://PhD//Research//EB_DCG//EC_DataProcessing//'&"multilayer_FP.grd",band=689)

FP1130 <- crop(FP1130,DEM_r)
FP1130 <- mask(FP1130,DEM_r)
FP1130 <- resample(FP1130,DEM_r) / 4
FP1130[FP1130<quantile(FP1130,0.9)] <- NA
FP1130 <- FP1130 / cellStats(FP1130,'sum') * 100

FP1030 <- crop(FP1030,DEM_r)
FP1030 <- mask(FP1030,DEM_r)
FP1030 <- resample(FP1030,DEM_r) / 4
FP1030[FP1030<quantile(FP1030,0.9)] <- NA
FP1030 <- FP1030 / cellStats(FP1030,'sum') * 100

FP1100 <- crop(FP1100,DEM_r)
FP1100 <- mask(FP1100,DEM_r)
FP1100 <- resample(FP1100,DEM_r) / 4
FP1100[FP1100<quantile(FP1100,0.9)] <- NA
FP1100 <- FP1100 / cellStats(FP1100,'sum') * 100

writeRaster(FP1130, path_data&'/'&'FP1130.tif', format = "GTiff",overwrite=TRUE)
writeRaster(FP1100, path_data&'/'&'FP1100.tif', format = "GTiff",overwrite=TRUE)
writeRaster(FP1030, path_data&'/'&'FP1030.tif', format = "GTiff",overwrite=TRUE)

FP_dates <- read.table('F://PhD//Research//EB_DCG//EC_DataProcessing//'&'FFP_dates.txt',header=F)
FP_flag <- read.table('F://PhD//Research//EB_DCG//EC_DataProcessing//'&'FFP_flag.txt',header=F)


ASWloc <- c(28.23967,85.55711)

d <- data.frame(lon=ASWloc[2], lat=ASWloc[1])
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4326") # WGS 84
coord <- spTransform(d,projec)

DEM_res <- aggregate(DEM_r,fact=1)
ID_Inside <- extract(DEM_res,GOutline,cellnumbers=T)
ID_AWS <- extract(DEM_r,coord,cellnumbers=T)[1]
DEM_dataframe <- cbind(seq(1,length(DEM_res[]),1),DEM_res[],coordinates(DEM_res))
DEM_dataframe_new <- DEM_dataframe[ID_Inside[[1]][,1],]
old <- Sys.time()
testFunc <- function(i){
#vf <- viewFactors_core(DEM_dataframe,DEM_r,res(DEM_res)[1],res(DEM_res)[1],NULL,NULL)
  vf <- viewFactors_core(DEM_dataframe_new[i:(i+1),],DEM_r,res(DEM_res)[1],res(DEM_res)[1],NULL,NULL)
#return(op)
}

# number of workers (leave one to keep system more responsive)
nworkers <- detectCores()-0

# run entire function parallelized
sfInit(parallel=T,cpu=nworkers, type="SOCK")                                                # initialize cluster
loadlib <- lapply(c('raster','truncnorm'),                                                  # load required packages in cluster
                  function(x) sfLibrary(x, character.only=T))
sfExportAll()                                                                               # transfer all variables to clusters
sfClusterSetupRNG()                                                                         # set up random number generator

tout.multi <- system.time(
  out <- sfLapply(c(2:length(DEM_dataframe_new[,1])-1)[1:3600], testFunc)
)
sfStop() 

write.table(unlist(out), path_data&'/tempout.txt', append = FALSE, sep = " ", dec = ".",col.names = TRUE)

# return view factors to raster

viewF_raster <- DEM_res * NA
viewF_raster[ID_Inside[[1]][,1]] <- unlist(out)

viewF_raster_highres <- disaggregate(viewF_raster, fact=1)
writeRaster(viewF_raster_highres, path_data&'/'&'debrisView_201610_LES.tif', format = "GTiff",overwrite=TRUE)
writeRaster(1-viewF_raster_highres, path_data&'/'&'skyView_201610_LES.tif', format = "GTiff",overwrite=TRUE)
writeRaster(viewF_raster_highres, path_data&'/'&'debrisView_201610_LES.tif', format = "GTiff",overwrite=TRUE)

SW_in <- 777.5 / (1 - viewF_raster_highres[ID_AWS]) * (1 - viewF_raster_highres)

# Derive Longwave from debris
# read gridded temperature Data
temp_r<-raster(path_data&'/'&temp_file)
projection(temp_r)<-projec
temp_r <- crop(temp_r,DEM_r) - 7.5
temp_r <- mask(temp_r,DEM_r)
temp_r <- resample(temp_r,DEM_r)

# read gridded emissivity Data
emis_r<-raster(path_data&'/'&emis_file)
projection(emis_r)<-projec
emis_r <- crop(emis_r,DEM_r)
emis_r <- mask(emis_r,DEM_r)
emis_r <- resample(emis_r,DEM_r)
emis_r[is.na(emis_r)] <- cellStats(emis_r,'mean')

#temp_res <- aggregate(temp_r,fact=4)
#emis_res <- aggregate(emis_r,fact=4)

meanTaround <- focal(temp_r, w=matrix(1,nrow=21,ncol=21), fun=mean)
meanearound <- focal(emis_r, w=matrix(1,nrow=21,ncol=21), fun=mean)

viewF_raster_highres <- raster(path_data&'/'&'debrisView_201610_LES.tif')
projection(viewF_raster_highres)<-projec
viewF_raster_highres <- crop(viewF_raster_highres,DEM_r)
viewF_raster_highres <- mask(viewF_raster_highres,DEM_r)
sigma <- 5.670367 * 10^(-8)
LW_grid <- viewF_raster_highres * meanearound * sigma * (meanTaround)^4
LW_grid_extra <- LW_grid - LW_grid[ID_AWS]
writeRaster(LW_grid, path_data&'/'&'LWdebris_201610_f3.tif', format = "GTiff",overwrite=TRUE)
writeRaster(LW_grid_extra, path_data&'/'&'LWdebris_extra_201610_f3.tif', format = "GTiff",overwrite=TRUE)

# make plots 

LW_grid_extra <- raster(path_data&'/'&'LWdebris_extra_201610_f3.tif')
LW_grid <- raster(path_data&'/'&'LWdebris_201610_f3.tif')

vFGrid <- raster(path_data&'/'&'skyView_201610_LES.tif')

pal <- colorRampPalette(c("white","red"))
pal2 <- colorRampPalette(c("white","black"))
png(file=path_figs&'\\SpatialRadiation.png', res = 160,width=2600,height=1300)
par(mar=c(5,7,2,6),cex.lab=1.5,cex.axis=1.5)
par(mfrow=c(1,2))
plot(vFGrid,ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list(text=expression('sky view [-]'), side=2, font=1, line=0.3, cex=1.5),col=pal2(100))
points(coord,col='white',lwd = 4,cex = 2.5)
grid(nx=NULL, ny=NULL)
plot(LW_grid_extra,ylab='',xlab='Longitude [m]',legend.args=list(text=expression('LW from debris [W '~ m^{-2} ~']'), side=2, font=1, line=0.3, cex=1.5),col=pal(100))
grid(nx=NULL, ny=NULL)
points(coord,col='white',lwd = 4,cex = 2.5)
dev.off()

# Make view factor Grids for all years

path_DEMs <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\DEMs\\reprocessed_lirung_data\\warped-dems_gradient-corrected'
listLanDEMs <- list.files(path = path_DEMs, pattern = '.tif', all.files = FALSE,
           full.names = TRUE)

# Call original code
source("F:\\PHD\\Research\\Collaborations\\PleunBonekamp\\LES\\viewFactors_core.R")
source("F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\viewFactors_core.R")
old <- Sys.time()
testFunc <- function(i){
  #vf <- viewFactors_core(DEM_dataframe,DEM_r,res(DEM_res)[1],res(DEM_res)[1],NULL,NULL)
  vf <- viewFactors_core(DEM_dataframe[i:(i+1),],DEMclip_r,res(DEMclip_r)[1],res(DEMclip_r)[1],NULL,NULL)
  #return(op)
}

for(k in 1:length(listLanDEMs)){
# read DEM
DEM_r<-raster(listLanDEMs[k])
projection(DEM_r)<-projec
DEMclip_r<-aggregate(DEM_r,fact = 5/res(DEM_r)[1])
DEM_dataframe <- cbind(seq(1,length(DEMclip_r[]),1),DEMclip_r[],coordinates(DEMclip_r))


# number of workers (leave one to keep system more responsive)
nworkers <- detectCores()-0

# run entire function parallelized
sfInit(parallel=T,cpu=nworkers, type="SOCK")                                                # initialize cluster
loadlib <- lapply(c('raster','truncnorm'),                                                  # load required packages in cluster
                  function(x) sfLibrary(x, character.only=T))
sfExportAll()                                                                               # transfer all variables to clusters
sfClusterSetupRNG()                                                                         # set up random number generator

tout.multi <- system.time(
  out <- sfLapply(which(!is.na(DEM_dataframe[,2])), testFunc)
)
sfStop() 

viewF_raster <- DEMclip_r * NA
viewF_raster[which(!is.na(DEM_dataframe[,2]))] <- unlist(out)

viewF_raster_highres <- disaggregate(viewF_raster, fact=5/res(DEM_r)[1])
writeRaster(viewF_raster_highres, paste('F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\debrisView_',k,'.tif'), format = "GTiff",overwrite=TRUE)

}
