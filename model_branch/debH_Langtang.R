################################################################################
# debH - extract point debris thickness measurements and aggregate to Raster 
# 
# debH_Langtang.R
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

projec<-'+proj=utm +zone=45N +datum=WGS84'
#projec<-'+proj=longlat +datum=WGS84'

path <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB'        
path_code <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code'
path_subcode <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\'
path_figs <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Figures'
path_data <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData'
path_dhdt <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt'
Lirung_DEMFolder <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt\\Lirung'
Langtang_DEMFolder <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt\\Langtang'

glacier_outlines <- 'debris_Lirung2006.shp'                       # Outlines for cliff features (write NULL if not applicable)

GOutline <- ogrInfo(path_data&'\\Outlines\\'&glacier_outlines)
GOutline<-readOGR(dsn=path_data&'\\Outlines\\'&glacier_outlines)
projection(GOutline)<-projec

deb_GPR <- read.csv(path_data&'\\DebrisThickness\\DebrisThickness_GPR.csv',header = T)
deb_Pits <- read.csv(path_data&'\\DebrisThickness\\DebrisThickness_Pits.csv',header = T)
deb_Pits_full <- read.csv(path_data&'\\DebrisThickness\\DebrisThickness_Pits_fullpitsonly.csv',header = T)

xy <- data.frame(ID = seq(1,length(deb_GPR$X),1), X = deb_GPR$X, Y = deb_GPR$Y)
coordinates(xy) <- cbind(xy$X,xy$Y)
projection(xy)<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
resGPR <- spTransform(xy, projec)

df <- as.data.frame(cbind(resGPR,deb_GPR$Debris.thickness..m.))
colnames(df)<-c('ID','WGSX','WGSY','d','X','Y')

e <- extent(GOutline)
rDepth <- raster(e,res=1)
projection(rDepth) <- projec
rDepth_GPR <- rasterize(df[, 5:6], rDepth, df[,4], fun=mean)
writeRaster(rDepth_GPR,path_data&'\\DebrisThickness\\GPRThickness_2015.tif','GTiff',overwrite=TRUE)

xy <- data.frame(ID = seq(1,length(deb_Pits$X),1), X = deb_Pits$X, Y = deb_Pits$Y)
coordinates(xy) <- cbind(xy$X,xy$Y)
projection(xy)<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
resPits <- spTransform(xy, projec)

df <- as.data.frame(cbind(resPits,deb_Pits$Debris.thickness..m.))
colnames(df)<-c('ID','WGSX','WGSY','d','X','Y')

e <- extent(GOutline)
rDepth <- raster(e,res=1)
projection(rDepth) <- projec
rDepth_Pits <- rasterize(df[, 5:6], rDepth, df[,4], fun=mean)
writeRaster(rDepth_Pits,path_data&'\\DebrisThickness\\PitThickness_2016.tif','GTiff',overwrite=TRUE)

xy <- data.frame(ID = seq(1,length(deb_Pits_full$X),1), X = deb_Pits_full$X, Y = deb_Pits_full$Y)
coordinates(xy) <- cbind(xy$X,xy$Y)
projection(xy)<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
resPits <- spTransform(xy, projec)

df <- as.data.frame(cbind(resPits,deb_Pits_full$Debris.thickness..m.))
colnames(df)<-c('ID','WGSX','WGSY','d','X','Y')

e <- extent(GOutline)
rDepth_full <- raster(e,res=1)
projection(rDepth_full) <- projec
rDepth_Pits_full <- rasterize(df[, 5:6], rDepth_full, df[,4], fun=mean)
writeRaster(rDepth_Pits_full,path_data&'\\DebrisThickness\\PitThickness_2016_full.tif','GTiff',overwrite=TRUE)

writeRaster(merge(rDepth_Pits_full,rDepth_GPR),path_data&'\\DebrisThickness\\MeasuredThickness_1m.tif','GTiff',overwrite=TRUE)


png(file=path_figs&'\\measuredDepth.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.5,cex.axis=1.5)
#layout(matrix(c(1,2), nrow = 2, ncol = 1,byrow=TRUE))
#par(xpd=FALSE)
plot(aggregate(rDepth_GPR,fact=10),legend.args=list(text=expression('thickness [m]')),xlab = 'Easting [m]', ylab ='Northing [m]',zlim=c(0.1,2.1))
plot(aggregate(rDepth_Pits,fact=10),add=T,zlim=c(0.1,2.1))
plot(GOutline,add=T)
grid(NULL,NULL)

dev.off()