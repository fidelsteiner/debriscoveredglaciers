################################################################################
# debEmergence - calculate distributed emergence velocity based on the flux gate approach
# 
# debEmer.R
#
# ReadMe:
# Based on Miles et al., 2018 (GRL)
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
library(pracma)

projec<-'+proj=utm +zone=45N +datum=WGS84'
#projec<-'+proj=longlat +datum=WGS84'

##################
# File Paths/File Names and basic Settings
##################
# Paths PC
path <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB'        
path_code <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code'
path_subcode <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Code\\subCodes\\'
path_figs <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Figures'
path_data <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData'

u_base_frac <- 0.3  # fraction of surface velocity for base sliding
distFGates <- 20    # Distanfce between Flux Gate locations
location <- c('Lirung')   # name of model location for file outputs
for(gl in 1:length(location)){

IceThickness <- 'ice.thickness_'&location[gl]&'.farinotti.tif'                  # Ice thickness raster
glacier_outlines <- location[gl]&'_2015.shp'                       # Outlines for cliff features (write NULL if not applicable)
velFolder <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\'&location[gl]&'\\Velocities'

emergenceFolder <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\emergenceVelocity\\Fluxgates'&location[gl]

for(flg in 1:length(list.files(emergenceFolder,pattern = "shp"))){
fluxgate <- list.files(emergenceFolder,pattern = "shp")[flg]

# Read Ice Thickness
HIce <- raster(path_data&'\\IceThickness\\'&IceThickness)
projection(HIce)<-projec

GOutline <- ogrInfo(path_data&'\\Outlines\\'&glacier_outlines)
GOutline<-readOGR(dsn=path_data&'\\Outlines\\'&glacier_outlines)
projection(GOutline)<-projec

# Make Flux Gate with segments 
FGate <- ogrInfo(emergenceFolder&'\\'&fluxgate)
FGate <-readOGR(dsn=emergenceFolder&'\\'&fluxgate)
projection(FGate)<-projec

num_len <- gLength(FGate)

int_len_seq <- seq(0, num_len, distFGates)

# Two functions below to separate the Gate line into evenly spaced segments
CreateSegment <- function(coords, from, to) {
  distance <- 0
  coordsOut <- c()
  biggerThanFrom <- F
  for (i in 1:(nrow(coords) - 1)) {
    d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i + 
                                                                             1, 2])^2)
    distance <- distance + d
    if (!biggerThanFrom && (distance > from)) {
      w <- 1 - (distance - from)/d
      x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
      y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
      coordsOut <- rbind(coordsOut, c(x, y))
      biggerThanFrom <- T
    }
    if (biggerThanFrom) {
      if (distance > to) {
        w <- 1 - (distance - to)/d
        x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
        y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
        coordsOut <- rbind(coordsOut, c(x, y))
        break
      }
      coordsOut <- rbind(coordsOut, c(coords[i + 1, 1], coords[i + 1, 
                                                               2]))
    }
  }
  return(coordsOut)
}

ls_crd_seg <- lapply(2:length(int_len_seq), function(i) {
  
  # extract coordinates of current line segment
  mat_segment <- CreateSegment(coordinates(FGate)[[1]][[1]], 
                               from = int_len_seq[i-1], to = int_len_seq[i])
  
  # end coordinate
  crd_out <- matrix(mat_segment[nrow(mat_segment), ], ncol = 2)
  
  # during the first iteration, also return the start coordinate
  if (i == 2) {
    crd_start <- matrix(mat_segment[1, ], ncol = 2)
    crd_out <- rbind(crd_start, crd_out)
  } 
  
  return(crd_out)
})

FGateCoords <- do.call("rbind", ls_crd_seg)

HIce <- crop(HIce,GOutline)
HIce <- mask(HIce,GOutline)

HIce <- mask(HIce,GOutline)

velUFiles <- list.files(velFolder&'\\etow-filtered')
velVFiles <- list.files(velFolder&'\\ntos-filtered')

# relevant TimeSteps
relSteps <- c(1,8,14,19,23,26)
gate.flux <- vector()
for(i in 1:8){
  velU <- raster(velFolder&'\\etow-filtered\\'&velUFiles[i])
  velV <- raster(velFolder&'\\ntos-filtered\\'&velVFiles[i])
  velMat <- sqrt(velU^2 + velV^2)


resHIce <- raster::extract(HIce,as.data.frame(SpatialPoints(velMat)))
velMatIce <- velMat
velMatIce[]<-resHIce

u_surf <- unlist(raster::extract(velMat,FGate))
thx <- unlist(raster::extract(velMatIce,FGate))
u_surf <- u_surf[!is.na(thx)]
thx <- thx[!is.na(thx)]

u_b <- u_base_frac * u_surf

meanU <- vector()
for (k in 1:length(u_surf)){

  thx[thx<1]<-1
  z <- seq(1,floor(thx[k]),1)
  meanU[k] <- mean(u_b[k]+u_surf[k]*(1-u_base_frac)*(1-((thx[k]-z)/thx[k])^(3+1)),na.rm=T);

}

Umean = velU*(mean(meanU,na.rm=T)/velMat);
Vmean = velV*(mean(meanU,na.rm=T)/velMat);

VelAcross <- meanU * thx * res(velMat)[1]

# Calculate Azimuth of Gate
gate.dXdY <- diff(FGateCoords[c(1,length(FGateCoords[,1])),])
gate.dist <- cbind(diff(FGateCoords[c(1,length(FGateCoords[,1])),])[1],diff(FGateCoords[c(1,length(FGateCoords[,1])),])[2])
gate.length <- sqrt(sum(gate.dist^2))

gate.dir <- gate.dist/gate.length
gate.fluxdir <- gate.dir%*%as.matrix(cbind(c(0,-1),c(1,0)))

gate.csx <- int_len_seq * gate.dir[1] + FGateCoords[,1]
gate.csy <- int_len_seq * gate.dir[2] + FGateCoords[,2]

gate.umean <- raster::extract(Umean,FGateCoords,method='bilinear')
gate.vmean <- raster::extract(Vmean,FGateCoords,method='bilinear')
gate.speed <- sqrt(gate.umean^2 + gate.vmean^2)
gate.magVel <- rbind(gate.fluxdir[1],gate.fluxdir[2])%*%gate.speed
gate.magSpeed <-  sqrt(gate.magVel[1,]^2+ gate.magVel[2,]^2)
#gate.magSpeed <- sqrt(gate.umean^2+ gate.vmean^2)
# Thickness below same points
gate.thx <- raster::extract(velMatIce,FGateCoords,method='bilinear')
gate.umean[is.na(gate.thx)] <- NA
gate.vmean[is.na(gate.thx)] <- NA

# Final Flux
gate.flux[i] <- sum(gate.magSpeed * gate.thx * 20,na.rm=T)

write.table(cbind(FGateCoords[,1],FGateCoords[,2], gate.umean,gate.vmean),path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg&'_TimeStep'&i&'_Velocity.txt')
}
write.table(gate.flux,path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg&'_Flux.txt')
write.table(gate.thx,path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg&'_Thickness.txt')

#write.table(thx,path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg&'_Thickness.txt')
}
}


##### Distribute emerging volume

# get emerging volume from first fluxgate
fluxgate1 <- list.files(emergenceFolder,pattern = "shp")[1]
FGate1 <- ogrInfo(emergenceFolder&'\\'&fluxgate1)
FGate1 <-readOGR(dsn=emergenceFolder&'\\'&fluxgate1)
projection(FGate1)<-projec

lpi_1 <- gIntersection(FGate1,GOutline)
blpi_1 <- gBuffer(lpi_1, width = 0.000001)
dpi_1 <- gDifference(GOutline, blpi_1)

lowerTongue <- SpatialPolygons(list(Polygons(list(dpi_1@polygons[[1]]@Polygons[[3]]), "1")))
projection(lowerTongue)<-projec
for(flg in 2 : length(list.files(emergenceFolder,pattern = "shp"))){
  fluxgate <- list.files(emergenceFolder,pattern = "shp")[flg]
  
  FGate_prev <- FGate
  # Make Flux Gate with segments 
  FGate <- ogrInfo(emergenceFolder&'\\'&fluxgate)
  FGate <-readOGR(dsn=emergenceFolder&'\\'&fluxgate)
  projection(FGate)<-projec
  
  lpi <- gIntersection(FGate,lowerTongue)
  blpi <- gBuffer(lpi, width = 0.000001)
  dpi <- gDifference(lowerTongue, blpi)
  
  flux1 <- unlist(read.table(path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg-1&'_Flux.txt'))
  flux2 <- unlist(read.table(path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg&'_Flux.txt'))
  emergenceArea <- SpatialPolygons(list(Polygons(list(dpi@polygons[[1]]@Polygons[[1]]), "1")))
  r_emer <- raster()
  extent(r_emer) <- extent(emergenceArea)
  res(r_emer) <- 1
  r_emer[] <- 1
  projection(r_emer)<-projec
  r_emer <- mask(r_emer,emergenceArea)

  
  for(i in 1:length(flux1)){
    r_emer[!is.na(r_emer[])] <- (flux1[i]-flux2[i])/area(emergenceArea)
    if(flg>=3){
      #browser()
      r_emer_old <- raster(path_data&'\\emergenceVelocity\\'&location[gl]&'emergenceField'&flg-1&'timestep'&i&'.tif')
      projection(r_emer_old) <- projec
      writeRaster(merge(r_emer,r_emer_old, tolerance = 0.4),path_data&'\\emergenceVelocity\\'&location[gl]&'emergenceField'&flg&'timestep'&i&'.tif','GTiff',overwrite=TRUE)
    }
    else if(flg<3){
    writeRaster(r_emer,path_data&'\\emergenceVelocity\\'&location[gl]&'emergenceField'&flg&'timestep'&i&'.tif','GTiff',overwrite=TRUE)
    }  
  }
  
  lowerTongue <- SpatialPolygons(list(Polygons(list(dpi@polygons[[1]]@Polygons[[2]]), "1")))
  projection(lowerTongue)<-projec
  
}

# last part of the tongue
r_emer <- raster()
extent(r_emer) <- extent(lowerTongue)
res(r_emer) <- 1
r_emer[] <- 1
projection(r_emer)<-projec
r_emer <- mask(r_emer,lowerTongue)

flux1 <- unlist(read.table(path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg-1&'_Flux.txt'))
flux2 <- unlist(read.table(path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg&'_Flux.txt'))
for(i in 1:length(flux1)){
  r_emer[!is.na(r_emer[])] <- (flux1[i])/area(lowerTongue)
  r_emer_old <- raster(path_data&'\\emergenceVelocity\\'&location[gl]&'emergenceField'&flg&'timestep'&i&'.tif')
  projection(r_emer_old) <- projec
  writeRaster(merge(r_emer,r_emer_old, tolerance = 0.4),path_data&'\\emergenceVelocity\\'&location[gl]&'FINALtimestep'&i&'.tif','GTiff',overwrite=TRUE)


  }

png(file=path_figs&'\\Emergence'&location[gl]&'.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
for(tim in 1:6){
Map <- raster(path_data&'\\emergenceVelocity\\'&location[gl]&'FINALtimestep'&tim&'.tif')
plot(Map)
for(flg in 1:length(list.files(emergenceFolder,pattern = "shp"))){
GateData <- read.table(path_data&'\\emergenceVelocity\\'&location[gl]&'Glacier_Fluxgate'&flg&'_TimeStep'&tim&'_Velocity.txt')
GateData[is.na(GateData)] <- 0
quiver(GateData[,1],GateData[,2], GateData[,3],GateData[,4], scale = 50, angle = 0, length = 0.1)
}
}
dev.off()

for(tim in 1:8){
Map <- raster(path_data&'\\emergenceVelocity\\'&location[gl]&'FINALtimestep'&tim&'.tif')
Map2 <- focal(Map, matrix(1, 21, 21), mean, pad = F, padValue = 0,na.rm=TRUE)
writeRaster(Map2,path_data&'\\emergenceVelocity\\'&location[gl]&'FINALtimestep'&tim&'_smudged.tif','GTiff',overwrite=TRUE)
}