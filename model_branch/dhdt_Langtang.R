################################################################################
# dhdt - make dhdt maps for Langtang catchment
# 
# dhdt_Langtang.R
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
path_emer <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\emergenceVelocity'
Lirung_DEMFolder <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt\\Lirung'
Langtang_DEMFolder <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt\\Langtang'

glacier_outlines <- 'Lirung_2015.shp'                       # Outlines for cliff features (write NULL if not applicable)

GOutline <- ogrInfo(path_data&'\\Outlines\\'&glacier_outlines)
GOutline<-readOGR(dsn=path_data&'\\Outlines\\'&glacier_outlines)
projection(GOutline)<-projec

  dhdt_Lirung_201305_201310 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[1]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[2])
  dhdt_Lirung_201310_201405 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[2]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[3])
  dhdt_Lirung_201405_201510 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[3]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[4])
  
  dhdt_Lirung_201310_201510 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[2]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[4])
  
  dhdt_Lirung_201510_201604 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[4]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[5])
  dhdt_Lirung_201604_201610 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[5]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[6])
  dhdt_Lirung_201610_201704 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[6]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[7])
  dhdt_Lirung_201704_201710 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[7]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[8])
  dhdt_Lirung_201710_201804 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[8]) - raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[9])
  
  projection(dhdt_Lirung_201305_201310)<-projec
  writeRaster(dhdt_Lirung_201305_201310, filename = path_dhdt&'\\2013Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  projection(dhdt_Lirung_201310_201405)<-projec
  writeRaster(dhdt_Lirung_201310_201405, filename = path_dhdt&'\\2013-2014Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  projection(dhdt_Lirung_201405_201510)<-projec
  writeRaster(dhdt_Lirung_201405_201510, filename = path_dhdt&'\\2014-2015Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  projection(dhdt_Lirung_201310_201510)<-projec
  writeRaster(dhdt_Lirung_201310_201510, filename = path_dhdt&'\\2013-2015Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  projection(dhdt_Lirung_201510_201604)<-projec
  writeRaster(dhdt_Lirung_201510_201604, filename = path_dhdt&'\\2015-2016Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  projection(dhdt_Lirung_201604_201610)<-projec
  writeRaster(dhdt_Lirung_201604_201610, filename = path_dhdt&'\\2016Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  projection(dhdt_Lirung_201610_201704)<-projec
  writeRaster(dhdt_Lirung_201610_201704, filename = path_dhdt&'\\2016-2017Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  projection(dhdt_Lirung_201704_201710)<-projec
  writeRaster(dhdt_Lirung_201704_201710, filename = path_dhdt&'\\2017Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  projection(dhdt_Lirung_201710_201804)<-projec
  writeRaster(dhdt_Lirung_201710_201804, filename = path_dhdt&'\\2017-2018Lirung_dhdt.tiff', 'GTiff',overwrite=TRUE)
  
  
  ## optionally correct dhdt map for emergence
  
  dhdt_emer <- function(emerg_map,dhdt_rast,filename,dayin,dayout){
  #emerg_map <- '2013Emergence_Lirung_smudged.tif'                   # raw emergence map
  emergence_Lirung_201305_201310 <- raster(path_emer&'\\EmergenceLirung\\'&emerg_map)
  projection(emergence_Lirung_201305_201310) <- projec
  projection(dhdt_rast)<-projec  
  emer_mask <- mask(emergence_Lirung_201305_201310,GOutline)
  dhdt_mask <- mask(dhdt_rast,GOutline)
  raster_domain <- raster()
  extent(raster_domain) <- extent(GOutline)
  res(raster_domain) <- 1
  raster_domain[] <- 0
  projection(raster_domain)<-projec
  raster_domain <- mask(raster_domain,GOutline)
  dhdtraw <- crop(dhdt_mask,raster_domain)
  dhdtraw <- mask(dhdtraw,GOutline)
  
  emergence <-crop(emer_mask,raster_domain)
  emergence <- mask(emer_mask,GOutline)
  emergence <- raster::resample(emergence,dhdtraw,'bilinear')
  emergence[is.na(emergence)] <- 0

  dhdtemer <- dhdtraw + emergence / (365) * (dayout-dayin)
  #dhdtEmer2 <- raster::extract(emer_mask,coordinates(dhdt_mask))
  #dhdt_emergence <- dhdt_mask
  #dhdt_emergence[which(!is.na(dhdtEmer2))] <- dhdt_emergence[which(!is.na(dhdtEmer2))] - (dhdtEmer2[which(!is.na(dhdtEmer2))]/365*(dayout-dayin))
  
  writeRaster(dhdtemer, filename = path_dhdt&'\\'&filename, 'GTiff',overwrite=TRUE)
  return(emer_mask/365*(dayout-dayin))
  }
  
  dayNums <- as.numeric(strftime(c('2013-05-18','2013-10-22','2014-05-01','2015-10-18','2016-04-30','2016-10-06','2017-04-20','2017-10-19','2018-04-28'), format = "%j"))


  emergence2013 <-  dhdt_emer('2013Emergence_Lirung_smudged.tif',dhdt_Lirung_201305_201310,'2013Lirung_dhdtEmergence.tiff',dayNums[1],dayNums[2])
  emergence2014 <-  dhdt_emer('2013-2014Emergence_Lirung_smudged.tif',dhdt_Lirung_201310_201405,'2013-2014Lirung_dhdtEmergence.tiff',dayNums[2],dayNums[3])
  emergence2015 <- dhdt_emer('2014-2015Emergence_Lirung_smudged.tif',dhdt_Lirung_201405_201510,'2014-2015Lirung_dhdtEmergence.tiff',dayNums[3],dayNums[4])
  emergencem2016 <- dhdt_emer('2015-2016Emergence_Lirung_smudged.tif',dhdt_Lirung_201510_201604,'2015-2016Lirung_dhdtEmergence.tiff',dayNums[4],dayNums[5])
  emergenceo2016 <- dhdt_emer('2016Emergence_Lirung_smudged.tif',dhdt_Lirung_201604_201610,'2016Lirung_dhdtEmergence.tiff',dayNums[5],dayNums[6])
  emergencem2017 <- dhdt_emer('2016-2017Emergence_Lirung_smudged.tif',dhdt_Lirung_201610_201704,'2016-2017Lirung_dhdtEmergence.tiff',dayNums[6],dayNums[7])
  emergenceo2017 <- dhdt_emer('2017Emergence_Lirung_smudged.tif',dhdt_Lirung_201704_201710,'2017Lirung_dhdtEmergence.tiff',dayNums[7],dayNums[8])
  emergencem2018 <- dhdt_emer('2017-2018Emergence_Lirung_smudged.tif',dhdt_Lirung_201710_201804,'2017-2018Lirung_dhdtEmergence.tiff',dayNums[8],dayNums[9])
  
  # Visualize mass loss with emergence corrected maps
  
  path_emergence_dhdt <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt\\Lirung\\dhdt_emergence'
  path_raw_dhdt <- 'F:\\PhD\\Research\\EB_DCG\\DistributedEB\\RawData\\dhdt\\Lirung\\dhdt_raw'

  rawList <- list.files(path = path_raw_dhdt, pattern = '.tif',full.names = T)
  
  glac_outline <- 'debris_Lirung2006.shp'                   # Outline of glacier tongue to be investigated
  pond_outlines_13 <- 'Ponds_Oct13.shp'                        # Outlines for pond features (write NULL if not applicable)
  cliff_outlines_13 <- 'Cliffs_Oct13.shp'                      # Outlines for cliff features (write NULL if not applicable)
  pond_outlines_14 <- 'Ponds_May14.shp'                        # Outlines for pond features (write NULL if not applicable)
  cliff_outlines_14 <- 'Cliffs_May14.shp'
  pond_outlines_15 <- 'Ponds_Oct15.shp'                        # Outlines for pond features (write NULL if not applicable)
  cliff_outlines_15 <- 'Cliffs_Oct15.shp'
  pond_outlines_m16 <- 'Ponds_Apr16.shp'                        # Outlines for pond features (write NULL if not applicable)
  cliff_outlines_m16 <- 'Cliffs_Apr16.shp'
  pond_outlines_o16 <- 'Ponds_Oct16.shp'                        # Outlines for pond features (write NULL if not applicable)
  cliff_outlines_o16 <- 'Cliffs_Oct16.shp'
  pond_outlines_m17 <- 'Ponds_Apr17.shp'                        # Outlines for pond features (write NULL if not applicable)
  cliff_outlines_m17 <- 'Cliffs_Apr17.shp'
  pond_outlines_o17 <- 'Ponds_Oct17.shp'                        # Outlines for pond features (write NULL if not applicable)
  cliff_outlines_o17 <- 'Cliffs_Oct17.shp'
  pond_outlines_m18 <- 'Ponds_Apr18.shp'                        # Outlines for pond features (write NULL if not applicable)
  cliff_outlines_m18 <- 'Cliffs_Apr18.shp'
  
  ogrInfo('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines\\'&glac_outline)
  glac_mask<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines\\'&glac_outline)
  
  # dhdt data
  dhdtraw_2013<-raster(rawList[1])
  dhdtraw_2014<-raster(rawList[2])
  dhdtraw_2015<-raster(rawList[9])
  dhdtraw_2015_2014<-raster(rawList[3])
  dhdtraw_m2016<-raster(rawList[4])
  dhdtraw_o2016<-raster(rawList[5])
  dhdtraw_m2017<-raster(rawList[6])
  dhdtraw_o2017<-raster(rawList[7])
  dhdtraw_m2018<-raster(rawList[8])
  ortho06 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\Figures\\Maps\\LargeOrthoLirung.tif')
  ortho13 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Orthos\\20131022_lirung_ortho_10cm_quakeshift.tif')
  ortho14 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Orthos\\20140501_lirung_ortho_10cm_quakeshift.tif')
  ortho15 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Orthos\\20151018_lirung_ortho_10cm.tif')
  orthom16 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Orthos\\20160430_lirung_ortho_10cm.tif')
  orthoo16 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Orthos\\20161006_lirung_ortho_10cm.tif')
  orthom17 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Orthos\\20170420_lirung_ortho_10cm.tif')
  orthoo17 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Orthos\\20171019_lirung_ortho_10cm.tif')  
  orthom18 <- brick('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Orthos\\20180428_lirung_ortho_10cm.tif') 
  
  ext_ortho <- extent(glac_mask) + c(0,500,-500,0)
  ext_ortho_plot <- ext_ortho + c(500,0,400,-1000)
  dhdtdataViz <- function(cliff_outlines,pond_outlines,ortho,dhdtmap,emergencemap, dayin,dayout,numyears){
  ogrInfo('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines)
  pond_mask<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines)
  ogrInfo('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines)
  cliff_mask<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines)
  
  projection(dhdtmap)<-projec

  raster_domain <- raster()
  extent(raster_domain) <- extent(glac_mask)
  res(raster_domain) <- 1
  raster_domain[] <- 0
  projection(raster_domain)<-projec
  raster_domain <- mask(raster_domain,glac_mask)
  dhdtraw <- crop(dhdtmap,raster_domain)
  dhdtraw <- mask(dhdtraw,glac_mask)
  emergencemap <-crop(emergencemap,raster_domain)
  emergencemap <- mask(emergencemap,glac_mask)
  emergencemap <- raster::resample(emergencemap,dhdtraw,'bilinear')
  emergencemap[is.na(emergencemap)] <- 0

  dhdtemer <- dhdtraw + emergencemap
  dhdtemer_cliff <- mask(dhdtemer,cliff_mask)
  dhdtemer_pond <- mask(dhdtemer,pond_mask)

  featuresID <- c(which(!is.na(dhdtemer_cliff[])),which(!is.na(dhdtemer_pond[])))
  dhdtemer_debris <- dhdtemer
  dhdtemer_debris[featuresID] <- NA
  #dhdtemer_debris[dhdtemer_debris>10]<- NA

  addeddays <- numyears * 365
  dhdtraw_debris <- crop(dhdtemer_debris,glac_mask) 
  dhdtemer_cliff <- crop(dhdtemer_cliff,glac_mask)
  dhdtemer_pond <- crop(dhdtemer_pond,glac_mask)

  dhdtraw_debris_ma <- raster::resample(dhdtemer_debris,raster_domain,method='bilinear') / (dayout + addeddays - dayin) * 365
  dhdtemer_cliff_ma <- raster::resample(dhdtemer_cliff,raster_domain,method='bilinear') / (dayout + addeddays - dayin) * 365
  dhdtemer_pond_ma <- raster::resample(dhdtemer_pond,raster_domain,method='bilinear') / (dayout + addeddays - dayin) * 365
  
  dhdtraw_debris <- raster::resample(dhdtemer_debris,raster_domain,method='bilinear')
  dhdtemer_cliff <- raster::resample(dhdtemer_cliff,raster_domain,method='bilinear')
  dhdtemer_pond <- raster::resample(dhdtemer_pond,raster_domain,method='bilinear')
  dhdtemer_all <- raster::resample(dhdtemer,raster_domain,method='bilinear')

  stackResult <- stack(dhdtraw_debris_ma,dhdtemer_cliff_ma,dhdtemer_pond,dhdtraw_debris,dhdtemer_cliff,dhdtemer_pond,dhdtemer_all)
    names(stackResult) <- c('debrismelt','cliffmelt','pondmelt','debrismelt_actual','cliffmelt_actual','pondmelt_actual','dhdtemer_all_actual')
  return(stackResult)
  }
  
  stack2013 <- dhdtdataViz(cliff_outlines_13,pond_outlines_13,ortho13,dhdtraw_2013,emergence2013,dayNums[1],dayNums[2],0)
  pond_mask_13<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_13)
  cliff_mask_13<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_13)
  #ortho13 <- crop(ortho13,ext_ortho)
  
  stack2014 <- dhdtdataViz(cliff_outlines_14,pond_outlines_14,ortho14,dhdtraw_2014,emergence2014,dayNums[2],dayNums[3],1)
  pond_mask_14<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_14)
  cliff_mask_14<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_14)
  #ortho14 <- crop(ortho14,ext_ortho)
  
  stack2015_2014 <- dhdtdataViz(cliff_outlines_15,pond_outlines_15,ortho15,dhdtraw_2015_2014,emergence2015,dayNums[3],dayNums[4],1)
  pond_mask_15<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_15)
  cliff_mask_15<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_15)
  #ortho15 <- crop(ortho15,ext_ortho)
  
  stack2015 <- dhdtdataViz(cliff_outlines_15,pond_outlines_15,ortho15,dhdtraw_2015,emergence2015,dayNums[2],dayNums[4],2)
  pond_mask_15<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_15)
  cliff_mask_15<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_15)
  #ortho15 <- crop(ortho15,ext_ortho)
  
  stackm2016 <- dhdtdataViz(cliff_outlines_m16,pond_outlines_m16,orthom16,dhdtraw_m2016,emergencem2016,dayNums[4],dayNums[5],1)
  pond_mask_m16<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_m16)
  cliff_mask_m16<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_m16)
  #orthom16 <- crop(orthom16,ext_ortho)
  
  stacko2016 <- dhdtdataViz(cliff_outlines_o16,pond_outlines_o16,orthoo16,dhdtraw_o2016,emergenceo2016,dayNums[5],dayNums[6],0)
  pond_mask_o16<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_o16)
  cliff_mask_o16<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_o16)
  #orthoo16 <- crop(orthoo16,ext_ortho)
  
  stackm2017 <- dhdtdataViz(cliff_outlines_m17,pond_outlines_m17,orthom17,dhdtraw_m2017,emergencem2017,dayNums[6],dayNums[7],1)
  pond_mask_m17<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_m17)
  cliff_mask_m17<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_m17)
  #orthom17 <- crop(orthom17,ext_ortho)
  
  stacko2017 <- dhdtdataViz(cliff_outlines_o17,pond_outlines_o17,orthoo17,dhdtraw_o2017,emergenceo2017,dayNums[7],dayNums[8],0)
  pond_mask_o17<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_o17)
  cliff_mask_o17<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_o17)
  #orthoo17 <- crop(orthoo17,ext_ortho)
  
  stackm2018 <- dhdtdataViz(cliff_outlines_m18,pond_outlines_m18,orthom18,dhdtraw_m2018,emergencem2018,dayNums[8],dayNums[9],1)
  pond_mask_m18<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&pond_outlines_m18)
  cliff_mask_m18<-readOGR(dsn='F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\Outlines_Features\\'&cliff_outlines_m18)
  #orthom18 <- crop(orthom18,ext_ortho)
  #orthom18[[1]][orthom18[[1]]==0] <- 255
  #orthom18[[2]][orthom18[[2]]==0] <- 255
  #orthom18[[3]][orthom18[[3]]==0] <- 255
  
  colorr <- brewer.pal(20, 'RdBu')
  
  png(file=path_figs&'\\dh_2014_problem.png', res = 160,width=3600,height=1800)
  par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
  #par(mfrow=c(3,1))
  layout(matrix(c(1,2,3,4,5,6), nrow = 1, ncol = 2, byrow = TRUE))
  plot(ext_ortho_plot,xlab = '',ylab='')
  plotRGB(ortho13,margins = T,add=T)
  axis(1,las = 1)
  axis(2,las = 3)
  plot(stack2014$debrismelt * (-1),add=T,zlim=c(-15,3),col=colorr[1:7],ylab='Latitude [m]',xlab='Longitude [m]',legend = F)
  #plot(cliff_mask_13,add=T,col='white')
  #plot(pond_mask_13,add=T,col='blue')
  plot(ext_ortho_plot,xlab = '',ylab='')
  plotRGB(orthom16,margins = T,add=T)
  axis(1,las = 1)
  axis(2,las = 3)
  plot(stack2015$debrismelt * (-1),add=T,zlim=c(-15,3),col=colorr[1:7],ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list(text=expression('mass loss [m/yr]'), side=2, font=1, line=0.3, cex=1))
  #plot(cliff_mask_m16,add=T,col='white')
  #plot(pond_mask_m16,add=T,col='blue')
  dev.off()
  
  
  stackwinter <- stack(stackm2016$debrismelt,stackm2017$debrismelt,stackm2018$debrismelt)
  stacksummer <- stack(stack2013$debrismelt,stacko2016$debrismelt,stacko2017$debrismelt)

  
  fun <- function(x) { mean(x,na.rm=T) }
  meanWinter <- calc(stackwinter, fun)
  meanSummer <- calc(stacksummer, fun)
  
  writeRaster(stackwinter[[1]], filename = path_dhdt&'\\2015-2016Lirung_dhdt_ma_winter.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stackwinter[[2]], filename = path_dhdt&'\\2016-2017Lirung_dhdt_ma_winter.tiff', 'GTiff',overwrite=TRUE)  
  writeRaster(stackwinter[[3]], filename = path_dhdt&'\\2017-2018Lirung_dhdt_ma_winter.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stacksummer[[1]], filename = path_dhdt&'\\2013Lirung_dhdt_ma_summer.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stacksummer[[2]], filename = path_dhdt&'\\2016Lirung_dhdt_ma_summer.tiff', 'GTiff',overwrite=TRUE)  
  writeRaster(stacksummer[[3]], filename = path_dhdt&'\\2017Lirung_dhdt_ma_summer.tiff', 'GTiff',overwrite=TRUE) 
  
  writeRaster(stack2013$dhdtemer_all_actual, filename = path_dhdt&'\\2013Lirung_dhdt_total.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stack2015$dhdtemer_all_actual, filename = path_dhdt&'\\2013-2015Lirung_dhdt_total.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stackm2016$dhdtemer_all_actual, filename = path_dhdt&'\\2015-2016Lirung_dhdt_total.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stacko2016$dhdtemer_all_actual, filename = path_dhdt&'\\2016Lirung_dhdt_total.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stackm2017$dhdtemer_all_actual, filename = path_dhdt&'\\2016-2017Lirung_dhdt_total.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stacko2017$dhdtemer_all_actual, filename = path_dhdt&'\\2017Lirung_dhdt_total.tiff', 'GTiff',overwrite=TRUE)
  writeRaster(stackm2018$dhdtemer_all_actual, filename = path_dhdt&'\\2017-2018Lirung_dhdt_total.tiff', 'GTiff',overwrite=TRUE)
  
  writeRaster(meanWinter, filename = path_dhdt&'\\meanwinter_dhdt_ma.tiff', 'GTiff',overwrite=TRUE) 
  writeRaster(meanSummer, filename = path_dhdt&'\\meansummer_dhdt_ma.tiff', 'GTiff',overwrite=TRUE)
  
  
  
  # Read ASTER GDEM for surrounding area
  DEM_NSIDC<-raster(path_data&'/DEMs/'&'langtang_NSIDCDEM.tif')
  projection(DEM_NSIDC) <- projec
  DEM_NSIDC_crop <- crop(DEM_NSIDC,ext_ortho_plot)
  ContoursLirung <- rasterToContour(DEM_NSIDC_crop)
  
  wd <- ext_ortho_plot[2]-ext_ortho_plot[1]
  hg <- ext_ortho_plot[4]-ext_ortho_plot[3]
  
  png(file=path_figs&'\\dh_2013-2018.png', res = 160,width = wd * 0.75,height=hg/2)
  par(mar=c(3,3,1,1),cex.lab=1.5,cex.axis=1.5)
  #par(mfrow=c(3,1))
  layout(matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, byrow = TRUE))
  plot(ext_ortho_plot,xlab = '',ylab='', xaxt='n')
  plotRGB(orthom16,margins = T,add=T)
  axis(1,las = 1,labels=F)
  axis(2,las = 3)
  plot(stack2013$debrismelt * (-1),add=T,zlim=c(-5,2),col=colorr[1:8],ylab='Latitude [m]',xlab='Longitude [m]',legend = F)
  plot(cliff_mask_13,add=T,col='white')
  plot(pond_mask_13,add=T,col='blue')
  plot(ContoursLirung,add=T)
  
  plot(ext_ortho_plot,xlab = '',ylab='', yaxt='n', xaxt='n')
  plotRGB(orthom16,margins = T,add=T)
  axis(1,las = 1,labels=F)
  axis(2,las = 3,labels=F)
  plot(stacko2016$debrismelt * (-1),add=T,zlim=c(-5,2),col=colorr[1:8],ylab='Latitude [m]',xlab='Longitude [m]',legend = F)
  plot(cliff_mask_o16,add=T,col='white')
  plot(pond_mask_o16,add=T,col='blue')
  plot(ContoursLirung,add=T)
  
  plot(ext_ortho_plot,xlab = '',ylab='', yaxt='n', xaxt='n')
  plotRGB(orthom16,margins = T,add=T)
  axis(1,las = 1,labels=F)
  axis(2,las = 3,labels=F)
  #legend.args=list(text=expression('mass loss [m/yr]'), side=2, font=1, line=0.3, cex=1)
  plot(stacko2017$debrismelt * (-1),add=T,zlim=c(-5,2),col=colorr[1:8],ylab='Latitude [m]',xlab='Longitude [m]', legend = F)
  plot(cliff_mask_o17,add=T,col='white')
  plot(pond_mask_o17,add=T,col='blue')
  plot(ContoursLirung,add=T)
  
  plot(ext_ortho_plot,xlab = '',ylab='')
  plotRGB(orthom16,margins = T,add=T)
  axis(1,las = 1)
  axis(2,las = 3)
  plot(stackm2016$debrismelt * (-1),add=T,zlim=c(-3,1),col=colorr[3:7],ylab='Latitude [m]',xlab='Longitude [m]',legend = F)
  plot(cliff_mask_m16,add=T,col='white')
  plot(pond_mask_m16,add=T,col='blue')
  plot(ContoursLirung,add=T)
  
  plot(ext_ortho_plot,xlab = '',ylab='',yaxt='n', xaxt='n')
  plotRGB(orthom16,margins = T,add=T)
  axis(1,las = 1)
  axis(2,las = 3,labels=F)
  plot(stackm2017$debrismelt * (-1),add=T,zlim=c(-3,1),col=colorr[3:7],ylab='Latitude [m]',xlab='Longitude [m]',legend = F)
  plot(cliff_mask_m17,add=T,col='white')
  plot(pond_mask_m17,add=T,col='blue')
  plot(ContoursLirung,add=T)
  
  plot(ext_ortho_plot,xlab = '',ylab='',yaxt='n', xaxt='n')
  plotRGB(orthom16,margins = T,add=T)
  axis(1,las = 1)
  axis(2,las = 3,labels=F)
  plot(stackm2018$debrismelt * (-1),add=T,zlim=c(-3,1),col=colorr[3:7],ylab='Latitude [m]',xlab='Longitude [m]',legend = F)
  plot(cliff_mask_m18,add=T,col='white')
  plot(pond_mask_m18,add=T,col='blue')
  plot(ContoursLirung,add=T)
  dev.off()  
  


  # Mass Loss according to elevation bands
  raster_domain <- raster()
  extent(raster_domain) <- extent(glac_mask)
  res(raster_domain) <- 1
  raster_domain[] <- 0
  projection(raster_domain)<-projec
  raster_domain <- mask(raster_domain,glac_mask)
  DEMNSIDC_resampled <- raster::resample(DEM_NSIDC_crop,raster_domain)
  

  histMassLossElevBand <- function(debrismelt,DEM,elevStep,MinVal,MaxVal){
  
  debrismelt[debrismelt<MinVal*(-1)|debrismelt>MaxVal] <- NA
    histMelt<-vector()
  k2<-1
elevStep<-20
  for(k in seq(4020, 4700, by=elevStep)){
    histMelt[k2] <- mean(debrismelt[DEM<k&DEM>k-elevStep],na.rm=T)
  k2 <- k2+1
  }
return(histMelt)
  }
  
  histo2013 <- histMassLossElevBand(stack2013$debrismelt,DEMNSIDC_resampled,20,2,5)
  histm2016 <- histMassLossElevBand(stackm2016$debrismelt,DEMNSIDC_resampled,20,2,5)
  histo2016 <- histMassLossElevBand(stacko2016$debrismelt,DEMNSIDC_resampled,20,2,5)
  histm2017 <- histMassLossElevBand(stackm2017$debrismelt,DEMNSIDC_resampled,20,2,5)
  histo2017 <- histMassLossElevBand(stacko2017$debrismelt,DEMNSIDC_resampled,20,2,5)
  histm2018 <- histMassLossElevBand(stackm2018$debrismelt,DEMNSIDC_resampled,20,2,5)
  
  ContoursLirung_20m <- rasterToContour(DEM_NSIDC_crop,levels=seq(4020, 4700, by=20))
  
  png(file=path_figs&'\\dh_summer_winter.png', res = 160,width = wd,height=hg)
  par(mar=c(3,3,1,4),cex.lab=1.5,cex.axis=1.5)
  #par(mfrow=c(3,1))
  layout(matrix(c(1,3,2,4), nrow = 2, ncol = 2, byrow = TRUE))
  plot(ext_ortho_plot,xlab = '',ylab='', xaxt='n')
  #plotRGB(ortho13,margins = T,add=T)
  axis(1,las = 1,labels=F)
  axis(2,las = 3)
  plot(meanSummer* (-1),add=T,zlim=c(-5,2),col=colorr[1:8],ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list(text='', side=2, font=1, line=0.3, cex=1))
  #plot(cliff_mask_13,add=T,col='white')
  #plot(pond_mask_13,add=T,col='blue')
  plot(ContoursLirung_20m,add=T)
  
  barplot(summerMassLoss,beside=T,horiz=T,xlim = c(-3,0.5))
  
  plot(ext_ortho_plot,xlab = '',ylab='', yaxt='n', xaxt='n')
  #plotRGB(ortho13,margins = T,add=T)
  axis(1,las = 1)
  axis(2,las = 3)
  plot(meanWinter *(-1),add=T,zlim=c(-3,1),col=colorr[3:7],ylab='Latitude [m]',xlab='Longitude [m]',legend.args=list( text = '',side=2, font=1, line=0.3, cex=1))
  #plot(cliff_mask_o16,add=T,col='white')
  #plot(pond_mask_o16,add=T,col='blue')
  plot(ContoursLirung_20m,add=T)
  
  barplot(winterMassLoss,beside=T,horiz=T,xlim = c(-3,0.5))
  
  dev.off()
  
#### Get mass loss properties
  library(spatialEco)
  
  DEMm2013 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[1])
  DEMm2013 <- crop(DEMm2013,raster_domain)
  DEMo2013 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[2])
  DEMo2013 <- crop(DEMo2013,raster_domain)
  DEMo2015 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[4])
  DEMo2015 <- crop(DEMo2015,raster_domain)
  DEMm2016 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[5])
  DEMm2016 <- crop(DEMm2016,raster_domain)
  DEMo2016 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[6])
  DEMo2016 <- crop(DEMo2016,raster_domain)
  DEMm2017 <- raster(Lirung_DEMFolder&'\\'&list.files(Lirung_DEMFolder)[7])
  DEMm2017 <- crop(DEMm2017,raster_domain)
  
  terrainprops <- function(DEMM,dhdtRaster){
    DEMM_coarse <- aggregate(DEMM,fact=10)
    
    curv <- curvature(DEMM_coarse,type='mcnab')
    curv <- disaggregate(curv,fact=10)
    curv <- mask(raster::resample(curv,dhdtRaster),glac_mask)   

    slop <- terrain(DEMM_coarse,opt='slope',unit='degrees')
    slop <- disaggregate(slop,fact=10)
    slop <- mask(raster::resample(slop,dhdtRaster),glac_mask)
    
    asp <- terrain(DEMM_coarse,opt='aspect',unit='degrees')
    asp <- disaggregate(asp,fact=10)
    asp <- mask(raster::resample(asp,dhdtRaster),glac_mask)
    
    prom <- terrain(DEMM_coarse,opt='tpi')
    prom <- disaggregate(prom,fact=10)
    prom <- mask(raster::resample(prom,dhdtRaster),glac_mask)
    
    dh <- dhdtRaster
    dh[dh>3] <- NA
    dh[dh< 0] <- NA
#browser()
    stackResult <- stack(curv,slop,asp,prom,dh)
    names(stackResult) <- c('curvature','slope','aspect','promotory','dh')
    return(stackResult)
  }
  
  terrainm2013 <- terrainprops(DEMm2013,stack2013$debrismelt_actual)
  #terraino2013 <- terrainprops(DEMo2013,stack2013$debrismelt_actual)
  #terraino2015 <- terrainprops(DEMo2015,stackm2016$debrismelt_actual)
  terrainm2016 <- terrainprops(DEMm2016,stacko2016$debrismelt_actual)
  #terraino2016 <- terrainprops(DEMo2016,stacko2016$debrismelt_actual)
  terrainm2017 <- terrainprops(DEMm2017,stacko2017$debrismelt_actual)
  
  library('openair')

  range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
  
  visTopovsMassLoss <- function(terrainStack,name){
  terrain_df <- as.data.frame(cbind(range01(terrainStack$curvature[]),terrainStack$slope[],terrainStack$aspect[],range01(terrainStack$promotory[]),terrainStack$dh[]))
  colnames(terrain_df) <- c('curvature','slope','aspect','promotory','dh')

  png(file=paste(path_figs&'\\polarPlot_',name,'_slope_massLoss.png',sep=""), res = 160,width = 1800,height= 1800)
  par(mar=c(3,3,1,4),cex.lab=1.5,cex.axis=1.5)
  polarPlot(terrain_df,x = 'slope',wd = 'aspect',pollutant='dh',force.positive = FALSE,exclude.missing=F,limits = c(0,3),statistic='median',upper=50,min.bin=10,mis.col='white',k=200,annotate = FALSE,key.header = "", key.footer = "",par.settings=list(fontsize=list(text=22),axis.line = list(col = "transparent")))
  dev.off()

  png(file=paste(path_figs&'\\polarPlot_',name,'_curvature_massLoss.png',sep=""), res = 160,width = 1800,height= 1800)
  par(mar=c(3,3,1,4),cex.lab=1.5,cex.axis=1.5)
  polarPlot(terrain_df,x = 'curvature',wd = 'aspect',pollutant='dh',force.positive = FALSE,exclude.missing=F,limits = c(0,3),statistic='median',upper=1,min.bin=10,mis.col='white',k=200,annotate = FALSE,key.header = "", key.footer = "",par.settings=list(fontsize=list(text=22),axis.line = list(col = "transparent")))
  dev.off()
  }
  
  visTopovsMassLoss(terrainm2013,'may2013')
  #visTopovsMassLoss(terraino2013,'oct2013')
  #visTopovsMassLoss(terraino2015,'oct2015')
  visTopovsMassLoss(terrainm2016,'may2016')
  #visTopovsMassLoss(terraino2016,'oct2016')
  visTopovsMassLoss(terrainm2017,'may2017')