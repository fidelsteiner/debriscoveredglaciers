################################################################################
# Convert mapped features on debris-covered glaciers to pixel values for energy balance model
# 
# featureRasterize_d2EB.R
#
# ReadMe:
# 
#
# Created:          2020/01/05
# Latest Revision:  
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
cl_thres <- 25                                            # threshold from which pixel is classified as cliff (%, relative area cover)
po_thres <- 25                                            # threshold from which pixel is classified as pond (%, relative area cover)
ModRes <- 15                                              # Spatial Resolution of Model [m]
DEM_domain <- '20160430_lirung_dem_20cm.tif'              # DEM for model domain (for lapsing)
cliffFileName_out <- '2018aprcliffsRasterized.tif'
pondFileName_out <- '2017aprpondsRasterized.tif'

############
# optional additional data
############
pond_outlines <- 'Ponds_Apr17.shp'                        # Outlines for pond features (write NULL if not applicable)
cliff_outlines <- 'Cliffs_Apr18.shp'                      # Outlines for cliff features (write NULL if not applicable)

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

# Read in Surface Features with different melt properties
if(!is.null(cliff_outlines)){
  ogrInfo(path_data&'\\Outlines_Features\\'&cliff_outlines)
  cliff_mask<-readOGR(dsn=path_data&'\\Outlines_Features\\'&cliff_outlines)
  cliff_mask@data[which(cliff_mask@data$Id != 1 & cliff_mask@data$Id != 2),'Id'] = 1     # If Id was neither 1 or 2, it is set to 1
  projection(cliff_mask)<-projec
  r_cliffmask_cover <- rasterize(cliff_mask, raster_domain, getCover=TRUE) # relative cover of cliff over pixel
  r_cliffmask <- rasterize(as(cliff_mask, "SpatialLines"), raster_domain,mask=F) # Individual cliff numbers (only pixels that lie on the polygon)
  r_cliffmask[is.na(r_cliffmask)] <- 0
  r_cliffmask_ins <- rasterize(cliff_mask, raster_domain,mask=F) # Individual cliff numbers (pixels inside the polygon)
  r_cliffmask_ins[r_cliffmask>0] <- 0
  r_cliffmask_ins[is.na(r_cliffmask_ins)] <- 0
  r_cliffmask <- r_cliffmask + r_cliffmask_ins                  # combine inside and outside
  r_cliffmask[r_cliffmask_cover < cl_thres/100] <- NA       # if cliff covers less than threshold of pixel, set to debris/ice cover
  r_cliffmask_Ids <- r_cliffmask # later filled with Ids
  if(!is.null(rasterToPolygons(r_cliffmask))){
    IdDF <- (intersect(cliff_mask,rasterToPolygons(r_cliffmask)))
    
    for(i in 1:max(IdDF$layer)){
      r_cliffmask_Ids[r_cliffmask==i] <- IdDF$Id[IdDF$layer==i][1]
    }
    
    cliffStack <- stack(r_cliffmask, r_cliffmask_Ids)   # raster stack for all cliffs where $indID is the number of the cliff and $typeID the identifier for cliff type 
    names(cliffStack) <- c('indID','typeID')
  } else {cliffStack <- NULL}
} else {cliffStack <- NULL}

writeRaster(cliffStack, filename=path_data&'\\Outlines_Features\\'&cliffFileName_out, options="INTERLEAVE=BAND", overwrite=TRUE)

if(!is.null(pond_outlines)){
  ogrInfo(path_data&'\\Outlines_Features\\'&pond_outlines)
  pond_mask<-readOGR(dsn=path_data&'\\Outlines_Features\\'&pond_outlines)
  projection(pond_mask)<-projec
  pond_mask@data[which(pond_mask@data$Id == 7),'Id'] = 3
  pond_mask@data[which(pond_mask@data$Id == 8),'Id'] = 4
  pond_mask@data[which(pond_mask@data$Id != 3 & pond_mask@data$Id != 4),'Id'] = 3     # If Id was neither 3 or 4, it is set to 3
  r_pondmask_cover <- rasterize(pond_mask, raster_domain, getCover=TRUE) # relative cover of cliff over pixel
  r_pondmask <- rasterize(as(pond_mask, "SpatialLines"), raster_domain,mask=F) # Individual pond numbers (only pixels that lie on the polygon)
  r_pondmask[is.na(r_pondmask)] <- 0
  r_pondmask_ins <- rasterize(pond_mask, raster_domain,mask=F) # Individual pond numbers (pixels inside the polygon)
  r_pondmask_ins[r_pondmask>0] <- 0
  r_pondmask_ins[is.na(r_pondmask_ins)] <- 0
  r_pondmask <- r_pondmask + r_pondmask_ins                  # combine inside and outside
  r_pondmask[r_pondmask_cover < po_thres/100] <- NA       # if pond covers less than threshold of pixel, set to debris/ice cover
  r_pondmask_Ids <- r_pondmask # later filled with Ids
  if(!is.null(rasterToPolygons(r_pondmask))){
    IdDF <- (intersect(pond_mask,rasterToPolygons(r_pondmask)))
    
    for(i in 1:max(IdDF$layer)){
      r_pondmask_Ids[r_pondmask==i] <- IdDF$Id[IdDF$layer==i][1]
    }
    
    pondStack <- stack(r_pondmask, r_pondmask_Ids)   # raster stack for all ponds where $indID is the number of the pond and $typeID the identifier for pond type 
    names(pondStack) <- c('indID','typeID')
  } else {pondStack <- NULL}
} else {pondStack <- NULL}

writeRaster(pondStack, filename=path_data&'\\Outlines_Features\\'&pondFileName_out, options="INTERLEAVE=BAND", overwrite=TRUE)
