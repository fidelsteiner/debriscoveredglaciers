################################################################################
# Calculting actual view factors
# 
# dependency of vEB.R
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
################################################################################

################################################################################
# 'HOR': CALCULATES HORIZON-FILE FOR EACH POINT          
#
# dfr:                dfr of structure: 'position' 'elevation' 'X'  'Y'
# dem_r:              digital elevation model (raster layer)
# cliff_r:            raster representation of ice cliff
# cliff_poly:         spatial polygon of ice cliff
# resol:		          ORIGINAL resolution of DEM
# newgridsize:        desired gridsize where local horizon has to be considered
#                      (if highresDEM = TRUE) [m]
# ID:                 identifying name, only shown in progress bar window
# simulation_folder:  name of simulation folder
# output:             single txt-file for each coordinate with elevation angle
#                      for every azimuth degree
viewFactors_core_cliff<-function(dfr,dem_r,cliff_r,cliff_poly,resol,
                       newgridsize,ID,simulation_folder,ID_cliffs){

			 endresult<-vector('list')

          for (i in 1:nrow(dfr)){

             # show progress
             cat('cell ' , i,'/',nrow(dfr),'     ',
                 round(resol,1),'m \n')

###########
             DFR_x<-dfr[i,3]
             DFR_y<-dfr[i,4]
             DFR_ele<-dfr[i,2]

             #newminx<-DFR_x-((newgridsize/resol)/2)*resol
			 #newmaxx<-DFR_x+((newgridsize/resol)/2)*resol
			 #newminy<-DFR_y-((newgridsize/resol)/2)*resol
			 #newmaxy<-DFR_y+((newgridsize/resol)/2)*resol
			 dem_cropped<-dem_r #crop(dem_r,extent(newminx,newmaxx,newminy,newmaxy))
             # newgrid DEM
             xy<-c(DFR_x,DFR_y)
             dist_r<-distanceFromPoints(dem_cropped,xy)
             # subtract elevation from point to create base dem
             z_r<-dem_cropped-DFR_ele
             # calculate tangent from point to each cell
             tan_r<-z_r/dist_r
             
             # avoid -Inf for central cell & replace with 0
             tan_r[is.infinite(tan_r)]<-0
             
             # projection has to be defined for terrain analysis
             projection(dist_r)<-'+proj=utm +zone=45N +datum=WGS84'
             # classification of 360 zones, one per degree of view
             zone_r<-terrain(dist_r,opt='aspect',unit='degrees',neighbors=8)
             # create raster stack
             stacked_r<-stack(tan_r,dist_r,z_r)
             names(stacked_r)<-c('tan','d','z')
             # get maximum value (= horizon) for each view direction
             # (older R-version: 'stat=' instead of 'fun=')
             newgridsize_hor<-data.frame(zonal(stacked_r$tan,zone_r,max,
                                               na.rm=TRUE))

          # cliff DEM
             # query polygon which contains specific i-pixel coordinate
             xy_pt<-SpatialPoints(cbind(xy[1],xy[2]),
                            proj4string=CRS('+proj=utm +zone=45N +datum=WGS84'))
             pol<-over(xy_pt,cliff_poly,returnList=FALSE) 
             # ->resulting variable should be a single integer 
             #   & is stored as a data frame
             
             if(typeof(pol) == 'integer'){
                pol<-which(names(ID_cliffs) == pol[1])# if stored as 'Named Int'
                }
             if(typeof(pol) == 'list'){
                pol<-which(ID_cliffs == pol$ID[1])# if stored as a data frame 
                # pol<-which(ID_cliffs == pol[1,1])# if stored as a data frame 
                }      
             pol<-cliff_poly[pol,]
             # make same extent
			 st_r<-crop(stacked_r,extent(cliff_r))
             # mask stacked raster with cliff polygon
             cl_r<-mask(st_r,pol)
             # classification of 360 zones, one per degree of view
             cl_zone_r<-terrain(st_r$d,opt='aspect',unit='degrees',neighbors=8)
             # get maximum value (= horizon) for each view direction
             # (older R-version: 'stat=' instead of 'fun=')
             cl_hor<-data.frame(zonal(cl_r$tan,cl_zone_r,max,na.rm=TRUE))
             # replace -Inf values with NA, otherwise conversion
             #  of -Inf to -90.0 further down while rad->deg converting
             is.na(cl_hor)<-do.call(cbind,lapply(cl_hor,is.infinite))
			 
          # both
             # list both horizons from newgrid and from cliff only
             # ->creating df with:
             #  'zone'  'horangle_newgrid'  'horangle_cliff'
             #    (many NA's in 'horangle_cliff' as some azimuth
             #     directions are not defined due to the lack of pixels)
             hor<-merge(newgridsize_hor,cl_hor,by='zone',all=TRUE)
             colnames(hor)<-c('zone','all_hor','cliff_hor')

             # converting radians to degrees
             # (older R-version: 'hor$max' instead of 'hor$value')
             hor$all_hor<-(atan(hor$all_hor)*180)/pi
             hor$cliff_hor<-(atan(hor$cliff_hor)*180)/pi
             hor<-round(hor,3)

             # use angle from cliff horizon if cliff_hor > all_hor
             #  (might happen rarely due to different zone-arrangement)
             idx<-which(hor$all_hor < hor$cliff_hor)
             hor[idx,2]<-hor[idx,3]

          # find where zones leap azimuth directions & interpolate there
             full.circ<-0:360
             idx_pres<-full.circ[full.circ %in% hor$zone]  #present azimuths
             #idx_miss<-full.circ[!full.circ %in% hor$zone]  #missing azimuths
             full.circ<-data.frame(full.circ)
             names(full.circ)<-'zone'
             # merge ideal circle with present data
             t1<-merge(full.circ,hor,by='zone',all=TRUE)
             # interpolate missing 'hor_all' linearly between present values
             t2<-na.approx(t1[min(idx_pres):max(idx_pres)+1,2])
             # create data frame with interpolated data
             #  (evt. start & end still missing)
             part.circ<-full.circ[min(idx_pres):max(idx_pres)+1,1]
             t3<-data.frame(cbind(part.circ,t2))
             names(t3)<-c('zone','all_hor')
             # merge again ideal circle with partly interpolated data
             #  & add also cliff horizons
             t4<-merge(full.circ,t3,by='zone',all=TRUE)
             t4<-data.frame(cbind(t4,t1$cliff_hor))
             # convert raster package directions (0°=south) 
             #  to 'normal' directions (0°=north)
             t5<-rbind(t4[180:360,],t4[1:179,])
             t6<-data.frame(cbind(1:360,t5[,2:3]))
             colnames(t6)<-c('az','el_all','el_ice')
             # interpolate still missing 'el_all' linearly 
			 #  between present values
			 #  (around 180°/S for northerly facing cliffs) 
             t7<-na.approx(t6$el_all)
             t7<-data.frame(cbind(t6$az,t7,t6$el_ice))
             colnames(t7)<-c('az','el_all','el_ice')
			 t8<-as.matrix(t7)
             # important: replace NA/Inf's with a number!
             t8[is.na(t8)]<- -9999
			 
###TESTING###
# plot(t7$el_all,type='n')
# lines(t7$el_all,lwd=5,col='orange')
# lines(t7$el_ice,lwd=5,col='red')
# plot(dem_cropped)
# plot(cliffs_p,add=T)
# points(x=xy[1],y=xy[2],pch=3,lwd=3)
# plot(atan(tan_r)*180/pi)
# plot(cliffs_p,add=T)
# points(x=xy[1],y=xy[2],pch=3,lwd=1)
#############
			 
             endresult[[i]]<-t8[,2:3]
             }
            close(pb)
            return(endresult)
          }
		  


################################################################################
viewFactors_core<-function(dfr,dem_r,resol,resol2,ID,simulation_folder){
  library(zoo)
			 endresult<-vector('list')
             for (i in 1:nrow(dfr)){
if(is.na(dfr[i,2])){endresult[[i]]<-NA
  
}
               else{
                   # show progress
                   cat('cell ' , i,'/',nrow(dfr),'     ',i/nrow(dfr)*100,'%, ',
                   round(resol2,1),'m \n')
                 


                   DFR_x<-dfr[i,3]
                   DFR_y<-dfr[i,4]

                   xy<-c(DFR_x,DFR_y)
                   dist_r<-distanceFromPoints(dem_r,xy)

                   # subtract elevation from point to create base dem
             # ?
             #####################
                   #z_r<-dem_r-dfr[,2]
                   z_r<-dem_r-dfr[i,2]
             #####################
             # ?
                   # calculate tangent from point to each cell
                   tan_r<-z_r/dist_r
                   
                   # avoid -Inf for central cell & replace with 0
                   tan_r[is.infinite(tan_r)]<-0
             
                   # projection has to be defined for terrain analysis
                   projection(dist_r)<-'+proj=utm +zone=45N +datum=WGS84'
                   # classification of 360 zones, one per degree of view
                   zone_r<-terrain(dist_r,opt='aspect',
                                      unit='degrees',neighbors=8)
                   # create raster stack
                   st_r<-stack(tan_r,dist_r,z_r)
                   names(st_r)<-c('tan','d','z')
                   # get maximum value (= horizon) for each view direction
                   # (older R-version: 'stat=' instead of 'fun=')
                   hor<-data.frame(zonal(st_r$tan,zone_r,max,na.rm=TRUE))
                   colnames(hor)<-c('zone','all_hor')
                   # converting radians to degrees
                   # (older R-version: 'hor$max' instead of 'hor$value')
                   hor$all_hor<-(atan(hor$all_hor)*180)/pi
                   hor<-round(hor,3)  

              # find where zones leap azimuth directions & interpolate there
                   full.circ<-0:360
                   #present azimuths
                   idx_pres<-full.circ[full.circ %in% hor$zone]
                   #missing azimuths
                   #idx_miss<-full.circ[!full.circ %in% hor$zone]
                   full.circ<-data.frame(full.circ)
                   names(full.circ)<-'zone'
                   # merge ideal circle with present data
                   t1<-merge(full.circ,hor,by='zone',all=TRUE)
                   # interpolate missing 'hor_all' linearly 
                   #  between present values
                   t2<-na.approx(t1[min(idx_pres):max(idx_pres)+1,2])
                   # create data frame with interpolated data
                   #  (evt. start & end still missing)
                   part.circ<-full.circ[min(idx_pres):max(idx_pres)+1,1]
                   t3<-data.frame(cbind(part.circ,t2))
                   names(t3)<-c('zone','all_hor')
                   # merge again ideal circle with partly interpolated data
                   t4<-merge(full.circ,t3,by='zone',all=TRUE)
                   # convert raster package directions (0°=south) 
                   #  to 'normal' directions (0°=north)
                   t5<-rbind(t4[180:360,],t4[1:179,])
                   t6<-data.frame(cbind(1:360,t5[,2]))
                   colnames(t6)<-c('az','el_all')
                   # interpolate still missing 'el_all' (around 180°/S) linearly
                   #  between present values
                   t7<-na.approx(t6$el_all)
                   t7<-data.frame(cbind(t6$az,t7))
                   colnames(t7)<-c('az','el_all')
							     t8<-as.matrix(t7)
							     #browser()

								 
                   # keep only horizon value & remove azimuth
							     
							     horizVals <- t8[,2]
							     horizVals[abs(horizVals) == 90] <- NA
							     horizVals[horizVals<0] <- 0
							     Vd <- sum((horizVals)/90,na.rm=T) / 360
							     Vs <- 1 - Vd
                   endresult[[i]]<- Vd
                   
                   
                   #browser()
                   
							     
             }
             }

              return(endresult)
              }

              