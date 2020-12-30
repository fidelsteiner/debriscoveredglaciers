################################################################################
# TSaggregate - Aggregate time series data to different time steps
# 
# TSaggregate.R
#
# ReadMe: :
# timStep:  30 - half hourly aggregation
#           60 - hourly aggregation
#           3600 - daily aggregation (needs to be removed!)
#           24 - daily aggregation
# aggrmode: mean,sigma or sum
#
# Created:          2017/02/05
# Latest Revision:  2018/11/17
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht 
# Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 
################################################################################

TSaggregate <- function(TS,TimStr,timStep,yearTS,aggrmode){
  
########### 
  # Aggregate to 30 min timesteps
  if(timStep == 30){
    minVal <- minute(TimStr)
    minVal[minVal<=30] <- 0
    minVal[minVal>=30] <- 30

    # Find unique lines
    timBind<-(cbind(year(TimStr),month(TimStr),lubridate::day(TimStr),hour(TimStr),minVal))
    timBinddf<-data.frame(timBind)
    df2<-unique(timBinddf)
    df2$ID <- 1:nrow(df2)
    timBinddf_unique<-merge(timBinddf,df2)
    DT <- data.table(timBinddf, key="V1,V2,V3,V4,minVal")
    DT[, Cluster_ID:=.GRP, by=key(DT)]
    DT<-as.data.frame(DT)

    # Aggregation of data
    timVec <- ISOdatetime(df2[,1],df2[,2],df2[,3],df2[,4],df2[,5],rep(0,dim(df2)[1]))
    valVec <- aggregate(TS,list(DT[,6]),aggrmode,na.rm=T)
    valVec_sd <- aggregate(TS,list(DT[,6]),sd,na.rm=T)
    valVec$x[is.nan(valVec$x)]<-NA
    TSaggregate <- cbind(timVec ,valVec[,2],valVec_sd[,2])
    
########### 
  # Aggregate to 60 min timesteps
  } else if(timStep == 60){
    minVal <- minute(TimStr)*0
  
    # Find unique lines
    timBind<-(cbind(year(TimStr),month(TimStr),lubridate::day(TimStr),hour(TimStr),minVal))
    timBinddf<-data.frame(timBind)
    df2<-unique(timBinddf)
    df2$ID <- 1:nrow(df2)
    timBinddf_unique<-merge(timBinddf,df2)
    DT <- data.table(timBinddf, key="V1,V2,V3,V4,minVal")
    DT[, Cluster_ID:=.GRP, by=key(DT)]
    DT<-as.data.frame(DT)
    
    # Aggregation of data
    timVec <- ISOdatetime(df2[,1],df2[,2],df2[,3],df2[,4],df2[,5],rep(0,dim(df2)[1]))
    valVec <- aggregate(TS,list(DT[,6]),aggrmode,na.rm=T)
    valVec_sd <- aggregate(TS,list(DT[,6]),sd,na.rm=T)
    valVec$x[is.nan(valVec$x)]<-NA
    TSaggregate <- cbind(timVec ,valVec[,2],valVec_sd[,2])  

###########     
    # Aggregate to 24 h timesteps
  } else if(timStep == 24){
    minVal <- minute(TimStr)*0
    hourVal <- hour(TimStr)*0
    
    # Find unique lines
    timBind<-(cbind(year(TimStr),month(TimStr),lubridate::day(TimStr),hourVal,minVal))
    timBinddf<-data.frame(timBind)
    df2<-unique(timBinddf)
    df2$ID <- 1:nrow(df2)
    timBinddf_unique<-merge(timBinddf,df2)
    DT <- data.table(timBinddf, key="V1,V2,V3,hourVal,minVal")
    DT[, Cluster_ID:=.GRP, by=key(DT)]
    DT<-as.data.frame(DT)
    
    # Aggregation of data
    timVec <- ISOdatetime(df2[,1],df2[,2],df2[,3],df2[,4],df2[,5],rep(0,dim(df2)[1]))
    valVec <- aggregate(TS,list(DT[,6]),aggrmode,na.rm=T)
    valVec_sd <- aggregate(TS,list(DT[,6]),sd,na.rm=T)
    valVec$x[is.nan(valVec$x)]<-NA
    TSaggregate <- cbind(timVec ,valVec[,2],valVec_sd[,2])  
    
###########
    ### Old 24 h Code (to be removed)
    # Aggregate to 24 h timesteps
    } else if(timStep == 3600){
    minVal <- minute(TimStr)*0
    hourVal <- hour(TimStr)*0
    timBind<-(cbind(year(TimStr),month(TimStr),lubridate::day(TimStr),hourVal,minVal))
    timBinddf<-data.frame(timBind)
    df2<-unique(timBinddf)
    df2$ID <- 1:nrow(df2)
    timBinddf_unique<-merge(timBinddf,df2)
    DT <- data.table(timBinddf, key="V1,V2,V3,hourVal,minVal")
    DT[, Cluster_ID:=.GRP, by=key(DT)]
    DT<-as.data.frame(DT)
    timVec <- ISOdatetime(df2[,1],df2[,2],df2[,3],df2[,4],df2[,5],rep(0,dim(df2)[1]))
    valVec <- aggregate(TS,list(DT[,6]),aggrmode,na.rm=T)
    valVec_sd <- aggregate(TS,list(DT[,6]),sd,na.rm=T)
    valVec$x[is.nan(valVec$x)]<-NA
    TSaggregate <- cbind(timVec ,valVec[,2],valVec_sd[,2])    
      }
}