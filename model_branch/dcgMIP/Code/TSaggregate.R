TSaggregate <- function(TS,TimStr,timStep,yearTS,aggrmode){

  
  ## CLEAN CODE!!
  if(timStep == 30){
minVal <- minute(TimStr)
minVal[minVal<=30] <- 0
minVal[minVal>=30] <- 30

timBind<-(cbind(month(TimStr),lubridate::day(TimStr),hour(TimStr),minVal))
timBinddf<-data.frame(timBind)
df2<-unique(timBinddf)
df2$ID <- 1:nrow(df2)
timBinddf_unique<-merge(timBinddf,df2)
DT <- data.table(timBinddf, key="V1,V2,V3,minVal")
DT[, Cluster_ID:=.GRP, by=key(DT)]
DT<-as.data.frame(DT)

timVec <- ISOdatetime(rep(yearTS,dim(df2)[1]),df2[,1],df2[,2],df2[,3],df2[,4],rep(0,dim(df2)[1]))
valVec <- aggregate(TS,list(DT[,5]),aggrmode,na.rm=T)
valVec_sd <- aggregate(TS,list(DT[,5]),sd,na.rm=T)
valVec$x[is.nan(valVec$x)]<-NA
#timVec[is.nan(timVec)] <- NA
#valVec[is.nan(valVec)] <- NA
#valVec_sd[is.nan(valVec_sd)] <- NA
TSaggregate <- cbind(timVec ,valVec[,2],valVec_sd[,2])

  } else if(timStep == 60){
    minVal <- minute(TimStr)*0
    
    timBind<-(cbind(month(TimStr),lubridate::day(TimStr),hour(TimStr),minVal))
    timBinddf<-data.frame(timBind)
    df2<-unique(timBinddf)
    df2$ID <- 1:nrow(df2)
    timBinddf_unique<-merge(timBinddf,df2)
    DT <- data.table(timBinddf, key="V1,V2,V3,minVal")
    DT[, Cluster_ID:=.GRP, by=key(DT)]
    DT<-as.data.frame(DT)
    
    timVec <- ISOdatetime(rep(yearTS,dim(df2)[1]),df2[,1],df2[,2],df2[,3],df2[,4],rep(0,dim(df2)[1]))
    valVec <- aggregate(TS,list(DT[,5]),aggrmode,na.rm=T)
    valVec_sd <- aggregate(TS,list(DT[,5]),sd,na.rm=T)
    valVec$x[is.nan(valVec$x)]<-NA
    #timVec[is.nan(timVec)] <- NA
    #valVec[is.nan(valVec)] <- NA
    #valVec_sd[is.nan(valVec_sd)] <- NA
    TSaggregate <- cbind(timVec ,valVec[,2],valVec_sd[,2])    
  } else if(timStep == 3600){
    minVal <- minute(TimStr)*0
    hourVal <- hour(TimStr)*0
    #browser()
    timBind<-(cbind(year(TimStr),month(TimStr),lubridate::day(TimStr),hourVal,minVal))
    timBinddf<-data.frame(timBind)
    df2<-unique(timBinddf)
    df2$ID <- 1:nrow(df2)
    timBinddf_unique<-merge(timBinddf,df2)
    DT <- data.table(timBinddf, key="V1,V2,V3,hourVal,minVal")
    DT[, Cluster_ID:=.GRP, by=key(DT)]
    DT<-as.data.frame(DT)
#browser()
    timVec <- ISOdatetime(df2[,1],df2[,2],df2[,3],df2[,4],df2[,5],rep(0,dim(df2)[1]))
    valVec <- aggregate(TS,list(DT[,6]),aggrmode,na.rm=T)
    valVec_sd <- aggregate(TS,list(DT[,6]),sd,na.rm=T)
    valVec$x[is.nan(valVec$x)]<-NA
    #timVec[is.nan(timVec)] <- NA
    #valVec[is.nan(valVec)] <- NA
    #valVec_sd[is.nan(valVec_sd)] <- NA
    TSaggregate <- cbind(timVec ,valVec[,2],valVec_sd[,2])    
    #browser()
  }
}

