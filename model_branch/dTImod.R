################################################################################
# point scale temperature index models for debris-covered glaciers
# 
# dTImod.R
#
# 
# Created:          2021/02/05
# Latest Revision:  2021/02/05
#
#
# Required Input:
# modType                 - 'dTI' or 'dETI' depending on model (simple TI or enhanced TI)
# Tdata                   - matrix with three columns of same length
#                             1 - numeric time
#                             2 - temperature data in Celsius
#                             3 - off-glacier lapse rate [C/m] (see Heynen et al. 2016)
#                             4 - on-glacier lapse rate [C/m] (see Steiner and Pellicciotti, 2016)
# timFrame                - vector with two values for numeric time of beginning and end of simulation 
#                           (start needs to be >lag time (TIparam[1]) after the beginning of the temperature data) 
# debThick                - debris thickness at pixel [m]
# ElevD                   - 3 value vector 
#                           in case of on-glacier data: c(elevation of measurement,elevation of pixel,0) 
#                           in case of off-glacier data: c(elevation of measurement, lowest elevation of glacier, elevation of pixel) 
# TIparam                 - 5 value vector: c(lag[h], TF1, TF2, SRF1,SRF2); SRF values can be left empty when dTI model is used
# Optional Input for the dETI model:
# SWdata                  - matrix with two columns of same length as Tdata
#                             1 - numeric time
#                             2 - incoming SW data
# alpha                   - one constant albedo value or variable vector of same length as SWdata
#
#
#
# Jakob F Steiner| ICIMOD | jakob.steiner@icimod.org | x-hydrolab.org 
################################################################################




dTImod <- function(modType,Tdata,timFrame,debThick,ElevD,TIparam,SWdata,alpha){

  T_thres <- 0        # temperature threshold for melt on clean ice, assumed for now, should become a variable that can be edited
  TF_clean <- 4.74 / 24    # Temperature factor created from data fro Yala (weighted based on the available days, using only TF>0), based on Litt et al. 2019

  idIN <- which(Tdata[,1]==timFrame[1])
  idOUT <- which(Tdata[,1]==timFrame[2])

  shiftTim <- seq(idIN - round(TIparam[1] *debThick), idOUT - round(TIparam[1] *debThick),1)
  

  #browser()
  # Catch locations with no debris and calculate clean ice melt
  if(debThick<=0.02){
    lapsedT <- Tdata[shiftTim,2] +
      (ElevD[3] - ElevD[1]) * Tdata[shiftTim,3]
    melt_mmwe <- TF_clean * (lapsedT - T_thres)}
    else{
      lapsedT <- Tdata[shiftTim,2] +
        (ElevD[2] - ElevD[1]) * Tdata[shiftTim,3] +
        (ElevD[3] - ElevD[2]) * Tdata[shiftTim,4]
  switch(modType,
         'dTI' = {
    melt_mmwe <- TIparam[2] * debThick^TIparam[3] * lapsedT
      },'dETI' = {
    shiftedSW <- SWdata[shiftTim,2]
    melt_mmwe <- TIparam[2] * debThick^TIparam[3] * lapsedT + TIparam[4] * exp(debThick*TIparam[5]) * shiftedSW * (1-alpha)
      })
    }
  #browser()
  melt_mmwe[melt_mmwe<0] <- 0
  return(cbind(Tdata[idIN:idOUT,1],melt_mmwe,shiftTim))  
}