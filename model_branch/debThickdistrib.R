################################################################################
# debThickdistrib - generate debris thickness values based on basic statistics
# 
# debThickdistrib.R
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

require(truncnorm)
require(EnvStats)


# Create Multiple Debris thickness vectors for different MC Setups

path_data <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData'

deb_thick_MC100 <- hist(rlnormTrunc(10000,meanlog=log(0.84),min=0,max=2.7),breaks=100)

n <- round(deb_thick_MC100$counts/200)
lp <- vector()
kl <-1
for(k in 1:length(deb_thick_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(deb_thick_MC100$mids[k], n[k]))}
}

rtruncnorm <- function(N, mean = 0, sd = 1, a = -Inf, b = Inf) {
  if (a > b) stop('Error: Truncation range is empty');
  U <- runif(N, pnorm(a, mean, sd), pnorm(b, mean, sd));
  rlnorm(U, mean, sd); }

sel_deb_thick_MC50 <- lp

n <- round(deb_thick_MC100$counts/100)
lp <- vector()
kl <-1
for(k in 1:length(deb_thick_MC100$counts)){
  if(n[k]>0){
  lp <- c(lp,rep(deb_thick_MC100$mids[k], n[k]))}
}

sel_deb_thick_MC100 <- lp

n <- round(deb_thick_MC100$counts/50)
lp <- vector()
kl <-1
for(k in 1:length(deb_thick_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(deb_thick_MC100$mids[k], n[k]))}
}

sel_deb_thick_MC200 <- lp

n <- round(deb_thick_MC100$counts/20)
lp <- vector()
kl <-1
for(k in 1:length(deb_thick_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(deb_thick_MC100$mids[k], n[k]))}
}

sel_deb_thick_MC500 <- lp

n <- round(deb_thick_MC100$counts/10)
lp <- vector()
kl <-1
for(k in 1:length(deb_thick_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(deb_thick_MC100$mids[k], n[k]))}
}

sel_deb_thick_MC1000 <- lp

sel_deb_thick_MC100 <- rlnormTrunc(100,meanlog=log(0.84),min=0,max=2.7)
sel_deb_thick_MC500 <- rlnormTrunc(500,meanlog=log(0.84),min=0,max=2.7)
sel_deb_thick_MC50 <- rlnormTrunc(50,meanlog=log(0.84),min=0,max=2.7)
sel_deb_thick_MC1000 <- rlnormTrunc(1000,meanlog=log(0.84),min=0,max=2.7)
sel_deb_thick_MC200 <- rlnormTrunc(200,meanlog=log(0.84),min=0,max=2.7)

write.csv(sel_deb_thick_MC50,paste(path_data,'\\DebrisThickness\\DebrisThickness_lognormal_n50.csv',sep=""))
write.csv(sel_deb_thick_MC100,paste(path_data,'\\DebrisThickness\\DebrisThickness_lognormal_n100.csv',sep=""))
write.csv(sel_deb_thick_MC200,paste(path_data,'\\DebrisThickness\\DebrisThickness_lognormal_n200.csv',sep=""))
write.csv(sel_deb_thick_MC500,paste(path_data,'\\DebrisThickness\\DebrisThickness_lognormal_n500.csv',sep=""))
write.csv(sel_deb_thick_MC1000,paste(path_data,'\\DebrisThickness\\DebrisThickness_lognormal_n1000.csv',sep=""))

