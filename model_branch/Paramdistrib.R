################################################################################
# Paramdistrib - generate parameter ranges for the d2EB Model
# 
# Paramdistrib.R
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




# Surface Roughness
z0_MC100 <- hist(rlnormTrunc(10000,meanlog=log(0.03),min=0.005,max=0.5),breaks=50)

n <- round(z0_MC100$counts/200)
lp <- vector()
kl <-1
for(k in 1:length(z0_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(z0_MC100$mids[k], n[k]))}
}

z0_MC50 <- lp

n <- round(z0_MC100$counts/100)
lp <- vector()
kl <-1
for(k in 1:length(z0_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(z0_MC100$mids[k], n[k]))}
}

z0_MC100_ <- lp

n <- round(z0_MC100$counts/50)
lp <- vector()
kl <-1
for(k in 1:length(z0_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(z0_MC100$mids[k], n[k]))}
}

z0_MC200 <- lp

n <- round(z0_MC100$counts/20)
lp <- vector()
kl <-1
for(k in 1:length(z0_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(z0_MC100$mids[k], n[k]))}
}

z0_MC500 <- lp

n <- round(z0_MC100$counts/10)
lp <- vector()
kl <-1
for(k in 1:length(z0_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(z0_MC100$mids[k], n[k]))}
}

z0_MC1000 <- lp

z0_MC50 <- rlnormTrunc(50,meanlog=log(0.03),min=0.005,max=0.5)
z0_MC100 <- rlnormTrunc(100,meanlog=log(0.03),min=0.005,max=0.5)
z0_MC200 <- rlnormTrunc(200,meanlog=log(0.03),min=0.005,max=0.5)
z0_MC500 <- rlnormTrunc(500,meanlog=log(0.03),min=0.005,max=0.5)
z0_MC1000 <- rlnormTrunc(1000,meanlog=log(0.03),min=0.005,max=0.5)
write.csv(sample(z0_MC50),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\z0_lognormal_n50.csv')
write.csv(sample(z0_MC100),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\z0_lognormal_n100.csv')
write.csv(sample(z0_MC200),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\z0_lognormal_n200.csv')
write.csv(sample(z0_MC500),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\z0_lognormal_n500.csv')
write.csv(sample(z0_MC1000),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\z0_lognormal_n1000.csv')

# Debris porosity

por_MC100 <- rtruncnorm(n = 100, a = 0.44-0.14, b = 0.44+0.14, mean = 0.44, sd = 0.07)
por_MC50 <- rtruncnorm(n = 50, a = 0.44-0.14, b = 0.44+0.14, mean = 0.44, sd = 0.07)
por_MC200 <- rtruncnorm(n = 200, a = 0.44-0.14, b = 0.44+0.14, mean = 0.44, sd = 0.07)
por_MC500 <- rtruncnorm(n = 500, a = 0.44-0.14, b = 0.44+0.14, mean = 0.44, sd = 0.07)
por_MC1000 <- rtruncnorm(n = 1000, a = 0.44-0.14, b = 0.44+0.14, mean = 0.44, sd = 0.07)

write.csv(sample(por_MC50),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\por_lognormal_n50.csv')
write.csv(sample(por_MC100),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\por_lognormal_n100.csv')
write.csv(sample(por_MC200),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\por_lognormal_n200.csv')
write.csv(sample(por_MC500),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\por_lognormal_n500.csv')
write.csv(sample(por_MC1000),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\por_lognormal_n1000.csv')

# Debris density
den_MC100 <- rtruncnorm(n = 100, a = 1588 - 175*2, b = 1588 + 175*2, mean = 1588, sd = 175)
den_MC50 <- rtruncnorm(n = 50, a = 1588 - 175*2, b = 1588 + 175*2, mean = 1588, sd = 175)
den_MC200 <- rtruncnorm(n = 200, a = 1588 - 175*2, b = 1588 + 175*2, mean = 1588, sd = 175)
den_MC500 <- rtruncnorm(n = 500, a = 1588 - 175*2, b = 1588 + 175*2, mean = 1588, sd = 175)
den_MC1000 <- rtruncnorm(n = 1000, a = 1588 - 175*2, b = 1588 + 175*2, mean = 1588, sd = 175)

write.csv(sample(den_MC50),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\den_lognormal_n50.csv')
write.csv(sample(den_MC100),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\den_lognormal_n100.csv')
write.csv(sample(den_MC200),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\den_lognormal_n200.csv')
write.csv(sample(den_MC500),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\den_lognormal_n500.csv')
write.csv(sample(den_MC1000),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\den_lognormal_n1000.csv')

# Debris conductivity
DebrisDataelse <- read.csv('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData\\DebrisData\\AggregatedDebrisDataRanges.csv')

meanBulk_rho <- mean(DebrisDataelse$bulk.density..kg.m3.,na.rm=T)
sdBulk_rho <- sd(DebrisDataelse$bulk.density..kg.m3.,na.rm=T)

meanmoist <- mean(DebrisDataelse$volumetric.soil.moisture..m3.m3.,na.rm=T)
sdmoist <- sd(DebrisDataelse$volumetric.soil.moisture..m3.m3.,na.rm=T)

meanpor <- mean(DebrisDataelse$calculated.porosity....,na.rm=T)
sdpor <- sd(DebrisDataelse$calculated.porosity....,na.rm=T)

# convert from diffusivity to conductivity
phi_deb <- meanpor
phiDeb_lit <- c(0.2,0.3,0.43,0.33,0.3)
c_deb <- 948
rho_d <- meanBulk_rho       # 19 samples 

rho_rock <- 2650# kg m-3
c_rock <- c(750,948,804,900,890,948) # J kg-1 K-1

c_w <- 4181.3 
c_air <- 1003.5
rho_w <<- 999.7
rho_air <<- 0.819

volum_heat_cap_void <- rho_w * c_w * meanmoist/0.2 + rho_air * c_air * (1-meanmoist/0.2)

volum_heat_cap_wet <- rho_w * c_w
volum_heat_cap_dry <- rho_air * c_air

k_diff_lit_mean <- c(0.3,0.38,0.3536667,0.95,0.689,0.6,0.9,0.82,0.89,0.82,0.6,0.59,0.42,0.61) / 10^6 # m2 s-1
k_diff_lit_sd <- c(0.05,0.02,0.01,0.01,0.1,0.1,0.01,0.01,0.01,0.01,0.03,0.06) / 10^6

cond_lit_dry <- c(0.585,0.637,1.29)
cond_lit_wet <- c(1.669,1.776)
cond_lit <- c(0.94,0.96,1.33,1.62,1.523,1.29,1.28,0.85)
cond_tot_distrib <- k_diff_lit_mean * (rho_rock * 890 * (1-phi_deb) + volum_heat_cap_void * phi_deb)
cond_tot_phi0.3 <- k_diff_lit_mean * (rho_rock * 890 * (1-0.3) + volum_heat_cap_void * 0.3)

cond_tot_Lirung_1 <- mean(0.94,0.74)/10^6 * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)
cond_tot_Lirung_2 <- mean(1.07,0.6)/10^6 * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)
cond_tot_Lirung_3 <- mean(0.51)/10^6 * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)

cond_tot_wet <- k_diff_lit_mean * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_wet * phi_deb)
cond_tot_dry <- k_diff_lit_mean * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_dry * phi_deb)
cond_tot <- k_diff_lit_mean * (rho_rock * 890* (1-phi_deb) + volum_heat_cap_void * phi_deb)

rgbeta <- function(n, mean, var, min = 0, max = 1)
{
  dmin <- mean - min
  dmax <- max - mean
  
  if (dmin <= 0 || dmax <= 0)
  {
    stop(paste("mean must be between min =", min, "and max =", max)) 
  }
  
  if (var >= dmin * dmax)
  {
    stop(paste("var must be less than (mean - min) * (max - mean) =", dmin * dmax))
  }
  
  # mean and variance of the standard beta distributed variable
  mx <- (mean - min) / (max - min)
  vx <- var / (max - min)^2
  
  # find the corresponding alpha-beta parameterization
  a <- ((1 - mx) / vx - 1 / mx) * mx^2
  b <- a * (1 / mx - 1)
  
  # generate standard beta observations and transform
  x <- rbeta(n, a, b)
  y <- (max - min) * x + min
  
  return(y)
}

num <- c(50,100,200,500,1000)
for(k in 1:5){
cond_distrib_range_wet <- rtruncnorm(num[k],a = min(cond_tot_wet), b = max(cond_tot_wet),mean = mean(cond_tot_wet),sd = sd(cond_tot_wet))
cond_distrib_range_dry <- rtruncnorm(num[k],mean = mean(cond_tot_dry),sd = sd(cond_tot_dry),a=min(cond_tot_dry),b=max(cond_tot_dry))
cond_distrib_range <- rtruncnorm(num[k],mean = mean(cond_tot),sd = sd(cond_tot),a=min(cond_tot),b=max(cond_tot))
cond_Lirung <- as.data.frame(cbind(cond_distrib_range_dry,cond_distrib_range_wet,cond_distrib_range))
colnames(cond_Lirung) <- c('conductivity_dry','conductivity_wet','conductivity_tot')
write.csv(cond_Lirung,paste('F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\debris_conductivity_n',as.character(num[k]),'.csv',sep=""))
}

# Debris Albedo
alpha_MC100 <- rtruncnorm(n = 100, a = 0.06, b = 0.36, mean = 0.13, sd = 0.03)
alpha_MC50 <- rtruncnorm(n = 50,  a = 0.06, b = 0.36, mean = 0.13, sd = 0.03)
alpha_MC200 <- rtruncnorm(n = 200,  a = 0.06, b = 0.36, mean = 0.13, sd = 0.03)
alpha_MC500 <- rtruncnorm(n = 500,  a = 0.06, b = 0.36, mean = 0.13, sd = 0.03)
alpha_MC1000 <- rtruncnorm(n = 1000,  a = 0.06, b = 0.36, mean = 0.13, sd = 0.03)

write.csv(sample(alpha_MC50),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\alpha_n50.csv')
write.csv(sample(alpha_MC100),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\alpha_n100.csv')
write.csv(sample(alpha_MC200),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\alpha_n200.csv')
write.csv(sample(alpha_MC500),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\alpha_n500.csv')
write.csv(sample(alpha_MC1000),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\alpha_n1000.csv')

# Climate uncertainty ranges
TAir_MC100 <- runif(100, -2, 2)
TAir_MC50 <- runif(50, -2, 2)
TAir_MC200 <- runif(200, -2, 2)
TAir_MC500 <- runif(500, -2, 2)
TAir_MC1000 <- runif(1000, -2, 2)

write.csv(sample(TAir_MC50),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\TAir_n50.csv')
write.csv(sample(TAir_MC100),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\TAir_n100.csv')
write.csv(sample(TAir_MC200),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\TAir_n200.csv')
write.csv(sample(TAir_MC500),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\TAir_n500.csv')
write.csv(sample(TAir_MC1000),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\TAir_n1000.csv')

RH_MC100 <- runif(100, -20, 20)
RH_MC50 <- runif(50, -20, 20)
RH_MC200 <- runif(200, -20, 20)
RH_MC500 <- runif(500, -20, 20)
RH_MC1000 <- runif(1000, -20, 20)

write.csv(sample(RH_MC50),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\RH_n50.csv')
write.csv(sample(RH_MC100),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\RH_n100.csv')
write.csv(sample(RH_MC200),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\RH_n200.csv')
write.csv(sample(RH_MC500),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\RH_n500.csv')
write.csv(sample(RH_MC1000),'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\RH_n1000.csv')

ws_MC100 <- runif(100, 1/2.5, 2.5)
ws_MC50 <- runif(50, 1/2.5, 2.5)
ws_MC200 <- runif(200, 1/2.5, 2.5)
ws_MC500 <- runif(500, 1/2.5, 2.5)
ws_MC1000 <- runif(1000, 1/2.5, 2.5)

write.csv(ws_MC50,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\ws_n50.csv')
write.csv(ws_MC100,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\ws_n100.csv')
write.csv(ws_MC200,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\ws_n200.csv')
write.csv(ws_MC500,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\ws_n500.csv')
write.csv(ws_MC1000,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\ws_n1000.csv')

LW_MC100 <- runif(100, -50, 50)
LW_MC50 <- runif(50, -50, 50)
LW_MC200 <- runif(200, -50, 50)
LW_MC500 <- runif(500, -50, 50)
LW_MC1000 <- runif(1000, -50, 50)

write.csv(LW_MC50,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\LW_n50.csv')
write.csv(LW_MC100,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\LW_n100.csv')
write.csv(LW_MC200,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\LW_n200.csv')
write.csv(LW_MC500,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\LW_n500.csv')
write.csv(LW_MC1000,'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\MCRanges\\LW_n1000.csv')