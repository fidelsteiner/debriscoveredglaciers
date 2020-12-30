require(truncnorm)
require(EnvStats)


# Create Multiple Debrius thickness vectors for different MC Setups

path_data <- 'F:\\PHD\\Research\\EB_DCG\\DistributedEB\\RawData'

deb_thick_MC100 <- hist(rlnormTrunc(10000,meanlog=log(0.84),min=0,max=2.7),breaks=50)

n <- round(deb_thick_MC100$counts/200)
lp <- vector()
kl <-1
for(k in 1:length(deb_thick_MC100$counts)){
  if(n[k]>0){
    lp <- c(lp,rep(deb_thick_MC100$mids[k], n[k]))}
}

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

write.csv(sel_deb_thick_MC50,)