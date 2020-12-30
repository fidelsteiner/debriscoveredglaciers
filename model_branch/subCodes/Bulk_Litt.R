################################################################################
# Bulk Flux Code based on Litt et al 2015
# Original Version coded by M Litt
# 
# Bulk_Litt.R
#
# ReadMe:
# Input:
#V windspeed in ms-1
#Ta air temperature in degC
#RH air relative humidity in %
#p air barometric pressure in mbars
#Zv,Zt,Zrh : height of the measurements (V,T,RH) in m 
#Ts: temperature of the surface in K
#z0,z0t,z0q: surface roughness lengths in m

# Created:          2017/01/12
# Latest Revision:  2017/03/23
#
# Jakob F Steiner | PhD candidate | Faculty of Geosciences | Universiteit Utrecht | 
# Heidelberglaan 2, 3584 CS Utrecht | W.C. van Unnik building | Room 124, Zonneveldvleugel | 
# j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org
#
################################################################################

Bulk_Litt <- function(V,Ta,RH,rH_s,Ts,p,Zv,Zt,Zrh,z0,z0t,z0q) { 

  # Constants
  g <- 9.81             # gravity acceleration ms-2
  k <- 0.41             # Karman
  LL1<-2830000          # Latent Heat of Sublimation
  LL2<-2.476*10^6          # Latent Heat of Evaporation
  M_a <- 0.02896  # Molar Mass of Air [kg mol-1]
  R <- 8.31       # Universal Gas Constant
  p_0 <- 101325   # air pressure at level [Pa]
  T_0 <- 288.15   # Temperature at Sea Level
  L <- 0.0065     # Temperature Lapse Rate
  #browser()
#choice of stability functions set  1=monin-obukhov, 2=brutsaert
func<-2

# Potential Temperatures
Tapot = (Ta + 273.15) * ((1013 / (p )) ^ 0.2869)-273.15
Tspot = (Ts + 273.15) * ((1013 / (p )) ^ 0.2869)-273.15

# Sublimation or Evaporation
LL <- LL1
LL[which(Ts >= 0)] <- LL2

# air density (Arrhenius substituded with rho at sea level)
#rho <- p * M_a / R / Ta * 100
#rho <- 1.29*p/1024
p <- p
#p <- p_0* ((1-(L*elev/T_0))^(g*M_a/(R*L))) / 100 # Air Pressure
rho <- (p * 100) * M_a / R / (Ta + 273.15)
rho <- 1.29 * p / (p_0/100)
#browser()
#RH extrapolated to the height of V measurements, assuming log profile (neutral) 
#alpha<-(RH-100)/(log(Zrh/z0q))
#RHz<-100+alpha*log(Zt/z0q)
if (!is.null(rH_s)) {
  rH_surf <- rH_s
} else {
  rH_surf <- 100
}

RH[is.na(RH)]<-100
RHz <- RH
#rH_surf <- 0;
# partial pressure of water in the air [mbar]
e <- (RHz/100)*6.1078*exp(17.27*(Ta)/(237.3+Ta))

# partial pressure of water at the surface [mbar]
es <- rH_surf/100*6.1078*exp(17.27*(Ts)/(237.3 + Ts))

q <- e * 0.622 / (p - (0.378 * e))    # specific humidity of air [kg/kg]
qs <- es * 0.622 / (p - (0.378 * es))  # specific humidity of surface [kg/kg]

# specific heat capacity of air
Cp = 1006 * (1 + 0.84 * q)

# Virtual Temperatures for Derivation of Monin Obukhov Length
Ta_v = (Ta + 273.15) * (1 + 0.61 * q)
Ts_v = (Ts + 273.15) * (1 + 0.61 * qs)

Tav_pot = (Ta_v) * ((1013 / p ) ^ 0.2869)-273.15
Tsv_pot = (Ts_v) * ((1013 / p ) ^ 0.2869)-273.15

# Flux Calculations
zsurL <- NA * V       #z/L stability Variable
if (Zv>z0 & Zt>z0t & Zrh>z0q & V>0 & is.na(V)==FALSE)
{flag<-1} else {flag<-0}

if (flag==1) {       
# Initial Values for z*,t*,q*,L*
 ustar_i <- k*V/log(Zv/z0)
 tstar_i <- k*(Ta-Ts)/log(Zt/z0t)
 qstar_i <- k*(q-qs)/log(Zrh/z0q)

 Lstar<-(Ta_v*ustar_i^2)/(k*g*(tstar_i+0.61*qstar_i*(Ta + 273.15))) *(-1)
 show(Lstar)
 zsurL <- Zv/Lstar} else { zsurL<- NA }

###### Core Calculations
if (flag==1 & is.na(zsurL)!=TRUE) {ite<-1
             delta<-1000
while (ite<51 & delta>0.01) {
  # Stability functions (Brutsaert, 1982)
  # 
  psi <- function(zsurL,tstab){    # Stability Correction for Momentum
    x <- (1-16*zsurL)^(1/4) 
    if (zsurL<0&is.na(zsurL)!=TRUE) {
      psimunstabl <- 2*log((1+x)/2)+log((1+x^2)/2)-2*atan(x)+pi/2
      psihunstabl <- 2*log((1+x^2)/2)
      if (tstab == 1){return(psimunstabl)} 
      if (tstab == 2){return(psihunstabl)}} 
    if (zsurL>0&zsurL<1&is.na(zsurL)!=TRUE) {
      psimstabl2<- -5*zsurL
      psihstabl2<- -5*zsurL
      if (tstab == 1){return(psimstabl2)} 
      if (tstab == 2){return(psihstabl2)}} 
    if (zsurL>1&is.na(zsurL)!=TRUE)  {
      psimstabl2<- -5*(log(zsurL)+1)
      psihstabl2<- -5*(log(zsurL)+1)
      if (tstab == 1){return(psimstabl2)} 
      if (tstab == 2){return(psihstabl2)}} 
  }

ustar <- k*(V)/(log(Zv/z0)-psi(Zv/Lstar,1)+psi(z0/Lstar,1))  
tstar <- k*(Ta-Ts)/(log(Zt/z0t)-psi(Zt/Lstar,2)+psi(z0t/Lstar,2))
qstar <- k*(q-qs)/(log(Zrh/z0q)-psi(Zrh/Lstar,2)+psi(z0q/Lstar,2))
Lstarold <- Lstar
Lstar <- (Ta_v*ustar^2) / (k*g*(tstar+0.61*qstar*(Ta + 273.15)))
delta <- abs((Lstar-Lstarold))
ite<-ite+1

}
            
# Flux Calculation
H <- rho * Cp * ustar * tstar
if(H<=-500){H<-NA}
LE <- LL * rho * ustar * qstar
if(LE<=-500){LE<-NA}
Lstar <- Lstar;
}
else { 
 # browser()
H<- NA
LE<- NA
ustar<- NA
tstar<- NA
qstar<- NA
zsurL<- NA
Lstar<-NA
}

 return(cbind(H,LE,ustar,tstar,qstar,zsurL,Lstar,rH_s,rH_surf,qs))
}
