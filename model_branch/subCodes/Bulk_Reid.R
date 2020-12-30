################################################################################
# Bulk Flux Code based on Reid et al 2010
# 
# Bulk_Reid.R
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

Bulk_Reid <- function(T_a, T_s, rH, rH_s, u, p_a_meas, z_0, z_q, z_insU,z_insT,elev,timeVec) {                  

  k   <- 0.41     # von Karman constant [-]
  p_0 <- 101325   # air pressure at level [Pa]
  T_0 <- 288.15   # Temperature at Sea Level
  g <- 9.81       # gravitational acceleration [m s-2]
  M_a <- 0.02896  # Molar Mass of Air [kg mol-1]
  R <- 8.31       # Universal Gas Constant
  L <- 0.0065     # Temperature Lapse Rate
  l_v <- 2.476*10^6 # latent heat of vaporization
  l_s <- 2830000    # Latent Heat of Sublimation
  
  # Sublimation or Evaporation
  l_vs <- l_v + T_s*0
  l_vs[which(T_s < 0)] <- l_s

  c_ad <- 1006                        # specific heat capacity of dry air
  
  #p_a <- p_0* ((1-(L*elev/T_0))^(g*M_a/(R*L))) # Air Pressure
  p_a <- p_a_meas;
  rho_air <- p_a  * M_a / R / (T_a + 273.15)      # air density

  z_0t <- z_q                         # surface roughness length of heat (assumed equal to z_0)
  z_0q <- z_q                         # surface roughness length of humidity (assumed equal to z_0)

  #browser()
  # Saturation Vapour Pressure (Teten's equation)
  e_s_air <- 0.61078 * exp(17.27 * T_a / (237.3 + T_a)) * 1000
  
  # Saturation Vapour Pressure Surface (Teten's equation)
  es_s <- 0.61078 * exp(17.27 * T_s / (237.3 + T_s)) * 1000
  
  #Actual Vapour Pressure
  e_a = rH * e_s_air / 100

  # Calculate specific humidity of air, assuming actual vapour pressure.
  q_a = 0.622 * e_a / (p_a - (0.378 * e_a));

  c_p <- c_ad * (1 + 0.84 * q_a)      # specific heat capacity of humid air
  
  # Must calculate q at surface from surface temp
  # Calculate surface saturation vapour pressure and hence specific humidity.
  #es_s <- 0.61078 * exp((17.27 * (T_s) / 237.3 + T_s)); #[Pa]

  if (!is.null(rH_s)) {
    rH_surf <- rH_s
  } else {
      rH_surf <- 100
  }
  es <- rH_surf * es_s / 100;
  q_s <- 0.622 * es / (p_a - (0.378 * es));

  T_m <- (T_a + T_s) / 2 + 273.15
  R_b <- g*(T_a - T_s)*(z_insT - z_0) / T_m / u^2
  R_b[abs(R_b)>100]<-NA
  phi <- R_b*NA
  R_b[R_b>=1/5] <- NA
  phi[R_b>=0&!is.na(R_b)] <- (1-5*R_b[R_b>=0&!is.na(R_b)])^2
  phi[R_b<0&!is.na(R_b)] <- (1-16*R_b[R_b<0&!is.na(R_b)])^0.75
  phi[is.na(R_b)]<-NA
  #phi[abs(R_b)<0.25] <- 1

  Q_H <- rho_air * c_p * k^2 * u * (T_a - T_s) / log(z_insU/z_0) / log(z_insT/z_0t) * phi
  #Q_H <- rho_air * c_p * k^2 * u * (T_a - T_s) / log(z_insU/z_0) / log(z_insT/z_0t)
  Q_H[abs(Q_H)>500]<-NA
  #browser()
  Q_LE <- rho_air * l_vs * k^2 * u * (q_a - q_s) / log(z_insU/z_0) / log(z_insT/z_0q) * phi
  Q_LE[abs(Q_LE)>500]<-NA
  return(cbind(timeVec,Q_H, Q_LE,R_b,phi))
}