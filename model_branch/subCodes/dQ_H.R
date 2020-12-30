################################################################################
# Sensible Heat Flux
# 
# Q_H.R
#
# ReadMe:
# Code to calculate Sensible Heat Flux
# 
# Requires .... Surface Temperature [C] as Input

# Created:          2017/03/17
# Latest Revision:  2017/03/27
#
# Jakob F Steiner | PhD candidate | Faculty of Geosciences | Universiteit Utrecht | 
# Heidelberglaan 2, 3584 CS Utrecht | W.C. van Unnik building | Room 124, Zonneveldvleugel | 
# j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org
#
################################################################################

dQ_H <- function(T_a,T_s,rH,rH_s,u,q_a,p_a,TurbMod){
  #browser()

           
           p_0 <- 101325   # air pressure at level [Pa]
           T_0 <- 288.15   # Temperature at Sea Level
           l_v <- 2.49*10^6 # latent heat of vaporization
           l_s <- 2830000    # Latent Heat of Sublimation
           
           # Sublimation or Evaporation
           l_vs <- l_v + T_a*0
           l_vs[which(T_s < 0)] <- l_s
           
           rho_air <- 1.29                   # air density [kg/m3]
           
           # Saturation Vapour Pressure (Teten's equation)
           e_s_air <- 0.61078 * exp(17.27 * (T_a - 273.15) / (237.3 + T_a - 273.15)) * 1000
           

           
           #Actual Vapour Pressure
           e_a = rH * e_s_air / 100
           
           p_a <- p_a  # p_a to mB
           # Calculate specific humidity of air, assuming actual vapour pressure.
           #browser()
           q_a = 0.622 * e_a / (p_a - (0.378 * e_a)); # all inputs in [Pa]
           

           
           z_0T <- z_0 * 0.05;
           
           A <- k_vk^2 / (log(z_a/z_0))/(log(z_a/z_0T))
           c_p <- c_ad * (1 + 0.84 * q_a)      # specific heat capacity of humid air
           #browser()
           dQ_H <- -1* rho_air * (p_a  / p0) * c_p * A * u
           return(dQ_H)
         }