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

Q_H <- function(T_a,T_s,rH,rH_s,u,q_a,p_a,TurbMod){
  #browser()
  switch(TurbMod,
         'nSC' = {
           
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
           
           # Saturation Vapour Pressure Surface (Teten's equation)
           es_s <- 0.61078 * exp(17.27 * (T_s - 273.15) / (237.3 + T_s - 273.15)) * 1000
           
           #Actual Vapour Pressure
           e_a = rH * e_s_air / 100
           
           p_a <- p_a * 100
           # Calculate specific humidity of air, assuming actual vapour pressure.
           q_a = 0.622 * e_a / (p_a - (0.378 * e_a));
           
           # Must calculate q at surface from surface temp
           # Calculate surface saturation vapour pressure and hence specific humidity.
           
           if (!is.null(rH_s)) {
             rH_surf <- rH_s
           } else {
             rH_surf <- 100
           }
           
           es <- rH_surf * es_s / 100;
           q_s <- 0.622 * es / (p_a - (0.378 * es));
           
           z_0T <- z_0 * 0.05;
           
           A <- k_vk^2 / (log(z_a/z_0))/(log(z_a/z_0T))
           c_p <- c_ad * (1 + 0.84 * q_a)      # specific heat capacity of humid air
           
           Q_H <- rho_air * (p_a / p0) * c_p * A * u * (T_a - T_s)
           
           Q_H[abs(Q_H)>500]<-500*(-1)

           return(Q_H)
         },
         
         'rSC' = {
           # Richardson Number
           R_b <- grav * (T_a - T_s_data) * (z_a - z_0) / (0.5*(T_a + T_s_data)) / u^2;
           R_b[is.infinite(R_b)&R_b<0] <- 0
           R_b[is.infinite(R_b)&R_b>0] <- 0
           R_b[is.na(R_b)] <- 0
           R_b[u=0] <- 0
           
           phi_stab <- vector()
           phi_stab[R_b>=0] <- (1 - 5*R_b[R_b>=0])^2
           phi_stab[R_b<0] <- (1 - 16*R_b[R_b<0])^0.75
           
           # Density of Air / Ideal Gas Law
           rho_a = p_a*Mair/Rgas/(T_a)
           
           z_0m <- z_0;
           z_0q <- z_0;
           
           c_a <- c_ad*(1+0.84*q_a)
           
           # Sensible Heat Flux
           Fl <- rho_a * c_a * k_vk^2 * u *(T_a - T_s_data) *phi_stab/ log(z_a/z_0m) / log(z_a/z_0q)# * phi_stab
           Fl[R_b=0] <- NA
           Fl[u=0] <- NA;
           
           Q_H <- Fl;
           
         })
  
  }
