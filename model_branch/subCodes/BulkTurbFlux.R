BulkTurbFlux <- function(rho_air,c_p,u,T_a,T_s,rH,z_ins,z_disp,z_0,z_0q,TimStr,Method){
  
  ################################################################################
  # Calculating 
  # 
  # BulkTurbFlux.R
  #
  # ReadMe:
  # Needed as Input:
  # rho_air         Air Density kg/m3
  # 
  # 
  # Created:          2017/01/12
  # Latest Revision:  2017/01/12
  #
  # Jakob F Steiner | PhD candidate | Faculty of Geosciences | Universiteit Utrecht | 
  # Heidelberglaan 2, 3584 CS Utrecht | W.C. van Unnik building | Room 124, Zonneveldvleugel | 
  # j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org
  #
  ################################################################################
  
  # Constant Parameters
  k   <- 0.41   # von Karman constant [-]
  p_a <- 1013   # air pressure at level
  g <- 9.81     # gravitational acceleration [m s-2]
  
  # Saturation Vapour Pressure (Teten's equation)
  e_s_air <- 0.6108 * exp(17.3 * T_a / (237.3 + T_a))
  
  # Saturation Vapour Pressure Surface (Teten's equation)
  e_s_surf <- 0.6108 * exp(17.3 * T_s / (237.3 + T_s))
  
  #Actual Vapour Pressure
  e_a = rH * e_s_air / 100
  
  q_a <- 0.622 * e_a / (p_a - (0.378 * e_a))
  q_s = 0
  e_s = 0.611;
  
  Reid <- function() {                  
    ##### Reid 2010 / All variables calculated according to the Publication

    c_ad <- 1005                        # specific heat capacity of dry air
    c_p <- c_ad * (1 + 0.84 * q_a)      # specific heat capacity of humid air
    
    p_a <- p_0 * ((1 - L * elev) / T_0) ^ (g*M_a/R/L) # Air Pressure
    z_0t <- z_0                         # surface roughness length of heat (assumed equal to z_0)
    z_0q <- z_0                         # surface roughness length of humidity (assumed equal to z_0)
    T_m <- (T_a + T_s) / 2
    R_b <- g*(T_a - T_s)*(z_ins - z_0) / T_m / u^2
    browser()
    psi <- R_b*NA
    psi[R_b>=0&is.na(R_b)==F] <- (1-5*R_b[R_b>=0&is.na(R_b)==F])^2
    psi[R_b<0&is.na(R_b)==F] <- (1-16*R_b[R_b<0&is.na(R_b)==F])^0.75
    psi[is.na(R_b)]<-NA
    Q_H <- rho_air * c_p * k^2 * u * (T_a - T_s) / log(z_ins/z_0) / log(z_ins/z_0t) / psi
    Q_LE <- rho_air * l_v * k^2 * u * (q_a - q_s) / log(z_ins/z_0) / log(z_ins/z_0q) / psi
    return(c(Q_LE, Q_H))
  }
  
  Ayala <-function() {                  
  ##### Ayala 2017
  Q_LE <- rho_air * l_v * k^2 * u * (e_a - e_s_surf) / (log((z_ins-z_disp)/z_0,base=10) - psi_m) / (log((z_ins-z_disp)/z_q,base=10) - psi_q)
  Q_H <- rho_air * c_p * k^2 * u * (T_air - T_surf) / (log((z_ins-z_disp)/z_0,base=10) - psi_m) / (log((z_ins-z_disp)/z_q,base=10) - psi_q)
  
  }
  
  Han <- function() {
  ##### Han 2010
  
  Q_LE <- 0.623 * rho_air / p_a * l_v * k^2 * u * (e_a - e_s_air) * 1000 / (log((z_ins-z_disp)/z_0,base=10) - psi_m) / (log((z_ins-z_disp)/z_q,base=10) - psi_q)
  
  }
  Q_Switch <- switch(Method,"aya" = Ayala(), "han" = Han(), "lit" = Litt(), "rei" = Reid())
  
  
  
  return(c(Q_LE, Q_H))
  
}