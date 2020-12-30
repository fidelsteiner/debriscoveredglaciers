BulkTurbFlux <- function(rho_air,p_a,c_p,l_v,u,T_air,T_surf,rH,z_ins,z_disp,z_0,z_q,psi_m,psi_q,TimStr,Method){
  
  ################################################################################
  # Calculating 
  # 
  # BulkTurbFlux.R
  #
  # ReadMe:
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
  k <- 0.41   # von Karman constant
  
  # Saturation Vapour Pressure (Teten's equation)
  e_s_air <- 0.6108 * exp(17.3 * T_air / (237.3 + T_air))
  
  # Saturation Vapour Pressure Surface (Teten's equation)
  e_s_surf <- 0.6108 * exp(17.3 * T_surf / (237.3 + T_surf))
  
  #Actual Vapour Pressure
  e_a = rH * e_s_air / 100
  
  q_a <- 0.622 * e_a / (p_a - (0.378 * e_a))
  q_s = 0
  e_s = 0.611;
  
  Q_Switch <- switch(Method,"aya" = Ayala(), "han" = Han(), "lit" = Litt())
                     
  
  Ayala <-function() {                  
  ##### Ayala 2017
  Q_LE <- rho_air * l_v * k^2 * u * (e_a - e_s_surf) / (log((z_ins-z_disp)/z_0,base=10) - psi_m) / (log((z_ins-z_disp)/z_q,base=10) - psi_q)
  Q_H <- rho_air * c_p * k^2 * u * (T_air - T_surf) / (log((z_ins-z_disp)/z_0,base=10) - psi_m) / (log((z_ins-z_disp)/z_q,base=10) - psi_q)
  
  }
  
  Han <- function() {
  ##### Han 2010
  
  
  Q_LE <- 0.623 * rho_air / p_a * l_v * k^2 * u * (e_a - e_s_air) * 1000 / (log((z_ins-z_disp)/z_0,base=10) - psi_m) / (log((z_ins-z_disp)/z_q,base=10) - psi_q)
  
  }
  Litt <-function() {                  
    ##### Litt 2015
    
  }
  
  
  return(c(Q_LE, Q_H))
  
}