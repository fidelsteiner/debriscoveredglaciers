BulkFlux_Sensible <- function(rho_air,c_p,l_v,u,T_air,T_surf,z_ins,z_disp,z_0,z_q,psi_m,psi_q,TimStr){
  
  # Constant Parameters
  k <- 0.41   # von Karman constant
    
    ##### Ayala 2017
    Q_H_Ayala2017 <- rho_air * c_p * k^2 * u * (T_air - T_surf) / (log((z_ins-z_disp)/z_0,base=10) - psi_m) / (log((z_ins-z_disp)/z_q,base=10) - psi_q)
  

    
    BulkFlux_Sensible <- cbind(Q_H_Ayala2017)
    
  
}