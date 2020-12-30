################################################################################
# Energy Balance Model Code
# 
# EBModCore.R
#
# ReadMe:
# Energy Balance Code for Clean and Debris Covered Glaciers

# Model needs: 
#         - S_downwelling
#         - Albedo and/or S_upwelling (if more than 20% of S_upwelling is NA Albedo is used)
# Part (1): Debris Covered Glaciers
#         : Distributed / Mask for Cliffs/Ponds
# 
# Created:          2017/03/17
# Latest Revision:  2017/03/27
#
# Jakob F Steiner | PhD candidate | Faculty of Geosciences | Universiteit Utrecht | 
# Heidelberglaan 2, 3584 CS Utrecht | W.C. van Unnik building | Room 124, Zonneveldvleugel | 
# j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org
#
################################################################################

EBModCore <- function(EBCoreInput){

  datasize <- length(EBCoreInput$SWin)

  library(lubridate)

  
  # Calculate saturation vapour pressure of air, using Teten's equation (uses air temp in Celsius).
  e_s_air <- vector();
  e_s_air <- 610.8 * exp(17.27 * (EBCoreInput$T_a - 273.15) / (237.3 + EBCoreInput$T_a - 273.15))
  
  # Calculate actual vapour pressure of air, from RH and saturation vapour pressure.
  e_a = vector();
  e_a = EBCoreInput$RH_a * e_s_air / 100;
  
  # Calculate specific humidity of air, assuming actual vapour pressure.
  q_a = 0.622 * e_a / (EBCoreInput$p_a - (0.378 * e_a));
  q_a[q_a<0] = mean(q_a,na.rm=T);
  q_a[is.na(q_a)] <- mean(q_a,na.rm=T);
  
  # Shortwave Radiation
  EBCoreInput$Snet <- EBCoreInput$SWin - EBCoreInput$SWout
# (1) Model DCG (based on original Model by Reid and Brock 2010)

  ##################################
  # DEFINE INITIAL CONDITIONS
  T_s <- vector()
  # Guess 'kickstart' surface temperature equal to initial air temperature
  T_s[1] <- EBCoreInput$T_a[1];
  # Set 'kickstart' debris temperature profile as linear between surface and ice bottom temperatures.
  # Bottom of debris is assumed to remain constant at freezing point T_f. There are N-1 internal ice layers.
  T_d <- matrix(0, nrow = datasize + 1, ncol = N - 1)
  
  for (y in 1:(N-1)){
  
  T_d[1,y] = T_s[1] + (T_f - T_s[1])*y/N;
  }
  k_d <- vector()
  for(t in 1:datasize){
#    if(t%%100==0){
#    print(paste("Progress:", floor(t/datasize*100),"%"))
#    }
    
  # debris conductivity as per season
    
    switch(seasonsDyn,
          on = {                       # Separate conductivity values for each Season
             if (as.numeric(strftime(EBCoreInput$timeline_str[t], format = "%j")) > as.numeric(strftime(paste(year(EBCoreInput$timeline_str[t]),'/',monin,sep=""), format = "%j")) 
                 && as.numeric(strftime(EBCoreInput$timeline_str[t], format = "%j")) < as.numeric(strftime(paste(year(EBCoreInput$timeline_str[t]),'/',monout,sep=""), format = "%j"))) 
               k_d[t] <- k_d_wet
             else
               k_d[t] <- k_d_dry},
          off = {
            k_d[t] <- k_d_wet
          }
  )
  # Make dummy array called Ts for the surface temperature calculation
  Ts <- vector()
  Lout <- vector()
  n <- 2;
  
  # Take first guess as equal to T_s 'kickstart' value
  Ts[n-1] <- T_s[t];
  # Make second value of Ts offset by 0.5
  Ts[n] = T_s[t] - 0.5;
  # 'While' loop runs until successive calculated values are less than 0.01 C apart, or until the limit of 100 iterations.

  while (abs(Ts[n]-Ts[n-1]) > 0.01 && n<100 && !is.na(Ts[n])){
    breakID <- 0
    stop = FALSE
  # Calculate ice temperature profile for the necessary values
  # NOTE: The function 'debristemp.m' must be saved in the same directory as this routine
  Td <- debristemp(Ts[n],T_s[t],T_d[t,],rho_d,c_d,h,timestep,N,k_d[t],T_f);
  Td_plus <- debristemp(Ts[n]+range,T_s[t],T_d[t,],rho_d,c_d,h,timestep,N,k_d[t],T_f);
  Td_minus <- debristemp(Ts[n]-range,T_s[t],T_d[t,],rho_d,c_d,h,timestep,N,k_d[t],T_f);
  
  #% Now do the surface temperature iteration (Newton-Raphson method) calling subfunctions for all fluxes
  
  Ts[n+1] = Ts[n] - (EBCoreInput$Snet[t] + EBCoreInput$LWin[t] - Lup(Ts[n],epsilon_d,sigma_sb) + 
                       Q_G(Ts[n],Td[1],k_d[t],h) + 
                       Q_H(EBCoreInput$T_a[t],Ts[n],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC') +
                       Q_LE(EBCoreInput$T_a[t],Ts[n] + range,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')) / 
                       #+ P(T_a(t),Ts(n),r(t))) 
                  ( (EBCoreInput$Snet[t] + EBCoreInput$LWin[t] - Lup(Ts[n] + range,epsilon_d,sigma_sb) + Q_G(Ts[n] + range,Td_plus[1],k_d[t],h) + 
                       Q_H(EBCoreInput$T_a[t],Ts[n] + range,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC') + 
                       Q_LE(EBCoreInput$T_a[t],Ts[n] + range,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')) #+ P(T_a(t),Ts(n)+range,r(t)))...
                  - (EBCoreInput$Snet[t] + EBCoreInput$LWin[t] - Lup(Ts[n] - range,epsilon_d,sigma_sb) + Q_G(Ts[n] - range,Td_minus[1],k_d[t],h) + 
                       Q_H(EBCoreInput$T_a[t],Ts[n] - range,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC') +
                       Q_LE(EBCoreInput$T_a[t],Ts[n] + range,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')) /
                    #+P(T_a(t),Ts(n)-range,r(t))) 
                    
                        (2*range) );
  
  # Exit loop if one of the fluxes is NA, making calculation of T impossible
  if(is.na(Ts[n+1])){
    breakID <- 1
    stop = TRUE
    break
  }
  if(stop){break}

  # Allow no steps bigger than 1 degree (protect from low-derivative problem) 

  if (Ts[n+1] - Ts[n] > 1){
    Ts[n+1] = Ts[n]+1;
  }
  if (Ts[n+1] - Ts[n] < -1){
    Ts[n+1] = Ts[n]-1;
  }
  n <- n + 1
  }
  
  # Now set the actual array T_s to the estimated value
  if(n == 100){T_s[t + 1] <- (Ts[99] + Ts[100]) / 2
  } else if(is.na(Ts[n])){
    T_s[t + 1] <- NA # This is solely for the case that one of the fluxes was NA and hence temperature calculation was not feasible.
    T_d[t+1,] <- NA
  } else {
    T_s[t + 1] <- Ts[n]}

  T_d[t+1,] <- debristemp(T_s[t+1],T_s[t],T_d[t,],rho_d,c_d,h,timestep,N,k_d[t],T_f);

  }

  # End of Model Loop
    
  # Delete first value in T_s/T_d array because it was just used to kickstart
  T_s <- T_s[-1];
  T_d <- T_d[-1,];
  # Extend debris temperature array to include debris surface and bottom temperatures
  T_f_array <- matrix(T_f, nrow = datasize, ncol = 1);
  T_d <- cbind(T_s,T_d,T_f_array);

  # (1b) Longwave Radiation
  ##################################
  Lup_v <- Lup(T_s, epsilon_d, sigma_sb);
  Lnet <- EBCoreInput$LWin - Lup_v;
  
  # (1c) Turbulent Fluxes
  ##################################
  Q_H_v <- Q_H(EBCoreInput$T_a,T_s,EBCoreInput$RH_a,NULL,EBCoreInput$u,q_a,EBCoreInput$p_a,'nSC')
  Q_LE_v <- Q_LE(EBCoreInput$T_a,T_s,EBCoreInput$RH_a,NULL,EBCoreInput$u,q_a,EBCoreInput$p_a,'nSC')  

  # (1d) Ground Heat Flux
  Q_G <- Q_G(T_s,T_d[,2],k_d,h)
  
  # (1d) Calculate Melt Rate
  ##################################
  #% Calculate melt via conductive flux G_i at bottom of debris (into the ice)
  G_i <- vector();
  melt <- vector();
  # G_i depends on the temperature gradient at the debris base
  G_i <- k_d * (T_d[,N] - T_f) / h; # [W m-2]

  # Convert to metres water eq.
  melt <- G_i * timestep / (rho_i * L_f);
  # Refreezing Potential
  refr <- melt
  refr[refr>0] <- 0

  melt[is.na(melt)] <- 0;
  melt[melt<0] <- 0;
  MeanDailyMelt <- mean(melt,na.rm=T)*(3600*24/timestep);
  
  cumM <- melt
  cumM[is.na(cumM)] <- 0
  CUMmelt <- cumsum(cumM*(3600/timestep));

# =====================
# Model Outputs
# =====================  
  
  EBModCore <- cbind(melt,CUMmelt,T_s,k_d,T_f,h,EBCoreInput$Snet,Lnet,Q_G,Q_H_v,Q_LE_v,T_d,refr)

}