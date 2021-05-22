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

#  library(lubridate)

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

  switch(EBCoreInput$mod,
         'debris' = {
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
          on = {                      # Separate conductivity values for each Season
             if (as.numeric(strftime(EBCoreInput$timeline_str[t], format = "%j")) > as.numeric(strftime(paste(year(EBCoreInput$timeline_str[t]),'/',monin,sep=""), format = "%j")) 
                 && as.numeric(strftime(EBCoreInput$timeline_str[t], format = "%j")) < as.numeric(strftime(paste(year(EBCoreInput$timeline_str[t]),'/',monout,sep=""), format = "%j"))) 
               k_d[t] <- k_d_wet
             else
               k_d[t] <- k_d_dry},
          off = {
            k_d[t] <- k_d_wet
          }
  )

    switch(EBModType,
           meas_TS = {T_s[t+1] <- EBCoreInput$T_s_data[t]
           #browser()
           },
           no_TS = {
             
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
    #browser()
  Td <- debristemp(Ts[n],T_s[t],T_d[t,],rho_d,c_d,h,timestep,N,k_d[t],T_f);
  Td_plus <- debristemp(Ts[n]+CDrange,T_s[t],T_d[t,],rho_d,c_d,h,timestep,N,k_d[t],T_f);
  Td_minus <- debristemp(Ts[n]-CDrange,T_s[t],T_d[t,],rho_d,c_d,h,timestep,N,k_d[t],T_f);
  
  #% Now do the surface temperature iteration (Newton-Raphson method) calling subfunctions for all fluxes
  #browser()
  Ts[n+1] = Ts[n] - (EBCoreInput$Snet[t] + EBCoreInput$LWin[t] - Lup(Ts[n],epsilon_d,sigma_sb) + 
                       Q_G(Ts[n],Td[1],k_d[t],h) + 
                       Q_H(EBCoreInput$T_a[t],Ts[n],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC') +
                       Q_LE(EBCoreInput$T_a[t],Ts[n],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')) / 
                       #+ P(T_a(t),Ts(n),r(t))) 
                  ( (EBCoreInput$Snet[t] + EBCoreInput$LWin[t] - Lup(Ts[n] + CDrange,epsilon_d,sigma_sb) + Q_G(Ts[n] + CDrange,Td_plus[1],k_d[t],h) + 
                       Q_H(EBCoreInput$T_a[t],Ts[n] + CDrange,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC') + 
                       Q_LE(EBCoreInput$T_a[t],Ts[n] + CDrange,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')) #+ P(T_a(t),Ts(n)+CDrange,r(t)))...
                  - (EBCoreInput$Snet[t] + EBCoreInput$LWin[t] - Lup(Ts[n] - CDrange,epsilon_d,sigma_sb) + Q_G(Ts[n] - CDrange,Td_minus[1],k_d[t],h) + 
                       Q_H(EBCoreInput$T_a[t],Ts[n] - CDrange,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC') +
                       Q_LE(EBCoreInput$T_a[t],Ts[n] - CDrange,EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')) /
                    #+P(T_a(t),Ts(n)-CDrange,r(t))) 
                    
                        (2 * CDrange) );
  
  # Exit loop if one of the fluxes is NA, making calculation of T impossible
  if(is.na(Ts[n+1])){
    breakID <- 1
    stop = TRUE
    break
  }
  if(stop){break}

  # Allow no steps bigger than 1 degree (protect from low-derivative problem) 

  if (Ts[n+1] - Ts[n] > 1){
    Ts[n+1] = Ts[n]+0.2;
  }
  if (Ts[n+1] - Ts[n] < -1){
    Ts[n+1] = Ts[n]-0.2;
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

  #if(T_s[t + 1] > 308.15){
  #browser()
  #  T_s[t + 1] <- 308.15
  #}
  #if(T_s[t + 1] > 293.15){
  #browser()
  #  T_s[t + 1] <- T_s[t + 1] * 0.95
  #}
})
    #browser()
  T_d[t+1,] <- debristemp(T_s[t+1],T_s[t],T_d[t,],rho_d,c_d,h,timestep,N,k_d[t],T_f);
  #browser()
  }

  # End of Model Loop
    #browser()
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
  #browser()
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
  #browser()
  # Convert to metres water eq.
  melt <- G_i * timestep / (rho_w * L_f);
  melt_ice <- G_i * timestep / (rho_i * L_f);
  # Refreezing Potential
  refr <- melt
  refr[refr>0] <- 0

  melt[is.na(melt)] <- 0;
  melt[melt<0] <- 0;
  melt[EBCoreInput$snow>0.05] <- 0
  
  melt_ice[is.na(melt_ice)] <- 0;
  melt_ice[melt_ice<0] <- 0;
  melt_ice[EBCoreInput$snow>0.05] <- 0
  
  MeanDailyMelt <- mean(melt,na.rm=T)*(3600*24/timestep);
  
  cumM <- melt
  cumM[is.na(cumM)] <- 0
  CUMmelt <- cumsum(cumM*(3600/timestep));

  },
  'clean' = {
  # (2) Model clean ice glacier
    
    # (2a) Calculate energy for melt 
    rH_s_cig <- EBCoreInput$RH_a * 0 + 100;   # set relative humidity to full saturation

    Q_M <- EBCoreInput$Snet + EBCoreInput$LWin - EBCoreInput$LWout + 
      Q_H(EBCoreInput$T_a,EBCoreInput$T_s_data,EBCoreInput$RH_a,rH_s_cig,EBCoreInput$u,q_a,EBCoreInput$p_a,'clean') + 
      Q_LE(EBCoreInput$T_a,EBCoreInput$T_s_data,EBCoreInput$RH_a,rH_s_cig,EBCoreInput$u,q_a,EBCoreInput$p_a,'clean')
    
    # (2b) Longwave Radiation
    ##################################
    Lnet <- EBCoreInput$LWin - EBCoreInput$LWout;
    
    # (2c) Turbulent Fluxes
    ##################################
    Q_H_v <- Q_H(EBCoreInput$T_a,EBCoreInput$T_s_data,EBCoreInput$RH_a,rH_s_cig,EBCoreInput$u,q_a,EBCoreInput$p_a,'clean')
    Q_LE_v <- Q_LE(EBCoreInput$T_a,EBCoreInput$T_s_data,EBCoreInput$RH_a,rH_s_cig,EBCoreInput$u,q_a,EBCoreInput$p_a,'clean')
    
    # (3) Convert to metres water eq.
    melt <- Q_M * timestep / (rho_i * L_f);
    # Refreezing Potential
    refr <- melt
    refr[refr>0] <- 0
    
    melt[is.na(melt)] <- 0;
    melt[melt<0] <- 0;
    MeanDailyMelt <- mean(melt,na.rm=T)*(3600*24/timestep);
    
    cumM <- melt
    cumM[is.na(cumM)] <- 0
    CUMmelt <- cumsum(cumM*(3600/timestep));
    
    # (4) Fill irrelevant columns
    k_d <- CUMmelt*NA
    T_s <- CUMmelt*NA
    T_f <- CUMmelt*NA
    h <- CUMmelt*NA
    Q_G <- CUMmelt*NA
    T_d <- CUMmelt*NA
    },
'pond' = {
# (3) Model pond (based on Miles et al., 2016)

# (3a) Calculate energy for melt available at the pond surface
    
    # Calculate solar elevation theta
    I_SC <- 1367                          # Solar Constant
    phi <- EBCoreInput$lat*pi/180;        # Latitude
    Gamma <- 2*pi*(as.numeric(strftime(EBCoreInput$timeline_str, format = "%j"))-79.6764-0.2422*(year(EBCoreInput$timeline_str)-1985)+floor((year(EBCoreInput$timeline_str)-1985)/4))/65.2422;
    E_0 <- 1.00011+0.034221*cos(Gamma)+0.00128*sin(Gamma)-0.000719*cos(2*Gamma)+0.000077*sin(2*Gamma);
    # Calculate Hour angle for the relevant time of day
    omega <- 15*(12-hour(EBCoreInput$timeline_str)+1)*pi/180;
    # Calculate declination angle % turned off by Reid
    delta <- 23.45*sin(((as.numeric(strftime(EBCoreInput$timeline_str, format = "%j"))+284)*360/365)*pi/180)*pi/180; # Simpler method! -- JS: Cooper PI 1968, SE
    # Calculate solar elevation angle (also based on Wikipedia!)
    theta_s <- asin(cos(omega)*cos(delta)*cos(phi)+sin(delta)*sin(phi));
    
    albedo_pond <- 0.78 * theta_s^(-0.45)
    
    # Calculate surface water temperature based on Miles et al (2018, GRL; see supplement info S3)
    Ta_movav <- ma(EBCoreInput$T_a,24);
    Tws_anom <- 0.17843 * (EBCoreInput$T_a - Ta_movav) - 0.0001605
    Tws24 <- 0.0841782 * Ta_movav + 0.58899
    
    Tws <- 1.5 + (EBCoreInput$T_a - Ta_movav) * mean(Tws_anom/(EBCoreInput$T_a - Ta_movav)) + Ta_movav * mean(Tws24/Ta_movav)
    
    # Shortwave Radiation
    S_net <- EBCoreInput$SWin * (1 - albedo_pond)      # so far ignored diffuse and terrestrial SW

    # Longwave Radiation
    LW_out <- 0.95 * sigma_sb * Tws^4
    
    # Turbulent Fluxes
    rH_s_pond <- EBCoreInput$RH_a * 0 + 100;          # full saturation of surface humidity
    QH_pond <- Q_H(EBCoreInput$T_a,Tws,EBCoreInput$RH_a,rH_s_pond,EBCoreInput$u,q_a,EBCoreInput$p_a,'clean')
    QLE_pond <- Q_LE(EBCoreInput$T_a,Tws,EBCoreInput$RH_a,rH_s_pond,EBCoreInput$u,q_a,EBCoreInput$p_a,'clean')
    
    
    Q_M <- S_net + EBCoreInput$LWin - LW_out + QH_pond + QLE_pond
    
    # (3b) Calculate energy flux at the base
    k_w <- 0.565      # conductivity for water, J K-1 m-1 s-1
    k_r <- 2          # conductivity for rock, J K-1 m-1 s-1
    por <- 0.3        # porosity
    k_d_pond <- k_r (1 - por) + k_w * por
    
    

# (2b) Longwave Radiation
##################################
Lnet <- EBCoreInput$LWin - EBCoreInput$LWout;

# (2c) Turbulent Fluxes
##################################
Q_H_v <- Q_H(EBCoreInput$T_a,EBCoreInput$T_s_data,EBCoreInput$RH_a,rH_s_cig,EBCoreInput$u,q_a,EBCoreInput$p_a,'clean')
Q_LE_v <- Q_LE(EBCoreInput$T_a,EBCoreInput$T_s_data,EBCoreInput$RH_a,rH_s_cig,EBCoreInput$u,q_a,EBCoreInput$p_a,'clean')

# (3) Convert to metres water eq.
melt <- Q_M * timestep / (rho_i * L_f);
# Refreezing Potential
refr <- melt
refr[refr>0] <- 0

melt[is.na(melt)] <- 0;
melt[melt<0] <- 0;
MeanDailyMelt <- mean(melt,na.rm=T)*(3600*24/timestep);

cumM <- melt
cumM[is.na(cumM)] <- 0
CUMmelt <- cumsum(cumM*(3600/timestep));

# (4) Fill irrelevant columns
k_d <- CUMmelt*NA
T_s <- CUMmelt*NA
T_f <- CUMmelt*NA
h <- CUMmelt*NA
Q_G <- CUMmelt*NA
T_d <- CUMmelt*NA
})
  ##################################
  # =====================
  # Model Outputs
  # =====================  
  EBModCore <- cbind(melt,CUMmelt,T_s,k_d,T_f,h,EBCoreInput$Snet,Lnet,Q_G,Q_H_v,Q_LE_v,T_d,melt_ice,refr)

}