EBModRounce <- function(EBCoreInput){
  
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



N <- N+1


Td <- matrix(0,N, datasize)
# Td is the temperature in the debris with row 1 being Ts
a_Crank <- matrix(0,N, datasize)
b_Crank <- matrix(0,N, datasize)
c_Crank<- matrix(0,N, datasize)
d_Crank <- matrix(0,N, datasize)
A_Crank <- matrix(0,N, datasize)
S_Crank <- matrix(0,N, datasize)
n_iterations <- matrix(0,1, datasize)
Ts_past <- matrix(0,1, datasize)
eS_Saturated <- matrix(0,1, datasize)
eS_dry <- matrix(0,1, datasize)
eS <- matrix(0,1, datasize)
eZ<- matrix(0,1, datasize)
LE_Benn<- matrix(0,1, datasize)
Rn <- matrix(0,1, datasize)
H_Benn <- matrix(0,1, datasize)
Qc <- matrix(0,1, datasize)
P_Flux <- matrix(0,1, datasize)
dLE_Benn <- matrix(0,1, datasize)
dRn <- matrix(0,1, datasize)
dH_Benn <- matrix(0,1, datasize)
dQc <- matrix(0,1, datasize)
dP_Flux <- matrix(0,1, datasize)
F_Ts_Benn <- matrix(0,1, datasize)
dF_Ts_Benn <- matrix(0,1, datasize)
Qc_ice <- matrix(0,1, datasize)
Melt <- matrix(0,1, datasize)

#k_d <- k_d_wet

n_iterations <- vector()
for(t in 1:datasize){

n_iterations[t] <- 0
Ts_past[t] <- 0
Td[N,t] <- 273.15

switch(seasonsDyn,
       on = {                      # Separate conductivity values for each Season
         if (as.numeric(strftime(EBCoreInput$timeline_str[t], format = "%j")) > as.numeric(strftime(paste(year(EBCoreInput$timeline_str[t]),'/',monin,sep=""), format = "%j")) 
             && as.numeric(strftime(EBCoreInput$timeline_str[t], format = "%j")) < as.numeric(strftime(paste(year(EBCoreInput$timeline_str[t]),'/',monout,sep=""), format = "%j"))) 
           k_d <- k_d_wet
         else
           k_d <- k_d_dry},
       off = {
         k_d <- k_d_wet
       }
)
k_d <- 0.8
C <- k_d*timestep/(2*rho_d*c_d*h^2)
# Initially assume Ts = Tair, for all other time steps assume it's equal to previous Ts
  if(t == 1){Td[1,t] <- EBCoreInput$T_a[t]
  } else {
    Td[1,t] <- Td[1,t-1]
    }

# Calculate temperature profile in the debris
# For t = 0, which is i = 1, assume initial condition of linear temperature profile in the debris
if(t == 1){
  Td_gradient <- (Td[1,1] - Td[N,1])/d
  for(j in 2:(N-1)){
  Td[j,1] <- Td[1,1] - (j*h)*Td_gradient
  }
} else {
#Perform Crank-Nicholson Scheme
for(j in 2:(N-1)){
  # Equations A8 in Reid and Brock (2010) 

  a_Crank[j,t] <- C
  b_Crank[j,t] <- 2*C+1
  c_Crank[j,t] <- C
  
  
  # Equations A9 in Reid and Brock (2010) 
  if(j == 2){
  d_Crank[j,t] <- C*Td[1,t] + C*Td[1,t-1] + (1-2*C)*Td[j,t-1] + C*Td[j+1,t-1]
  } else if(j < (N-1)){
  d_Crank[j,t] <- C*Td[j-1,t-1] + (1-2*C)*Td[j,t-1] + C*Td[j+1,t-1]
  } else if(j == (N-1)){
  d_Crank[j,t] <- 2*C*Td[N,t] + C*Td[N-2,t-1] + (1-2*C)*Td[N-1,t-1];
  }
  # note notation:
    # "i-1" refers to the past
  # "j-1" refers to the cell above it
  # "j+1" refers to the cell below it          
  
  # Equations A10 and A11 in Reid and Brock (2010)
  if(j == 2){
  A_Crank[j,t] <- b_Crank[j,t];
  S_Crank[j,t] <- d_Crank[j,t];
  } else {
    A_Crank[j,t] <- b_Crank[j,t] - a_Crank[j,t]/A_Crank[j-1,t]*c_Crank[j-1,t];
  S_Crank[j,t] <- d_Crank[j,t] + a_Crank[j,t]/A_Crank[j-1,t]*S_Crank[j-1,t];
  }

}
  
  # Equations A12 in Reid and Brock (2010)
  for(j in (N-1):2){
  if(j == (N-1)){
  Td[j,t] <- S_Crank[j,t]/A_Crank[j,t]
  } else {
    Td[j,t] <- 1/A_Crank[j,t]*(S_Crank[j,t]+c_Crank[j,t]*Td[j+1,t]);
  }
  }

}


Qc <- k_d*(Td[2,t] - Td[1,t])/h;
Rn <- EBCoreInput$Snet[t] + EBCoreInput$LWin[t] - Lup(Td[1,t],epsilon_d,sigma_sb)
QH <- Q_H(EBCoreInput$T_a[t],Td[1,t],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')
QLE <- Q_H(EBCoreInput$T_a[t],Td[1,t],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')
F_Ts_Benn <- Rn + QH + QLE + Qc;

dF_Ts_Benn <- -4*epsilon_d*(5.67*10^-8)*Td[1,t]^3 + (-k_d/h) +dQ_H(EBCoreInput$T_a[t],Td[1,t],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC') + dQ_LE(EBCoreInput$T_a[t],Td[1,t],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')



# Newton-Raphson method to solve for surface temperature
while(abs(Td[1,t]-Ts_past[t]) > 0.01 && n_iterations[t] < 100){
  

n_iterations[t] <- n_iterations[t] + 1;
Ts_past[t] <- Td[1,t];
# max step size is 1 degree C
Td[1,t] <- Ts_past[t] - F_Ts_Benn/dF_Ts_Benn;
if(Td[1,t] - Ts_past[t] >= 1){
Td[1,t] <- Ts_past[t] + 1
} else if(Td[1,t] - Ts_past[t] <= -1){
Td[1,t] <- Ts_past[t] - 1;
}

if(t == 1){
  Td_gradient <- (Td[1,1] - Td[N,1])/d
  for(j in 2:(N-1)){
    Td[j,1] <- Td[1,1] - (j*h)*Td_gradient
  }
} else {
  #Perform Crank-Nicholson Scheme
  for(j in 2:(N-1)){
    # Equations A8 in Reid and Brock (2010) 
    
    a_Crank[j,t] <- C
    b_Crank[j,t] <- 2*C+1
    c_Crank[j,t] <- C
    
    
    # Equations A9 in Reid and Brock (2010) 
    if(j == 2){
      d_Crank[j,t] <- C*Td[1,t] + C*Td[1,t-1] + (1-2*C)*Td[j,t-1] + C*Td[j+1,t-1]
    } else if(j < (N-1)){
      d_Crank[j,t] <- C*Td[j-1,t-1] + (1-2*C)*Td[j,t-1] + C*Td[j+1,t-1]
    } else if(j == (N-1)){
      d_Crank[j,t] <- 2*C*Td[N,t] + C*Td[N-2,t-1] + (1-2*C)*Td[N-1,t-1];
    }
    # note notation:
    # "i-1" refers to the past
    # "j-1" refers to the cell above it
    # "j+1" refers to the cell below it          
    
    # Equations A10 and A11 in Reid and Brock (2010)
    if(j == 2){
      A_Crank[j,t] <- b_Crank[j,t];
      S_Crank[j,t] <- d_Crank[j,t];
    } else {
      A_Crank[j,t] <- b_Crank[j,t] - a_Crank[j,t]/A_Crank[j-1,t]*c_Crank[j-1,t];
      S_Crank[j,t] <- d_Crank[j,t] + a_Crank[j,t]/A_Crank[j-1,t]*S_Crank[j-1,t];
    }
    
  }
  
  # Equations A12 in Reid and Brock (2010)
  for(j in (N-1):2){
    if(j == (N-1)){
      Td[j,t] <- S_Crank[j,t]/A_Crank[j,t]
    } else {
      Td[j,t] <- 1/A_Crank[j,t]*(S_Crank[j,t]+c_Crank[j,t]*Td[j+1,t]);
    }
  }
  
}

Qc <- k_d*(Td[2,t] - Td[1,t])/h;
Rn <- EBCoreInput$Snet[t] + EBCoreInput$LWin[t] - Lup(Td[1,t],epsilon_d,sigma_sb)
QH <- Q_H(EBCoreInput$T_a[t],Td[1,t],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')
QLE <- Q_H(EBCoreInput$T_a[t],Td[1,t],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')
F_Ts_Benn <- Rn + QH  +QLE + Qc;

dF_Ts_Benn <- -4*epsilon_d*(5.67*10^-8)*Td[1,t]^3 + (-k_d/h) +dQ_H(EBCoreInput$T_a[t],Td[1,t],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC') + dQ_LE(EBCoreInput$T_a[t],Td[1,t],EBCoreInput$RH_a[t],0,EBCoreInput$u[t],q_a[t],EBCoreInput$p_a[t],'nSC')

if(n_iterations == 100){
  Td[1,t] <- (Td[1,t] + Ts_past[t]) / 2;
}

}




if(n_iterations[t] == 0){browser()}
}
browser()
EBModRounce <- Td[1,]
# (1b) Longwave Radiation
##################################
Lup_v <- Lup(Td[1,], epsilon_d, sigma_sb);
Lnet <- EBCoreInput$LWin - Lup_v;

# (1c) Turbulent Fluxes
##################################
#browser()
Q_H_v <- Q_H(EBCoreInput$T_a,Td[1,],EBCoreInput$RH_a,NULL,EBCoreInput$u,q_a,EBCoreInput$p_a,'nSC')
Q_LE_v <- Q_LE(EBCoreInput$T_a,Td[1,],EBCoreInput$RH_a,NULL,EBCoreInput$u,q_a,EBCoreInput$p_a,'nSC')  

# (1d) Ground Heat Flux
Q_G <- Q_G(Td[1,],Td[2,],k_d,h)

# (1d) Calculate Melt Rate
##################################
#% Calculate melt via conductive flux G_i at bottom of debris (into the ice)
G_i <- vector();
melt <- vector();
# G_i depends on the temperature gradient at the debris base
G_i <- k_d * (Td[N-1,] - T_f) / h; # [W m-2]
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
T_s<-Td[1,];

#EBModRounce <- cbind(melt,CUMmelt,T_s,k_d,T_f,h,EBCoreInput$Snet,Lnet,Q_G,Q_H_v,Q_LE_v,melt_ice,refr)
}

                                          