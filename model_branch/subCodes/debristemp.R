################################################################################
# Debris Temperatur Profile
# 
# debristemp.R
#
# ReadMe:
# Code to calculate the Debris Temperature Profile for the Energy Balance
# 
# Requires .... Surface Temperature [C] as Input
#% Requires these inputs: 
#  % New surface temp (T_s(t+1), here called Ts_t1)
#% Surface temp on last timestep (T_s(t), here called Ts_t)
#% Debris temp profile on last timestep, (T_d(t,:), here called Td)
# Created:          2017/03/17
# Latest Revision:  2017/03/27
#
# Jakob F Steiner | PhD candidate | Faculty of Geosciences | Universiteit Utrecht | 
# Heidelberglaan 2, 3584 CS Utrecht | W.C. van Unnik building | Room 124, Zonneveldvleugel | 
# j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org
#
################################################################################

debristemp <- function(Ts_t1,Ts_t,Td,rho_d,c_d,h,timestep,N,k_d,T_f){
#function [debristemp] = debristemp(Ts_t1,Ts_t,Td)
#browser()
#global rho_d c_d h timestep N k_d T_f
# Set size of debristemp output array
dt <- rep(0, N-1);
# C, A and S are intermediate variables for debris temp calculation.
C <- k_d * timestep / (2 * rho_d * c_d * h^2);
A <- rep(0, N-1);
S <- rep(0, N-1);
# Find value of A for top layer.
A[1] <- 2 * C + 1;
# Find values of A for every other layer.
for (i in 2:(N-1)){    
A[i] <- 2*C+1 - C^2 / A[i-1]
}
# Find value of S for top layer.
S[1] <- C*Ts_t1 + C*Ts_t + (1-2*C)*Td[1] + C*Td[2];
# Find values of S for all the middle layers.
for (y in 2:(N-2)){
S[y] <- C*Td[y-1] + (1-2*C)*Td[y] + C*Td[y+1] + C*S[y-1]/A[y-1]}
# Find value of S for bottom layer.
S[N-1] <- 2*C*T_f + C*Td[N-2] + (1-2*C)*Td[N-1] + C*S[N-2]/A[N-2];
# Find new ice temperature in bottom layer.
dt[N-1] <- S[N-1] / A[N-1];
# Find new ice temperature at all other layers.
for (y in 1:(N-2)){
dt[N-1-y] <- (S[N-1-y] + C*dt[N-1-y+1]) / A[N-1-y]}

debristemp <- dt
}