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

debristemp_loop <- function(Ts_t1,Ts_t,Td,rho_d,c_d,h,timestep,N,k_d,T_f){


  a_Crank <- matrix(0,N, 1)
  b_Crank <- matrix(0,N, 1)
  c_Crank<- matrix(0,N, 1)
  d_Crank <- matrix(0,N, 1)
  A_Crank <- matrix(0,N, 1)
  S_Crank <- matrix(0,N, 1)
  #N <- N + 1
  #browser()
  Td_old <- c(Td,T_f)
  Td <- c(Td,T_f)

  #Td[N,t] <- 273.15
  
  C <- k_d*timestep/(2*rho_d*c_d*h^2)
  
    #Perform Crank-Nicholson Scheme
    for(j in 2:(N-1)){
      # Equations A8 in Reid and Brock (2010) 
      
      a_Crank[j] <- C
      b_Crank[j] <- 2*C+1
      c_Crank[j] <- C
      
      
      # Equations A9 in Reid and Brock (2010) 
      if(j == 2){
        d_Crank[j] <- C*Td[1] + C*Td_old[1] + (1-2*C)*Td_old[j] + C*Td_old[j+1]
      } else if(j < (N-1)){
        d_Crank[j] <- C*Td_old[j-1] + (1-2*C)*Td_old[j] + C*Td_old[j+1]
      } else if(j == (N-1)){
        d_Crank[j] <- 2*C*Td[N] + C*Td_old[N-2] + (1-2*C)*Td_old[N-1];
      }
      # note notation:
      # "i-1" refers to the past
      # "j-1" refers to the cell above it
      # "j+1" refers to the cell below it          
      
      # Equations A10 and A11 in Reid and Brock (2010)
      if(j == 2){
        A_Crank[j] <- b_Crank[j];
        S_Crank[j] <- d_Crank[j];
      } else {
        A_Crank[j] <- b_Crank[j] - a_Crank[j]/A_Crank[j-1]*c_Crank[j-1];
        S_Crank[j] <- d_Crank[j] + a_Crank[j]/A_Crank[j-1]*S_Crank[j-1];
      }
      #browser()
    }
    
    # Equations A12 in Reid and Brock (2010)
    for(j in (N-1):2){
      if(j == (N-1)){
        Td[j] <- S_Crank[j]/A_Crank[j]
      } else {
        Td[j] <- 1/A_Crank[j]*(S_Crank[j]+c_Crank[j]*Td[j+1]);
      }
    }
    #browser()
  return(Td[1:N-1])
  }  