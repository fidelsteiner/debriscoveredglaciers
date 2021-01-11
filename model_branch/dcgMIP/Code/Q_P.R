################################################################################
# Precipitation Heat Flux
# 
# Q_P.R
#
# ReadMe:
# Code to calculate Precipitation Heat Flux
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

Q_P <- function(Ta,Ts,rho_w,rain,c_w){
  Q_P <- rho_w * rain * c_w * (Ta - Ts);
}