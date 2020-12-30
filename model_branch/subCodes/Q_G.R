################################################################################
# Ground Heat Flux
# 
# Q_G.R
#
# ReadMe:
# Code to calculate Ground/Conductive Heat Flux
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

Q_G <- function(Ts,Td1,k_d,h){

Q_G <- -k_d*(Ts-Td1)/h;
}