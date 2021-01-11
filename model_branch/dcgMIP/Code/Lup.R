################################################################################
# Upwelling Longwave Radiation
# 
# Lup.R
#
# ReadMe:
# Code to calculate upwelling longwave radiation from surface temperature
# 
# Requires Surface Temperature [C] as Input

# Created:          2017/03/17
# Latest Revision:  2017/03/27
#
# Jakob F Steiner | PhD candidate | Faculty of Geosciences | Universiteit Utrecht | 
# Heidelberglaan 2, 3584 CS Utrecht | W.C. van Unnik building | Room 124, Zonneveldvleugel | 
# j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org
#
################################################################################

Lup <- function(Tsfc,epsilon_d,sigma){

L <- epsilon_d*sigma*(Tsfc)^4;

# Gap Filling???
L[is.na(L)]<-mean(L,na.rm=T)
Lup <- L;
}
