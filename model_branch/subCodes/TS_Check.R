TS_Check <- function(TS,maxV,minV){
  
  library(tsoutliers)
  m_TS <- TS
  
  m_TS[m_TS>=maxV]<-NA
  m_TS[m_TS<=minV]<-NA
  #browser()
  for(k in 2:length(m_TS)){
if(!is.na(m_TS[k])&!is.na(m_TS[k-1])){
    if (abs(m_TS[k])>200&abs(m_TS[k]-m_TS[k-1])>5*abs(m_TS[k-1])){
      #browser()
      m_TS[k] <- NA}}
  }
  
  TS_Check <- data.frame(cbind(m_TS))
  TS_Check <- TS_Check[,1]
  
}