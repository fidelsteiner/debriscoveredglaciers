################################################################################
diu_cycle<-function(dataframe){
  
  # CALCULATE DIURNAL CYCLE (+SD) PER DATA-VECTOR ->TS x VARIABLES
  #
  # dataframe:   Dataframe with timesteps (numeric or POSIX-format) 
  #               & values for each timestep (each column one variable)
  #
  # e.g.:
  #  Date                  L_s   L_d   H      Q_m   melt
  #  2013-05-19 00:00:00   220.3 127.5 6.793  48.44 0.0006
  #  2013-05-19 01:00:00   216.2 126.2 3.954  40.19 0.0005
  #  2013-05-19 02:00:00   210.8 122.5 2.172  29.27 0.0004
  #  2013-05-19 03:00:00   208.3 119.5 2.453  24.14 0.0003
  #  2013-05-19 04:00:00   210.6 120.2 3.912  28.48 0.0004
  #  2013-05-19 05:00:00   212.3 121.6 4.608  32.30 0.0004
  colnames(dataframe)[1]<-'Date' 
  
  # extract only hours per timestep
  date_h<-as.POSIXlt(dataframe$Date,origin='1970-01-01')$hour
  
  # remove 'Date'-column
  dataframe<-dataframe[!colnames(dataframe) %in% 'Date']
  
  # aggregate flux values to diurnal cycle (and compute sd)
  Flux_dc_ls<-list()
  Flux_sd_ls<-list()
  for(v in 1:ncol(dataframe)){          
    Flux_dc_ls[[v]]<-aggregate(dataframe[[v]],list(date_h),mean)[,2]
    Flux_sd_ls[[v]]<-aggregate(dataframe[[v]],list(date_h),sd)[,2]
  }
  
  # dirunal cycle                                       
  DC_df<-data.frame(1:24,do.call(cbind,Flux_dc_ls))
  colnames(DC_df)<-c('hours',colnames(dataframe)) 
  
  # standard deviation
  SD_df<-data.frame(1:24,do.call(cbind,Flux_sd_ls))
  colnames(SD_df)<-c('hours',colnames(dataframe)) 
  
  output_ls<-list(DC_df,SD_df)
  names(output_ls)<-c('DC','SD')
  return(output_ls)
}