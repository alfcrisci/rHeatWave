inch2Millimeter <- function(val, ...) {

  val_new <- val * 25.4
  val_new <- round(val_new, ...)

  return(val_new)
}

knots2mps <- function(val, ...) {

  val_new <- val * 0.514444
  val_new <- round(val_new, ...)

  return(val_new)
}



knots2ms <- function(val, ...) {

  kmperhour = val*1.8535
  val_new = kmperhour*(1000/3600)
  val_new <- round(val_new, ...)

  return(val_new)
}



retrieve_GSOD=function(usaf,WBAN="99999",start_year = NA, end_year = NA, dsn = ".")
{

  fls_gz <- sapply(start_year:end_year, function(year) {
    dlbase <- paste0(as.character(usaf), "-", WBAN, "-", year, ".op.gz")
    dlurl <- paste0("ftp://ftp.ncdc.noaa.gov/pub/data/gsod/",year, "/", dlbase)
    dlfile <- paste0(dsn, "/", dlbase)

    if (file.exists(dlfile)) {
      cat("File", dlfile, "already exists. Proceeding to next file... \n")
    }
    else {
      try(download.file(dlurl, dlfile), silent = FALSE)
    }
    return(dlfile)  })

}


is.leapyear=function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}




findindexes=function(x,window=15) {
  j=window+x
  indexdayj=c(c(365-window):365,c(1:365),c(1:window))
  indexdayj[(j-window+1):(j+window-1)]
}



####################################################
#' selected_italian_weather
#'
#' @description  Main features of a selected GSOD italian weather station.
#'
#' @references   Data are obtained from GSOD NOAA site.  \url{https://data.noaa.gov/dataset/global-surface-summary-of-the-day-gsod}
#' @format data.frame
#'
#' @source \url{https://data.noaa.gov/dataset/global-surface-summary-of-the-day-gsod}
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it}
#'

####################################################
#' geo_selected_italian_weather
#'
#' @description  Geographical information of a selected GSOD italian weather station.
#'
#' @references   Data are obtained from GSOD NOAA site. \url{https://data.noaa.gov/dataset/global-surface-summary-of-the-day-gsod}
#' @format SpatialPointDataframe
#'
#' @source \url{https://data.noaa.gov/dataset/global-surface-summary-of-the-day-gsod}
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it}
#' @keywords  meteo spatial
#'

####################################################
#' calculate_HWdays
#'
#' Calculate HWdays days and Critical ones adding some columns at data.frame.
#' Mean air temperature, mean relative humidity and extreme thermal values are required.
#'
#'
#' @param data data.frame Daily time series with colums' names: date(YYYY-MM-DD),tmed,tmax,tmin,rhum,vmed.
#' @param clim_array data.frame Name of location or the time series analized.
#' @param method character Method to assess critical and heatwave days. Default is "EuroHeat" and other are "TMAX","TMED","TMIN".
#' @param index character Biometerological index to calculate apparent temperature with ATI or UTCI. Default is "ATI" corresponding to the shaded version of ATI Steadman (2004).
#' @param wind logical  If TRUE the wind in biometeorological index are considered.
#' @param months numeric The month her HW days are calculated. Default are c(5:9) corresponding to May-Dec.
#' @param P_intesity numeric If one of the "TMAX","TMED","TMIN" methods was selected is the level of intensity in term of percentiles.Default are c(95),
#' @param  duration numeric Minimum spell of critical days to be considered for heatave definition. Default is 2.
#' @param time_reference numeric Temporal format of the referenced climatology: DAY (daily) or MON (monthly).Default is "DAY".
#' @return Return a list of data.frame respectively with monthly and daily values.
#'
#'
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso Crisci \email{a.crisci@@ibimet.cnr.it}
#' @keywords heat wave
#'
#' @export
#'
#'
#'
#'
#'
#'


calculate_HWdays=function(data,
                          clim_array,
                          method="EuroHeat",
                          index="ATI",
                          wind=T,
                          P_intesity=c(95),
                          duration=2,
                          time_reference="DAY"
) {

  if (class(data) != "data.frame") {stop(paste("data argument must be a daily data.frame object with column's names equal to:",
                                               "date( %Y-%m-%d ),tmed,tmax,tmin,rhum,vmed"))
  }

  if (time_reference=="DAY") {if  (nrow(clim_array)!=365) {stop(paste("Valid climatology is needed!")) } }
  if (time_reference=="MONTH") {if  (nrow(clim_array)!=12) {stop(paste("Valid climatology is needed!")) } }
  if (wind==T ) { if (length(grep("vmed",names(data)))==0)  {stop(paste("Wind data are needed!")) } }
  if (length(grep(P_intesity,names(clim_array)))<1) {stop(paste("Valid P_intesity Percentile intensity value is needed!"))}

  temp_mat=data;
  temp_mat$yday=strptime(as.Date(data$date), format = "%Y-%m-%d")$yday+1
  temp_mat$yday365=temp_mat$yday
  temp_mat$yday365[grep("-02-29",temp_mat$date)]=28
  temp_mat$yday365[which(temp_mat$yday==366)]=365
  temp_mat$CritDay=0
  temp_mat$HWDay=0
  temp_matTappmax=NA
  temp_mat$Tappmed=NA
  temp_mat$mese=as.numeric(format(as.Date(temp_mat$date),"%m"))

  if ( index=="ATI")
  {
    temp_mat$Tappmax=round(steadman_indoor(temp_mat$tmax,temp_mat$rhum),1)
    temp_mat$Tappmed=round(steadman_indoor(temp_mat$tmed,temp_mat$rhum),1)

    if ( wind == T)
    {  temp_mat$Tappmax=round(steadman_outdoor_shade(temp_mat$tmax,temp_mat$rhum,temp_mat$vmed),1)
    temp_mat$Tappmed=round(steadman_outdoor_shade(temp_mat$tmed,temp_mat$rhum,temp_mat$vmed),1)
    }
  }

  if ( index=="UTCI")

  {
    temp_mat$Tappmax=round(UTCI(temp_mat$tmax,temp_mat$rhum,rep(0.1,nrow(temp_mat)),temp_mat$tmax),1)
    temp_mat$Tappmed=round(UTCI(temp_mat$tmed,temp_mat$rhum,rep(0.1,nrow(temp_mat)),temp_mat$tmed),1)


    if ( wind == T)
    { temp_mat$Tappmax=round(UTCI(temp_mat$tmax,temp_mat$rhum,temp_mat$vmed,temp_mat$tmax),1)
    temp_mat$Tappmed=round(UTCI(temp_mat$tmed,temp_mat$rhum,temp_mat$vmed,temp_mat$tmed),1)
    }
  }

  #####################################################################################################################
  if ( time_reference == "DAY") {

    if (method == "EuroHeat"  & wind==T)

    {
      for (i in 1:nrow(temp_mat))
      { iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
      if (is.na(temp_mat$Tappmax[i])) {temp_mat$CritDay[i]=0;next}
      if (is.na(temp_mat$tmin[i])) {temp_mat$CritDay[i]=0;next}
      a=ifelse(temp_mat$Tappmax[i] >= clim_array$TmaxApp_V_P90_DAY[temp_mat$yday365[i]],1,0)
      b=ifelse(temp_mat$tmin[i] >= clim_array$Tmin_P90_DAY[temp_mat$yday365[i]] && temp_mat$Tappmax[i] >= clim_array$TmaxApp_V_P50_DAY[temp_mat$yday365[i]],1,0)
      if(a==1 || b==1) {temp_mat$CritDay[i]=1}
      if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }

    if (method == "EuroHeat"  & wind==F)

    {
      for (i in 1:nrow(temp_mat))
      { iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
      if (is.na(temp_mat$Tappmax[i])) {temp_mat$CritDay[i]=0;next}
      if (is.na(temp_mat$tmin[i])) {temp_mat$CritDay[i]=0;next}
      a=ifelse(temp_mat$Tappmax[i] >= clim_array$TmaxApp_P90_DAY[temp_mat$yday365[i]],1,0)
      b=ifelse(temp_mat$tmin[i] >= clim_array$Tmin_P90_DAY[temp_mat$yday365[i]] && temp_mat$Tappmax[i] >= clim_array$TmaxApp_P50_DAY[temp_mat$yday365[i]],1,0)
      if(a==1 || b==1) {temp_mat$CritDay[i]=1}
      if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }

    if (method == "TMED" )

    {
      for (i in 1:nrow(temp_mat))
      {
        iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
        if (is.na(temp_mat$tmed[i])) {temp_mat$CritDay[i]=0;next}
        a=ifelse(temp_mat$tmed[i] >= clim_array[temp_mat$yday365[i],paste0("Tmax_P",P_intesity,"_DAY")],1,0)
        if( a==1)  {temp_mat$CritDay[i]=1}
        if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }

    if (method == "TMAX" )

    {
      for (i in 1:nrow(temp_mat))
      {
        iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
        if (is.na(temp_mat$tmax[i])) {temp_mat$CritDay[i]=0;next}
        a=ifelse(temp_mat$tmax[i] >= clim_array[temp_mat$yday365[i],paste0("Tmax_P",P_intesity,"_DAY")],1,0)
        if( a==1)  {temp_mat$CritDay[i]=1}
        if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }
    if (method == "TMIN" )

    {
      for (i in 1:nrow(temp_mat))
      {
        iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
        if (is.na(temp_mat$tmin[i])) {temp_mat$CritDay[i]=0;next}
        a=ifelse(temp_mat$tmin[i] >= clim_array[temp_mat$yday365[i],paste0("Tmin_P",P_intesity,"_DAY")],1,0)
        if( a==1)  {temp_mat$CritDay[i]=1}
        if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }

    ########################################################################################################################

  } # day loop
  ########################################################################################################################

  if ( time_reference == "MON") {

    if (method == "EuroHeat"  & wind==T)

    {
      for (i in 1:nrow(temp_mat))
      { iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
      if (is.na(temp_mat$Tappmax[i])) {temp_mat$CritDay[i]=0;next}
      if (is.na(temp_mat$tmin[i])) {temp_mat$CritDay[i]=0;next}
      a=ifelse(temp_mat$Tappmax[i] >= clim_array$TmaxApp_V_P90_MON[temp_mat$mese[i]],1,0)
      b=ifelse(temp_mat$tmin[i] >= clim_array$Tmin_P90_MON[temp_mat$mese[i]] && temp_mat$Tappmax[i] >= clim_array$TmaxApp_V_P50_MON[temp_mat$mese[i]],1,0)
      if(a==1 || b==1) {temp_mat$CritDay[i]=1}
      if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }

    if (method == "EuroHeat"  & wind==F)

    {
      for (i in 1:nrow(temp_mat))
      { iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
      if (is.na(temp_mat$Tappmax[i])) {temp_mat$CritDay[i]=0;next}
      if (is.na(temp_mat$tmin[i])) {temp_mat$CritDay[i]=0;next}
      a=ifelse(temp_mat$Tappmax[i] >= clim_array$TmaxApp_P90_MON[temp_mat$mese[i]],1,0)
      b=ifelse(temp_mat$tmin[i] >= clim_array$Tmin_P90_MON[temp_mat$mese[i]] && temp_mat$Tappmax[i] >= clim_array$TmaxApp_P50_MON[temp_mat$mese[i]],1,0)
      if(a==1 || b==1) {temp_mat$CritDay[i]=1}
      if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }

    if (method == "TMED" )

    {
      for (i in 1:nrow(temp_mat))
      {
        iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
        if (is.na(temp_mat$tmed[i])) {temp_mat$CritDay[i]=0;next}
        a=ifelse(temp_mat$tmed[i] >= clim_array[temp_mat$mese[i],paste0("Tmax_P",P_intesity,"_MON")],1,0)
        if( a==1)  {temp_mat$CritDay[i]=1}
        if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }

    if (method == "TMAX" )

    {
      for (i in 1:nrow(temp_mat))
      {
        iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
        if (is.na(temp_mat$tmax[i])) {temp_mat$CritDay[i]=0;next}
        a=ifelse(temp_mat$tmax[i] >= clim_array[temp_mat$mese[i],paste0("Tmax_P",P_intesity,"_MON")],1,0)
        if( a==1)  {temp_mat$CritDay[i]=1}
        if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }

    if (method == "TMIN" )
    {
      for (i in 1:nrow(temp_mat))
      {
        iHW=ifelse(i-(duration-1)>0,i-(duration-1),1)
        if (is.na(temp_mat$tmin[i])) {temp_mat$CritDay[i]=0;next}
        a=ifelse(temp_mat$tmin[i] >= clim_array[temp_mat$mese[i],paste0("Tmin_P",P_intesity,"_MON")],1,0)
        if( a==1)  {temp_mat$CritDay[i]=1}
        if (sum(temp_mat$CritDay[iHW:i]) == duration) {temp_mat$HWDay[i]=1}
      }
    }
  }
  # MON Loop

  #######################################################################################################################
                          
  res_HW=temp_mat["date","mese","yday365","CritDay","HWDay"]
   
  return(res_HW)

}


#' create_clim_param_HW
#'
#' Calculate the parameters to assess heat-wave daily status from a climatic time series monthly.
#' Mean air temperature , mean relative humidity and extreme thermal values are required.
#' Parameters
#'
#'
#' @param data data.frame Daily time series with colums' names: date(YYYY-MM-DD),tmed,tmax,tmin,rhum,vmed.
#' @param name character Name of location or ID of time series.
#' @param index character Biometerological index to calculate apparent temperature.
#' @param iniclim character Initial year of climatological reference time period, Default is "1981".
#' @param yfinclim character Final year of climatological reference time period, Default is "2010".
#' @param PT numeric  Three levels of percentile used as thresholds for intensity. Default are c(90, 95, 98).
#' @param window_days numeric Half Day span to calculate daily climatology. Default is 15.
#' @param save_clim logical If True save file daily and montly in RDS and CSV format. Default is F.
#' @return Return a list of data.frame respectively with monthly and daily values.
#'
#'
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso Crisci \email{a.crisci@@ibimet.cnr.it}
#' @keywords heat wave
#'
#' @export
#'
#'
#'
#'


create_clim_param_HW=function(data,
                              name,
                              index="ATI",
                              yiniclim="1981",
                              yfinclim="2010",
                              PT=c(90,95,98),
                              window_days=15,
                              save_clim=F) {

  if (class(data) != "data.frame") {stop(paste("data argument must be a daily data.frame object with column's names equal to:",
                                               "date( %Y-%m-%d ),tmed,tmax,tmin,rhum,vmed"))
  }
  if (length(PT) != 3 ) {stop(paste("Please  three values of PT parameter are required!"))}

  dates=seq(as.Date(paste0(yiniclim,"-01-01")),as.Date(paste0(yfinclim,"-12-31")),1)
  indexday=which(as.character(data$date) %in% as.character(dates) ==T)
  anni=unique(as.numeric(format(as.Date(data$date),"%Y")))
  if (length(which(anni %in% yiniclim:yfinclim==T)) < length(yiniclim:yfinclim) ) {stop(paste("Please refine year climate range!"))}

  ##########################################################################################################################################################################



  PTD=PT/100
  temp_mat=data;
  temp_mat$yday=strptime(as.Date(data$date), format = "%Y-%m-%d")$yday+1
  temp_mat$yday365=temp_mat$yday
  temp_mat$yday365[grep("-02-29",temp_mat$date)]=28
  temp_mat$yday365[which(temp_mat$yday==366)]=365

  ##########################################################################################################################################################################

  temp_mat$Tappmax=NA
  temp_mat$Tappmed=NA
  temp_mat$Tappmax_v=NA
  temp_mat$Tappmed_v=NA

  if ( index=="ATI")
  {
    temp_mat$Tappmax=round(steadman_indoor(temp_mat$tmax,temp_mat$rhum),1)
    temp_mat$Tappmed=round(steadman_indoor(temp_mat$tmed,temp_mat$rhum),1)

    if ( length(grep("vmed",names(temp_mat)))>0)
    {  temp_mat$Tappmax_v=round(steadman_outdoor_shade(temp_mat$tmax,temp_mat$rhum,temp_mat$vmed),1)
    temp_mat$Tappmed_v=round(steadman_outdoor_shade(temp_mat$tmed,temp_mat$rhum,temp_mat$vmed),1)
    }
  }

  if ( index=="UTCI")
  {
    temp_mat$Tappmax=round(UTCI(temp_mat$tmax,temp_mat$rhum,rep(0.1,nrow(temp_mat)),temp_mat$tmax),1)
    temp_mat$Tappmed=round(UTCI(temp_mat$tmed,temp_mat$rhum,rep(0.1,nrow(temp_mat)),temp_mat$tmed),1)


    if ( length(grep("vmed",names(temp_mat)))>0)
    { temp_mat$Tappmax_v=round(UTCI(temp_mat$tmax,temp_mat$rhum,temp_mat$vmed,temp_mat$tmax),1)
      temp_mat$Tappmed_v=round(UTCI(temp_mat$tmed,temp_mat$rhum,temp_mat$vmed,temp_mat$tmed),1)
    }
  }

  ##########################################################################################################################################################################


  temp_mat_climdata=temp_mat[indexday,]
  temp_mat_climdata$mese=as.numeric(format(as.Date(temp_mat_climdata$date),"%m"))

  ##########################################################################################################################################################################

  temp_clim=data.frame(Tmax_P50=NA,Tmax_PT1=NA,Tmax_PT2=NA,Tmax_PT3=NA,
                       TmaxApp_P50=NA,TmaxApp_PT1=NA,TmaxApp_PT2=NA,TmaxApp_PT3=NA,
                       TmaxApp_V_P50=NA, TmaxApp_V_PT1=NA,TmaxApp_V_PT2=NA,TmaxApp_V_PT3=NA,
                       Tmin_P50=NA,Tmin_PT1=NA,Tmin_PT2=NA,Tmin_PT3=NA,
                       Tmed_P50=NA,Tmed_PT1=NA,Tmed_PT2=NA,Tmed_PT3=NA)

  temp_clim_daily=cbind(id=1:365,temp_clim)
  temp_clim_monthly=cbind(id=1:12,temp_clim)

  #######################################################################################################################################################à

  for ( i in 1:365) {

    temp_clim_daily$Tmax_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmax,na.rm=T))
    temp_clim_daily$Tmax_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmax,PTD[1],na.rm=T)[1])
    temp_clim_daily$Tmax_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmax,PTD[2],na.rm=T)[1])
    temp_clim_daily$Tmax_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmax,PTD[3],na.rm=T)[1])

    temp_clim_daily$TmaxApp_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$Tappmax,na.rm=T))
    temp_clim_daily$TmaxApp_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$Tappmax,PTD[1],na.rm=T)[1])
    temp_clim_daily$TmaxApp_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$Tappmax,PTD[2],na.rm=T)[1])
    temp_clim_daily$TmaxApp_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$Tappmax,PTD[3],na.rm=T)[1])

    temp_clim_daily$TmaxApp_V_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$Tappmax_v,na.rm=T))
    temp_clim_daily$TmaxApp_V_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$Tappmax_v,PTD[1],na.rm=T)[1])
    temp_clim_daily$TmaxApp_V_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$Tappmax_v,PTD[2],na.rm=T)[1])
    temp_clim_daily$TmaxApp_V_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$Tappmax_v,PTD[3],na.rm=T)[1])


    temp_clim_daily$Tmed_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmed,na.rm=T))
    temp_clim_daily$Tmed_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmed,PTD[1],na.rm=T)[1])
    temp_clim_daily$Tmed_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmed,PTD[2],na.rm=T)[1])
    temp_clim_daily$Tmed_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmed,PTD[3],na.rm=T)[1])

    temp_clim_daily$Tmin_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmin,na.rm=T))
    temp_clim_daily$Tmin_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmin,PTD[1],na.rm=T)[1])
    temp_clim_daily$Tmin_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmin,PTD[2],na.rm=T)[1])
    temp_clim_daily$Tmin_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$yday365 %in% findindexes(i,window_days))==TRUE),]$tmin,PTD[3],na.rm=T)[1])

  }



  for ( i in 1:12) {

    temp_clim_monthly$Tmax_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmax,na.rm=T))
    temp_clim_monthly$Tmax_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmax,PTD[1],na.rm=T)[1])
    temp_clim_monthly$Tmax_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmax,PTD[2],na.rm=T)[1])
    temp_clim_monthly$Tmax_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmax,PTD[3],na.rm=T)[1])

    temp_clim_monthly$TmaxApp_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$Tappmax,na.rm=T))
    temp_clim_monthly$TmaxApp_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$Tappmax,PTD[1],na.rm=T)[1])
    temp_clim_monthly$TmaxApp_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmax,PTD[2],na.rm=T)[1])
    temp_clim_monthly$TmaxApp_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmax,PTD[3],na.rm=T)[1])

    temp_clim_monthly$TmaxApp_V_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$Tappmax_v,na.rm=T))
    temp_clim_monthly$TmaxApp_V_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$Tappmax_v,PTD[1],na.rm=T)[1])
    temp_clim_monthly$TmaxApp_V_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$Tappmax_v,PTD[2],na.rm=T)[1])
    temp_clim_monthly$TmaxApp_V_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$Tappmax_v,PTD[3],na.rm=T)[1])


    temp_clim_monthly$Tmed_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmed,na.rm=T))
    temp_clim_monthly$Tmed_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmed,PTD[1],na.rm=T)[1])
    temp_clim_monthly$Tmed_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmed,PTD[2],na.rm=T)[1])
    temp_clim_monthly$Tmed_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmed,PTD[3],na.rm=T)[1])

    temp_clim_monthly$Tmin_P50[i]=as.numeric(median(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmin,na.rm=T))
    temp_clim_monthly$Tmin_PT1[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmin,PTD[1],na.rm=T)[1])
    temp_clim_monthly$Tmin_PT2[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmin,PTD[2],na.rm=T)[1])
    temp_clim_monthly$Tmin_PT3[i]=as.numeric(quantile(temp_mat_climdata[which((temp_mat_climdata$mese == i)),]$tmin,PTD[3],na.rm=T)[1])

  }


  #######################################################################################################################################################################

  names(temp_clim_daily)=gsub("T1",PT[1],names(temp_clim_daily))
  names(temp_clim_daily)=gsub("T2",PT[2],names(temp_clim_daily))
  names(temp_clim_daily)=gsub("T3",PT[3],names(temp_clim_daily))
  names(temp_clim_monthly)=gsub("T1",PT[1],names(temp_clim_monthly))
  names(temp_clim_monthly)=gsub("T2",PT[2],names(temp_clim_monthly))
  names(temp_clim_monthly)=gsub("T3",PT[3],names(temp_clim_monthly))

  names(temp_clim_daily)=paste0(names(temp_clim_daily),"_DAY")
  names(temp_clim_monthly)=paste0(names(temp_clim_monthly),"_MON")

  if ( save_clim == T) {
    saveRDS(temp_clim_daily,paste0(name,"_DAYCLIM.rds"))
    write.csv(temp_clim_daily,paste0(name,"_DAYCLIM.csv"))
    saveRDS(temp_clim_monthly,paste0(name,"_MONCLIM.rds"))
    write.csv(temp_clim_monthly,paste0(name,"_MONCLIM.csv"))

  }



  return(list(MONCLIM=temp_clim_monthly,
              DAYCLIM=temp_clim_daily))
}







#############################################################################################################################################
#' manage_GSOD
#'
#' Download and build a specific data.frame for heatwave analisys from a GSOD NOAA repository station based on its unique USAF code.
#'
#'
#' @param id character USAF code of weather station.
#' @param name character Name of location or ID of time series.
#' @param ini_year character Initial year of data time period, Default is "1979".
#' @param last_year character Final year of climatological reference time period, Default is "2016".
#' @param save_files logical If True save GSOD data file in RDS and CSV format. Default is F.
#' @return  Return  data.frame ready to use for heat-wave analitics.
#'
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso Crisci \email{a.crisci@@ibimet.cnr.it}
#' @keywords GSOD
#' @importFrom GSODTools gzGsodStations
#' @importFrom lubridate ymd
#' @export
#'
#'
#'

manage_GSOD=function(id,name,ini_year = 1979, last_year = 2015, missing=NA,save_files=T) {
  retrieveGSOD(usaf=id,start_year = ini_year, end_year = last_year)
  temp_daily<- gzGsodStations(x,start_year = ini_year, end_year = last_year)
  temp_daily$Date=as.Date(ymd(temp_daily$YEARMODA))
  temp_daily$YEARMODA=NULL
  temp_daily$tmed <- toCelsius(temp_daily$TEMP, digits = 1)
  temp_daily$tmax <- toCelsius(temp_daily$MAX, digits = 1)
  temp_daily$tmin <- toCelsius(temp_daily$MIN, digits = 1)
  temp_daily$prec <-inch2Millimeter(temp_daily$PRCP, digits = 1)
  temp_daily$DEWP <- toCelsius(temp_daily$DEWP, digits = 1)
  temp_daily$rhum <- 100*(exp((17.625*temp_daily$DEWP)/(243.04+temp_daily$DEWP))/exp((17.625*temp_daily$tmed)/(243.04+temp_daily$tmed)))
  temp_daily$slp=temp_daily$SLP
  temp_daily$prec[which(temp_daily$PRCPFLAG=="I")]=missing
  temp_daily$vmed=knots2ms(temp_daily$WDSP,1)
  temp_daily$vmax=knots2ms(temp_daily$MXSPD,1)
  temp_daily=temp_daily[c("date","tmed","tmax","tmin","rhum","DEWP","vmed","vmax","prec","slp")]
  if ( save_files == T) {
                         saveRDS(temp_daily,paste0(name,"_DAYCLIM.rds"))
                         write.csv(temp_daily,paste0(name,"_DAYCLIM.csv"))

                        }
  return(temp_daily)
}








