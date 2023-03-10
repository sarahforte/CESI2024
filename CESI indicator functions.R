# functions for CESI - to calculate indicators
#  1) get_ann_mean_flow
#  2) get_flood_metrics
#  3) drought_days

#--------------------------------------------------
#1) Get_ann_mean_flow
#'
#' @description
#' Get_ann_mean_flow calculates the annual mean flow based on a station and a year
#' vector indicating which year the mean flow will be calculated.
#'
#' @param station a character indicating the station number to calculate the annual mean flow;
#' @param year a character indicating the list of years where the annual mean flow is being
#' calculated, or by default equals to 'all' indicating all recorded years at the station
#'
#' @return
#' a data frame of the format:
#' || station || Year || annual mean flow||
#'
#' @export

get_ann_mean_flow <- function(station, year='all'){
  mean_flow = c()
  flow.daily <- hy_daily_flows(station)
  flow.daily$Year <- as.numeric(format(flow.daily$Date, "%Y"))
  if (year[1]=="all") {
    valid_year=unique(flow.daily$Year)
    absent = c()
  }else if(any(!year %in% flow.daily$Year)){
    valid_year=year[(year %in% unique(flow.daily$Year))]
    absent = year[which(!year %in% flow.daily$Year)]
  }else{
    valid_year=year
    absent = c()
  }
  for (i in valid_year){
    flow.dates = flow.daily[flow.daily$Year==i,]
    if (nrow(flow.dates %>% filter(!is.na(Value)))<150){
      ann_mean_flow = NA
    }else{
      if (nrow(flow.dates%>%filter(!is.na(Value))) > 0.9 * ifelse(flow.dates$Year[1]%%4==0, 366, 365)){
        ann_mean_flow <- flow.dates %>% mutate(mean_flow = mean(Value, na.rm = TRUE))
        ann_mean_flow <- ann_mean_flow[1, "mean_flow"] %>% round(2)
        ann_mean_flow <- ann_mean_flow$mean_flow
      } else {
        Mode <- flow.daily  %>% group_by(Year) %>% filter(!is.na(Value)) %>%
          summarise(MODE=names(sort(table(Value), decreasing = TRUE)[1]))
        Mode.multi <- mean(as.numeric(Mode$MODE))
        if (Mode.multi <0.01){
          ann_mean_flow <- flow.dates %>% filter(!is.na(Value)) %>% group_by(Year) %>%
            summarise(TOTAL=sum(Value)/n())
          ann_mean_flow <- round(ann_mean_flow$TOTAL,2)
        } else {
          ann_mean_flow <- NA
        }
      }
    }
    mean_flow <- c(mean_flow, ann_mean_flow)
  }
  
  ret = data.frame(STATION_NUMBER=rep(station, (length(valid_year)+length(absent))),
                   Year = c(valid_year, absent),
                   ann_mean_flow = c(mean_flow, rep(NA, length(absent))))
  return(ret)
}

#--------------------------------------------------
#2) get_flood_metrics
#'
#' @description
#' Get_flood_metrics takes in a character representing the station number within 1014
#' RHBN_U stations and a vector of numerics or a character, calculates four flood metrics:
#' pot_threshold: the upper threshold that was calculated in Threshold_build function from above;
#' pot_days:days over the threshold within a year;
#' pot_events:independent flood events over a year;
#' pot_ma_dur: the maximum consecutive duration over the threshold within a year;
#' The four metrics will then be recorded in a table and returned.
#'
#' @details
#' Note that there are requirements regarding the independence conditions for pot_event metric,
#' we implement a check for the correlation between flood events.
#' Referring to Lang etal (1999), we set two conditions:
#' 1) R > 5+log(A/1.609^2);
#' 2) Xmin < 0.75 * min(Xi, Xj);
#' where R is the time in days between two peaks, A is the watershed area in squared kilometers;
#' Xmin is the minimum flow between the adjacent flow peaks Xi and Xj.
#' In order for two seperate peaks to be considered as independent events of drought or flood,
#' both conditions need to be met. If either of the conditions fail, we consider two events
#' as the same one. R is defined within the iteration.
#'
#' @references
#' Slater et al.(2016). On the impact of gaps on trend detection un extreme
#' streamflow time series. Int.J.climatol. 37:3976-3983(2017). doi: 10.1002/joc.4954
#'
#' @param station a character representing a station we want the flood metrics to be calculated at;
#' @param year a vector of numeric representing the list of years we want thew flood metrics to be calculated at,
#' defaulted to 'all' indicating all the years at the station;
#'
#' @return
#' A dataframe of the format:
#' || STATION_NUMBER || Year || pot_days || pot_events || pot_threshold || pot_max_dur ||
#
#' @export
get_flood_metrics <- function(date, value, threshold, R){
  # initialize output data frame
  output <- data.frame(pot_days=NA, pot_events=NA, pot_max_dur=NA, X1_day_max=NA)
  
  flow.data <- data.frame(Date=date, Value=value)
  
  output$X1_day_max <- max(flow.data$Value, na.rm=TRUE)

  if (output$X1_day_max >= threshold){
    factor <- case_when( threshold >= 1000 ~ 1000,
                         threshold >= 100 ~ 100,
                         threshold >= 10 ~ 10,
                         threshold >= 0 ~ 1)
    pot2 <- possibly(evir::pot, otherwise = NA)
    pot <- pot2(flow.data$Value[!is.na(flow.data$Value)]/factor, threshold/factor)
    if (!is.na(pot[1])) {
      output$pot_days <- length(pot$data)
    
      pot.dates <- data.frame(dates = attr(pot$data, "times"))
      pot.dates$lag <- pot.dates$dates - lag(pot.dates$dates, 1)
      seq_blocks <- split(pot.dates$dates, cumsum(c(TRUE, diff(pot.dates$dates)!=1)))
      # We iterate within seq_blocks to find maximum flow within each events, before testing
      # their independence
      events_count <- c()
      for (a in 1:length(seq_blocks)){
        events_count <- append(events_count, length(seq_blocks[[as.character(a)]]))
      }
      events_array <- c()
      for (b in 1:length(events_count)){
        c <- events_count[b]
        while (c > 0){
          events_array<-append(events_array,b)
          c <- c-1
        }
      }
      pot_dat <- data.frame(data = pot$data, dates=attr(pot$data, "times"),event_num=events_array)
      max_peak <- pot_dat %>% group_by(event_num) %>% summarise(max=max(data))
      max_peak_central<-c()
      for (d in 1:length(max_peak$max)){
        max_peak_dates<-(pot_dat %>% filter(data==max_peak$max[d],dates>=seq_blocks[[d]][1],
                              dates<=seq_blocks[[d]][length(seq_blocks[[d]])]))$dates
        if (length(max_peak_dates)%%2==1){
          max_peak_dates=max_peak_dates[(length(max_peak_dates)+1)/2]
        } else {
          max_peak_dates=(max_peak_dates[length(max_peak_dates)/2]+
                          max_peak_dates[length(max_peak_dates)/2+1])/2
        }
        max_peak_central=append(max_peak_central, max_peak_dates)
        max_peak_dates=c()
      }
      max_peak$central_dates<-max_peak_central
      R_vec<-c()
      for (t in 1:length(max_peak$max)){
        R_vec<-append(R_vec, R)
      }
      max_peak$R_vec<-R_vec
      Xmin_vec<-c()
      if (length(max_peak$max)==1){
        Xmin_vec<-append(Xmin_vec, NA)
      }else{ 
        for( e in 1:(length(max_peak$max)-1)){
          Xmin_vec<-append(Xmin_vec, 0.75*min(max_peak$max[e],max_peak$max[e+1])*factor)
        }
        Xmin_vec<-append(Xmin_vec, NA)
      }
      max_peak$Xmin_vec<-Xmin_vec
      ind_events<-c()
      if (length(seq_blocks)>1){
        for (j in 1:(length(seq_blocks)-1)){
          range_flow <- flow.data %>% filter(Date>=as.Date(max_peak$central_dates[j], origin=paste0(yr, "-01-01")),
                                            Date<=as.Date(max_peak$central_dates[j+1],origin=paste0(yr, "-01-01")))
          if (all(min(range_flow$Value, na.rm=TRUE)<Xmin_vec[j],
                ((max_peak$central_dates[j+1] - max_peak$central_dates[j])>R))){
            ind_events <- append(ind_events, 1)
          }else{
            ind_events <- append(ind_events, 0)
          }
        }
        ind_events <- append(ind_events, 1)
      } else { 
        ind_events <- c(1)
      }
      max_peak$ind_events<-ind_events
      output$pot_events <- sum(ind_events)
    
      max_dur_count <- c()
      count<-0
      for (f in 1:length(ind_events)){
        if (ind_events[f]==1 | is.na(ind_events[f])){
          max_dur_count = append(max_dur_count, count + events_count[f])
          count = 0
        } else { # i.e. in_events[i] == 0
          count <- count + events_count[f]}
        }
      output$pot_max_dur <- max(max_dur_count)
  
    } else {
      output[1:3] <- 0
    }
  } # close loop when for when it is possible to calculate peaks over threshold

  return(output)
}

#--------------------------------------------------
### 3) Function to calculate hydrological drought events using lfstats package
### that is an implementation of methods detailed in: World Meteorological 
### Organization (2008) Manual on low-flow estimation and prediction
##### Input:  Y = numeric value of year
#####         2 vectors of same length
#####             M_D = days in format MM-DD
#####             flow = flow for each of the days in M-D
#####         dr_thresh = numeric value of drought threshold
#####         stn = alphanumeric station name
##### Output: dataframe of the format: STATION_NUMBER || type || year || dr_events || 
#####                                  dr_max_dur || dr_days

drought_days <- function(Y, M_D, flow, dr_thresh, stn) {
  # initialize output data frame
  output <- data.frame(STATION_NUMBER=stn, type="perennial", year=Y, dr_events=NA, 
                       dr_max_dur=NA, dr_days=NA)
  
  # transform the data into object to use lfstat package
  data <- data.frame(day=substr(M_D,4,5), month=substr(M_D,1,2), 
                       year=rep(Y,length(M_D)), flow)
  data.lf <- createlfobj(data, hyearstart=1)
  flowunit(data.lf) <- "m^3/s"
  
  # lfstat function
  dr <- find_droughts(data.lf, threshold=dr_thresh)
  if (max(dr$event.no)>0) { # if there are drought events
    ## pool adjacent events with less than 5 days between and 
    ## Vabove / min(Vbelowi, Vbelowj)>=0.1 where, Vabove is the volume
    ## exceeding the threshold, Vbelowi and Vbelowj are the volumes 
    ## below the threshold during events i and j.
    if (max(dr$event.no, na.rm=TRUE)>1) { dr <- pool_ic(dr, tmin=5, ratio=0.1) }
    
    summary2 <- possibly(base::summary, otherwise = data.frame(event.no=NULL))
    dur<-summary2(dr) #summary function that eliminates droughts less than 5 days
    
    
    ## write messages about pooling/elimination
    if (!is.null(attributes(dur)$deficit$n.pooled)) {
      if (attributes(dur)$deficit$n.pooled>0) {
        message(paste0(stn," has ", attributes(dur)$deficit$n.pooled,
                       " pooled event(s) in ",Y))
      }
    }
    if (attributes(dur)$deficit$omitted>0) {
      message(paste0(stn," has ", attributes(dur)$deficit$omitted,
                     " eliminated event(s) in ", Y))
    }
    
    ## if all droughts were eliminated, note that there were none
    if (nrow(dur)==0) { 
      output[1,4:6] <- 0
    } else { # note results
      output[1,4:6] <- c(nrow(dur), max(dur$duration), sum(dur$duration))
    }
  } else { ## for those years without drought
    output[1,4:6] <- 0
  }
  
  return(output)
}
