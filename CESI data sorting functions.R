# functions for CESI - to sort hydrological data
#  1) ice.free
#  2) find.median.freshet
#  3) after.freshet
#  4) ret.freq.curve
#  5) Threshold_Build
#  6) percentile.ranked
#  7) bday_calculation
#  8) stn.locations

#-------------------------------------------------------------------------------
### 1) Function to calculate ice-free period for a single year of data using "B"
### flag in hydat database which indicates ice influence
#####  Input: 2 vectors of same length
#####             date = in format MM-DD
#####             symbol = site condition code
##### Output: vector with same length as input with TRUE during ice-free period
#####           and FALSE in ice affected period

ice.free <- function(date, symbol) {
  
  output <- rep(TRUE, length(date))
  
  month <- as.numeric(substr(date, 6 ,7))
  data <- data.frame(Date=date, Month=month, Symbol=symbol)
  
  # First B day
  firstbday <- data %>% filter(Month > 8) %>%  filter(Symbol == "B")  %>% head(1)
  firstbday.date <- ifelse(nrow(firstbday)>0, firstbday$Date, NA)
  if (!is.na(firstbday.date)) {
    FU <- which(as.Date(data$Date)==as.Date(firstbday.date))
    output[FU:length(date)] <- FALSE
  }

  # Last B day
  lastbday <- data %>% filter(Month < 6) %>% filter(Symbol == "B") %>% tail(1)
  lastbday.date <- ifelse(nrow(lastbday)>0, lastbday$Date, NA)
  if (!is.na(lastbday.date)) {
    BU <- which(as.Date(data$Date)==as.Date(lastbday.date))
    output[1:BU] <- FALSE
  }
  
  return(output)
}

#-------------------------------------------------------------------------------
### 2) Function to calculate median freshet date from median hydrograph during
### reference period as the date with maximum flow between 1 march and 31 july
#####  Input: stn = alphanumeric station number
#####         start.yr = number, first year of reference period
#####         end.yr = number, last year of reference period
##### Output: date of median freshet in format MM-DD
find.median.freshet <- function(stn, start.yr, end.yr) {
  flow <- hy_daily_flows(stations[j,1], start_date=paste0(start.yr,"-01-01"),
                         end_date=paste0(end.yr,"-12-31"))
  flow$M_D <- substr(flow$Date, 6, 10)
  
  median.flow <- data.frame(Date=format(seq.Date(as.Date("03-01",format="%m-%d"),
                        as.Date("07-31",format="%m-%d"), by=1), format="%m-%d"),
                        flow=NA)

  for (b in 1:nrow(median.flow)) {
    median.flow$flow[b] <- median(flow$Value[flow$M_D==median.flow$Date[b]],
                                na.rm = TRUE)
  }
  
  output <- median.flow %>% filter(flow==max(median.flow$flow, na.rm=TRUE)) %>%
                    pull(Date)
  return(output)
}

#-------------------------------------------------------------------------------
### 3) Function to calculate date of freshet for specific year as date of max flow 
### in 6 week window centered on date of median freshet for that station
#####  Input: 2 vectors of same length
#####             date = in format MM-DD
#####             flow = flow value
#####         date.med.fresh = median freshet date in format MM-DD
##### Output: vector with same length as input with TRUE during period after
#####           freshet and FALSE prior to it

after.freshet <- function(date, flow, date.med.fresh) {

  output <- rep(TRUE, length(date))
  
  window.start <- which(date==format((as.Date(date.med.fresh, format="%m-%d")-21),
                    format="%m-%d"))
  window.end <- which(date==format((as.Date(date.med.fresh, format="%m-%d")+21),
                  format="%m-%d"))
  
  if (length(window.start)>0 & length(window.end)>0) {
    FM <- which(flow==max(flow[window.start:window.end], na.rm=TRUE))
    if (!is.na(FM[1])) { output[1:FM[1]] <- FALSE }
  }
  
  return(output)
}

#-------------------------------------------------------------------------------
### 4) Function to fit a curve to the lower end of a non-exceedence distribution
### to find a return frequency flow
##### Input:  ret.freq = numeric value of drought return frequency
#####         2 vectors of same length
#####             years = numeric value of years
#####             annual.minima = numeric value of annual minima
#####         ref.years = vector of numeric values of years from which to calculate return frequency flow
##### Output: dataframe of the format: censored || freq0 || best.dist || R2 ||
#####           param1 || param2 || param3|| ret.freq.flow
#####         These represent:  censored - True/False, whether there are 0 values
#####                           freq0 - proportion of the values that are 0
#####                           best.dist - which of the 7 distributions fits best
#####                           R2 - how well distribution fits data (1 is perfect)
#####                           param1,2,3 - parameters of best fit curve
#####                           ret.freq.flow - return frequency flow

ret.freq.curve <- function(ret.freq, years, annual.minima, ref.years) {
  # initialize output dataframe
  output <- data.frame(censored=NA, freq0=NA, best.dist=NA, R2=NA, param1=NA, 
                       param2=NA, param3=NA, ret.freq.flow=NA)
  
  # Curves for 7 different distribution types are tested to keep the best fit
  distributions <- c("pe3", # Pearson type III
                     "wei", # Weibull
                     "ln3", # lognormal
                     "nor", # normal
                     "exp", # exponential
                     "gevR", # generalized extreme value
                     "gum") # Gumbell

  ann.min.ref <- annual.minima[which(years %in% ref.years)] %>% na.omit()
  # only calculate curve if more than 75% of annual minima values are non-zero
  if ((sum(ann.min.ref!=0)>(0.75*length(ann.min.ref))) & (length(ann.min.ref)>=10)) { 
    fit <- suppressWarnings(evquantile(evfit(ann.min.ref, distribution = distributions,
             zeta=0, extreme = "minimum"), ret.freq))
    
    # note in output if series is censored and frequency of zero values
    output$censored <- fit$is.censored
    output$freq0 <- fit$freq.zeros
    #determine best fitting function based on largest R^2 value
    ## estimates are only made on censored values (zero values are removed), so 
    ## create data frame with censored values and estimates
    cens <- data.frame(values=fit$values[fit$values!=0],fit$estimates)
    cens <- cens[order(cens$values),]
    ## since we are only interested in the fit of the smaller values (lowest half),
    ## find what the half point is
    half <- ceiling(nrow(cens)/2)
    ## calculate R2
    R2 <- cor(cens$values[1:half],cens[1:half,2:8])^2
    #### This will be empty when the lowest half values are identical, in which
    #### a curve should not be calculated
    if (any(!is.na(R2))) { #only for those data for which there is at least one curve
      best.dist <- which(R2==max(R2)) #find which distribution curve has the highest R2
      # if best distribution gives return frequency <0.001, when there are no 0
      # values in the minimum data, try next best option
      if (fit$is.censored==FALSE & fit$T_Years_Event[best.dist]<0.001) {
        best.dist2 <- which(R2==sort(R2,partial=6)[6])
        if (fit$T_Years_Event[best.dist2]>0.001) {
          best.dist <- best.dist2
        }
      }
      # note for output details of curve
      output$best.dist <- distributions[best.dist]
      output$R2 <- R2[best.dist]
      output$param1 <- unlist(fit$parameters[best.dist])[1]
      output$param2 <- unlist(fit$parameters[best.dist])[2]
      output$param3 <- unlist(fit$parameters[best.dist])[3]
     # set to 0 negative return frequency flows, an artifact of calculations and not possible
      output$ret.freq.flow <- ifelse(fit$T_Years_Event[best.dist]<0,  0,
                                     fit$T_Years_Event[best.dist])
    }
  } else { #for stations with too many zeros
    output$censored <- TRUE
    output$freq0 <- ">25%"
  }

  return(output)
}

#-------------------------------------------------------------------------------
#5) Threshold_Build
#'
#' @description
#' Threshold_Building takes into a list of stations, year_from and year_to representing the beginning
#' and the end of the reference period, and the two ratio numbers hi and low representing the
#' percentile(high and low). bday_calculation function will be needed to calculate the threshold.
#'
#'
#' @param stations a vector of characters indicating the list of station we want to calculate the
#' threshold on;
#' @param yfrom a numeric indicating the year threshold calculation starts, defaulted to 1981;
#' @param yto a numeric indicating the year threshold calculation ends, defaulted to 2010;
#' @param hi a numeric indicating the upper confident threshold, defaulted to 0.95;
#' @param low a numeric indicating the lower confident threshold, defaulted to 0.05
#' @param write_csv a boolean value indicating whether ot not the output will be written
#' to csv, defaulted to FALSE;
#'
#' @return
#' A dataframe of the format:
#' || station ||UpperThreshold ||LowerThreshold ||
#'
#' @export

Threshold_Build <- function(stations, yfrom=1981, yto=2010, hi=0.95, low=0.05, write_csv=FALSE){
  Upper = c()
  Lower = c()
  for (i in 1:length(stations)){
    stn.id <- stations[i]
#    print(stn.id)
    flow.daily <- hy_daily_flows(stn.id)
    flow.daily$Year = as.numeric(substr(flow.daily$Date, 1, 4))
    ref <- flow.daily %>% filter(Year>=yfrom, Year<=yto)
    summ = ref %>% group_by(Year) %>% summarise(count=length(Value[!is.na(Value)]))
    if (length(summ$Year)<20 | length(summ$count[summ$count>=150])<length(summ$count)){
      Upper=append(Upper, NA)
      Lower=append(Lower, NA)
#      print(paste0(stn.id, "---not enough data"))
    }else{
      Upper_All=ref$Value
      for (i in 1:length(Upper_All)){
        if (is.na(Upper_All[i])){
          Upper_All[i]=0
        }
      }
      yr=summ$Year
      num_leap_yrs=length(yr[as.numeric(yr)%%4==0])
      total_length=length(yr)*365+num_leap_yrs
      Upper_All=append(Upper_All, rep.int(0, total_length-length(Upper_All)))
      Lower_All=c()
      # First B day
      flow.daily$Month = as.numeric(format(flow.daily$Date, "%m"))
      for (j in 1:length(summ$Year)){
        year.op=summ$Year[j]
        fbday=bday_calculation(stn.id, year.op)
        lastbday = (fbday %>% filter(Year==year.op))[["LastBday"]]
        firstbday = (fbday %>% filter(Year==year.op))[["FirstBday"]]
        Value_within=(ref %>% filter(as.Date(Date)>=as.Date(lastbday)&as.Date(Date)<=as.Date(firstbday)))$Value
        Value_within=Value_within[Value_within!=0]
        Lower_All=append(Lower_All, Value_within)
      }
      Upper=append(Upper,quantile(Upper_All, probs = hi, names = FALSE, na.rm = TRUE))
      Lower=append(Lower,quantile(Lower_All, probs = low, names = FALSE, na.rm = TRUE))
#      print(stn.id)
    }
  }
  Threshold = data.frame(STATION_NUMBER = stations, Q_Upper=Upper, Q_Lower=Lower)
  if (write_csv){
    write.csv(Threshold, file = "../Dependencies/Threshold_New.csv", row.names = FALSE)
  }
  return(Threshold)
}

#-------------------------------------------------------------------------------
#6) percentile.ranked
#' @description
#' Calculates the pencentage of values in a vector smaller than value
#' @param a.vector a vector of any comparable data type, ideally numeric;
#'
#' @param value a comparable data type to the components in a.vector, ideally numeric
#'
#' @return the percentage of values within the vector that's smaller than the value itself.
#' @export
percentile.ranked <- function(a.vector, value) {
  numerator <- length(sort(a.vector)[a.vector < value])
  denominator <- length(a.vector)
  round(numerator/denominator,2)*100
}

#-------------------------------------------------------------------------------
#7) bday_calculation
#'
#' @description
#' bday_calculation calculates the last and first bday for some years with
#' record at a station.
#'
#' @param station a character representing the station number that exists in hydat
#' @param year a vectoer of numeric representing the list of years used for calculating
#' the threshold
#'
#' @return a dataframe of the format || station || year || lastbday || firstbday ||
#' @export
bday_calculation <- function(station, year="all"){
  flow.daily <- hy_daily_flows(station)
  flow.daily$Year <- as.numeric(format(flow.daily$Date, "%Y"))
  flow.daily$Month<- as.numeric(format(flow.daily$Date, "%m"))
  if (year=="all"){
    valid_year=unique(flow.daily$Year)
  }else if(any(!year %in% flow.daily$Year)){
    valid_year=year[(year %in% unique(flow.daily$Year))]
    return("Some of the year provided are not in HYDAT")
  }else{
    valid_year=year
  }
  first_vec = c()
  last_vec = c()
  for (i in valid_year){
#    print(i)
    # First B day
    flow.late <- flow.daily %>% filter(Year==i & Month > 8)
    firstbday <- flow.late %>% filter(Symbol == "B", Month >= 10) %>% head(1)
    firstbday.date <- ifelse(nrow(firstbday)>0, substring(firstbday$Date,1,10), NA)
    
    firstbday.long <- ifelse(is.na(firstbday.date),paste0(i,"-10-01"),
                             ifelse(firstbday.date > as.Date(paste0(i, "-10-01")),
                                    paste0(i, "-10-01"),firstbday.date))
    first_vec = c(first_vec, firstbday.long)
    # Last B day
    flow.early <- flow.daily %>% filter(Year == i & Month < 6)
    lastbday <- flow.early %>% filter(Symbol == "B", Month < 6) %>% tail(1)
    lastbday.date <- ifelse(nrow(lastbday)>0, substring(lastbday$Date,1,10), NA)
    
    lastbday.long<- ifelse(is.na(lastbday.date),paste0(i,"-05-01"),
                           ifelse(lastbday.date < as.Date(paste0(i, "-05-01")),
                                  paste0(i,"-05-01"),lastbday.date))
    last_vec = c(last_vec, lastbday.long)
  }
  ret = data.frame(Station=rep(station, length(valid_year)),
                   Year = valid_year,
                   LastBday = last_vec,
                   FirstBday = first_vec)
  return(ret)
}

#-------------------------------------------------------------------------------
### 8) Function to find the province or territory of a list of stations, and if
### stations are in the United States, they are assigned to the province or
### territory to which they are hydrolically connected
##### Input:  stn_list = list of alphanumeric station names
##### Output: data frame with two columns: station list that was input and 
#####         abbreviation for province or territory where station is located

stn.locations <- function(stn_list) {
  PT <- hy_stations(stn_list) %>% select(STATION_NUMBER, PROV_TERR_STATE_LOC)
  colnames(PT) <- c("STATION_NUMBER", "PT")
  ## Stations located in the United States are counted in the adjacent territory or 
  ## province.
  AKtoYT <- "09ED001"
  AKtoBC <- c("08AB002", "08BB005", "08CF003")
  PT$PT[PT$STATION_NUMBER %in% AKtoYT] <- "YT"
  PT$PT[PT$STATION_NUMBER %in% AKtoBC] <- "BC"
  PT$PT[PT$PT=="WA"] <- "BC"
  PT$PT[PT$PT=="ID"] <- "BC"
  MTtoBC <- c("08NG039", "08NP001")
  MTtoAB <- c("05AD030", "05AD029", "05AD031", "05AD032", "05AD033", "05AE029",
              "05AE030", "05AE031","05AE032","05AE033", "05AE036","05AE028",
              "05AE034", "05AE035","11AA030", "11AA032", "11AA033", "11AA031")
  MTtoSK <- c("11AB105", "11AB107", "11AD001", "11AE007", "11AE004", "11AE005", 
              "11AE006", "11AE008", "11AF004", "11AE009")
  PT$PT[PT$STATION_NUMBER %in% MTtoBC] <- "BC"
  PT$PT[PT$STATION_NUMBER %in% MTtoAB] <- "AB"
  PT$PT[PT$STATION_NUMBER %in% MTtoSK] <- "SK"
  NDtoSK <- c("05NB026", "05NB027", "05ND007")
  NDtoMB <- c("05PC018", "05OA005", "05OB031", "05OB033", "05OC004")
  PT$PT[PT$STATION_NUMBER %in% NDtoSK] <- "SK"
  PT$PT[PT$STATION_NUMBER %in% NDtoMB] <- "MB"
  MNtoMB <- c("05OD030", "05OD032", "05OD031", "05PD001")
  MNtoON <- c("05NF012", "02AA001")
  PT$PT[PT$STATION_NUMBER %in% MNtoMB] <- "MB"
  PT$PT[PT$STATION_NUMBER %in% MNtoON] <- "ON"
  PT$PT[PT$PT=="MI"] <- "ON"
  PT$PT[PT$PT=="ME"] <- "NB"
  
  return(PT)
}