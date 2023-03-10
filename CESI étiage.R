#######################################################################
#
#    -###- CANADIAN ENVIRONMENTAL SUSTAINABILITY INDICATORS -###-
#              -- Water Quantity Indicator Calculator --
#
#  This project contains scripts used to automate the calculation
#          of CESI hydrological drought indicator.
#
#  Environment and Climate Change Canada
#  Created: November 2022
#  Updated: February 2023
#
#######################################################################

# This script calculates the number of hydrological drought days during the 
# summer at each hydrometric station based on a 30 year reference period.
# Different methodologies are used for perennial and intermittent
# watercourses, but both involve defining a station-specific threshold
# and comparing the yearly summer data to the threshold.
# Trends over 50 years are calculated and the results are displayed as a kriged map.
# The specific steps are:
#    1) Define variables
#    2) Create tables to compile results
#    3) Retrieve hydrometric data from hydat database
#    4) Check if there is sufficient data for calculations
#    5) Identify summer period
#    6A) For perennial watercourses:
#       i. calculate rolling average
#      ii. define threshold by fitting curve to 30-year reference data
#     iii. find droughts using threshold while pooling adjacent events and
#            eliminating short events
#    6B) For intermittent watercourses:
#       i. calculate dry spell durations
#      ii. define threshold from cumulative distribution of dry spell duration
#     iii. evaluate each year against threshold
#    7) Calculate 50-year trends
#    8) Identify outlier trends and eliminate them
#    9) Save output as csv files
#    10) Krige results and produce map

# Assuming there are three sub-folders in working directory: 
# 1. "Dependencies" - with following files
#             country boundaries [CanadaBound.shp]
#             outlines of big lakes [CANwaterbodies.shp]
#             CESI north arrow [esri6north.png]
#             arrows to interpret trends [up_level_down_arrows.png]
#             station list [RHBN_U.csv] (optional, there are other ways of getting station list)
#             kriging zones [5 merged ecozones.shp] (only if regional kriging to be used)
# 2. "Output" - in which to save the output
# 3. "R" - in which 5 script files are located:
#             [CESI étiage.R-THIS ONE]
#             [CESI data sorting functions.R]
#             [CESI indicator functions.R]
#             [CESI trend functions.R]
#             [CESI mapping sorting functions.R]

# Assuming hydat database has been downloaded and is available. If not, run:
# download_hydat()

#######################################################################

# load necessary libraries, the order in which you import them is important not
# to mask certain functions (specifically select from dplyr)
library(readxl)    # for importing excel file
library(tidyhydat) # for interacting with hydat database
library(lfstat)    # for fitting curves to low flow statistics
library(purrr)     # function 'possibly' is easiest way to deal with errors
library(MASS)      # for function in Negative Binomial & hurdle tests
library(countreg)  # for hurdle test
library(trend)     # for wald-wolfowitz stationarity test
library(EnvStats)  # for outlier test
library(sf)        # to create spatial objects
library(terra)     # to work with raster data
library(corrplot)  # for colourbar legend
library(automap)   # to model variogram
library(gstat)     # for Kriging
library(tidyverse) # to merge dataframes, includes library dplyr
library(png)       # to import image


source('./R/CESI data sorting functions.R')
source('./R/CESI indicator functions.R')
source('./R/CESI trend functions.R')
source('./R/CESI mapping functions.R')

######### STEP 1) Define variables #########
# Three options for creating hydrometric station list for calculations. Only
# RHBN stations are used.
# 1-simplest command which should work if hydat was up to date
# stations<- hy_stations() %>% filter(RHBN)  %>% select(STATION_NUMBER)
# 2-also simple but creates dependency
# stations <- read.csv("./Dependencies/RHBN_U.csv", header = TRUE)
# 3-a bit more complex, but with no dependency
url <- "https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/RHBN/RHBN_Metadata.xlsx"
destfile <- "./Output/RHBN_Metadata.xlsx"
download.file(url, destfile, method="curl")
stations <- read_xlsx(destfile, range="A3:P1286")
stations <- stations[-(1),]
stations <- filter(stations, Evaluation_Year==2020, DATA_TYPE=="Q")

yrs.of.int <- c(1970:2021) #period of interest
yrs.of.ref <- c(1991:2020) #30-year reference period
n.day <- 7  #number of days for rolling average
return.freq <- 5 #return frequency in years for n.day average
perc <- 0.9 #percentile for threshold dry spell duration
summer.beg <- "05-01" #first day of summer period
summer.end <- "10-01" #last day of summer period


######### STEP 2) Create tables to compile results #########
threshold.per <- data.frame(STATION_NUMBER=NA, censored=NA, freq0=NA,
                    best.dist=NA, R2=NA,param1=NA, param2=NA, param3=NA,
                    ret.freq.flow=NA)
threshold.int <- data.frame(STATION_NUMBER=NA, duration=NA)
drought <- data.frame(STATION_NUMBER=NA, type= NA, year=NA, dr_events=NA,
                      dr_max_dur=NA, dr_days=NA)
trends <- data.frame(STATION_NUMBER=stations$STATION_NUMBER, type= NA,
                     hurdlechk=NA, hurdlemes=NA, negbinchk=NA, negbinmes=NA, 
                     slope=NA, intercept=NA, CATTrend=NA, years.for.trend=NA, 
                     test=NA, mapslope=NA)

######### STEP 3) Retrieve hydrometric data #########
for (j in 1:length(stations$STATION_NUMBER)) { #start loop that cycles through each station
  print(paste(j, "-", stations$STATION_NUMBER[j]))
  flow.daily <- hy_daily_flows(stations[j,1], start_date=paste0(min(yrs.of.int),"-01-01"),
                               end_date=paste0(max(yrs.of.int),"-12-31"))
  flow.daily$Year <- as.numeric(format(flow.daily$Date, "%Y")) #re-arrange year to faciliate sorting
  valid_year <- unique(flow.daily$Year) #note which years have data
  len.v.year <- length(valid_year)
  flow.daily$M_D <- substr(flow.daily$Date, 6, 10) #also reformat months and days for use later


######### STEP 4) Check if there is sufficient data #########
  # requirements for calculating thresholds are at least 150 days of data in at
  # least 20 years during the reference period
  ref <- flow.daily %>% filter(Year %in% yrs.of.ref)
  summ <- ref %>% group_by(Year) %>% summarise(count=length(Value[!is.na(Value)]))
  if (sum(summ$count>=150)>=20) { #use if loop to stop calculations in case of insufficient data
    ref_years <-  summ %>% filter(summ$count>=150) %>% pull(Year) #note which of the reference years are available for that station

    
######### STEP 5) Identify summer period #########
    # 3 criteria are used to identify summer:
    #     1) from May 1 to October 1
    #     2) ice-free period during this window
    #     3) period after freshet during this window
    # Start by creating tables to compile summer flow and rolling averages
    # first table is for summer flows each year
    days <- format(seq.Date(as.Date("05-01",format="%m-%d"), as.Date("10-01",
              format="%m-%d"), by = 1), format="%m-%d") #create vector with all possible days of summer (154)
    flow.summer <- data.frame(matrix(ncol=155,nrow=len.v.year))
    colnames(flow.summer) <- c("Year", days)
    flow.summer[,1] <- valid_year
    
    # next table is for average flows, with size of table dependant on length of 
    # averaging window
    max_w_per_year <- as.numeric(as.Date(summer.end, format="%m-%d") -
                        as.Date(summer.beg, format="%m-%d")) - n.day + 2
    windows <- format(seq.Date(as.Date(summer.beg, format="%m-%d")+(n.day-1),
                 as.Date(summer.end, format="%m-%d"), by = 1), format = "%m-%d")
    averages <- data.frame(matrix(ncol=max_w_per_year+2,nrow=len.v.year))
    colnames(averages) <- c("Year", windows, "ann.min")
    averages[,1] <- valid_year

    # Calculate median freshet date from 30-year median hydrograph as the date
    # with max flow between 1 march and 31 july
    med.freshet <- find.median.freshet(stations$STATION_NUMBER[j], min(yrs.of.ref), max(yrs.of.ref))

    # Find summer period for all years
    for (k in 1:len.v.year) { #cycle through all the years for that station
      #isolate data from all summer days (May 1-Oct 1)
      summer <- merge(days, flow.daily[(flow.daily$Year==valid_year[k] &
                  flow.daily$M_D>="05-01" & flow.daily$M_D<="10-01"),c(7,2,4,5)],
                  by.x="x", by.y="M_D", all.x=TRUE)
      #only keep flow values for ice-free period (as identified by "B" symbol in Hydat)
      summer$Value[!ice.free(summer$Date, summer$Symbol)] <- NA
      #only keep flow values after freshet
      summer$Value[!after.freshet(summer$x, summer$Value, med.freshet)] <- NA
      # compile the k year data into the larger table
      flow.summer[k,2:155] <- summer$Value


######### STEP 6A) For perennial watercourses: #########
#----- i. calculate rolling average for summer windows  -----
      for (l in 1:max_w_per_year) {
        if (all(!is.na(flow.summer[k,l:(l+n.day-1)]))) { # only calculate averages when there are values for all the days in the window
          averages[k,l+1] <- mean(as.numeric(flow.summer[k,l:(l+n.day-1)]))
        }
      }
      #annual minimum window
      if (sum(!is.na(averages[k,-1]))>0) { # only calculate for those years with at least one average value during summer
        averages$ann.min[k] <- min(averages[k,-1],na.rm=TRUE)
      }
    } #this closes the loop cycling through all the years for that station

#----- ii. define threshold by fitting curve to 30-year reference data -----
    curve.data <- ret.freq.curve(return.freq, averages$Year, averages$ann.min, ref_years)
    threshold.per[nrow(threshold.per)+1,] <- c(stations$STATION_NUMBER[j], curve.data)

#----- iii. use function from lfstat library to calculate drought events -----
    # This uses methodology described in World Meteorological Organization (2008)
    # Manual on low-flow estimation and prediction
    # if there is a threshold above detection limit, identify stations as perennial
    if ((curve.data[1]==FALSE | curve.data[1]==0) & !is.na(curve.data[8]) & curve.data[8]>=0.001) { 
      trends$type[j] <- "perennial"
      for (k in 1:len.v.year) { #cycle through all the years for that station
        if (sum(!is.na(flow.summer[k,-1]))>60) { #only do years with at least 60 calculated averages
          y_dr <- drought_days(Y=valid_year[k], M_D=days,  
                               flow=as.numeric(flow.summer[k,-1]),
                               dr_thresh=as.numeric(curve.data[8]),
                               stn=stations$STATION_NUMBER[j])
        } else { # for those years with insufficient data
          y_dr <- data.frame(STATION_NUMBER=stations$STATION_NUMBER[j],
                             type="perennial", year=valid_year[k], dr_events=NA, 
                             dr_max_dur=NA, dr_days=NA)
        }
        
        # add results for that station to table
        drought <- rbind(drought, y_dr)
      } # close loop for all years at that station
    } # close loop for perennial stations (those for which a curve could be calculated)
      

######### STEP 6B) For intermittent watercourses: #########
    # This uses methodology described in van Huijgevoort et al. (2012) A generic 
    # method for hydrological drought identification across different climate regions
    # Identify hydrometric stations on intermittent watercourses as any of those that:
    #     1) have more than 25% annual n-day minima as 0 during 30 reference
    #          years (censored=TRUE);
    #     2) have a return flow frequency of less than 0.001 m3/s (detection limit); or
    #     3) do not have a calculated return frequency flow.
    if (curve.data[1]==TRUE | curve.data[1]==1 | curve.data[8]<0.001 | is.na(curve.data[8])) {
      trends$type[j] <- "intermittent" # flag station as intermittent
      #----- i. calculate dry spell durations  -----
      no.flow.spells <- data.frame(matrix(ncol=79,nrow=len.v.year)) #number of columns is max possible number of dry spells (total summer days/2) + 1
      no.flow.spells[,1] <- valid_year

      for (k in 1:len.v.year) { #cycle through all the years for that station
        # the counter counts the number of consecutive days with flow of 0
        counter <- ifelse(!is.na(flow.summer$'05-01'[k]) &
                      flow.summer$'05-01'[k]==0, 1, 0)
        # column is to know where in the table the next total from counter goes
        column <- 2
        if (sum(!is.na(flow.summer[k,]))>1) { #only calculate for years with data
          for (t in 2:ncol(flow.summer)) { #cycle through each of the columns with flow data
            if (!is.na(flow.summer[k,t])){ #for those columns with a value
              if (flow.summer[k,t]==0){ #if the flow value is zero, increase the counter by one
                counter <- counter+1
              } else { #still for cases when column has flow value
                if (flow.summer[k,t]>0 & counter!=0) { #if value is above zero and counter is not zero, that denotes the end of a dry spell
                  no.flow.spells[k,column] <- counter #note the length of the dry spell
                  column <- column+1  #increase the column number not to overwrite for the next noting
                  counter <- 0 #set the counter back to zero
                }
              }
              # case when summer ends in drought or NA
              if ((flow.summer[k,ncol(flow.summer)]==0 |
                   is.na(flow.summer[k,ncol(flow.summer)])) & counter>0) {
                no.flow.spells[k,column] <- counter
              }
            }
          }
        }
      }

#----- ii. define threshold from cumulative distribution of dry spell duration  -----
      # create variable with all the length of no flow spells during the reference period
      dry.spells <- no.flow.spells %>% filter(no.flow.spells$X1 %in% ref_years &
                     !is.na(no.flow.spells$X2)) %>% dplyr::select(2:79) %>%
                     apply(1,function(x) na.omit(x)) %>%  unlist()
      # There needs to be a minimum of 11 dry spells during the 30 year reference
      # period to calculate a threshold
      if (length(dry.spells)>10) {
        Th <- quantile(dry.spells, probs=perc, na.rm=TRUE)
        threshold.int[nrow(threshold.int)+1,] <- c(stations$STATION_NUMBER[j], Th)

#----- iii. evaluate each year against threshold  -----
        # events are defined as no flow spells longer than the threshold, and this
        # compiles how many there are each year
        events <- rowSums(no.flow.spells[,-1]>Th, na.rm=TRUE)
        most <- max(events)
        if (most == 0) { # if there were no events the durations were 0
          days.above <- rep(0, length(events))
          longest.above <- rep(0, length(events))
        } else {
          # initialize data frame to compile duration of no flow events
          dry <- data.frame(matrix(nrow=length(events), ncol=most))
          for (i in 1:length(events)) {
            if (events[i]>0) { # for those years where there are events
              # note the number of days exceeding the threshold for each event
              dry[i,] <- no.flow.spells[i,which(no.flow.spells[i,-1]>Th)+1]-Th
              if (events[i]==1 & most>1) {dry[i,2:most] <- NA}
              if (events[i]==2 & most>2) {dry[i,3:most] <- NA}
            }
          }
          days.above <- round(rowSums(dry, na.rm=TRUE))
          longest.above <- round(do.call(pmax, c(dry, na.rm=TRUE)))
        }
        
        # add results for that station to table
        drought <- rbind(drought, data.frame(STATION_NUMBER=rep(stations$STATION_NUMBER[j],
                      len.v.year), type=rep("intermittent", len.v.year), year=valid_year,
                      dr_events=events, dr_max_dur=longest.above, dr_days=days.above))
      } #close loop for intermittent stations with threshold
    } # close loop for all intermittent stations

    
######### STEP 7) Calculate 50-year trends #########
    # Look for trends using 3 tests, when possible a hurdle test for those series
    # with 3 or more zeros, a negative binomial test otherwise, and finally a
    # stationarity check using all zeros or Wald-Wolfowitz test
    ## Minimum data requirements: some data 1970-1975, >=30 points & no gap over 10 years
    data <- drought %>% filter(drought$STATION_NUMBER==stations$STATION_NUMBER[j]) %>%
                dplyr::select(year, dr_days) %>% na.omit()
    goodyears <- data$year #find which years have data
    gap.check <- na.omit(goodyears - lag(goodyears)) #check for holes in data continuity
    if (all(any(goodyears %in% c(1970:1975)), (length(goodyears) >= 30),
            (max(gap.check) <= 11))){
      # calculate trends
      trends[j,3:12] <- identify_trends(x_var=data$year, y_var=data$dr_days)
    }
  } #close loop for calculations on those stations which have sufficient data
} #close loop cycling through all the stations


######### STEP 8) Test for outliers and eliminate them #########
# using Rosner's and Pettitt tests
out_filename <- "outliers_étiage"
outliers <- find_outliers(tr.stn=trends$STATION_NUMBER, tr.slope=trends$slope,
                          tr.int=trends$intercept, tr.test=trends$test, 
                          data.stn=drought$STATION_NUMBER, data.yr=drought$year, 
                          data.val=drought$dr_days, 
                          start.yr=1970, end.yr=2021,
                          file.name=out_filename)
# set the mapslope of outliers to NA so they are not plotted
for (m in 1:length(outliers)) { trends$mapslope[outliers[m]] <- NA }


######### STEP 9) Save output as csv files #########
# delete the first row of the threshold tables and drought tables as they will be NA
threshold.per <- threshold.per[-1,]
threshold.int <- threshold.int[-1,]
drought <- drought[-1,]
# harmonize TRUE/FALSE entries
threshold.per$censored[threshold.per$censored==0] <- FALSE
threshold.per$censored[threshold.per$censored==1] <- TRUE

# then write all the files, omitting the row numbers
write.csv(threshold.per, "./Output/CESIétiage-perennial thresholds.csv", row.names = FALSE)
write.csv(threshold.int, "./Output/CESIétiage-intermittent thresholds.csv", row.names = FALSE)
write.csv(drought, "./Output/CESIétiage-drought.csv", row.names = FALSE)
write.csv(trends, "./Output/CESIétiage-trends.csv", row.names = FALSE)


######### STEP 10) Krige results and produce map #########
surface <- krige_data_1zone(stns=trends$STATION_NUMBER, values=trends$mapslope,
                            file.name="variogram_1Z_étiage")

map_results_surface(input=surface, echelle=c(-0.3,0.3), divisions=6, high.col="orange", 
            leg.title="Trends in number of low flow days", map.title=NA, 
            file.name="CESIétiage - map of trends", type="PNG")
