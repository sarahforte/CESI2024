#######################################################################
#
#    -###- CANADIAN ENVIRONMENTAL SUSTAINABILITY INDICATORS -###-
#              -- Water Quantity Indicator Calculator --
#
#  This project contains scripts used to automate the calculation
#          of CESI flood/very high flow indicator.
#
#  Environment and Climate Change Canada
#  Created: February 2023
#
#######################################################################

# This script calculates the very high flow days during the ice-free period at 
# hydrometric stations based on thresholds defined from a 30 year reference period.
# Trends over 50 years are calculated and the results are displayed as a kriged map.
# The specific steps are:
#    1) Define variables
#    2) Create tables to compile results
#    3) Retrieve hydrometric data from hydat database
#    4) Calculate thresholds
#    5) Identify ice-free period
#    6) Calculate flood metrics
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
#             [CESI crue.R-THIS ONE]
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
library(tidyverse) # to merge dataframes
library(zoo)       # For working with time series
library(evir)      # For POT
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
# RHBN stations are used
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


######### STEP 2) Create tables to compile results #########
flood <- data.frame(STATION_NUMBER=NA, Year=NA, pot_threshold=NA, pot_days=NA,
                    pot_events=NA, pot_max_dur=NA, X1_day_max=NA)
trends <- data.frame(STATION_NUMBER=stations$STATION_NUMBER, hurdlechk=NA, 
                     hurdlemes=NA, negbinchk=NA, negbinmes=NA, slope=NA, 
                     intercept=NA, CATTrend=NA, years.for.trend=NA, test=NA, 
                     mapslope=NA)


######### STEP 3) Retrieve hydrometric data #########
for (j in 1:length(stations$STATION_NUMBER)) {
  station <- stations$STATION_NUMBER[j]
  print(paste(j, "-", station))
  flow.daily <- hy_daily_flows(station, start_date=paste0(min(yrs.of.int),"-01-01"),
                               end_date=paste0(max(yrs.of.int),"-12-31"))
  flow.daily$Year <- as.numeric(format(flow.daily$Date, "%Y")) #re-arrange year to faciliate sorting
  valid_year <- unique(flow.daily$Year) #note which years have data

######### STEP 4) Calculate thresholds #########
# requirements for calculating thresholds are at least 150 days of data in at
# least 20 years during the reference period
  ref <- flow.daily %>% filter(Year %in% yrs.of.ref)
  summ <- ref %>% group_by(Year) %>% summarise(count=length(Value[!is.na(Value)]))
  if (sum(summ$count>=150)>=20) { #use if loop to stop calculations in case of insufficient data
    thresholds <- Threshold_Build(station, yfrom=min(yrs.of.ref), 
                                  yto=max(yrs.of.ref), hi=0.95, low=0.05)
    # confirm a threshold was calculated to continue
    if (!is.na(thresholds$Q_Upper)) {
      # also calculate statistic to check for the correlation between flood events
      stn_area <- hy_stations(station) %>% pull(DRAINAGE_AREA_GROSS)
      R <- 5 + log(stn_area/(1.609^2)) %>% round(1)
      # confirm a threshold was calculated to continue, sometimes no station area is availble
      if (!is.na(R)) {
        
######### STEP 5) Identify ice-free period #########
        for (yr in valid_year) { #cycle through all the years for that station
          flow <- flow.daily %>% filter(Year==yr)
          #only keep flow values for ice-free period (as identified by "B" symbol in Hydat)
          flow$Value[!ice.free(flow$Date, flow$Symbol)] <- NA

######### STEP 6) Calculate flood metrics #########
# calculate using peaks over threshold method. For each year of interest:
#     - maximum flow,
#     - a threshold at 95th percentile of flow in reference period,
#     - number of days over threshold flow,
#     - number of events, and
#     - maximum duration of events
          flood[nrow(flood)+1,] <- c(station, yr, thresholds$Q_Upper, 
                  get_flood_metrics(flow$Date, flow$Value, thresholds$Q_Upper, R))
        } # close loop cycling through years

######### STEP 7) Calculate 50-year trends #########
# Look for trends using 3 tests, when possible a hurdle test for those series
# with 3 or more zeros, a negative binomial test otherwise, and finally a
# stationarity check using all zeros or Wald-Wolfowitz test
## Minimum data requirements: some data 1970-1975, >=30 points & no gap over 10 years
        data <- flood %>% filter(flood$STATION_NUMBER==stations$STATION_NUMBER[j]) %>%
                    dplyr::select(Year, pot_days) %>% na.omit()
        goodyears <- data$Year #find which years have data
        gap.check <- na.omit(goodyears - lag(goodyears)) #check for holes in data continuity
        if (all(any(goodyears %in% c(1970:1975)), (length(goodyears) >= 30),
              (max(gap.check) <= 11))){
          # calculate trends
          trends[j,2:11] <- identify_trends(x_var=data$Year, y_var=data$pot_days)
        }
      } # close loop for all stations with calculated Q_Upper threshold
    } # close loop for all stations with calculated R threshold
  } #close loop for calculations on those stations which have sufficient data
} #close loop cycling through all the stations


######### STEP 8) Test for outliers and eliminate them #########
# using Rosner's and Pettitt tests
out_filename <- "outliers_crue"
outliers <- find_outliers(tr.stn=trends$STATION_NUMBER, tr.slope=trends$mapslope,
                          tr.int=trends$intercept, tr.test=trends$test, 
                          data.stn=flood$STATION_NUMBER, data.yr=flood$Year, 
                          data.val=flood$pot_days, 
                          start.yr=1970, end.yr=2021,
                          file.name=out_filename)
# set the mapslope of outliers to NA so they are not plotted
for (m in 1:length(outliers)) { trends$mapslope[outliers[m]] <- NA }


######### STEP 9) Save output as csv files #########
# delete the first row of the flood table as it will be NA
flood <- flood[-1,]

# then write all the files, omitting the row numbers
write.csv(flood, "./Output/CESIcrue-flood.csv", row.names = FALSE)
write.csv(trends, "./Output/CESIcrue-trends.csv", row.names = FALSE)


######### STEP 10) Krige results and produce map #########
surface <- krige_data_1zone(stns=trends$STATION_NUMBER, values=trends$mapslope,
                            file.name="variogram_1Z_crue")

map_results_surface(input=surface, echelle=c(-0.3,0.3), divisions=6, high.col="blue",
            leg.title="Trends in number of high flow days", map.title=NA, 
            file.name="CESIcrue - map of trends", type="PNG")
