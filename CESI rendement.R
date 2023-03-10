#######################################################################
#
#    -###- CANADIAN ENVIRONMENTAL SUSTAINABILITY INDICATORS -###-
#              -- Water Quantity Indicator Calculator --
#
#  This project contains scripts used to automate the calculation
#          of CESI yield (annual  water quantity) indicator.
#
#  Environment and Climate Change Canada
#  Created: February 2023
#
#######################################################################

# This script calculates the yield during the ice-free period at hydrometric 
# stations. Three different products are created with the yield. The first is
# done for all stations with sufficient data and the second two only for 
# stations included in the RHBN network
# I. Map where yield at each station is classified as high-normal-low when
#     compared to normal from 30 year reference period. Map represents single 
#     year.OUTPUT step #9
# II. Maps of yield at each station for a year is ranked as percentile compared
#       to normal from 30 year reference period and gridded to be plotted as a
#       surface map. A series of maps for years of interest is created to be 
#       assembled as a video.OUTPUT step #10
# III. Gridded surface map of trend in yield at each station over 50 years. 
#         OUTPUT step #11
#
# The specific steps are:
#    1) Define variables
#    2) Create tables to compile results
#    3) Calculate annual mean flow and yield
#    4) Calculate thresholds
#    5) Find flow classification and percentile ranking of yield
#    6) Calculate 50-year trends
#    7) Identify outlier trends and eliminate them
#    8) Save output as csv files
#    9) Make station map for high-normal-low classification 
#    10) Krige percentile results and produce mapsKrige results and produce map
#    11) Krige trend results and produce map
#
# Assuming there are three sub-folders in working directory: 
# 1. "Dependencies" - with following files
#             country boundaries [CanadaBound.shp]
#             outlines of big lakes [CANwaterbodies.shp]
#             CESI north arrow [esri6north.png]
#             curly brackets to interpret ranks [curly_brackets.png]
#             arrows to interpret trends [up_level_down_arrows.png]
#             kriging zones [5 merged ecozones.shp] (only if regional kriging to be used)
# 2. "Output" - in which to save the output
# 3. "R" - in which 5 script files are located:
#             [CESI rendement.R-THIS ONE]
#             [CESI data sorting functions.R]
#             [CESI indicator functions.R]
#             [CESI trend functions.R]
#             [CESI mapping functions.R]

# Assuming hydat database has been downloaded and is available. If not, run:
# download_hydat()

#######################################################################

# load necessary libraries, the order in which you import them is important not
# to mask certain functions
library(readxl)    # for importing excel file
library(tidyhydat) # for interacting with hydat database
library(tidyverse) # to merge dataframes
library(zoo)       # For working with time series
library(evir)      # For POT
library(purrr)     # function 'possibly' is easiest way to deal with errors
library(zyp)       # for function in Mann-Kendall test
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

yrs.of.ref <- c(1991:2020) #30-year reference period
map.year <- 2021 #year for water quantity at monitoring stations map
yrs.for.class <- c(2000:2021) #years in which to classify yield
yrs.for.trend <- c(1970:2021) #period for calculating trends
smallest.min <- min(min(yrs.of.ref), map.year, min(yrs.for.class),min(yrs.for.trend))
biggest.max <- max(max(yrs.of.ref), map.year, max(yrs.for.class),max(yrs.for.trend))
yrs.range <- c(smallest.min:biggest.max)
time <-  365*24*3600 #number of seconds in a year
thres.high <- 0.85
thres.low <- 0.15

# Create list of hydrometric station based on years of interest, and note operation
# schedule
stations <- hy_stn_data_range() %>% filter(DATA_TYPE=="Q", Year_to>=smallest.min, 
                                           RECORD_LENGTH>=25)
op <- hy_stn_data_coll(stations$STATION_NUMBER) %>% filter(DATA_TYPE=="Flow") %>% 
  group_by(STATION_NUMBER) %>% summarize(OPERATION=last(OPERATION))
stations <- merge(stations, op, by="STATION_NUMBER")

# import list of RHBN stations
url <- "https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/RHBN/RHBN_Metadata.xlsx"
destfile <- "./Output/RHBN_Metadata.xlsx"
download.file(url, destfile, method="curl")
RHBN <- read_xlsx(destfile, range="A3:P1286")
RHBN <- RHBN[-(1),]
RHBN <- filter(RHBN, Evaluation_Year==2020, DATA_TYPE=="Q")

######### STEP 2) Create tables to compile results #########
yield <- data.frame(STATION_NUMBER=NA, Year=NA, ann_mean_flow=NA, 
                    ann_mean_yield=NA, flow_status=NA, rank=NA)

thresholds <- data.frame(STATION_NUMBER=stations$STATION_NUMBER,  Q_Upper= NA, 
                         Q_Lower=NA)

trends <- data.frame(STATION_NUMBER=RHBN$STATION_NUMBER,  test= NA, slope=NA, 
                     intercept=NA, CATTrend=NA, years.for.trend=NA)


for (j in 1:length(stations$STATION_NUMBER)) {
  station <- stations$STATION_NUMBER[j]
  print(paste(j, "-", station))
  
  # get station area for yield calculation
  stn_area <- hy_stations(station) %>% pull(DRAINAGE_AREA_GROSS)
  # check if there is a station area before continuing
  if (!is.na(stn_area)) {
    
    
######### STEP 3) Calculate annual mean flow and yield #########
    y_rows <- nrow(yield)
    flow <- get_ann_mean_flow(station, year=yrs.range)
    yield[(y_rows+1):(y_rows+nrow(flow)),1:3] <- flow
    yield$ann_mean_yield[(y_rows+1):(y_rows+nrow(flow))] <- 
      yield$ann_mean_flow[(y_rows+1):(y_rows+nrow(flow))]*time/(stn_area*10^3)


######### STEP 4) Calculate thresholds #########
# requirements for calculating thresholds are at least 20 years with sufficient
# data during the reference period, which is 174 days for seasonal stations and
# 292 days for continuous stations
    # 
    # if there is sufficient flow data
    if (sum(!is.na(yield$ann_mean_flow[yield$STATION_NUMBER== station & 
                yield$Year %in% yrs.of.ref]))>20) {
      # retrieve and sort data
      cutoff <- case_when(stations$OPERATION[j]=="Continuous" ~ 292,
                          stations$OPERATION[j]=="Seasonal" ~ 174) 
      flow.daily <- hy_daily_flows(station, start_date=paste0(min(yrs.of.ref),"-01-01"),
                        end_date=paste0(max(yrs.of.ref),"-12-31"))
      flow.daily$Year <- as.numeric(substr(flow.daily$Date,1,4))
      summ = flow.daily %>% group_by(Year) %>% summarise(count=length(Value[!is.na(Value)]))
      # determine if there is sufficient data in each year to proceed
      if (sum(summ$count>=cutoff) >= 20) {
        ref_period <- yield %>% filter(STATION_NUMBER==station & Year %in% yrs.of.ref) %>%
          dplyr::select(Year, ann_mean_yield)
        # remove years with insufficient data
        if (sum(summ$count<cutoff) > 0) {
          list <- summ$Year[summ$count<cutoff]
          ref_period <- ref_period[-which(flow.daily$Year %in% list),]
        }
        # calculate threshold
        thresholds[j,2:3] <- quantile(ref_period$ann_mean_yield, probs=c(thres.high,
                                      thres.low), na.rm=TRUE)
      
        
######### STEP 5) Find flow classification and percentile ranking of yield #########
        for (yr in yrs.for.class) {
          # find line of table corresponding to year
          line <- which(yield$STATION_NUMBER==station & yield$Year==yr)
          # check if there was a flow calculated for that year
          if (!is.na(yield$ann_mean_flow[line])) {
            # classify flow
            yield$flow_status[line] <- case_when(
                yield$ann_mean_yield[line] > thresholds$Q_Upper[j] ~ "high",
                (yield$ann_mean_yield[line] <= thresholds$Q_Upper[j]) &
                  (yield$ann_mean_yield[line] >= thresholds$Q_Lower[j]) ~ "normal",
                yield$ann_mean_yield[line] < thresholds$Q_Lower[j] ~ "low")
        
            # calculate rank for RHBN stations
            if (station %in% RHBN$STATION_NUMBER) {
              yield$rank[line] <- percentile.ranked(ref_period$ann_mean_yield, 
                                                  yield$ann_mean_yield[line])
            }
          } # close loop for years with a calculated flow
        } # close loop cycling through all years


######### STEP 6) Calculate 50-year trends for RHBN stations#########
# Look for trends using 2 tests; start with Mann-Kendall test and if results are
# inconclusive do stationarity check using all zeros or Wald-Wolfowitz test
## Minimum data requirements: some data 1970-1975, >=30 points & no gap over 10 years
        if (station %in% RHBN$STATION_NUMBER) {
          jj <- which(trends$STATION_NUMBER==station)
          data <- yield %>% filter(STATION_NUMBER==station) %>% 
                    dplyr::select(Year, ann_mean_yield) %>% na.omit()
          goodyears <- data$Year #find which years have data
          gap.check <- na.omit(goodyears - lag(goodyears)) #check for holes in data continuity
          if (all(any(goodyears %in% c(1970:1975)), (length(goodyears) >= 30),
                (max(gap.check) <= 11))){
            # calculate trends - Mann-Kendall test
            print("Mann-Kendall test for trend")
            trends[jj,2:6] <- c("M-K", mk_test(data$ann_mean_yield, data$Year))
            # if results are uncertain, check for stationarity
            if (is.na(trends$CATTrend[jj]) | (!is.na(trends$CATTrend[jj]) & 
                    trends$CATTrend[jj]=="Uncertain")) {
              print("stationarity check")
              trends[jj,2:6] <- c("W-W", stationarity_test(data$ann_mean_yield)[1:4])
            }
          } # close loop of stations with sufficient data for trend calculations
        } # close loop of RHBN stations
      } # close loop for stations where threshold can be calculated
    } # close loop for stations with sufficient flow data
  } # close loop for stations with surface area available
} # close loop for all stations
  
######### STEP 7) Test for outliers in trends and eliminate them #########
# using Rosner's and Pettitt tests
trends$mapslope <- trends$slope
out_filename <- "outliers_rendement"
outliers <- find_outliers(tr.stn=trends$STATION_NUMBER, tr.slope=trends$mapslope,
                          tr.int=trends$intercept, tr.test=trends$test, 
                          data.stn=yield$STATION_NUMBER, data.yr=yield$Year, 
                          data.val=yield$ann_mean_yield, 
                          start.yr=1970, end.yr=2021,
                          file.name=out_filename)
# set the mapslope of outliers to NA so they are not plotted
for (m in 1:length(outliers)) { trends$mapslope[outliers[m]] <- NA }


######### STEP 8) Save output as csv files #########
# delete the first row of the yield table as it will be NA
yield <- yield[-1,]

# then write all the files, omitting the row numbers
write.csv(yield, "./Output/CESIrendement-yield.csv", row.names = FALSE)
write.csv(trends, "./Output/CESIrendement-trends.csv", row.names = FALSE)

    
######### STEP 9) Make station map for high-normal-low classification #########
# only make for single year
map_results_HNL(stns=yield$STATION_NUMBER[yield$Year==map.year], 
                status=yield$flow_status[yield$Year==map.year], map.title=NA, 
                file.name=paste0("CESIrendement-flow status-",map.year),
                type="PNG")


######### STEP 10) Krige percentile results and produce maps #########
if (!dir.exists("./Output/yield rank series")) {dir.create("./Output/yield rank series")}
for (yr in yrs.for.class) {
  message(paste("Creating map for", yr))
  data <- yield %>% filter(Year==yr) %>% dplyr::select(STATION_NUMBER, rank) %>% na.omit()
  surface <- krige_data_1zone(stns=data$STATION_NUMBER, values=data$rank,
                  file.name=paste0("yield rank series/variogram_1Z_rendement_rank_",yr))
  map_results_contour(input=surface, high.col="blue", map.title=yr,
                      file.name=paste0("yield rank series/CESIrendement-rank of yields-", yr),
                      type="PNG", legend_interp=TRUE, north_a=TRUE, scale_bar=TRUE)
}


######### STEP 11) Krige trend results and produce map #########
surface <- krige_data_1zone(stns=trends$STATION_NUMBER, values=trends$mapslope,
                            file.name="variogram_1Z_rendement_trend")
#surface5 <- krige_data_5zones(stns=trends$STATION_NUMBER, values=trends$mapslope,
#                            file.name="variogram_5Z_rendement_trend")
#surfaceIDW <- idw_data(stns=trends$STATION_NUMBER, values=trends$mapslope,
#                       power = 5)

map_results_surface(input=surface, echelle=c(-2,2), divisions=10, high.col="blue",
                    leg.title="Trends in millimeters per year", map.title=NA, 
                    file.name="CESIrendement - map of trends1Z", type="PNG")
