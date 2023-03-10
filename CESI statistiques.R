#######################################################################
#
#    -###- CANADIAN ENVIRONMENTAL SUSTAINABILITY INDICATORS -###-
#              -- Summary statistics Calculator --
#
#  This project contains scripts used to automate the calculation
#                   of CESI summary statistics.
#
#  Environment and Climate Change Canada
#  Created: February 2023
#
#######################################################################

# This script calculates the summary statistics that are used in the CESI
# narrative. It relies on the output of three other CESI scripts: CESIétiage, 
# CESIcrue and CESIrendement. The output includes a text file with statistics
# included in the narrative text, a figure of the classification of yield over
# time at stations, 5 csv files for each of the tables in the document and a
# figure with station locations.


# The specific steps are:
#    1) Identify files to import
#    2) Define variables
#    3) Find percentage of flow classifications each year
#    4) Find percentage of trends
#    5) Location maps of stations
#    6) Trend stats for example stations
#
# Assuming there are two sub-folders in working directory: 
# 1. "Dependencies" - with following files
#             country boundaries [CanadaBound.shp]
#             outlines of big lakes [CANwaterbodies.shp]
#             CESI north arrow [esri6north.png]
# 2. "Output" - in which to save the output and find the output from previous scripts
# 3. "R" - in which 4 script files are located:
#             [CESI statistiques.R-THIS ONE]
#             [CESI data sorting functions.R]
#             [CESI trend functions.R]
#             [CESI mapping functions.R]

# Assuming hydat database has been downloaded and is available. If not, run:
# download_hydat()

#######################################################################

# load necessary libraries, the order in which you import them is important not
# to mask certain functions
library(tidyhydat) # for interacting with hydat database
library(sf)        # to create spatial objects
library(terra)     # to work with raster data
library(tidyverse) # to merge dataframes, includes library dplyr
library(png)       # to create image

source('./R/CESI data sorting functions.R')
source('./R/CESI trend functions.R')
source('./R/CESI mapping functions.R')

######### STEP 1) Identify files to import #########
yield <- read.csv("./Output/CESIrendement-yield.csv", header=TRUE) #yield classifications
trendsY <- read.csv("./Output/CESIrendement-trends.csv", header=TRUE) #trends in yield
HFdays <- read.csv("./Output/CESIcrue-flood.csv", header=TRUE) #trends in high flows
trendsHF <- read.csv("./Output/CESIcrue-trends.csv", header=TRUE) #trends in high flows
LFdays <- read.csv("./Output/CESIétiage-drought.csv", header=TRUE) #trends in low flows
trendsLF <- read.csv("./Output/CESIétiage-trends.csv", header=TRUE) #trends in low flows


######### STEP 2) Define variables #########
map.year <- 2021
yrs.for.class <- c(2001:2021) #years in which to classify yield
orange <- rgb(161,83,34, maxColorValue=255)
green <- rgb(133,161,66, maxColorValue = 255)
blue <- rgb(15,76,106, maxColorValue=255)
P_Tlist <- c("NL", "PE", "NS", "NB", "QC", "ON", "MB", "SK", "AB", "BC", "YT", 
             "NT", "NU")
y.ex.stns <- c("08GD004", "07JD002") # list of stations to calculate decadal trend stats
HF.ex.stns <- c("05BG006", "02PJ007") # list of stations to calculate decadal trend stats
LF.ex.stns <- c("02NE011", "05OB016") # list of stations to calculate decadal trend stats

# prepare text file to compile results
txt.file <- file("./Output/CESIstatistiques.txt")
writeLines(c("CESI summary statistics", format(Sys.Date(), format="%B %d %Y"), "-------"), txt.file)
close(txt.file)
txt.file <- "./Output/CESIstatistiques.txt"


######### STEP 3) Find percentage of flow classifications each year #########
# create table to compile H-N-L results per year
yieldnoNA <- yield %>% filter(!is.na(flow_status))
H_N_L <- data.frame(Year=yrs.for.class, num_stns=NA, High=NA, Normal=NA, Low=NA)
for (j in 1:length(yrs.for.class)) {
  y <- yieldnoNA %>% filter(Year==yrs.for.class[j])
  H_N_L$num_stns[j] <- nrow(y)
  H_N_L$High[j] <- round(sum(y$flow_status=="high")/nrow(y),2)*100
  H_N_L$Normal[j] <- round(sum(y$flow_status=="normal")/nrow(y),2)*100
  H_N_L$Low[j] <- round(sum(y$flow_status=="low")/nrow(y),2)*100
}
write.csv(H_N_L, "./Output/T_A1-water quantity at monitoring stations.csv", row.names=FALSE)

# create a graph with the results
png(file="./Output/H-N-L percentage of stations.png", width=5, height=3, units="in", res=600)
par(mar=c(2, 2, 2, 4))
plot.new()
plot.window(xlim=c(min(yrs.for.class),max(yrs.for.class)), ylim=c(0,100))
axis(side=1, at=seq(min(yrs.for.class),max(yrs.for.class),2), tck=FALSE, 
     cex.axis=0.5, mgp=c(1,0,0))
axis(side=2, at=seq(0,100,20), tck=FALSE, las=1, lty=0, cex.axis=0.5,
     mgp=c(1, 0.5, 0))
abline(h=seq(20,100,20), col="gray")
points(H_N_L$Year, H_N_L$High, col=blue , type='o', pch=16, lwd=1.5, cex=0.9)
points(H_N_L$Year, H_N_L$Normal, col=green , type='o', pch=16, lwd=1.5, cex=0.9)
points(H_N_L$Year, H_N_L$Low, col=orange , type='o', pch=16, lwd=1.5, cex=0.9)
loc <- par("usr")
text(loc[1], loc[4], xpd=TRUE, "                Percentage of stations", pos = 3, cex=0.5)
legend("right", inset=c(-0.2,0), xpd=TRUE, cex=0.5, legend=c("High", "Normal", "Low"), 
       bty="n", pch=16, col=c(blue, green, orange), pt.cex=0.9, lty=1, lwd=1.5)
dev.off()

# note the values for the map year in the text file
write(c(paste("National water quantity in", map.year),
  paste(H_N_L$num_stns[H_N_L$Year==map.year], "stations were used for", map.year, "statistics"),
  paste("higher than normal at", H_N_L$High[H_N_L$Year==map.year], "% of stations"),
  paste("normal at", H_N_L$Normal[H_N_L$Year==map.year], "% of stations"),
  paste("lower than normal at", H_N_L$Low[H_N_L$Year==map.year], "% of stations"),
  "-------"), file=txt.file, append=TRUE)

# create table to compile H-N-L results by province/territory in map year and save csv
y <- yield %>% filter(Year==map.year, !is.na(flow_status))
PT <- stn.locations(y$STATION_NUMBER)
out <- PT %>% group_by(PT) %>% count()
write.csv(out, paste0("./Output/T_1-number water quantity monitoring stations-", 
            map.year, ".csv"), row.names=FALSE)


######### STEP 4) Find percentage of trends #########
trends <- list(trendsY, trendsHF, trendsLF)
trends.txt <- c("trends in yield", "trends in high flow", "trends in low flow")

# loop through each of the trend files
for (t in 1:length(trends)) {
  tr_perc <-  data.frame(P_T=P_Tlist, num_stns=rep(0, length(P_Tlist)), 
                increase=rep(0, length(P_Tlist)), stable=rep(0, length(P_Tlist)),
                decrease=rep(0, length(P_Tlist)), uncertain=rep(0, length(P_Tlist)))
  # filter for only stations with sufficient data to test for trends and find location
  tr <-  trends[[t]] %>% filter(!is.na(test))
  tr <- merge(stn.locations(tr$STATION_NUMBER), tr, by="STATION_NUMBER")
  # loop through each province and territory
  for (pt in P_Tlist) {
    n_stns <- ifelse((nrow(filter(tr,PT==pt))>0),nrow(filter(tr,PT==pt)),0)
    if (n_stns>0) {
      tr_perc$num_stns[tr_perc$P_T==pt] <- n_stns
      tr_perc$increase[tr_perc$P_T==pt] <- round((nrow(filter(tr, PT==pt, mapslope>0))/n_stns)*100)
      tr_perc$stable[tr_perc$P_T==pt] <- round((nrow(filter(tr, PT==pt, mapslope==0))/n_stns)*100)
      tr_perc$decrease[tr_perc$P_T==pt] <- round((nrow(filter(tr, PT==pt, mapslope<0))/n_stns)*100)
      tr_perc$uncertain[tr_perc$P_T==pt] <- round((nrow(filter(tr, PT==pt, is.na(mapslope)))/n_stns)*100)
    }
  }
  # note output
  write(paste(nrow(tr),"stations have calculated", trends.txt[t]),
        file=txt.file, append=TRUE)
  write.csv(tr_perc, paste0("./Output/T_A", t+1, "-", trends.txt[t], 
            " percentages.csv"), row.names=FALSE)
}


######### STEP 5) Location map of stations #########
# map combining location of stations used for national statistics and trend
# calculations

# call on a mapping function to make station location figure for national statistics
map_stn_loc_natsts(stns=yield$STATION_NUMBER[yield$Year==map.year & !is.na(yield$flow_status)], 
                  file.name=paste("national stats station locations", map.year), 
                  type="PNG", north_a=TRUE, scale_bar=TRUE)

# make a list of stations that have been used for trends
tY <- trendsY %>% filter(!is.na(mapslope)) %>% select(STATION_NUMBER)
tHF <- trendsHF %>% filter(!is.na(mapslope)) %>% select(STATION_NUMBER)
tLF <- trendsHF %>% filter(!is.na(mapslope)) %>% select(STATION_NUMBER)
tSTNS <- list(tY, tHF, tLF) %>% reduce(full_join, by="STATION_NUMBER")     
# call on a mapping function to make station location figure for trends
map_stn_loc_trs(stns=tSTNS$STATION_NUMBER, file.name="trends station locations", 
                type="PNG", north_a=TRUE, scale_bar=TRUE)


######### STEP 6) Trend stats for example stations #########
# open pdf file in which to save results
pdf(file="./Output/individual station trend stats.pdf", width=8.5, height=11)
par(mfrow = c(4,3))

# plot stations from each trend
ex_stn_trend_stats(stn.list=y.ex.stns, tr.type="yield", tr.stn=trendsY$STATION_NUMBER, 
                tr.slope=trendsY$mapslope, tr.int=trendsY$intercept, 
                data.stn=yield$STATION_NUMBER, data.yr=yield$Year, 
                data.val=yield$ann_mean_yield)
ex_stn_trend_stats(stn.list=HF.ex.stns, tr.type="high flow days", 
                tr.stn=trendsHF$STATION_NUMBER, tr.slope=trendsHF$mapslope, 
                tr.int=trendsHF$intercept, data.stn=HFdays$STATION_NUMBER, 
                data.yr=HFdays$Year, data.val=HFdays$pot_days)
ex_stn_trend_stats(stn.list=LF.ex.stns, tr.type="low flow days", 
                tr.stn=trendsLF$STATION_NUMBER, tr.slope=trendsLF$mapslope, 
                tr.int=trendsLF$intercept, data.stn=LFdays$STATION_NUMBER, 
                data.yr=LFdays$year, data.val=LFdays$dr_days)

# close file
dev.off()