# functions for CESI - to create surfaces and maps
#  1) krige_data_5zones
#  2) krige_data_1zone
#  3) map_results_surface
#  4) map_results_HNL
#  5) idw_data
#  6) map_results_contour
#  7) map_stn_loc_natsts
#  8) map_stn_loc_trs

# All functions require some files, expected in the folder "./Dependencies"
# These are:  CanadaBound.shp (1,2,3,4,5,6) > boundary of Canada, territories and provinces
#             5 merged ecozones.shp (1) > merged ecozone boundaries for kriging
#             CANwaterbodies.sho (3,4,6) > main waterbodies in Canada to add to map
#             up_level_down_arrows.png (3) > arrows to explain trends in legend
#             esri6north.png (3,4,6) > north arrow used by IID to put on maps
#             curly_brackets.png (6) > curly brackets to explain percentiles

#--------------------------------------------------
### 1) Function to Krige data using 5 isotropic variograms for 5 combined ecozones
##### Input:  stns = list of station numbers
#####         values = list of same length as stns with values to be kriged
#####         file = name of output file
##### Output: pdf file with specified name showing the 5 kriging zones and the
#####           variograms calculated for each of these zones
#####         grid of kriged values

krige_data_5zones <- function(stns, values, file.name) {
  # create data frame
  data <- data.frame(STATION_NUMBER=stns, value=values) %>% na.omit()
  
  # read in shapefile with canadian boundary to set canvas boundary
  Canada <- st_read(dsn="./Dependencies", "CanadaBound")
  map_crs <- st_crs(Canada)
  
  # read in shapefiles with kriging zones which are grouped ecozones
  Kzones <- st_read(dsn="./Dependencies", "5 merged ecozones")
  Kzones$ZONE_ID <- 1:nrow(Kzones)
  Kzones <- st_transform(Kzones, map_crs)
  
  # get station coordinates
  stn_coor <- as.data.frame(hy_stations(data$STATION_NUMBER) %>% dplyr::select(STATION_NUMBER, LATITUDE, LONGITUDE))
  data <- merge(data, stn_coor, by="STATION_NUMBER") #combine data with station coordinates
  data <- st_as_sf(data, coords=c("LONGITUDE", "LATITUDE"), crs="+proj=longlat +datum=WGS84")
  data <- st_transform(data, map_crs)
  
  # create grid with 10km cells covering the whole country onto which data
  # should be interpolated
  Cangrid <- st_make_grid(Canada, cellsize=10000, crs=map_crs, what="centers")
  grids <- list()

  # variogram model curve types
  Mtype <- c("Sph", "Exp")
  Mcol <- c("red", "green")

  # create pdf file to populate with variograms
  pdf(file = paste0("./Output/", file.name, ".pdf"), width = 8.5, height = 11)
  par(mfrow = c(5,4))
  
  # add map of zones to show which are for which variogram
  # transform zone outline and create labels
  par(mar = c(0, 0, 0, 0))
  Zoutline <- Kzones$geometry
  labels <- st_centroid(Zoutline)
#  labels@coords[1,2] <- 100000
  ID <- 1:length(labels)
  plot(Zoutline, col=rainbow(5)[ID], border="black", lwd = 0.4, lty=2)
  plot(Canada$geometry, add=TRUE, col=NA, border = "grey45")
  text(st_coordinates(labels)[,1],st_coordinates(labels)[,2], ID, cex=1)
  title(main="Grouped ecozones", line=-1)
  

  # work on each kriging zone separately
  for (z in 1:length(Zoutline)) {
    message(paste("kriging zone", z))
    ## create buffer around zone
    Szone <- Kzones[Kzones$ZONE_ID==z,]
    Bzone <- st_buffer(Szone, dist=300000, nQuadSegs=100) # buffer zone is 300 km
    Mzone <- st_buffer(Szone, dist=25000, nQuadSegs=100) # overlap zone is 50 km, 25 on each side
    ## extract stations in buffered zone
    Zstns <- st_intersection(Bzone, data)

    ## find average distance between stations in zone
    dist <- as.data.frame(st_distance(Zstns))
    dist[upper.tri(dist, diag=TRUE)] <- NA 
    MDstns <- mean(apply(dist[-1,],1,FUN=min, na.rm=TRUE))
    
    ## create variogram
    var <- variogram(Zstns$value~1, Zstns, width=(MDstns/3), cutoff=500000) #calculates sample variogram values with 500km max range and lag distance proportional to the mean min distance between points
    
    ## two options for fitting variogram model, try automap version first, and if
    ## it creates output with nugget of 0, try gstat version.
    ### automap fitting - entirely independent including creation of new variogram,
    ### including determining which model of exponential or spherical fits best
    out <- autofitVariogram(Zstns$value~1, as(Zstns, "Spatial"), model=c("Exp", "Sph"), width= MDstns, cutoff=500000)# the best fitting model from exponential or spherical is chosen
    v_mod <- out$var_model
    if (v_mod$psill[1]==0) {
      ### gstat fitting - works best with start values, so calculate these firt
      ### create new variogram
      v <- variogram(Zstns$value~1, Zstns)
      ### note data from variogram as start values 
      nugget0 <- min(v$gamma)
      psill0 <- ((max(v$gamma)+median(v$gamma))/2)-nugget0
      range0 <- 0.1*sqrt((st_bbox(Zstns)[1] - st_bbox(Zstns)[3])^2+
                         (st_bbox(Zstns)[2] - st_bbox(Zstns)[4])^2)
      ### find best model between exponential and spherical
      vmods <- list()
      SSErr <- rep(NA, length(Mtype))
      for (k in 1:length(Mtype)) {
        vmods[[k]] <- fit.variogram(v, vgm(psill0, model=Mtype[k], range0, nugget0))
        SSErr[k] <- attr(vmods[[k]], "SSErr")
      }
      k <- which(SSErr==min(SSErr))
      ### if alternate model is better, use it
      if (vmods[[k]]$psill[1] > model$psill[1]) { 
        message("'gstat' variogram fitting has better results that 'automap' fitting")
        v_mod <- vmods[[k]]
      }
    }
  
    ## sort variogram data to plot, and change units from m to km
    dist <- var$dist/1000
    vario <- var$gamma
    ## plot variogram first
    par(mar = c(4, 4, 2, 1))
    plot(dist, vario, xlab="", ylab="semivariance", xlim=c(0,500), ylim=c(0,
      (ceiling(max(var$gamma, na.rm=TRUE)*100))/100), las=1, col="blue")
    title(main=paste0("Variogram zone ",z), line=0.5)
    title(xlab="distance [km]", line=2)
    ## sort model data to plot
    dist0 <- seq(0, max(dist), length.out=100)
    nugget <- v_mod$psill[1]
    sill <- sum(v_mod$psill[1:2])
    range <- v_mod$range[2]/1000
    k <- case_when(v_mod$model[2]=="Exp" ~ 2, #exponential
                   v_mod$model[2]=="Sph" ~ 1) #spherical
    MV <- case_when(k==2 ~ (sill-nugget)*(1-exp(-dist0/range))+nugget, #exponential
                    k==1 ~ ifelse(dist0>range,sill,(sill-nugget)*(1.5*dist0/range-0.5*(dist0/range)^3)+nugget)) #spherical
    rangeC <- case_when(k==2 ~ 3*range, #exponential
                        k==1 ~ range) #spherical
    ## plot model and write its parameters
    lines(dist0, MV, col=Mcol[k])
    mtext(paste0("model: ", v_mod$model[2], "   "), side=1, line=-4, adj=1, cex=0.4, col=Mcol[k])
    mtext(paste0("nugget: ", round(nugget, 3), "   "), side=1, line=-3, adj=1, cex=0.4, col=Mcol[k])
    mtext(paste0("sill: ", round(sill, 3), "   "), side=1, line=-2, adj=1, cex=0.4, col=Mcol[k])
    mtext(paste0("range: ", round(rangeC, 0), "   "), side=1, line=-1, adj=1, cex=0.4, col=Mcol[k])
  
  
    ## krige
    OK <- krige(formula=value~1, locations=Zstns, newdata=Cangrid, model=v_mod, 
                maxdist=500000, nmin=8)

    ## transform results into  a raster
    OKraster <- rast(data.frame(st_coordinates(OK),OK$var1.pred), crs=map_crs$wkt)
    grids[[z]] <- mask(OKraster, vect(Mzone))
  }
  
  # combine grids from all the zones, using the mean for the 50 km (25 km on each side) overlap
  OKmosaic <- do.call(mosaic, grids)
  output <- mask(OKmosaic, vect(Canada))
  
  # note limits
  zmax <- minmax(output)[2]
  zmin <- minmax(output)[1]
  par(mar = c(0, 0, 0, 0))
  plot(Canada$geometry, col=NA, border = NA)
  mtext(paste0("data min: ", round(min(values, na.rm=NA),2) , " max: ", round(max(values, na.rm=NA),2)), 1, -12, adj=0, cex=0.8)
  mtext(paste0("grid min: ", round(zmin,2) , " max: ", round(zmax,2)), 1, -9, adj=0, cex=0.8)

  #close pdf file
  dev.off()
  
  # end
  return(output)
}


#--------------------------------------------------
### 2) Function to Krige data using a single isotropic variogram
##### Input:  stns = list of station numbers
#####         values = list of same length as stns with values to be kriged
#####         file.name = name of output file
##### Output: pdf file with specified name with calculated variogram
#####         grid of kriged values

krige_data_1zone <- function(stns, values, file.name) {
  # create data frame
  data <- data.frame(STATION_NUMBER=stns, value=values) %>% na.omit() 

  # read in shapefile with canadian boundary to set canvas boundary
  Canada <- st_read(dsn="./Dependencies", "CanadaBound")
  map_crs <- st_crs(Canada)

  # get station coordinates
  stn_coor <- as.data.frame(hy_stations(data$STATION_NUMBER) %>% dplyr::select(STATION_NUMBER, LATITUDE, LONGITUDE))
  data <- merge(data, stn_coor, by="STATION_NUMBER") #combine data with station coordinates
  data <- st_as_sf(data, coords=c("LONGITUDE", "LATITUDE"), crs="+proj=longlat +datum=WGS84")
  data <- st_transform(data, map_crs)
  
  # create grid with 10km cells covering the whole country onto which data
  # should be interpolated
  Cangrid <- st_make_grid(Canada, cellsize=10000, crs=map_crs, what="centers")
  grids <- list()
  
  # variogram model curve types
  Mtype <- c("Sph", "Exp")
  Mcol <- c("red", "green")
  
  # create pdf file to populate with variograms
  pdf(file = paste0("./Output/", file.name, ".pdf"), width = 8.5, height = 11)
  par(mfrow = c(5,4))

  ## find average distance between stations in zone
  dist <- as.data.frame(st_distance(data))
  dist[upper.tri(dist, diag=TRUE)] <- NA 
  MDstns <- mean(apply(dist[-1,],1,FUN=min, na.rm=TRUE))
    
  ## create variogram
  var <- variogram(data$value~1, data, width=(MDstns/3), cutoff=500000) #calculates sample variogram values with 500km max range and lag distance proportional to the mean min distance between points
    
  ## two options for fitting variogram model, try automap version first, and if
  ## it creates output with nugget of 0, try gstat version.
  ### automap fitting - entirely independent including creation of new variogram,
  ### including determining which model of exponential or spherical fits best
  out <- autofitVariogram(data$value~1, as(data, "Spatial"), model=c("Exp", "Sph"), width= MDstns, cutoff=500000)# the best fitting model from exponential or spherical is chosen
  v_mod <- out$var_model
  if (v_mod$psill[1]==0) {
    ### gstat fitting - works best with start values, so calculate these firt
    ### create new variogram
    v <- variogram(data$value~1, data)
    ### note data from variogram as start values 
    nugget0 <- min(v$gamma)
    psill0 <- ((max(v$gamma)+median(v$gamma))/2)-nugget0
    range0 <- 0.1*sqrt((st_bbox(data)[1] - st_bbox(data)[3])^2+
                         (st_bbox(data)[2] - st_bbox(data)[4])^2)
    ### find best model between exponential and spherical
    vmods <- list()
    SSErr <- rep(NA, length(Mtype))
    for (k in 1:length(Mtype)) {
      vmods[[k]] <- fit.variogram(v, vgm(psill0, model=Mtype[k], range0, nugget0))
      SSErr[k] <- attr(vmods[[k]], "SSErr")
    }
    k <- which(SSErr==min(SSErr))
    ### if alternate model is better, use it
    if (vmods[[k]]$psill[1] > model$psill[1]) { 
      message("'gstat' variogram fitting has better results that 'automap' fitting")
      v_mod <- vmods[[k]]
    }
  }
    
  ## sort variogram data to plot, and change units from m to km
  dist <- var$dist/1000
  vario <- var$gamma
  ## plot variogram first
  par(mar = c(4, 4, 2, 1))
  plot(dist, vario, xlab="", ylab="semivariance", xlim=c(0,500), ylim=c(0,
          (ceiling(max(var$gamma, na.rm=TRUE)*100))/100), las=1, col="blue")
  title(main=paste0("Variogram"), line=0.5)
  title(xlab="distance [km]", line=2)
  ## sort model data to plot
  dist0 <- seq(0, max(dist), length.out=100)
  nugget <- v_mod$psill[1]
  sill <- sum(v_mod$psill[1:2])
  range <- v_mod$range[2]/1000
  k <- case_when(v_mod$model[2]=="Exp" ~ 2, #exponential
                 v_mod$model[2]=="Sph" ~ 1) #spherical
  MV <- case_when(k==2 ~ (sill-nugget)*(1-exp(-dist0/range))+nugget, #exponential
                  k==1 ~ ifelse(dist0>range,sill,(sill-nugget)*(1.5*dist0/range-0.5*(dist0/range)^3)+nugget)) #spherical
  rangeC <- case_when(k==2 ~ 3*range, #exponential
                      k==1 ~ range) #spherical
  ## plot model and write its parameters
  lines(dist0, MV, col=Mcol[k])
  mtext(paste0("model: ", v_mod$model[2], "   "), side=1, line=-4, adj=1, cex=0.4, col=Mcol[k])
  mtext(paste0("nugget: ", round(nugget, 3), "   "), side=1, line=-3, adj=1, cex=0.4, col=Mcol[k])
  mtext(paste0("sill: ", round(sill, 3), "   "), side=1, line=-2, adj=1, cex=0.4, col=Mcol[k])
  mtext(paste0("range: ", round(rangeC, 0), "   "), side=1, line=-1, adj=1, cex=0.4, col=Mcol[k])
  
    
  ## krige
  OK <- krige(formula=value~1, locations=data, newdata=Cangrid,
              model=v_mod, maxdist=500000, nmin=8)

  ## transform results into  a raster
  OKraster <- rast(data.frame(st_coordinates(OK),OK$var1.pred), crs=map_crs$wkt)
  output <- mask(OKraster, vect(Canada))

  # note limits
  zmax <- minmax(output)[2]
  zmin <- minmax(output)[1]
  par(mar = c(0, 0, 0, 0))
  plot(Canada$geometry, col=NA, border = NA)
  mtext(paste0("data min: ", round(min(values, na.rm=NA),2) , " max: ", round(max(values, na.rm=NA),2)), 1, -12, adj=0, cex=0.8)
  mtext(paste0("grid min: ", round(zmin,2) , " max: ", round(zmax,2)), 1, -9, adj=0, cex=0.8)#  substring("variogramLF_all_elim",13,nchar("variogramLF_all_elim"))
  
  #close pdf file
  dev.off()
  
  # end
  return(output)
}

#--------------------------------------------------
### 3) Function to create surface map using CESI colours and accoutrements
##### Input:  input = gridded data to map
#####         echelle = list of two numbers with min and max of colour scale
#####         divisions = number of intervals on  colour scale
#####         high.col = either "blue" or "orange", describing the order of the colour bar
#####         leg.title = character string with legend title
#####         file.name = character string of name of output file
#####         type = character string for type of output file, either "PNG" or "PDF"
#####         north_a = logical, whether to include north arrow
#####         ISD_a = logical, whether to include increasing, stable and decreasing arrows
#####         scale_bar = logical, whether to include a scale_bar
##### Output: png or pdf file with map

map_results_surface <- function(input, echelle, divisions, high.col="blue", 
                                leg.title, map.title=NA, file.name, type, 
                                north_a=TRUE, ISD_a=TRUE, scale_bar=TRUE) {
  
  # import objects to map
  Canada <- st_read(dsn="./Dependencies", "CanadaBound")
  map_crs <- st_crs(Canada)
  Lakes <- st_read(dsn="./Dependencies", "CANwaterbodies", crs=map_crs)
  
  # prepare colour pallet - CESI colours
  orange <- rgb(161,83,34, maxColorValue=255)
  white <- rgb(247,247,247, maxColorValue=255)
  blue <- rgb(15,76,106, maxColorValue=255)
  blue_lakes <- rgb(150,198,250, maxColorValue=255)
  grey_background <- rgb(164,163,173, maxColorValue=255)
  grey_borders <- rgb(130,130,130, maxColorValue=255)
  if (high.col=="orange") {
    colours <- colorRampPalette(c(blue,white,orange))
    high.col <- orange
    low.col <- blue
  } else {
    colours <- colorRampPalette(c(orange, white, blue))
    high.col <- blue
    low.col <- orange
  }
  
  # find raster min max
  zmax <- minmax(input)[2]
  zmin <- minmax(input)[1]
  
  # find legend placement based on plot limits
  leg.col.ext.x = c(ext(input)[1]+((ext(input)[2]-ext(input)[1])*0.75),
                    ext(input)[1]+((ext(input)[2]-ext(input)[1])*0.8))
  leg.col.ext.y = c(ext(input)[3]+((ext(input)[4]-ext(input)[3])*0.59),
                    ext(input)[3]+((ext(input)[4]-ext(input)[3])*0.89))
  leg.col.tit = c(ext(input)[1]+((ext(input)[2]-ext(input)[1])*0.77),
                  ext(input)[3]+((ext(input)[4]-ext(input)[3])*0.92))
  
  # open file depending on type
  if (type=="PNG") {
    png(file=paste0("./Output/", file.name, ".png"), width=5, height=5,
        units="in", res=600)
  } else {
    pdf(file=paste0("./Output/", file.name, ".pdf"), width=5, height=5)
  }
  par(mar = c(0, 0, 0, 0))

  #plot
  plot(Canada$geometry, col=grey_background, border = NA)
  colorlegend(colbar=colours(101), 
              labels=seq(echelle[1],echelle[2], ((echelle[2]-echelle[1])/divisions)), 
              xlim=leg.col.ext.x, ylim=leg.col.ext.y, 
              align="l", cex=0.6,
              lim.segment=c(0.43,0.7))
  text(leg.col.tit[1], leg.col.tit[2], leg.title, adj=0.5, cex=0.6)
  if (echelle[2] < zmax) { plot(input, add=TRUE, col=high.col, range=c(echelle[2],zmax), legend=FALSE, axes=FALSE) }
  if (echelle[1] > zmin) { plot(input, add=TRUE, col=low.col, range=c(zmin,echelle[1]), legend=FALSE, axes=FALSE) }
  plot(input, add=TRUE, col=colours(1000), range=echelle, legend=FALSE, axes=FALSE)
  plot(Canada$geometry, add=TRUE, col=NA, border=grey_borders, lwd = 0.5)
  plot(Lakes$geometry, add=TRUE, col=blue_lakes, border=grey_borders, lwd=0.1)
  legend("topright", inset = c(0.155, 0.43), cex = 0.6, bty="n",
         border= grey_borders, fill = grey_background, legend="No data")
  if (!is.na(map.title)) { title(main=map.title, line=-1) }
  
# add ISD arrows if requested
  if(ISD_a) {
    # increasing, stable, decreasing arrows
    ISD_arrows <- readPNG("./Dependencies/up_level_down_arrows.png")
    centre_ISD_arrows <- c(2300000,mean(leg.col.ext.y))
    ISD_text <- "Increasing\n\nStable\n\nDecreasing"
    ISD_text_loc <- c(2450000,mean(leg.col.ext.y))
    # calculations to find corners of ISD arrows
    size_ISD_arrows <- c(250000, dim(ISD_arrows)[1]*250000/dim(ISD_arrows)[2])
    x_left <- centre_ISD_arrows[1]-size_ISD_arrows[1]/2
    x_right <- centre_ISD_arrows[1]+size_ISD_arrows[1]/2
    y_bottom <- centre_ISD_arrows[2]-size_ISD_arrows[2]/2
    y_top <- centre_ISD_arrows[2]+size_ISD_arrows[2]/2
    # plot ISD arrows
    rasterImage(ISD_arrows, x_left, y_bottom, x_right, y_top)
    # add accompanying text
    text(ISD_text_loc[1], ISD_text_loc[2], ISD_text, adj=0, cex=0.6)
  }
  
# add north arrow if requested
  if (north_a) {
    # north arrow and position
    n_arrow <- readPNG("./Dependencies/esri6north.png")
    centre_n_arrow <- c(2056967,-417531)
    # calculations to find angle of north arrow
    dim_can_lat <- c(41.91, 83.11)
    dim_can_fig <- ext(Canada)[3:4]
    loc_n_fig <- ((dim_can_fig[2]-dim_can_fig[1])/(dim_can_lat[2]-dim_can_lat[1])*(90-dim_can_lat[2]))+dim_can_fig[2]
    angle <- atan((centre_n_arrow[1]-mean(ext(Canada)[1:2]))/(loc_n_fig-centre_n_arrow[2]))*pi
    # calculations to find corners of north arrow
    size_n_arrow <- c(150000, dim(n_arrow)[1]*150000/dim(n_arrow)[2])
    x_left <- centre_n_arrow[1]-size_n_arrow[1]/2
    x_right <- centre_n_arrow[1]+size_n_arrow[1]/2
    y_bottom <- centre_n_arrow[2]-size_n_arrow[2]/2
    y_top <- centre_n_arrow[2]+size_n_arrow[2]/2
    # plot north arrow
    rasterImage(n_arrow, x_left, y_bottom, x_right, y_top, angle=angle)
  }
  
# add scale bar if requested
  if (scale_bar) {
    # scale bar position
    loc_scale=c(2465000, -500000)
    # plot scale
    sbar(d=500000, xy=loc_scale, type="line", divs=2, lonlat=FALSE,
         label=c(0,250,500), adj=c(0.5,-0.75), ticks=TRUE, cex=0.5, lwd=0.5)
    mtext("km ", side=1, line=-4.8, adj=1, cex=0.5)
  }

  #close file
  dev.off()
}

#--------------------------------------------------
### 4) Function to create high-normal-low station map using CESI colours and accoutrements
##### Input:  2 vectors of the same length:
#####             stns = list of stations
#####             status = list of high-normal-low status for each station
#####         file.name = character string of name of output file
#####         type = character string for type of output file, either "PNG" or "PDF"
#####         north_a = logical, whether to include north arrow
#####         scale_bar = logical, whether to include a scale_bar
##### Output: png or pdf file with map
map_results_HNL <- function(stns, status, map.title=NA, file.name, type, 
                                north_a=TRUE, scale_bar=TRUE) {
  
  # import objects to map
  Canada <- st_read(dsn="./Dependencies", "CanadaBound")
  map_crs <- st_crs(Canada)
  Lakes <- st_read(dsn="./Dependencies", "CANwaterbodies", crs=map_crs)
  
  # prepare colour pallet - CESI colours
  orange <- rgb(161,83,34, maxColorValue=255)
  green <- rgb(133,161,66, maxColorValue = 255)
  white <- rgb(247,247,247, maxColorValue=255)
  blue <- rgb(15,76,106, maxColorValue=255)
  blue_lakes <- rgb(150,198,250, maxColorValue=255)
  grey_background <- rgb(164,163,173, maxColorValue=255)
  grey_borders <- rgb(130,130,130, maxColorValue=255)
  
  # get station coordinates and information
  data <- data.frame(STATION_NUMBER=stns, status=status) %>% na.omit()
  op <- hy_stn_data_coll(data$STATION_NUMBER) %>% filter(DATA_TYPE=="Flow") %>% 
          group_by(STATION_NUMBER) %>% summarize(OPERATION=last(OPERATION))
  reg <- hy_stn_regulation(data$STATION_NUMBER) %>% select(STATION_NUMBER, REGULATED)
  coor <- hy_stations(data$STATION_NUMBER) %>% select(STATION_NUMBER, LATITUDE, LONGITUDE)
  df_list <- list(data, op, reg, coor)      
  #merge all data frames together
  data <- df_list %>% reduce(full_join, by="STATION_NUMBER")
  data <- st_as_sf(data, coords=c("LONGITUDE", "LATITUDE"), crs="+proj=longlat +datum=WGS84")
  data <- st_transform(data, map_crs)

  # create symbols for each station
  # shapes
  sym.shp <- rep(NA, nrow(data))
  sym.shp <- case_when(data$REGULATED==FALSE ~ 21,
                       data$REGULATED==TRUE ~ 24)
  # colours
  sym.col <- rep(NA, nrow(data))
  sym.col <- case_when(data$status=="high" ~ blue,
                       data$status=="normal" ~ green,
                       data$status=="low" ~ orange)
  sym.bg <- sym.col
  sym.bg[data$OPERATION=="Seasonal"] <- grey_background

  # Legend Items
  local.items<- c(as.expression(bquote(bold("Yearly"))), "Low", "Normal", "High",
                  as.expression(bquote(bold("Seasonal"))), "Low", "Normal", "High")
  local.leg.col <- c(NA, orange, green, blue, NA, orange, green, blue)
  local.leg.pt.bg<- c(NA,orange, green, blue, NA,rep(grey_background, 3))
  
  # open file depending on type
  if (type=="PNG") {
    png(file=paste0("./Output/", file.name, ".png"), width=5, height=5,
        units="in", res=600)
  } else {
    pdf(file=paste0("./Output/", file.name, ".pdf"), width=5, height=5)
  }
  par(mar = c(0, 0, 0, 0))
  
  #plot
  plot(Canada$geometry, col=grey_background, border = grey_borders, lwd = 0.5)
  plot(Lakes$geometry, add=TRUE, col=blue_lakes, border=grey_borders, lwd=0.1)
  plot(data$geometry, add=TRUE, pch=sym.shp, col=sym.col, bg=sym.bg, cex=0.4)
  # legend is plotted in 4 parts to get output looking like what IID do
  legend("topright", inset=c(0.2, 0.1), cex=0.6, bty="n", 
         title="Natural", title.cex=0.7, title.adj=0, legend=local.items[1:4],
         pch=21, col=local.leg.col[1:4], pt.bg=local.leg.pt.bg[1:4])
  legend("topright", inset = c(0.05, 0.1), cex=0.6, bty="n", 
         title=NA, title.cex=0.7, legend=local.items[5:8],
         pch=21, col=local.leg.col[5:8], pt.bg=local.leg.pt.bg[5:8])
  legend("topright", inset = c(0.2, 0.24), cex=0.6, bty="n", 
         title="Regulated", title.cex=0.7, title.adj=0, legend=local.items[1:4],
         pch=24, col=local.leg.col[1:4], pt.bg=local.leg.pt.bg[1:4])
  legend("topright", inset = c(0.05, 0.24), cex=0.6, bty="n", 
         title=NA, title.cex=0.7, legend=local.items[5:8],
         pch=24, col=local.leg.col[5:8], pt.bg=local.leg.pt.bg[5:8])
  if (!is.na(map.title)) { title(main=map.title, line=-1) }
  
  
  # add north arrow if requested
  if (north_a) {
    # north arrow and position
    n_arrow <- readPNG("./Dependencies/esri6north.png")
    centre_n_arrow <- c(2056967,-417531)
    # calculations to find angle of north arrow
    dim_can_lat <- c(41.91, 83.11)
    dim_can_fig <- ext(Canada)[3:4]
    loc_n_fig <- ((dim_can_fig[2]-dim_can_fig[1])/(dim_can_lat[2]-dim_can_lat[1])*(90-dim_can_lat[2]))+dim_can_fig[2]
    angle <- atan((centre_n_arrow[1]-mean(ext(Canada)[1:2]))/(loc_n_fig-centre_n_arrow[2]))*pi
    # calculations to find corners of north arrow
    size_n_arrow <- c(150000, dim(n_arrow)[1]*150000/dim(n_arrow)[2])
    x_left <- centre_n_arrow[1]-size_n_arrow[1]/2
    x_right <- centre_n_arrow[1]+size_n_arrow[1]/2
    y_bottom <- centre_n_arrow[2]-size_n_arrow[2]/2
    y_top <- centre_n_arrow[2]+size_n_arrow[2]/2
    # plot north arrow
    rasterImage(n_arrow, x_left, y_bottom, x_right, y_top, angle=angle)
  }
  
  # add scale bar if requested
  if (scale_bar) {
    # scale bar position
    loc_scale=c(2465000, -500000)
    # plot scale
    sbar(d=500000, xy=loc_scale, type="line", divs=2, lonlat=FALSE,
         label=c(0,250,500), adj=c(0.5,-0.75), ticks=TRUE, cex=0.5, lwd=0.5)
    mtext("km ", side=1, line=-4.8, adj=1, cex=0.5)
  }
  
  #close file
  dev.off()
}

#--------------------------------------------------
### 5) Function to do inverse distance weighted interpolation
##### Input:  stns = list of station numbers
#####         values = list of same length as stns with values to be kriged
#####         power = number, power of inverse distance weighting, typically
#####                     between 1 and 5. Higher numbers have more smoothing
##### Output: grid of inverse distance weighted interpolated values

idw_data <- function(stns, values, power=5) {
  # create data frame
  data <- data.frame(STATION_NUMBER=stns, value=values) %>% na.omit() 
  
  # read in shapefile with canadian boundary to set canvas boundary
  Canada <- st_read(dsn="./Dependencies", "CanadaBound")
  map_crs <- st_crs(Canada)
  
  # get station coordinates
  stn_coor <- as.data.frame(hy_stations(data$STATION_NUMBER) %>% dplyr::select(STATION_NUMBER, LATITUDE, LONGITUDE))
  data <- merge(data, stn_coor, by="STATION_NUMBER") #combine data with station coordinates
  data <- st_as_sf(data, coords=c("LONGITUDE", "LATITUDE"), crs="+proj=longlat +datum=WGS84")
  data <- st_transform(data, map_crs)
  
  # create grid with 10km cells covering the whole country onto which data
  # should be interpolated
  Cangrid <- st_make_grid(Canada, cellsize=10000, crs=map_crs, what="centers")
  
  # interpolate
  IDW <- idw(formula=value~1, locations=data, newdata=Cangrid,
             maxdist=500000, idp=power)
  
  # transform results into a raster and mask the raster to the country outline
  IDWraster <- rast(data.frame(st_coordinates(IDW),IDW$var1.pred), crs=map_crs$wkt)
  output <- mask(IDWraster, vect(Canada))
  
  # end
  return(output)
}

#--------------------------------------------------
### 6) Function to create contour map using CESI colours and accoutrements
###     Assumes that the scale will be 0 to 100, and that scalebar should highlight
###     only values that would be identified as high or low flows (> 85 or <15
###     percentiles)
##### Input:  input = gridded data to map
#####         high.col = either "blue" or "orange", describing the order of the colour bar
#####         file.name = character string of name of output file
#####         type = character string for type of output file, either "PNG" or "PDF"
#####         legend_interp = logical, whether to include high, normal, low classifications
#####         north_a = logical, whether to include north arrow
#####         scale_bar = logical, whether to include a scale_bar
##### Output: png or pdf file with map

map_results_contour <- function(input, high.col="blue", map.title=NA, file.name, 
                                type, legend_interp=TRUE, north_a=TRUE, 
                                scale_bar=TRUE) {
  
  # import objects to map
  Canada <- st_read(dsn="./Dependencies", "CanadaBound")
  map_crs <- st_crs(Canada)
  Lakes <- st_read(dsn="./Dependencies", "CANwaterbodies", crs=map_crs)
  
  # prepare colour pallet - CESI colours
  orange <- rgb(161,83,34, maxColorValue=255)
  white <- rgb(247,247,247, maxColorValue=255)
  blue <- rgb(15,76,106, maxColorValue=255)
  blue_lakes <- rgb(150,198,250, maxColorValue=255)
  grey_background <- rgb(164,163,173, maxColorValue=255)
  grey_borders <- rgb(130,130,130, maxColorValue=255)
  if (high.col=="orange") {
    colours <- colorRampPalette(c(blue,rep(white,6),orange))
  } else {
    colours <- colorRampPalette(c(orange,rep(white,6),blue))
  }
  col.classes <- data.frame(from=c(95,90,85,15,10,5,0), to=c(100,95,90,85,15,10,5),
        colour=c(colours(20)[20], colours(20)[19], colours(20)[18], colours(20)[10],
        colours(20)[3], colours(20)[2], colours(20)[1]))
  
  # populate legend items
  leg.items <- c("95-100", "90-95", "85-90", "15-85", "10-15", "5-10", "0-5", "No data")
  leg.col <- c(col.classes$colour, grey_background)  
  
  # open file depending on type
  if (type=="PNG") {
    png(file=paste0("./Output/", file.name, ".png"), width=5, height=5,
        units="in", res=600)
  } else {
    pdf(file=paste0("./Output/", file.name, ".pdf"), width=5, height=5)
  }
  par(mar = c(0, 0, 0, 0))
  
  #plot
  plot(Canada$geometry, col=grey_background, border = NA)
  plot(input, add=TRUE, type="classes", col=col.classes)
  plot(Canada$geometry, add=TRUE, col=NA, border=grey_borders, lwd = 0.5)
  plot(Lakes$geometry, add=TRUE, col=blue_lakes, border=grey_borders, lwd=0.1)
  legend("topright", inset = c(0.16, 0.15), cex = 0.6, bty="n", title="Percentile",
         border=grey_borders, fill=leg.col, legend=leg.items)
  if (!is.na(map.title)) { title(main=map.title, line=-1) }
  
  # add legend interpretation for High-Normal-Low if requested
  if (legend_interp) {
    # curly brackets
    curly <- readPNG("./Dependencies/curly_brackets.png")
    centre_curly <- c(2300000,2880000)
    curly_text <- "High\n\nNormal\n\nLow"
    curly_text_loc <- c(2400000,2880000)
    # calculations to find corners of curly brackets
    size_curly <- c(140000, dim(curly)[1]*140000/dim(curly)[2])
    x_left <- centre_curly[1]-size_curly[1]/2
    x_right <- centre_curly[1]+size_curly[1]/2
    y_bottom <- centre_curly[2]-size_curly[2]/2
    y_top <- centre_curly[2]+size_curly[2]/2
    # plot curly brackets
    rasterImage(curly, x_left, y_bottom, x_right, y_top)
    # add accompanying text
    text(curly_text_loc[1], curly_text_loc[2], curly_text, adj=0, cex=0.6)
  }

  # add north arrow if requested
  if (north_a) {
    # north arrow and position
    n_arrow <- readPNG("./Dependencies/esri6north.png")
    centre_n_arrow <- c(2056967,-417531)
    # calculations to find angle of north arrow
    dim_can_lat <- c(41.91, 83.11)
    dim_can_fig <- ext(Canada)[3:4]
    loc_n_fig <- ((dim_can_fig[2]-dim_can_fig[1])/(dim_can_lat[2]-dim_can_lat[1])*(90-dim_can_lat[2]))+dim_can_fig[2]
    angle <- atan((centre_n_arrow[1]-mean(ext(Canada)[1:2]))/(loc_n_fig-centre_n_arrow[2]))*pi
    # calculations to find corners of north arrow
    size_n_arrow <- c(150000, dim(n_arrow)[1]*150000/dim(n_arrow)[2])
    x_left <- centre_n_arrow[1]-size_n_arrow[1]/2
    x_right <- centre_n_arrow[1]+size_n_arrow[1]/2
    y_bottom <- centre_n_arrow[2]-size_n_arrow[2]/2
    y_top <- centre_n_arrow[2]+size_n_arrow[2]/2
    # plot north arrow
    rasterImage(n_arrow, x_left, y_bottom, x_right, y_top, angle=angle)
  }
  
  # add scale bar if requested
  if (scale_bar) {
    # scale bar position
    loc_scale=c(2465000, -500000)
    # plot scale
    sbar(d=500000, xy=loc_scale, type="line", divs=2, lonlat=FALSE,
         label=c(0,250,500), adj=c(0.5,-0.75), ticks=TRUE, cex=0.5, lwd=0.5)
    mtext("km ", side=1, line=-4.8, adj=1, cex=0.5)
  }
  
  #close file
  dev.off()
}


#--------------------------------------------------
### 7) Function to create station location map for national statistics
##### Input:  stns = list of alphanumeric station names
#####         file.name = character string of name of output file
#####         type = character string for type of output file, either "PNG" or "PDF"
#####         north_a = logical, whether to include north arrow
#####         scale_bar = logical, whether to include a scale_bar
##### Output: png or pdf file with map
map_stn_loc_natsts <- function(stns, file.name, type, north_a=TRUE, scale_bar=TRUE) {
  
  # import objects to map
  Canada <- st_read(dsn="./Dependencies", "CanadaBound")
  map_crs <- st_crs(Canada)
  Lakes <- st_read(dsn="./Dependencies", "CANwaterbodies", crs=map_crs)
  
  # prepare colour pallet - CESI colours
  accent6 <- rgb(2,147,164, maxColorValue=255)
  accent12 <- rgb(242,144,0, maxColorValue=255)
  blue_lakes <- rgb(150,198,250, maxColorValue=255)
  grey_background <- rgb(164,163,173, maxColorValue=255)
  grey_borders <- rgb(130,130,130, maxColorValue=255)
  
  # get station coordinates and information
  data <- data.frame(STATION_NUMBER=stns)
  op <- hy_stn_data_coll(data$STATION_NUMBER) %>% filter(DATA_TYPE=="Flow") %>% 
    group_by(STATION_NUMBER) %>% summarize(OPERATION=last(OPERATION))
  reg <- hy_stn_regulation(data$STATION_NUMBER) %>% select(STATION_NUMBER, REGULATED)
  coor <- hy_stations(data$STATION_NUMBER) %>% select(STATION_NUMBER, LATITUDE, LONGITUDE)
  df_list <- list(data, op, reg, coor)      
  #merge all data frames together
  data <- df_list %>% reduce(full_join, by="STATION_NUMBER")
  data <- st_as_sf(data, coords=c("LONGITUDE", "LATITUDE"), crs="+proj=longlat +datum=WGS84")
  data <- st_transform(data, map_crs)

  # create symbols for each station
  # shapes
  sym.shp <- rep(NA, nrow(data))
  sym.shp <- case_when(data$REGULATED==FALSE ~ 21,
                       data$REGULATED==TRUE ~ 24)
  # colours
  sym.col <- rep(NA, nrow(data))
  sym.col <- case_when(data$REGULATED==FALSE ~ accent6,
                       data$REGULATED==TRUE ~ accent12)
  sym.bg <- sym.col
  sym.bg[data$OPERATION=="Seasonal"] <- grey_background
  
  # Legend Items
  local.items<- c("Yearly", "Seasonal")

  # open file depending on type
  if (type=="PNG") {
    png(file=paste0("./Output/", file.name, ".png"), width=5, height=5,
        units="in", res=600)
  } else {
    pdf(file=paste0("./Output/", file.name, ".pdf"), width=5, height=5)
  }
  par(mar = c(0, 0, 0, 0))
  
  #plot
  plot(Canada$geometry, col=grey_background, border = grey_borders, lwd = 0.5)
  plot(Lakes$geometry, add=TRUE, col=blue_lakes, border=grey_borders, lwd=0.1)
  plot(data$geometry, add=TRUE, pch=sym.shp, col=sym.col, bg=sym.bg, cex=0.4)
  legend("topright", inset=c(0.2, 0.15), cex=0.6, bty="n", 
         title="Natural", title.cex=0.6, title.adj=0, legend=local.items,
         pch=21, col=accent6, pt.bg=c(accent6, grey_background))
  legend("topright", inset = c(0.2, 0.24), cex=0.6, bty="n", 
         title="Regulated", title.cex=0.6, title.adj=0, legend=local.items,
         pch=24, col=accent12, pt.bg=c(accent12, grey_background))

  # add north arrow if requested
  if (north_a) {
    # north arrow and position
    n_arrow <- readPNG("./Dependencies/esri6north.png")
    centre_n_arrow <- c(2056967,-417531)
    # calculations to find angle of north arrow
    dim_can_lat <- c(41.91, 83.11)
    dim_can_fig <- ext(Canada)[3:4]
    loc_n_fig <- ((dim_can_fig[2]-dim_can_fig[1])/(dim_can_lat[2]-dim_can_lat[1])*(90-dim_can_lat[2]))+dim_can_fig[2]
    angle <- atan((centre_n_arrow[1]-mean(ext(Canada)[1:2]))/(loc_n_fig-centre_n_arrow[2]))*pi
    # calculations to find corners of north arrow
    size_n_arrow <- c(150000, dim(n_arrow)[1]*150000/dim(n_arrow)[2])
    x_left <- centre_n_arrow[1]-size_n_arrow[1]/2
    x_right <- centre_n_arrow[1]+size_n_arrow[1]/2
    y_bottom <- centre_n_arrow[2]-size_n_arrow[2]/2
    y_top <- centre_n_arrow[2]+size_n_arrow[2]/2
    # plot north arrow
    rasterImage(n_arrow, x_left, y_bottom, x_right, y_top, angle=angle)
  }
  
  # add scale bar if requested
  if (scale_bar) {
    # scale bar position
    loc_scale=c(2465000, -500000)
    # plot scale
    sbar(d=500000, xy=loc_scale, type="line", divs=2, lonlat=FALSE,
         label=c(0,250,500), adj=c(0.5,-0.75), ticks=TRUE, cex=0.5, lwd=0.5)
    mtext("km ", side=1, line=-4.8, adj=1, cex=0.5)
  }
  
  #close file
  dev.off()
}

#--------------------------------------------------
### 8) Function to create station location map for trends
##### Input:  stns = list of alphanumeric station names
#####         file.name = character string of name of output file
#####         type = character string for type of output file, either "PNG" or "PDF"
#####         north_a = logical, whether to include north arrow
#####         scale_bar = logical, whether to include a scale_bar
##### Output: png or pdf file with map
map_stn_loc_trs <- function(stns, file.name, type, north_a=TRUE, scale_bar=TRUE) {
  
  # import objects to map
  Canada <- st_read(dsn="./Dependencies", "CanadaBound")
  map_crs <- st_crs(Canada)
  Lakes <- st_read(dsn="./Dependencies", "CANwaterbodies", crs=map_crs)
  
  # prepare colour pallet - CESI colours
  accent6 <- rgb(2,147,164, maxColorValue=255)
  blue_lakes <- rgb(150,198,250, maxColorValue=255)
  grey_background <- rgb(164,163,173, maxColorValue=255)
  grey_borders <- rgb(130,130,130, maxColorValue=255)
  
  # get station coordinates and information
  data <- hy_stations(stns) %>% select(STATION_NUMBER, LATITUDE, LONGITUDE)
  data <- st_as_sf(data, coords=c("LONGITUDE", "LATITUDE"), crs="+proj=longlat +datum=WGS84")
  data <- st_transform(data, map_crs)
  
  # open file depending on type
  if (type=="PNG") {
    png(file=paste0("./Output/", file.name, ".png"), width=5, height=5,
        units="in", res=600)
  } else {
    pdf(file=paste0("./Output/", file.name, ".pdf"), width=5, height=5)
  }
  par(mar = c(0, 0, 0, 0))
  
  #plot
  plot(Canada$geometry, col=grey_background, border = grey_borders, lwd = 0.5)
  plot(Lakes$geometry, add=TRUE, col=blue_lakes, border=grey_borders, lwd=0.1)
  plot(data$geometry, add=TRUE, pch=19, col=accent6, cex=0.4)

  # add north arrow if requested
  if (north_a) {
    # north arrow and position
    n_arrow <- readPNG("./Dependencies/esri6north.png")
    centre_n_arrow <- c(2056967,-417531)
    # calculations to find angle of north arrow
    dim_can_lat <- c(41.91, 83.11)
    dim_can_fig <- ext(Canada)[3:4]
    loc_n_fig <- ((dim_can_fig[2]-dim_can_fig[1])/(dim_can_lat[2]-dim_can_lat[1])*(90-dim_can_lat[2]))+dim_can_fig[2]
    angle <- atan((centre_n_arrow[1]-mean(ext(Canada)[1:2]))/(loc_n_fig-centre_n_arrow[2]))*pi
    # calculations to find corners of north arrow
    size_n_arrow <- c(150000, dim(n_arrow)[1]*150000/dim(n_arrow)[2])
    x_left <- centre_n_arrow[1]-size_n_arrow[1]/2
    x_right <- centre_n_arrow[1]+size_n_arrow[1]/2
    y_bottom <- centre_n_arrow[2]-size_n_arrow[2]/2
    y_top <- centre_n_arrow[2]+size_n_arrow[2]/2
    # plot north arrow
    rasterImage(n_arrow, x_left, y_bottom, x_right, y_top, angle=angle)
  }
  
  # add scale bar if requested
  if (scale_bar) {
    # scale bar position
    loc_scale=c(2465000, -500000)
    # plot scale
    sbar(d=500000, xy=loc_scale, type="line", divs=2, lonlat=FALSE,
         label=c(0,250,500), adj=c(0.5,-0.75), ticks=TRUE, cex=0.5, lwd=0.5)
    mtext("km ", side=1, line=-4.8, adj=1, cex=0.5)
  }
  
  #close file
  dev.off()
}