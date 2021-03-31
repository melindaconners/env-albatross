############################################################################
# Create Sci Pol figure:
# Albatross tracks and SST contours El Nino vs La Nina vs 2100

# SST Months: Dec 2009, Dec 2010
# Add Forecasted SST for 2100
# Location: start at 16N, 58N and -179 to -120
# Product: 
  # - starting with monthly composites from JPL NASA Multi-scale Ultra-high Resolution (MUR) SST Analysis fv04.1, Global, 0.01°, 2002-present, Monthly: https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.html
  # - 2100 sst calculated from Bruno et al trend_yearmean_ensemble_tos_RCP85.n

#############################################################################
# Set up Environment
#############################################################################
rm(list = ls())

library(tidyverse)
library(plotdap)
library(sf)
library(velox)
library(lubridate)
library(sp)
library(viridis) # Perceptually correct palettes
library(raster) # Dealing with rasters
library(rasterVis) # Plotting rasters 
library(ncdf4) # Required - viewing and working with climate data in netCDF format
library(rgdal) # Working with shapefiles
library(rgeos) # Geometry engine
library(gridExtra) # Grid arranged graphics
library(grid) # Grid arrangement
library(rerddap)
library(rerddapXtracto)
library(oceanmap)
library(fields)
data('cmap')


dir<- '~/Dropbox/Academia/03a Science Policy/Environmental/'
years<-2009
# years<-2010

#############################################################################
# Extract, save, plot SST grids for 2009 Incubation for central North Pacific
#############################################################################
# requires downloading two grids (on either side of dateline) and then merging E and W

# Western Chunk:
xpos <- c(150, 180)  # larger window for now
ypos <- c(20, 60)
tpos <- c(paste0(years,"-12-16"), paste0(years,"-12-16"))
SSTInfo <- rerddap::info('jplMURSST41mday')
xsst_w <- rxtracto_3D(SSTInfo, parameter = 'sst', xcoord = xpos, ycoord = ypos, tcoord = tpos)
save(xsst_w,file=paste0(dir,"SST_December_Final/SST_jplMURSST41mday_Dec162009_WEST.Rdata"))
plotBBox(xsst_w, plotColor = 'temperature')

# Eastern Chunk:
xpos <- c(-179.99, -120)  # larger window for now
ypos <- c(20, 60)
tpos <- c(paste0(years,"-12-16"), paste0(years,"-12-16"))
SSTInfo <- rerddap::info('jplMURSST41mday')
xsst_e <- rxtracto_3D(SSTInfo, parameter = 'sst', xcoord = xpos, ycoord = ypos, tcoord = tpos)
save(xsst_e,file=paste0(dir,"SST_December_Final/SST_jplMURSST41mday_Dec162009_EAST.Rdata"))
plotBBox(xsst_e, plotColor = 'temperature')


#############################################################################
# Prep SST data for merging and solve rounding error of MUR dataset
# Save raster of SST for full central North Pacific 
#############################################################################

# Wrap longitudes to 360 for merging in new extent
xsst_e$longitude <- Lon180to360(xsst_e$longitude)
xsst_w$longitude <- Lon180to360(xsst_w$longitude)

      ################################################################
      # MUR SST has rounding error, so R not recognizing it as a regular grid
      # As a result, it won't create a raster from the grid. 
      # A work around (below) is to redefine your grid and then repopulate it
      # with your extracted SST values.
      
      # Before doing this, just see if xyzToRaster() works on your grid! 
      # If it does, then jump straight to rasterToContour()
      
      # Troubleshooting advice here from discussions with Michael Sumner and Roy Mendellsohn
      
      ################################################################
      #reset extent to solve rounding error of MUR dataset and save as raster
      x_spacing <- 0.01/2
      y_spacing <- 0.01/2
      ex_east <- raster::extent(min(xsst_e$longitude) - x_spacing,
                                max(xsst_e$longitude) + x_spacing, 
                                min(xsst_e$latitude) - y_spacing, 
                                max(xsst_e$latitude) + y_spacing)
      
      r_east <- setExtent(raster::brick(xsst_e$sst[, ncol(xsst_e$sst):1, ,  drop = FALSE],
                                        transpose = TRUE,
                                        crs = "+proj=longlat +datum=WGS84"), ex_east)
      
      ex_west <- raster::extent(min(xsst_w$longitude) - x_spacing,
                                max(xsst_w$longitude) + x_spacing, 
                                min(xsst_w$latitude) - y_spacing, 
                                max(xsst_w$latitude) + y_spacing)
      
      r_west <- setExtent(raster::brick(xsst_w$sst[, ncol(xsst_w$sst):1, ,  drop = FALSE],
                                        transpose = TRUE,
                                        crs = "+proj=longlat +datum=WGS84"), ex_west)

# Merge east and west
rPacific<-merge(r_west, r_east)
plot(rPacific)
title("2009")
save(rPacific,file=paste0(dir,"SST_December_Final/SST_jplMURSST41mday_Dec16_2009_FULLPACIFIC.Rdata"))


#############################################################################
# Create and save SST Isotherms of interest 
#############################################################################

cr2009_18<- rasterToContour(rPacific, levels=18) #TZCF proxy 2009
cr2009_12<- rasterToContour(rPacific, levels=12) #LAAL proxy 2009

cr2010_18<- rasterToContour(rPacific, levels=18) #TZCF proxy 2010
cr2010_12<- rasterToContour(rPacific, levels=12) #LAAL proxy 2010

# plot(cr2010_18)
# plot(cr2010_12, add=TRUE)
# title("2010 TZCF")


##########################################################
# Calculate 2100 TZCF Proxy from FORECASTED SST  (RCP85)
#########################################################

# Import RCP85 warming trend from Bruno et al 
sst85_dir <- '~/Dropbox/Academia/03a Science Policy/Environmental/2100_sst/BrunoEtAl/trend_yearmean_ensemble_tos_RCP85.nc'
sst85 <- raster(sst85_dir) # Read in raster

# crop to Hawaii Raster
extent(sst85) <- extent(rPacific) # Change extent
image(sst85)

# Check CRS
sst85@crs
rPacific@crs

# Change resolution of SST2010 Raster .... 
# rPacific
# class       : RasterLayer 
# dimensions  : 4001, 9001, 36013001  (nrow, ncol, ncell)
# resolution  : 0.01, 0.01  (x, y)
# extent      : 150, 240, 20, 60.01  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
# data source : in memory
# names       : layer 
# values      : -1.489, 28  (min, max)

# needs to be
# sst85
# class       : RasterLayer 
# dimensions  : 180, 360, 64800  (nrow, ncol, ncell)
# resolution  : 0.25, 0.2223  (x, y)
# extent      : 150, 240, 20, 60.01  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
# data source : ~/Dropbox/Academia/03a Science Policy/Environmental/2100_sst/BrunoEtAl/trend_yearmean_ensemble_tos_RCP85.nc 
# names       : tos 
# z-value     : 1121.46826239546 
# zvar        : tos

# Create new extent and resolution and resample SST grids into new grid
e <- extent(sst85)
ra = raster(e) 
res(ra)=c(0.25, 0.2223)
crs(ra) <- crs(sst85)
rPac_2009_re <-resample(rPacific,ra,method='bilinear')
rPac_2010_re <-resample(rPacific,ra,method='bilinear')

# Create new contours for full Pacific
cr2009_18<- rasterToContour(rPac_2009_re, levels=18) #TZCF proxy 2009
cr2009_12<- rasterToContour(rPac_2009_re, levels=12) #LAAL proxy 2009

cr2010_18<- rasterToContour(rPac_2010_re, levels=18) #TZCF proxy 2010
cr2010_12<- rasterToContour(rPac_2010_re, levels=12) #LAAL proxy 2010

# plot(rPac_2009_re)
# plot(cr2009_18, add=TRUE)
# plot(cr2009_12, add=TRUE)
# title("2009 TZCF")
# 
# plot(rPac_2010_re)
# plot(cr2010_18, add=TRUE)
# plot(cr2010_12, add=TRUE)
# title("2010 TZCF")


###########################################################
# Raster Calculation for SST 2100 based off of warming trend
############################################################
nyears <- 2100-2010
sstDelta <- sst85*nyears #(#degreeChangePerYear * NumberOfYears)
sstStack <- stack(rPac_2010_re,sstDelta) # add sstDelta to 2010sst (s)
sst2100<-sum(sstStack)

###########################################################
# 2100 SST Contours: 
########################################################
ctr_18_2100 <- rasterToContour(sst2100, levels=18) #TZCF proxy 2100
ctr_12_2100 <- rasterToContour(sst2100, levels=12) #LAAL pref proxy 2100

###########################################################
# Save all SST Contours: 
########################################################
save(rPac_2009_re, rPac_2010_re, sst2100, cr2010_18, cr2010_12, cr2009_18, cr2009_12,  ctr_18_2100, ctr_12_2100, file='/Dropbox/Academia/03a Science Policy/Environmental/all_stt_contours/TZCF_SST_Contour_Files_2019-07-28.Rdata')

##########################################################
# Quick Plots of all contours
##########################################################
data("wrld_simpl")

#TZCF 18C
image(rPac_2010_re,col=alpha(cmap$sst, 0.65), main=NA)
plot(cr2009_18, add=TRUE)
plot(cr2010_18, add=TRUE)
plot(ctr_18_2100, add=TRUE)
title("December SST 2010 with 18C Proxy for TZCF: 2009, 2010, 2100")
plot(wrld_simpl, add=TRUE, col="black")
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE, col="black")

#LAAL 12C
image(rPac_2010_re,col=alpha(cmap$sst, 0.65), main=NA)
plot(cr2009_12, add=TRUE)
plot(cr2010_12, add=TRUE)
plot(ctr_12_2100, add=TRUE)
title("December SST 2010 with 12C Proxy for LAAL preference: 2009, 2010, 2100")
plot(wrld_simpl, add=TRUE, col="black")
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE, col="black")


image(rPac_2010_re,col=alpha(cmap$sst, 0.65), main=NA)
plot(cr2009_12, add=TRUE)
plot(cr2009_18, add=TRUE)

plot(ctr_12_2100, add=TRUE, lty='dotted')
plot(ctr_18_2100, add=TRUE, lty='dotted')
plot(wrld_simpl, add=TRUE, col="black")
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE, col="black")
title("TZCF 2009 vs 2100")



#################################################################################
# Quick Plots of TZCF BAND Dec 2009-10 and Dec 1010-11 w/ incubating LAAL tracks.
#################################################################################

# load contours ... 
load('~/Dropbox/Academia/03a Science Policy/Environmental/SST_December_Final/TZCF_SST_Contour_Files_2019-07-23.Rdata')

# Import bird data ....
bird_m <- read.csv('~/Dropbox/Academia/03a Science Policy/Albatross/location_data/TEIS_All_GPS_PTT_Combined/teis_all_albatross_gps_ptt_0506_1112.csv')
bird_m$lon <- Lon180to360(bird_m$lon)

# Isolate incubating LAAL in each year.....
bird_m$ptime<-as.POSIXct(bird_m$datetime)
bird_m$month<-month(as.Date(bird_m$ptime))
bird_m$year <- year(as.Date(bird_m$ptime))
spp_num <- substr(bird_m$animID,1,2)
bird_m$species <- ifelse(spp_num=="28", "BFAL", "LAAL")
# laal Dec 2009
m1<-bird_m %>% filter(species=="LAAL", breeding_phase=="incubation", yearspan=="20092010")
# laal Dec 2010
m2<-bird_m %>% filter(species=="LAAL", breeding_phase=="incubation", yearspan=="20102011")

image(rPac_2009_re,col=alpha(cmap$sst, 0.65), main=NA)
plot(cr2009_12, add=TRUE, lty='longdash')
plot(cr2009_18, add=TRUE, lty='longdash')
lines(m1$lon, m1$lat,lwd=1)
plot(wrld_simpl, add=TRUE, col="black")
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE, col="black")
title("TZCF 2009 - El Nino")

image(rPac_2010_re,col=alpha(cmap$sst, 0.65), main=NA)
plot(cr2010_12, add=TRUE, lty='longdash')
plot(cr2010_18, add=TRUE, lty='longdash')
lines(keepbirds$lon, keepbirds$lat,lwd=1)
plot(wrld_simpl, add=TRUE, col="black")
plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE, col="black")
title("TZCF 2010 - La Nina")




