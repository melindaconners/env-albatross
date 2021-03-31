
######################################
# Script to Plot Hawaii Winds
######################################
require(sf)
require(tidyverse)
require(raster)
library(nngeo)
library(lubridate)
library(rgdal)
library(stars)
library(ncmeta)
library(ncdf4)
library(lattice)
library(rCAT) #rad2deg()
library(rasterVis) #vectorplot
library(colorRamps) #matlab.like
library(viridisLite)
library(colorspace)
library(DescTools) #closest
library(imputeTS) #na.interpolation
library(sp)
library(maptools)
wrap360 = function(lon) {lon360<-ifelse(lon<0,lon+360,lon);return(lon360)}


########################################################################
# 1. Import Wind Data (Downloaded as NetCDF grids from ECMWF ERA5 hourly)
########################################################################
# Stack all wind datasets.
# Each stack will be associated with a date and time (24 per day)

nc_dir<-'~/Dropbox/Academia/04 SUNY/01 Project Components/11-Global_Env_Data/Wind/era5_hourly/Jan_Feb2019/'
setwd(nc_dir)

# Read wind netcdf file ....................
nc1 <- read_ncdf("era5_uv_1000_jan312019.nc") 
nc2 <- read_ncdf("era5_uv_1000_Feb1_Feb10_2019.nc") 
nc_t <- c(nc1,nc2) # all wind stacked together

# Set extent and crop ....................
xmin360=wrap360(160)
xmax360=wrap360(-172)
bb <- st_bbox(c(xmin=xmin360, xmax=xmax360, ymin=20, ymax=47), crs=st_crs(nc1))
bb

# Crop stars object to desired extent .........................................
nc_hawaii<-nc_t[bb]
names_hc <- st_get_dimension_values(nc_hawaii, "time")
nc_hawaii_u<-as(nc_hawaii[1,,,], "Raster")
nc_hawaii_v<-as(nc_hawaii[2,,,], "Raster")
names(nc_hawaii_u) <- names_hc
names(nc_hawaii_v) <- names_hc

#get the date from the names of the layers and extract the month ...................
date_indices <- format(as.Date(names(nc_hawaii_u), format = "X%Y.%m.%d"), format = "%m")
named_months <- format(as.Date(names(nc_hawaii_u), format = "X%Y.%m.%d"), format = "%b")
indices <- as.numeric(date_indices)

#sum layers ..........................
MonthWindU<- stackApply(nc_hawaii_u, indices, fun = mean)
names(MonthWindU) <- unique(named_months)

MonthWindV<- stackApply(nc_hawaii_v, indices, fun = mean)
names(MonthWindV) <- unique(named_months)

# Take Feb only for plotting .... 
feb_v <- MonthWindV$Feb
feb_u <- MonthWindU$Feb
feb_uv_brick <- brick(feb_u,feb_v)
wind_slope <- sqrt(feb_u^2 + feb_v^2)

vectorplot(feb_uv_brick*1.5,isField="dXY",region=wind_slope,margin=FALSE,
           par.settings=rasterTheme(region = matlab.like(10)),narrows=2000,lwd.arrows=.8,length=.02)

# Plot Bird tracks and NWHI coastline on top of wind field.......................... 
mat<-read.csv('~/Dropbox/Academia/04 SUNY/01 Project Components/4-MIAT_2019/MIAT_Albatross_Rproj/MIAT_albatross/02 summarized data/gps01_init_quality_check_and_interp/1hr_int/miat2019_alltracks_int1hr_3kmbuffer.csv')
mat$species <- substr(mat$species, 1, 4)
la<-mat[mat$species=="LAAL",]
bf<-mat[mat$species=="BFAL",]
# Convert tracks to Spatial Lines (function at end of this script)
tracks_sp<-points_to_line(data=mat, long="lon", lat="lat", id_field="animID")  # all  
tracks_la<-points_to_line(data=la, long="lon", lat="lat")  # LAAL
tracks_bf<-points_to_line(data=bf, long="lon", lat="lat")  # BFAL

nwhi_sf <- st_read("~/Dropbox/Academia/10 R Programming Universal/earth_data/Coastline__Northwest_Hawaiian_Islands_NWHI/Coastline__Northwest_Hawaiian_Islands_NWHI.shp") 
# http://geoportal.hawaii.gov/datasets/coastline-northwest-hawaiian-islands-nwhi
nwhi_sp <- as(nwhi_sf, "Spatial")
nwhi_sp360 <- sp::spTransform(nwhi_sp,"+proj=longlat +datum=WGS84 +lon_wrap=360")

tiff("~/Desktop/Wind_ERA5_Hawaii.tiff", units="in", width=6, height=5, res=300)
vectorplot(feb_uv_brick*1, isField = "dXY", region = wind_slope, margin = FALSE, 
           par.settings = rasterTheme(region =  rev(brewer.pal(5,"Blues"))), 
           narrows = 700, lwd.arrows=.7, length=.04) + 
  layer(sp.lines(tracks_la,lwd=2, col="white")) + 
  layer(sp.lines(tracks_bf,lwd=2, col="black")) + 
  layer(sp.polygons(nwhi_sp360, col="grey")) 
dev.off()


