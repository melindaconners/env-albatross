######################################
# Set Environment
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
wrap360 = function(lon) {lon360<-ifelse(lon<0,lon+360,lon);return(lon360)}


########################################################################
# 1. Import Wind Data (Downloaded as NetCDF grids from ECMWF)
########################################################################
# Stack all wind datasets.
# Each stack will be associated with a date and hourly timestamp (24 per day)
nc_dir<-'~/Dropbox/Academia/04 SUNY/01 Project Components/11-Global_Env_Data/Wind/era5_hourly/Jan_Feb2019/'
setwd(nc_dir)

# Read wind netcdf file ....................
nc1 <- read_ncdf("era5_uv_1000_jan312019.nc")  # Dates for January days, Hawaii deployments
nc2 <- read_ncdf("era5_uv_1000_Feb1_Feb10_2019.nc") # Dates for February days, Hawaii deployments
nc_t <- c(nc1,nc2) # all hawaii wind stacked together
# NOTE: February 2021: Just found blog for scripting batch downloads of ERA5 data: https://retostauffer.org/code/Download-ERA5/. Might be useful. 

nc_t
####################################################################################################
# stars object with 3 dimensions and 2 attributes
# attribute(s), summary of first 1e+05 cells:
#   u [m/s]           v [m/s]        
# Min.   :-19.484   Min.   :-21.0782  
# 1st Qu.: -1.032   1st Qu.: -3.1835  
# Median :  2.257   Median :  0.6856  
# Mean   :  2.738   Mean   :  0.5616  
# 3rd Qu.:  6.472   3rd Qu.:  3.8007  
# Max.   : 20.838   Max.   : 20.3578  
# dimension(s):
#   from   to                  offset   delta                       refsys point values    
# longitude    1 1440                  -0.125    0.25 +proj=longlat +datum=WGS8...    NA   NULL [x]
# latitude     1  721                  90.125   -0.25 +proj=longlat +datum=WGS8...    NA   NULL [y]
# time         1  312 2019-01-30 23:30:00 UTC 1 hours                      POSIXct    NA   NULL    
####################################################################################################


####################################################################################################
# 2. Read in bird locations, rediscretized to one hour, and gather wind data
####################################################################################################
m <- read.csv('~/Dropbox/Academia/04 SUNY/01 Project Components/4-MIAT_2019/MIAT_Albatross_Rproj/MIAT_albatross/02 summarized data/gps01_init_quality_check_and_interp/1hr_int/miat2019_alltracks_int1hr_3kmbuffer.csv')


# subset wind data to local Hawaii region -----------------------------------------------
bb = st_bbox(c(xmin = min(m$lon) + 0.50,
               ymin = min(m$lat) - 0.50,
               xmax = max(m$lon) + 0.50,
               ymax = max(m$lat) - 0.50), crs = st_crs(nc_t))

# chop chop -----------------------------------------------------------------------------
nc_hawaii<-st_crop(nc_t, bb) #Crop wind to Hawaii region
# dim(nc_hawaii[[1]]) #[1] 1440  721  312

# Create vector of times available in nc_hawaii -------------------------------------------
time_j <- vector(mode="numeric", length=dim(nc_hawaii[[1]])[3])
for (j in 1:dim(nc_hawaii[[1]])[3]) {
  
  if (j==1) {
    time_j[j]<-as.numeric(attr(nc_hawaii[,,,2], "dimensions")[[3]]$offset) 
  }else{
    time_j[j]<-as.numeric(attr(nc_hawaii[,,,2], "dimensions")[[3]]$offset + (j-1)*3600)
  }
  
}
rm(j)

# ---------------------------------------------------------------------------------------
# Extract Wind Along Hourly Tracks
# ---------------------------------------------------------------------------------------

for (j in 1:length(m)) {
  
  # At time j
  # Find closest wind stack timestamp
  bird_time<-as.numeric(as.POSIXct(as.character(m$datetime[j]), tz="UTC"))
  time_match_ix <- which(abs(time_j-bird_time)==min(abs(time_j-bird_time)))
  
  # Isolate that layer
  nc_hawaii_timej<-nc_hawaii[,,,time_match_ix]
  
  # Overlay that data point to that stack ------------------------------------------------
  # Extract wind from grid
  bird.locs <- st_as_sf(x = m, 
                        coords = c("lon", "lat"),
                        crs = "+proj=longlat +datum=WGS84")
  
  uwind_at_locs <- raster_extract(nc_hawaii_timej[1,,,], bird.locs)
  vwind_at_locs <- raster_extract(nc_hawaii_timej[2,,,], bird.locs)
  m$u<-uwind_at_locs
  m$v<-vwind_at_locs
}

# ---------------------------------------------------------------------------------------
# Match Times and Interpolate for High-Res GPS Data (1 min locs)
# ---------------------------------------------------------------------------------------
# Add winds to higher resolution bird data, then interpolate wind (u and v components)

# Import 60 sec GPS tracks..........
m_hi <- read.csv('~/Dropbox/Academia/04 SUNY/01 Project Components/4-MIAT_2019/MIAT_Albatross_Rproj/MIAT_albatross/02 summarized data/gps01_init_quality_check_and_interp/1min_int/tracks_int_buff3km_ALLCOMBINED/all_miat2019.csv')
# dataframe has: datetime, species, ID, tripID, lat, lon

# Calculate Speed, Bearing, Distance between successive locations..........
# (Calculate bearing between two locations https://www.rdocumentation.org/packages/swfscMisc/versions/1.2/topics/bearing)
m_hi$ground_speed_kmHr <- NA
m_hi$distanceij_km <- NA
m_hi$bearingij <- NA

trips<- unique(m_hi$tripID)
for (i in 1:length(trips)) {
  tripi<-m_hi[m_hi$tripID==trips[i],]
  
  for (j in 1:length(tripi$X)-1) {
    hour_int<- 1/60 # 1 minute is 1 60th an hour
    tripi$distanceij_km[j]<-swfscMisc::distance(tripi$lat[j],tripi$lon[j],tripi$lat[j+1],tripi$lon[j+1], units = "km", method = "haversine")
    tripi$ground_speed_kmHr[j]<-tripi$distanceij_km[j]/hour_int
    tripi$bearingij[j]<-as.numeric(bearing(tripi$lat[j],tripi$lon[j],tripi$lat[j+1],tripi$lon[j+1])[1])
  }
  
  m[m$tripID==trips[i],]<-tripi
  
}


# placeholders for wind components
m_hi$u <- NA
m_hi$v <- NA

# Loop through m (lo-res) and match closest time in m_hi (high res)
# Add wind
hi_res_timematch<-as.numeric(as.POSIXct(as.character(m_hi$datetime)))

for (j in 1:length(m$species)) {
  j_match<-as.numeric(as.POSIXct(as.character(m$datetime[j])))
  match_ix<-Closest( hi_res_timematch, j_match, which=TRUE) # Which in m_hi is nearest j_match (m$datetime[j])
  m_hi$u[match_ix] <- m$u[j]
  m_hi$v[match_ix] <- m$v[j]
}

# -----------------------------------------------------------------------------
# Interpolate u and v vectors to hi resolution.
# -----------------------------------------------------------------------------
# NOTE: Here, used linear interpolation for prelim anal, but will explore other options: such as Guassian processes (recommended by Petar Djuric) >> ask Levi and Andy
m_hi$uint<-na_interpolation(m_hi$u, option = "linear")
m_hi$vint<-na_interpolation(m_hi$v, option = "linear")


# -----------------------------------------------------------------------------
# Get Wind Direction and Wind Velocity from U and V Components:
# -----------------------------------------------------------------------------
# Requires mathematical to meteorological adjustment:
# Important: our wind vectors were given in mathematical notation, not in meteorological convention. 
# Need to adjust by 270* - this is what uv2ddff does.
# http://colaweb.gmu.edu/dev/clim301/lectures/wind/wind-uv
# from uv2ddff: Wind direction is either in meteorological degrees (0 from North, from 90 East,
# 180 from South, and 270 from West) or in mathematical radiant if input \code{rad = TRUE}.

### I double-checked this formula and it was validated by another wind direction on the web.
ddff <- uv2ddff(m_hi$uint,m_hi$vint)
m_hi$wind_vel <- ddff$ff
m_hi$wind_dir360<-ddff$dd


##########################################################
#  Bird-Wind Angle
##########################################################
# note!!! Bearing is on same compass as wind direction (0-360) but it's going IN the direction (while wind is coming FROM that direction). In other words, wind with a direction of 2 degrees is coming from the north. Whereas a bird with a bearing of 2 degrees is going TO the north. 

m_hi$BWA<-abs(Lon360to180(m_hi$bearingij-m_hi$wind_dir360)) # bird bearing 0-180, wind 0-360
# m_hi$BWA360<-abs(m_hi$bearingij-m_hi$wind_dir360) # Use the above calculation

# Definitions from Spear and Ainley 1997
# Headwind  = abs(bird-wind): 0-59 
# Crosswind = abs(bird-wind): 60-119 
# Tailwind =  abs(bird-wind): 120-180 

# After discussion with LT, I used slightly larger windows for sidewinds than Spear and Ainley
headx <- which(m$BWA < 50)
crossx <- which(m$BWA >= 50 & m$BWA <= 130)
tailx <- which(m$BWA > 130)

m$bwaClass[headx] <- "Head"
m$bwaClass[crossx] <- "Cross"
m$bwaClass[tailx] <- "Tail"
m$bwaClass <- as.factor(m$bwaClass)



