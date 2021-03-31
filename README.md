# env-albatross
Contains scripts for importing and visualizing environmental data in relation to albatross movements

#### plot_TZCF-contours_albatross.R
This script was used to create figure "TZCF_Albatross_SciencyPolicy.pdf"
- Extracts monthly MUR SST data and creates SST contours/ TZCF proxies for El Nino and La Nina Years
- Uses RCP85 ensemble data to calculate forcasted TZCF proxy for 2100
- Overlays albatross tracks with SST contours

#### plot_WindField_tracks.R
This script uses ERA5 wind data to plot wind fields overlaid with albatross tracks

#### extract-ERA5wind_calc_bird-wind-angle.R
This script stacks rasters of hourly wind data used to extract wind along GPS tracks of albatross.
It also calculates the angle of the bird relative to wind direction (head, side, tail winds)
