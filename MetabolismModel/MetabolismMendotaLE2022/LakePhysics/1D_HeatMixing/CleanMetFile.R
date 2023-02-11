# Script for converting met file
library(tidyverse)
# Switches
PrecipMultiplier = 1000
BP = 98220
# FileToConvert = rep('bc/LakeEnsemblR_meteo_standard.csv',3)
MeteoFiles = rep('bc/NLDAS2_Mendota_1979_2016_cell_5_GLMReady_cut_timezonecorr.csv',3)
# Save results here
CleanedFile = "bc/NLDAS2_Mendota_1979_2016_cell_5_GLMReady_cut_timezonecorr_CLEANED.csv"
# CleanedFile = "bc/NLDAS2_Mendota_2010_2016_cell_5_GLMReady_cut_timezonecorr_CLEANED.csv"
# Subet by dates if true
SubsetByDates = TRUE
StartDate = "1995-01-01 UTC" #Earliest start date
EndDate   = "2015-12-31 UTC" #Latest end date
StartDate = "2000-01-01 UTC" 
EndDate   = "2007-12-31 UTC" 
RunPhysTest = FALSE

# Parameters applied after saving
windfactor = 0.8

# Load the bad file
meteo = read_csv(MeteoFiles[1])
#head(meteo)
daily_meteo = meteo
# Subset usable range
# daily_meteo = daily_meteo[14:3000,]
daily_meteo = na.omit(daily_meteo)
# Convert date format
daily_meteo$datetime = as.POSIXct(as.character(daily_meteo$datetime),format='%m/%d/%y %H:%M')
daily_meteo$Precipitation_millimeterPerDay = daily_meteo$Precipitation_millimeterPerDay * PrecipMultiplier
# Subset to date range
if (SubsetByDates){
  d = daily_meteo$datetime
  d2 = format(d,format='%Y-%m-%d')
  iSimDates = which(d2>=as.POSIXct(StartDate,format = '%Y-%m-%d') & d2<=as.POSIXct(EndDate,format = '%Y-%m-%d'))
  daily_meteo = daily_meteo[iSimDates,]
}
# Save date as char
daily_meteo$datetime = as.character(daily_meteo$datetime)
# if not NA then update BP
if (!is.na(BP)){
  daily_meteo$Surface_Level_Barometric_Pressure_pascal = BP
}
# Save it
write_csv(daily_meteo,CleanedFile)

if (RunPhysTest){
  # Reload the saved file and test it in the physics
  test = read_csv(CleanedFile)
  head(test)
  #head(daily_meteo)
  daily_meteo = test
  daily_meteo$date = daily_meteo$datetime
  daily_meteo = na.omit(daily_meteo)
  
  daily_meteo$Cloud_Cover <- gotmtools::calc_cc(date = as.POSIXct(daily_meteo$date),
                                                airt = daily_meteo$Air_Temperature_celsius,
                                                relh = daily_meteo$Relative_Humidity_percent,
                                                swr = daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared,
                                                lat = 43, lon = -89.41,
                                                elev = 258)
  daily_meteo$dt <- as.POSIXct(daily_meteo$date) - (as.POSIXct(daily_meteo$date)[1]) + 1
  daily_meteo$ea <- (daily_meteo$Relative_Humidity_percent * (4.596 * exp((17.27*(daily_meteo$Air_Temperature_celsius))/
                                                                            (237.3 + (daily_meteo$Air_Temperature_celsius) )))/100)
  daily_meteo$ea <- (101.325 * exp(13.3185 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15))) -
                                     1.976 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**2 -
                                     0.6445 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**3 -
                                     0.1229 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**4)) *daily_meteo$Relative_Humidity_percent/100
  daily_meteo$ea <- (daily_meteo$Relative_Humidity_percent/100) * 10^(9.28603523 - 2322.37885/(daily_meteo$Air_Temperature_celsius + 273.15))
  startDate <- daily_meteo$datetime[1]
  
  ## calibration parameters
  daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared <-
    daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared 
  daily_meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond <-
    daily_meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond * windfactor# wind speed multiplier
  
}

