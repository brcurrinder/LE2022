# Script for converting met file

# Save results here
CleanedFile = "bc/NLDAS2_Mendota_1979_2016_cell_5_GLMReady_cut_timezonecorr_CLEANED.csv"

windfactor = WindFactor
#MeteoFiles = rep('bc/LakeEnsemblR_meteo_standard.csv',3)
MeteoFiles = rep('bc/NLDAS2_Mendota_1979_2016_cell_5_GLMReady_cut_timezonecorr.csv',3)

# Load the bad file
meteo = read_csv(MeteoFiles[1])
#head(meteo)
daily_meteo = meteo
# Subset usable range
daily_meteo = daily_meteo[14:3000,]
# Convert date format
daily_meteo$datetime = as.POSIXct(as.character(daily_meteo$datetime),format='%m/%d/%y %H:%M')
# Save it
write_csv(daily_meteo,CleanedFile)

# Reload the saved file and test it in the physics
test = read_csv(CleanedFile)
head(test)
#head(daily_meteo)
daily_meteo = test
daily_meteo$date = daily_meteo$datetime


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

