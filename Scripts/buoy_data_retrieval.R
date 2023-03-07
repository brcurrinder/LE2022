## Script to download data from the limnological buoy on Lake Mendota
## and clean it for use in GLEON Lake Expedition 2022 analyses

# Cleaning is based on metadata and communications with LTER data manager. See
# files "Data/Buoy/Metadata.docx", "Data/Buoy/Buoy details from Mark
# Gahler.docx", "Data/Buoy/me_cp_compare_2019.csv", "mendota_buoy_sensors_v3.csv".

# This script outputs daily and hourly data files.

library(tidyverse)
library(lubridate)


# Explore differences between old and new sensors -------------------------------------

# Load data
mendota_buoy_sensors_v3 <- read_csv("Data/Buoy/mendota_buoy_sensors_v3.csv")
me_cp = read_csv("Data/Buoy/me_cp_compare_2019.csv")

# View comparison between Turner and YSI sensors when they were both deployed (2019)
me_cp %>% 
  filter(!is.na(chlor_turner),
         !is.na(chlor_ysi),
         chlor_ysi < 20) %>% 
ggplot()+
  geom_point(aes(x = chlor_turner/1000, y = chlor_ysi))

# View distribution of observations from the Turner sensor
hist(me_cp$chlor_turner/1000)

#View timeseries of both sensors in 2019
me_cp%>% 
  filter(!is.na(chlor_turner),
         !is.na(chlor_ysi),
         chlor_ysi < 20) %>% 
  mutate(sampledate = as_date(sampledate, format = "%m/%d/%Y"),
         dt = make_datetime(year = 2019, month = month(sampledate), day = day(sampledate), 
                            hour = hour(sampletime), min = minute(sampletime), 
                            sec = second(sampletime))) %>% 
  # uncomment next line to zoom in on one month
  # filter(month(dt) == 4) %>% 
  ggplot(aes(x = dt))+
  geom_point(aes(y = chlor_turner/1000), color = "blue")+
  geom_point(aes(y = chlor_ysi))+
  facet_wrap(~month(sampledate), scales= "free")

# Make clean daily data ---------------------------------------------------

# Load daily data
daily_raw = read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-ntl.129.34&entityid=cba9ed12834b8f315d6b10675bb60c5a")

# make sure data values are NA when flagged, preserving as many points as possible
daily_long = daily_raw %>% 
  pivot_longer(cols = starts_with("avg"), names_prefix = "avg_", names_to = "variable", values_to = "avg") %>% 
  pivot_longer(cols = starts_with("flag"), names_prefix = "flag_avg_", names_to = "var_flag", values_to = "flag") %>% 
  filter(variable == var_flag) %>% 
  mutate(avg = case_when(flag != "" ~ as.numeric(NA),
                         T ~ avg)) %>%  
  select(-var_flag, -flag) %>% 
  filter(!is.na(avg))

#visualize data to check
ggplot(daily_long)+
  geom_point(aes(x = sampledate, y = avg))+
  facet_wrap(~variable, scales = "free")

#specifically look at chlor/phyco
daily_long %>%  filter(variable %in% c("chlor_rfu", "phyco_rfu")) %>% 
  ggplot()+
  geom_point(aes(x = sampledate, y = avg, color = variable))+
  facet_wrap(~year(sampledate), scales = "free")
#ok - none of these look like outrageous outliers

#recreate wide table
daily_buoy = daily_long %>% 
  pivot_wider(names_from = variable, values_from = avg, names_prefix = "avg_") 

#scaling chlor and phyco RFU values
daily_buoy = daily_buoy %>%
  left_join(daily_buoy %>% group_by(year4) %>%
# new transformation: standard normal transform
              summarize(mean_chlor = mean(avg_chlor_rfu, na.rm = T),
                        mean_phyco = mean(avg_phyco_rfu, na.rm = T),
                        sd_chlor = sd(avg_chlor_rfu, na.rm = T),
                        sd_phyco = sd(avg_phyco_rfu, na.rm = T))) %>% 
  mutate(scaled_chlor_rfu = (avg_chlor_rfu - mean_chlor)/sd_chlor,
         scaled_phyco_rfu = (avg_phyco_rfu - mean_phyco)/sd_phyco) %>% 
  select(-c(mean_chlor:sd_phyco)) %>% 
# old transformation  

  #             summarise(max_chlor = max(avg_chlor_rfu, na.rm = T),
  #                       min_chlor = min(avg_chlor_rfu, na.rm = T),
  #                       max_phyco = max(avg_phyco_rfu, na.rm = T),
  #                       min_phyco = min(avg_phyco_rfu, na.rm = T))) %>% 
  # mutate(scaled_chlor_rfu = (avg_chlor_rfu- min_chlor)/max_chlor,
  #        scaled_phyco_rfu = (avg_phyco_rfu - min_phyco)/max_phyco) %>% 
  # select(-c(max_chlor:min_phyco)) %>% 
  # create trimmed / sensor adjusted values
  mutate(avg_chlor_rfu = case_when(year(sampledate) > 2018 ~ avg_chlor_rfu*1000,
                                   year(sampledate) < 2008 ~ as.numeric(NA),
                                   T ~ avg_chlor_rfu),
         avg_phyco_rfu = case_when(year(sampledate) > 2018 ~ avg_phyco_rfu*1000,
                                   year(sampledate) < 2008 ~ as.numeric(NA),
                                   T ~ avg_phyco_rfu))

#visualize data to check
daily_buoy %>% 
  pivot_longer(avg_chlor_rfu:scaled_phyco_rfu) %>% 
  ggplot()+
  geom_point(aes(x = sampledate, y = value))+
  facet_wrap(~name, scales = "free")

## NOTES: some extremely high variables of phyco (and to less degree chl) in 2020

write_csv(daily_buoy, "Data/Buoy/daily_buoy.csv")

rm(daily_buoy, daily_long, daily_raw)
# Make clean hourly data --------------------------------------------------

# Load hourly data
hourly_raw = read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-ntl.129.34&entityid=72494d432fe1e977f5326100a733cece")


# make sure data values are NA when flagged, preserving as many points as possible
hourly_long = hourly_raw %>% 
  pivot_longer(cols = starts_with("avg"), names_prefix = "avg_", names_to = "variable", values_to = "avg") %>% 
  pivot_longer(cols = starts_with("flag"), names_prefix = "flag_avg_", names_to = "var_flag", values_to = "flag") %>% 
  filter(variable == var_flag) %>% 
  mutate(avg = case_when(flag != "" ~ as.numeric(NA),
                         T ~ avg)) %>%  
  select(-var_flag, -flag) %>% 
  filter(!is.na(avg))

#visualize data to check
ggplot(hourly_long)+
  geom_point(aes(x = sampledate, y = avg))+
  facet_wrap(~variable, scales = "free")

#recreate wide table
hourly_buoy = hourly_long %>% 
  pivot_wider(names_from = variable, values_from = avg, names_prefix = "avg_") 

#scaling chlor and phyco RFU values
hourly_buoy = hourly_buoy %>%
  left_join(hourly_buoy %>% group_by(year4) %>% 
# new transformation: standard normal transform
              summarize(mean_chlor = mean(avg_chlor_rfu, na.rm = T),
                        mean_phyco = mean(avg_phyco_rfu, na.rm = T),
                        sd_chlor = sd(avg_chlor_rfu, na.rm = T),
                        sd_phyco = sd(avg_phyco_rfu, na.rm = T))) %>% 
  mutate(scaled_chlor_rfu = (avg_chlor_rfu - mean_chlor)/sd_chlor,
         scaled_phyco_rfu = (avg_phyco_rfu - mean_phyco)/sd_phyco) %>% 
  select(-c(mean_chlor:sd_phyco)) %>% 
  # old transformation                
  #             summarise(max_chlor = max(avg_chlor_rfu, na.rm = T),
  #                       min_chlor = min(avg_chlor_rfu, na.rm = T),
  #                       max_phyco = max(avg_phyco_rfu, na.rm = T),
  #                       min_phyco = min(avg_phyco_rfu, na.rm = T))) %>% 
  # mutate(scaled_chlor_rfu = (avg_chlor_rfu- min_chlor)/max_chlor,
  #        scaled_phyco_rfu = (avg_phyco_rfu - min_phyco)/max_phyco) %>% 
  # select(-c(max_chlor:min_phyco)) %>% 
  # create trimmed / sensor adjusted values
  mutate(avg_chlor_rfu = case_when(year(sampledate) > 2018 ~ avg_chlor_rfu*1000,
                                   year(sampledate) < 2008 ~ as.numeric(NA),
                                   T ~ avg_chlor_rfu),
         avg_phyco_rfu = case_when(year(sampledate) > 2018 ~ avg_phyco_rfu*1000,
                                   year(sampledate) < 2008 ~ as.numeric(NA),
                                   T ~ avg_phyco_rfu))

#visualize data to check
hourly_buoy %>% 
  pivot_longer(avg_chlor_rfu:scaled_phyco_rfu) %>% 
  ggplot()+
  geom_point(aes(x = sampledate, y = value))+
  facet_wrap(~name, scales = "free")

## NOTES: some extremely high variables of phyco (and to less degree chl) in 2020

write_csv(hourly_buoy, "Data/Buoy/hourly_buoy.csv")

