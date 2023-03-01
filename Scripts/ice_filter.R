library(tidyverse)
library(lubridate)

ice <- read.csv('Data/IceCover/ntl33_v7.csv') %>% 
  filter(lakeid == "ME",
         year4 >= 1983) %>%
  mutate(ice_off = mdy(ice_off),
         ice_on = mdy(ice_on),
         year = year4) %>%
  select(ice_on, ice_off, season, year)

df_list <- list()
for(i in 1983:2019){
  dates = ice %>% 
    filter(year == i)
  
  df = data.frame(Date = seq(dates$ice_on[1],dates$ice_off[1],by = 1),
                  Ice = 'Yes')
  
  df_list[[i-1982]] = df
}

ice_dates <- do.call(rbind,df_list)  

dates_shell <- data.frame(Date = seq(ice_dates$Date[1],ice_dates$Date[nrow(ice_dates)],by = 1))

ice_cover <- merge(x = dates_shell, y = ice_dates,
                   by.x='Date', by.y = 'Date',
                   all = TRUE) %>% 
  mutate(Ice = case_when(is.na(Ice) ~ 'No',
                         !is.na(Ice) ~ 'Yes'))

add_ons <- data.frame(Date = seq(ymd('2020-01-12'),ymd('2022-12-31'),by=1),
                      Ice = NA) %>% 
  mutate(Ice = case_when(Date >= ymd('2020-01-12') & Date <= ymd('2020-03-22') |
                   Date >= ymd('2021-01-03') & Date <= ymd('2021-03-20') |
                   Date >= ymd('2022-01-07') & Date <= ymd('2022-04-02') |
                   Date >= ymd('2022-12-25') ~ 'Yes')) %>% 
  mutate(Ice = case_when(is.na(Ice) ~ 'No',
                         !is.na(Ice) ~ 'Yes'))

full_ice_cover <- rbind(ice_cover,add_ons)

write.csv(full_ice_cover,file='Data/IceCover/cleaned_ice_cover.csv')
