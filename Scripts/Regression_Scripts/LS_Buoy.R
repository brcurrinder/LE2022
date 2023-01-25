library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)
library(broom)

buoy <- read.csv('C:/Users/maxgl/Documents/RPI/GLEON/LakeExpedition/MetabolismData/mendota_buoy/daily_buoy.csv') %>% 
  mutate(sampledate = ymd(sampledate)) %>% 
  filter(flag_avg_air_temp!='C'&flag_avg_rel_hum!='C'&flag_avg_wind_speed!='C'&flag_avg_wind_dir!='C'&flag_avg_chlor_rfu!='C'&flag_avg_phyco_rfu!='C'&flag_avg_par!='C'&flag_avg_par_below!='C'&flag_avg_do_wtemp!='C'&flag_avg_do_sat!='C'&flag_avg_do_raw!='C'&flag_avg_pco2_ppm!='C'&flag_avg_ph!='C'&flag_avg_fdom!='C'&flag_avg_turbidity!='C'&flag_avg_spec_cond!='C') %>% 
  select(-contains('flag'),-year4)

buoy_long <- buoy %>% 
  pivot_longer(cols=2:17,names_to = 'variable',values_to = 'value')

ggplot(buoy_long,aes(sampledate,value))+
  geom_point()+
  facet_wrap(~variable,scales='free')

ls8 <- read.csv('C:/Users/maxgl/Documents/RPI/GLEON/LakeExpedition/rs_Data/LS8_Buoy100m.csv') %>% 
  mutate(DATE_ACQUIRED = ymd(DATE_ACQUIRED)) %>% 
  rename('Aerosol' = SR_B1,
         'Blue' = SR_B2,
         'Green' = SR_B3,
         'Red' = SR_B4,
         'NIR' = SR_B5,
         'SWIR1' = SR_B6,
         'SWIR2' = SR_B7,
         'Thermal' = ST_B10) %>% 
  filter(Aerosol>0,Blue>0,Green>0,Red>0,NIR>0,SWIR1>0,SWIR2>0) %>% 
  filter(pixelCount == 32)

band_ts <- ls8 %>% 
  pivot_longer(cols=2:9,names_to='band',values_to='value')

ggplot(band_ts,aes(DATE_ACQUIRED,value))+
  geom_point()+
  geom_line()+
  facet_wrap(~band,scales='free')

combined <- merge(x = buoy, y = ls8,
                  by.x = 'sampledate', by.y = 'DATE_ACQUIRED',
                  all.x = F, all.y = T)

band_combos <- combined %>% 
  select(sampledate,18:25)

names <- c('Aerosol','Blue','Green','Red','NIR','SWIR1','SWIR2','Thermal')

all_ratios <- calc.rapports(band_combos, noms=names, log = FALSE, isoler = FALSE )

ar_merge <- merge(x = buoy, y = all_ratios,
                  by.x = 'sampledate', by.y = 'sampledate',
                  all.x = F, all.y = T)

allratioreg <- ar_merge %>% 
  pivot_longer(cols = 18:45, names_to = 'ratio', values_to = 'value')

corrplot(cor(ar_merge[,-1],use='pairwise'),type='lower')

cors <- cor(ar_merge[,-1],use='pairwise')

#Chl--------------------------------

ggplot(allratioreg%>% filter(!is.na(avg_chlor_rfu)),aes(avg_chlor_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')


#Phyco------------------------------

ggplot(allratioreg%>% filter(!is.na(avg_phyco_rfu)),aes(avg_phyco_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#fdom------------------------------

ggplot(allratioreg%>% filter(!is.na(avg_fdom)),aes(avg_fdom,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#turbidity------------------------------

ggplot(allratioreg%>% filter(!is.na(avg_turbidity)),aes(avg_turbidity,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')
