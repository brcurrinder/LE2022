library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)
library(broom)

#used these instructions https://community.rstudio.com/t/how-to-read-csv-file-from-googledrive/135922 to get the data directly from our google drive so that it will work on anyone's computer

id = "1fDTIKP_XCfe9OobnxjTjiJsRw1kEJDhh"
buoy = read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id)) %>%  
  mutate(sampledate = ymd(sampledate)) %>% 
  select(-year4)
rm(id) 

# make sure data values are NA when flagged, preserving as many points as possible
buoy_long = buoy %>% 
  pivot_longer(cols = starts_with("avg"), names_prefix = "avg_", names_to = "variable", values_to = "avg") %>% 
  pivot_longer(cols = starts_with("flag"), names_prefix = "flag_avg_", names_to = "var_flag", values_to = "flag") %>% 
  filter(variable == var_flag) %>% 
  mutate(avg = case_when(flag != "" ~ as.numeric(NA),
                         T ~ avg)) %>% 
  select(-var_flag, -flag)

#visualize buoy data
ggplot(na.omit(buoy_long),aes(sampledate,avg))+
  geom_point()+
  facet_wrap(~variable,scales='free')

#recreate wide table
buoy = buoy_long %>% 
  pivot_wider(names_from = variable, values_from = avg, names_prefix = "avg_")

# read LS data

id = "1gDUSm9cgG1x8ails5ZkCtVsECnjXX2Vh"
ls8 = read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id)) %>% 
#ls8 <- read.csv('C:/Users/maxgl/Documents/RPI/GLEON/LakeExpedition/rs_Data/LS8_Buoy100m.csv') %>% 
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
rm(id)

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
  pivot_longer(cols = 18:53, names_to = 'ratio', values_to = 'value')

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
allratioreg%>% filter(!is.na(avg_turbidity)) %>% 
  filter(avg_turbidity < 7) %>% #remove 1 outlier?
ggplot(aes(avg_turbidity,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')
