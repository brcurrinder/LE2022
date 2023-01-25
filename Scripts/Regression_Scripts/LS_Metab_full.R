library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)
library(broom)

metab <- read.csv('C:/Users/maxgl/Documents/RPI/GLEON/LakeExpedition/MetabolismData/SimResultsMatrix_MetabData_run12jan23.csv') %>% 
  mutate(SimDate = ymd(SimDate))

ggplot(metab,aes(SimDate,EpiNPP_mgC_L))+
  geom_point()

ls8_full <- read.csv('C:/Users/maxgl/Documents/RPI/GLEON/LakeExpedition/rs_Data/LS8_FullLake.csv') %>% 
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
  filter(pixelCount >40000)

band_ts <- ls8_full %>% 
  pivot_longer(cols=2:9,names_to='band',values_to='value')

ggplot(band_ts,aes(DATE_ACQUIRED,value))+
  geom_point()+
  geom_line()+
  facet_wrap(~band,scales='free')

combined <- merge(x = metab, y = ls8_full,
                  by.x = 'SimDate', by.y = 'DATE_ACQUIRED',
                  all.x = F, all.y = T) %>% 
  filter(!is.na(EpiPOC_mgC_L)) %>% 
  filter(month(SimDate) >= 4 & month(SimDate) < 12)

band_combos <- combined %>% 
  select(SimDate,11:18)

names <- c('Aerosol','Blue','Green','Red','NIR','SWIR1','SWIR2','Thermal')

all_ratios <- calc.rapports(band_combos, noms=names, log = FALSE, isoler = FALSE )

ar_merge <- merge(x = metab, y = all_ratios,
                  by.x = 'SimDate', by.y = 'SimDate',
                  all.x = F, all.y = T)

allratioreg <- ar_merge %>% 
  pivot_longer(cols = 18:45, names_to = 'ratio', values_to = 'value')

corrplot(cor(ar_merge[,-1],use='pairwise'),type='lower')

cors <- cor(ar_merge[,-1],use='pairwise')

#NPP--------------------------------

ggplot(allratioreg,aes(EpiNPP_mgC_L,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth()

#R--------------------------

ggplot(allratioreg,aes(EpiR_mgC_L,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#DO--------------------------

ggplot(allratioreg,aes(EpiDO_mgO2_L,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')


#Secchi------------------------------

ggplot(allratioreg,aes(log(Secchi_m),value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth()

#DOC------------------------------

ggplot(allratioreg,aes(EpiDOC_mgC_L,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth()

#Temp------------------------------

ggplot(allratioreg %>% filter(ratio == 'Thermal'),aes(EpiT_C,value))+
  geom_point()+
  geom_smooth(method='lm')

