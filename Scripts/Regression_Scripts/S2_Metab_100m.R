library(tidyverse)
library(lubridate)

metab <- read.csv('C:/Users/maxgl/Documents/RPI/GLEON/LakeExpedition/MetabolismData/SimResultsMatrix_MetabData_run12jan23.csv') %>% 
  mutate(SimDate = ymd(SimDate))

ggplot(metab,aes(SimDate,EpiNPP_mgC_L))+
  geom_point()

s2 <- read.csv('C:/Users/maxgl/Documents/RPI/GLEON/LakeExpedition/rs_Data/S2_Buoy100mMean.csv') %>% 
  mutate(Date = substr(DATATAKE_IDENTIFIER, 6, 13),
         Year = substr(Date,0,4),
         Month = substr(Date,5,6),
         Day = substr(Date,7,8),
         Date = ymd(paste(Year,Month,Day,sep='-'))) %>% 
  select(-DATATAKE_IDENTIFIER,-Year,-Month,-Day) %>% 
  relocate(Date) %>% 
  rename('Aerosol' = B1,
         'Blue' = B2,
         'Green' = B3,
         'Red' = B4,
         'RedEdge1' = B5,
         'RedEdge2' = B6,
         'RedEdge3' = B7,
         'NIR' = B8,
         'RedEdge4' = B8A,
         'SWIR1' = B11,
         'SWIR2' = B12) %>% 
  filter(Aerosol>0,Blue>0,Green>0,Red>0,NIR>0,SWIR1>0,SWIR2>0) %>% 
  filter(pixelCount >300) %>% 
  filter(MGRS_TILE == '15TYH')

band_ts <- s2 %>% 
  pivot_longer(cols=4:14,names_to='band',values_to='value')

ggplot(band_ts,aes(Date,value))+
  geom_point()+
  geom_line()+
  facet_wrap(~band,scales='free')

combined <- merge(x = metab, y = s2,
                  by.x = 'SimDate', by.y = 'Date',
                  all.x = F, all.y = T) %>% 
  filter(!is.na(EpiPOC_mgC_L)) %>% 
  filter(month(SimDate) >= 4 & month(SimDate) < 12)


band_combos <- combined %>% 
  select(SimDate,13:23)

names <- c('Aerosol','Blue','Green','Red','RedEdge1','RedEdge2','RedEdge3','NIR','RedEdge4','SWIR1','SWIR2')

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
  geom_smooth(method='lm')

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

