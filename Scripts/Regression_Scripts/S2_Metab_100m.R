library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)

#read in metabolism data
metab <- read.csv('Data/Metabolism_Model_Output/SimResultsMatrix_MetabData_run12jan23.csv') %>% 
  mutate(SimDate = ymd(SimDate))

metab_long <- metab %>% 
  pivot_longer(cols=2:10,names_to = 'variable',values_to='value')

ggplot(metab_long,aes(SimDate,value))+
  geom_point()+
  facet_wrap(~variable,scales='free')

#filter NPP outliers, it seems like above 0.5 is outlier

metab <- metab %>% 
  filter(EpiNPP_mgC_L < 0.5)

#read in satellite data
s2 <- read.csv('Data/Sentinel2/S2_Buoy100mMean.csv') %>% 
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

#merge datasets
combined <- merge(x = metab, y = s2,
                  by.x = 'SimDate', by.y = 'Date',
                  all.x = F, all.y = T) %>% 
  filter(!is.na(EpiPOC_mgC_L)) %>% 
  filter(month(SimDate) >= 4 & month(SimDate) < 12)

#create band ratios
band_combos <- combined %>% 
  select(SimDate,13:23)

names <- c('Aerosol','Blue','Green','Red','RedEdge1','RedEdge2','RedEdge3','NIR','RedEdge4','SWIR1','SWIR2')

all_ratios <- calc.rapports(band_combos, noms=names, log = FALSE, isoler = FALSE)

ar_merge <- merge(x = metab, y = all_ratios,
                  by.x = 'SimDate', by.y = 'SimDate',
                  all.x = F, all.y = T)

#filter out days with ice cover
ice <- read.csv('Data/IceCover/cleaned_ice_cover.csv') %>% 
  select(-X) %>% 
  mutate(Date = ymd(Date))

ice_filt <- unique(merge(x = ar_merge, y = ice,
                         by.x = 'SimDate', by.y = 'Date',
                         all.x = T, all.y = F)) %>% 
  filter(Ice == 'No')
#removes 0 rows

allratioreg <- ice_filt %>% 
  pivot_longer(cols = 11:76, names_to = 'ratio', values_to = 'value')

corrplot(cor(ar_merge[,-1],use='pairwise'),type='lower')

cors <- cor(ar_merge[,-1],use='pairwise')

#NPP--------------------------------

ggplot(allratioreg,aes(EpiNPP_mgC_L,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
npp_reg <- ice_filt %>% 
  select(EpiNPP_mgC_L,11:76) %>% 
  filter(!is.na(EpiNPP_mgC_L))

bands <- colnames(npp_reg)
npp_lm <- data.frame(band = rep(NA,66),
                     r2 = rep(NA,66),
                     p = rep(NA,66),
                     slope = rep(NA,66),
                     int = rep(NA,66))

for(k in 2:length(bands)){
  df = npp_reg %>% 
    select(EpiNPP_mgC_L,bands[k])
  
  npp_lm$band[k-1] = bands[k]
  npp_lm$r2[k-1] = summary(lm(EpiNPP_mgC_L~.,data=df))$adj.r.squared
  npp_lm$p[k-1] = summary(lm(EpiNPP_mgC_L~.,data=df))$coefficients[2,4]
  npp_lm$slope[k-1] = summary(lm(EpiNPP_mgC_L~.,data=df))$coefficients[2]
  npp_lm$int[k-1] = summary(lm(EpiNPP_mgC_L~.,data=df))$coefficients[1]
}

sig_npp_lm <- npp_lm %>% 
  filter(p <= 0.05)

#does not seem to have a meaningful relationship with any band ratios
#greatest R2 is 0.13

#R--------------------------

ggplot(allratioreg,aes(EpiR_mgC_L,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
r_reg <- ice_filt %>% 
  select(EpiR_mgC_L,11:76) %>% 
  filter(!is.na(EpiR_mgC_L))

r_lm <- data.frame(band = rep(NA,66),
                     r2 = rep(NA,66),
                     p = rep(NA,66),
                     slope = rep(NA,66),
                     int = rep(NA,66))

for(k in 2:length(bands)){
  df = r_reg %>% 
    select(EpiR_mgC_L,bands[k])
  
  r_lm$band[k-1] = bands[k]
  r_lm$r2[k-1] = summary(lm(EpiR_mgC_L~.,data=df))$adj.r.squared
  r_lm$p[k-1] = summary(lm(EpiR_mgC_L~.,data=df))$coefficients[2,4]
  r_lm$slope[k-1] = summary(lm(EpiR_mgC_L~.,data=df))$coefficients[2]
  r_lm$int[k-1] = summary(lm(EpiR_mgC_L~.,data=df))$coefficients[1]
}

sig_r_lm <- r_lm %>% 
  filter(p <= 0.05)

#pretty similar to npp

#DO--------------------------

ggplot(allratioreg,aes(EpiDO_mgO2_L,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
do_reg <- ice_filt %>% 
  select(EpiDO_mgO2_L,11:76) %>% 
  filter(!is.na(EpiDO_mgO2_L))

do_lm <- data.frame(band = rep(NA,66),
                   r2 = rep(NA,66),
                   p = rep(NA,66),
                   slope = rep(NA,66),
                   int = rep(NA,66))

for(k in 2:length(bands)){
  df = do_reg %>% 
    select(EpiDO_mgO2_L,bands[k])
  
  do_lm$band[k-1] = bands[k]
  do_lm$r2[k-1] = summary(lm(EpiDO_mgO2_L~.,data=df))$adj.r.squared
  do_lm$p[k-1] = summary(lm(EpiDO_mgO2_L~.,data=df))$coefficients[2,4]
  do_lm$slope[k-1] = summary(lm(EpiDO_mgO2_L~.,data=df))$coefficients[2]
  do_lm$int[k-1] = summary(lm(EpiDO_mgO2_L~.,data=df))$coefficients[1]
}

sig_do_lm <- do_lm %>% 
  filter(p <= 0.05)


#Secchi------------------------------

ggplot(allratioreg,aes(Secchi_m,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
secchi_reg <- ice_filt %>% 
  select(Secchi_m,11:76) %>% 
  filter(!is.na(Secchi_m))

secchi_lm <- data.frame(band = rep(NA,66),
                    r2 = rep(NA,66),
                    p = rep(NA,66),
                    slope = rep(NA,66),
                    int = rep(NA,66))

for(k in 2:length(bands)){
  df = secchi_reg %>% 
    select(Secchi_m,bands[k])
  
  secchi_lm$band[k-1] = bands[k]
  secchi_lm$r2[k-1] = summary(lm(Secchi_m~.,data=df))$adj.r.squared
  secchi_lm$p[k-1] = summary(lm(Secchi_m~.,data=df))$coefficients[2,4]
  secchi_lm$slope[k-1] = summary(lm(Secchi_m~.,data=df))$coefficients[2]
  secchi_lm$int[k-1] = summary(lm(Secchi_m~.,data=df))$coefficients[1]
}

sig_secchi_lm <- secchi_lm %>% 
  filter(p <= 0.05)

#DOC------------------------------

ggplot(allratioreg,aes(EpiDOC_mgC_L,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
doc_reg <- ice_filt %>% 
  select(EpiDOC_mgC_L,11:76) %>% 
  filter(!is.na(EpiDOC_mgC_L))

doc_lm <- data.frame(band = rep(NA,66),
                        r2 = rep(NA,66),
                        p = rep(NA,66),
                        slope = rep(NA,66),
                        int = rep(NA,66))

for(k in 2:length(bands)){
  df = doc_reg %>% 
    select(EpiDOC_mgC_L,bands[k])
  
  doc_lm$band[k-1] = bands[k]
  doc_lm$r2[k-1] = summary(lm(EpiDOC_mgC_L~.,data=df))$adj.r.squared
  doc_lm$p[k-1] = summary(lm(EpiDOC_mgC_L~.,data=df))$coefficients[2,4]
  doc_lm$slope[k-1] = summary(lm(EpiDOC_mgC_L~.,data=df))$coefficients[2]
  doc_lm$int[k-1] = summary(lm(EpiDOC_mgC_L~.,data=df))$coefficients[1]
}

sig_doc_lm <- doc_lm %>% 
  filter(p <= 0.05)

#Temp--------------------------

ggplot(allratioreg,aes(EpiT_C,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
temp_reg <- ice_filt %>% 
  select(EpiT_C,11:76) %>% 
  filter(!is.na(EpiT_C))

temp_lm <- data.frame(band = rep(NA,66),
                     r2 = rep(NA,66),
                     p = rep(NA,66),
                     slope = rep(NA,66),
                     int = rep(NA,66))

for(k in 2:length(bands)){
  df = temp_reg %>% 
    select(EpiT_C,bands[k])
  
  temp_lm$band[k-1] = bands[k]
  temp_lm$r2[k-1] = summary(lm(EpiT_C~.,data=df))$adj.r.squared
  temp_lm$p[k-1] = summary(lm(EpiT_C~.,data=df))$coefficients[2,4]
  temp_lm$slope[k-1] = summary(lm(EpiT_C~.,data=df))$coefficients[2]
  temp_lm$int[k-1] = summary(lm(EpiT_C~.,data=df))$coefficients[1]
}

sig_temp_lm <- temp_lm %>% 
  filter(p <= 0.05)








