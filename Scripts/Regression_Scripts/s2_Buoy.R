library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)

#import buoy data
buoy <- read.csv("Data/Buoy/daily_buoy.csv") %>%  
  mutate(sampledate = ymd(sampledate)) %>% 
  select(-year4)

#visualize buoy databuoy_long <- buoy %>% 
buoy_long <- buoy %>% 
  pivot_longer(cols=2:19,names_to = 'variable',values_to = 'value')

ggplot(na.omit(buoy_long),aes(sampledate,value))+
  geom_point()+
  facet_wrap(~variable,scales='free')

#import sentinel data

s2 <- read.csv("Data/Sentinel2/S2_Buoy100mMean.csv") %>% 
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
combined <- merge(x = buoy, y = s2,
                  by.x = 'sampledate', by.y = 'Date',
                  all.x = F, all.y = T)

#seasonal distribution of available data
combined %>%
  filter(!is.na(avg_chlor_rfu)) %>%
  ggplot()+
  geom_histogram(aes(x=month(sampledate)))

#create all unique band ratios
band_combos <- combined %>% 
  select(sampledate,22:32)

names <- c('Aerosol','Blue','Green','Red','RedEdge1','RedEdge2','RedEdge3','NIR','RedEdge4','SWIR1','SWIR2')

all_ratios <- calc.rapports(band_combos, noms=names, log = FALSE, isoler = FALSE )

ar_merge <- merge(x = buoy, y = all_ratios,
                  by.x = 'sampledate', by.y = 'sampledate',
                  all.x = F, all.y = T)

#filter out days with ice cover
ice <- read.csv('Data/IceCover/cleaned_ice_cover.csv') %>% 
  select(-X) %>% 
  mutate(Date = ymd(Date))

ice_filt <- unique(merge(x = ar_merge, y = ice,
                  by.x = 'sampledate', by.y = 'Date',
                  all.x = T, all.y = F)) %>% 
  filter(Ice == 'No')
#removes 37 rows

allratioreg <- ice_filt %>% 
  pivot_longer(cols = 20:85, names_to = 'ratio', values_to = 'value')

corrplot(cor(ar_merge[,-1],use='pairwise'),type='lower')

cors <- cor(ar_merge[,-1],use='pairwise')


#Chl--------------------------------

#make corrplot for chl-a
chlcors = cor(ar_merge[,c(6,18:53)], use = 'pairwise')
corrplot(chlcors, type = 'lower', tl.cex = 0.7)

as_tibble(chlcors, rownames = NA) %>% 
  rownames_to_column() %>% 
  select(1:2) %>% 
  mutate(abscor = abs(avg_chlor_rfu)) %>% 
  arrange(-abscor) %>% 
  head(10) %>% 
  knitr::kable()

#plot regressions

#unscaled
ggplot(allratioreg%>% filter(!is.na(avg_chlor_rfu)),aes(log(avg_chlor_rfu),value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#scaled
# ggplot(allratioreg%>% filter(!is.na(scaled_chlor_rfu)),aes(scaled_chlor_rfu,value))+
#   geom_point()+
#   facet_wrap(~ratio,scales='free')+
#   geom_smooth(method='lm')

#run linear regressions
chl_reg <- ice_filt %>% 
  select(avg_chlor_rfu,20:85) %>% 
  filter(!is.na(avg_chlor_rfu))

bands <- colnames(chl_reg)
chl_lm <- data.frame(band = rep(NA,66),
                   r2 = rep(NA,66),
                   p = rep(NA,66),
                   slope = rep(NA,66),
                   int = rep(NA,66))

for(k in 2:length(bands)){
  df = chl_reg %>% 
    select(avg_chlor_rfu,bands[k])
  
  chl_lm$band[k-1] = bands[k]
  chl_lm$r2[k-1] = summary(lm(avg_chlor_rfu~.,data=df))$adj.r.squared
  chl_lm$p[k-1] = summary(lm(avg_chlor_rfu~.,data=df))$coefficients[2,4]
  chl_lm$slope[k-1] = summary(lm(avg_chlor_rfu~.,data=df))$coefficients[2]
  chl_lm$int[k-1] = summary(lm(avg_chlor_rfu~.,data=df))$coefficients[1]
}

sig_chl_lm <- chl_lm %>% 
  filter(p <= 0.09)

#chl does not seem to have a meaningful relationship with any band ratios
#greatest R2 is 0.13

#Phyco------------------------------

#plot regressions

#unscaled
ggplot(allratioreg%>% filter(!is.na(avg_phyco_rfu)),aes(avg_phyco_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#scaled
# ggplot(allratioreg%>% filter(!is.na(scaled_phyco_rfu)),aes(scaled_phyco_rfu,value))+
#   geom_point()+
#   facet_wrap(~ratio,scales='free')+
#   geom_smooth(method='lm')

#run linear regressions
phyco_reg <- ice_filt %>% 
  select(avg_phyco_rfu,20:85) %>% 
  filter(!is.na(avg_phyco_rfu))

phyco_lm <- data.frame(band = rep(NA,66),
                     r2 = rep(NA,66),
                     p = rep(NA,66),
                     slope = rep(NA,66),
                     int = rep(NA,66))

for(k in 2:length(bands)){
  df = phyco_reg %>% 
    select(avg_phyco_rfu,bands[k])
  
  phyco_lm$band[k-1] = bands[k]
  phyco_lm$r2[k-1] = summary(lm(avg_phyco_rfu~.,data=df))$adj.r.squared
  phyco_lm$p[k-1] = summary(lm(avg_phyco_rfu~.,data=df))$coefficients[2,4]
  phyco_lm$slope[k-1] = summary(lm(avg_phyco_rfu~.,data=df))$coefficients[2]
  phyco_lm$int[k-1] = summary(lm(avg_phyco_rfu~.,data=df))$coefficients[1]
}

sig_phyco_lm <- phyco_lm %>% 
  filter(p <= 0.05)

#the best predictor is red/red edge 1, also blue/green is solid, others are good too
#greatest R2 is 0.45
#fdom------------------------------


ggplot(allratioreg%>% filter(!is.na(avg_fdom)),aes(avg_fdom,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
fdom_reg <- ice_filt %>% 
  select(avg_fdom,20:85) %>% 
  filter(!is.na(avg_fdom))

fdom_lm <- data.frame(band = rep(NA,66),
                       r2 = rep(NA,66),
                       p = rep(NA,66),
                       slope = rep(NA,66),
                       int = rep(NA,66))

for(k in 2:length(bands)){
  df = fdom_reg %>% 
    select(avg_fdom,bands[k])
  
  fdom_lm$band[k-1] = bands[k]
  fdom_lm$r2[k-1] = summary(lm(avg_fdom~.,data=df))$adj.r.squared
  fdom_lm$p[k-1] = summary(lm(avg_fdom~.,data=df))$coefficients[2,4]
  fdom_lm$slope[k-1] = summary(lm(avg_fdom~.,data=df))$coefficients[2]
  fdom_lm$int[k-1] = summary(lm(avg_fdom~.,data=df))$coefficients[1]
}

sig_fdom_lm <- fdom_lm %>% 
  filter(p <= 0.05)


#turbidity------------------------------

ggplot(allratioreg%>% filter(!is.na(avg_turbidity)),aes(avg_turbidity,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#run linear regressions
turb_reg <- ice_filt %>% 
  select(avg_turbidity,20:85) %>% 
  filter(!is.na(avg_turbidity))

turb_lm <- data.frame(band = rep(NA,66),
                      r2 = rep(NA,66),
                      p = rep(NA,66),
                      slope = rep(NA,66),
                      int = rep(NA,66))

for(k in 2:length(bands)){
  df = turb_reg %>% 
    select(avg_turbidity,bands[k])
  
  turb_lm$band[k-1] = bands[k]
  turb_lm$r2[k-1] = summary(lm(avg_turbidity~.,data=df))$adj.r.squared
  turb_lm$p[k-1] = summary(lm(avg_turbidity~.,data=df))$coefficients[2,4]
  turb_lm$slope[k-1] = summary(lm(avg_turbidity~.,data=df))$coefficients[2]
  turb_lm$int[k-1] = summary(lm(avg_turbidity~.,data=df))$coefficients[1]
}

sig_turb_lm <- turb_lm %>% 
  filter(p <= 0.05)
#very solid regressions
