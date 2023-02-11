library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)
library(broom)


#import buoy data
buoy = read.csv("Data/Buoy/daily_buoy.csv") %>%  
  mutate(sampledate = ymd(sampledate)) %>% 
  select(-year4)

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

#import sentinel data

s2 <- read.csv("Data/Sentinel2/S2_FullLake.csv") %>% 
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
  filter(pixelCount >300000) %>% 
  filter(MGRS_TILE == '15TYH')

band_ts <- s2 %>% 
  pivot_longer(cols=4:14,names_to='band',values_to='value')

ggplot(band_ts,aes(Date,value))+
  geom_point()+
  geom_line()+
  facet_wrap(~band,scales='free')

# ggplot(s2,aes(Date,Blue/Green))+
#   geom_point()+
#   geom_line()
# 
# ggplot(buoy %>%
#          filter(sampledate > as.Date('2019-01-01')),aes(sampledate,avg_phyco_rfu))+
#   geom_point()+
#   geom_line()

combined <- merge(x = buoy, y = s2,
                  by.x = 'sampledate', by.y = 'Date',
                  all.x = F, all.y = T)

#seasonal distribution of available data
combined %>%
  filter(!is.na(avg_chlor_rfu)) %>%
  ggplot()+
  geom_histogram(aes(x=month(sampledate)))

band_combos <- combined %>% 
  select(sampledate,20:30)

names <- c('Aerosol','Blue','Green','Red','RedEdge1','RedEdge2','RedEdge3','NIR','RedEdge4','SWIR1','SWIR2')

all_ratios <- calc.rapports(band_combos, noms=names, log = FALSE, isoler = FALSE )

ar_merge <- merge(x = buoy, y = all_ratios,
                  by.x = 'sampledate', by.y = 'sampledate',
                  all.x = F, all.y = T)

allratioreg <- ar_merge %>% 
  pivot_longer(cols = 18:83, names_to = 'ratio', values_to = 'value')

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
ggplot(allratioreg%>% filter(!is.na(avg_chlor_rfu)),aes(avg_chlor_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

#Phyco------------------------------

ggplot(allratioreg%>% filter(!is.na(avg_phyco_rfu)),aes(avg_phyco_rfu,value))+
  geom_point()+
  facet_wrap(~ratio,scales='free')+
  geom_smooth(method='lm')

phyco_reg <- ar_merge %>% 
  select(avg_phyco_rfu,18:83) %>% 
  filter(!is.na(avg_phyco_rfu))

summary(lm(avg_phyco_rfu~Red.RedEdge1.r,data=phyco_reg))

ind_lm_ <- phyco_reg %>% 
  pivot_longer(cols = 2:67, names_to = 'band', values_to = 'value') %>% 
  group_by(band) %>% 
  nest() %>% 
  mutate(trend = map(data,~lm(avg_phyco_rfu~value,data=.x))) %>% 
  mutate(slope = map(trend,~tidy(.x))) %>% 
  unnest(slope)

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

