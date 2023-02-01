library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)
library(broom)

s2 <- read.csv('./Data/Sentinel2/S2_FullLake.csv') %>% 
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

band_combos <- s2 %>% 
  select(Date,4:14)

names <- c('Aerosol','Blue','Green','Red','RedEdge1','RedEdge2','RedEdge3','NIR','RedEdge4','SWIR1','SWIR2')

all_ratios <- calc.rapports(band_combos, noms=names, log = FALSE, isoler = FALSE )

S2_Red_RedEdge1 <- all_ratios %>%
  select(Date,Red.RedEdge1.r)

JulianDay <- S2_Red_RedEdge1 %>%
  select(Date) %>%
  format("%j") %>%
  rename("DOY" = "Date")

S2_Red_RedEdge1 <- c(JulianDay,S2_Red_RedEdge1)
