library(tidyverse)
library(lubridate)
library(corrplot)
library(SARP.compo)
library(broom)


# read buoy data ----------------------------------------------------------

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

# read LS data ------------------------------------------------------------

ls8 = read.csv("Data/Landsat/LS8_Buoy100m.csv") %>% 
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
  filter(pixelCount == 32) %>% 
  select(-pixelCount)

band_ts <- ls8 %>% 
  pivot_longer(cols=2:9,names_to='band',values_to='value')

ggplot(band_ts,aes(DATE_ACQUIRED,value))+
  geom_point()+
  geom_line()+
  facet_wrap(~band,scales='free')


# describe algorithms -----------------------------------------------------

#https://oceancolor.gsfc.nasa.gov/atbd/chlor_a/
# CI = G - (B + (Yg-Yb)/(Yr-Yb)*(R-B))
# chlorA = 10^(-0.4287 + 230.47*CI)

# OCX from oreilly 1998
# MCP: Chl = 10^(a0 + a1*R + a2*R^2 + a3*R^3) + a4
# OC4
#R = log((Rrs443 > Rrs490 > Rrs510)/Rrs555)
# R = log(pmax(Aerosol, Blue)/Green)
#OC4_chl = 10^(a0 + a1*R + a2*R^2 + a3*R^3) + a4
# a0 = 0.4708; a1 = -3.8469; a2 = 4.5338; a3 = -2.4434; a4 = -0.0414


# https://www.mdpi.com/1660-4601/12/9/10391
# NIR/aerosol, NIR/blue, NIR/green, NIR/red

# https://www.mdpi.com/2072-4292/13/22/4607 Dallosch and Creed 2021
# choose algorithms based on OWT (mendota is likely A, BC, or F)
# A: B/R or (R/G)*N
# B: B/G
# C: (B/G)*(R/G) or B*N
# F: (G/R)*N or G*(B+G+R)
# Global: (R/B)*(R/N)

#https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/eap.1708
# SABI (N-R)/(B+G)
# KIVU (B-R)/G
# NDVI (N-R)/(N+R)
# 2BDA N/R
# Kab1 1.67 − 3.94*ln(B) + 3.78*ln(G)
# Kab2 6.92274 − 5.7581*(ln(A)/ln(G))

ls8_algs = ls8 %>% 
  mutate(SABI = (NIR-Red)/(Blue+Green),
         KIVU = (Blue-Red)/Green,
         NDVI = (NIR-Red)*(NIR+Red),
         Kab1 = 1.67-3.94*log(Blue) + 3.78*log(Green),
         Kab2 = 6.92274 - 5.7581*(log(Aerosol)/log(Green)),
         B_R = Blue/Red,
         R_G.N = (Red/Green)*NIR,
         B_G = Blue/Green,
         B_G.R_G = (Blue/Green)*(Red/Green),
         B.N = Blue*NIR,
         G_R.N = (Green/Red)*NIR,
         G.BGR = Green*(Blue+Green+NIR),
         R_B.R_N = (Red/Blue)*(Red/NIR),   
         N_A = NIR/Aerosol,
         N_B = NIR/Blue,
         N_G = NIR/Green,
         N_R = NIR/Red,
         CI_nasa = Green - (Blue +(560-480)/(655-480)*(Red-Blue)),
         chla_nasa = 10^(-0.4287 + 230.47*CI_nasa),
         R_oreilly = log(pmax(Aerosol, Blue)/Green),
         OC4 = 10^(0.4708 - 3.8469*R_oreilly + 4.5338*R_oreilly^2 - 2.4434*R_oreilly^3)-0.0414
         )

# combine data ------------------------------------------------------------

buoy_sub = buoy %>% 
  select(sampledate, avg_air_temp, avg_chlor_rfu, avg_phyco_rfu, 
         avg_do_wtemp, avg_pco2_ppm, avg_fdom, avg_turbidity)


# chl-a -------------------------------------------------------------------


# check out different combinations against chlorophyll
chl_LS <- merge(x = buoy_sub[,c(1,3)], y = ls8_algs,
                  by.x = 'sampledate', by.y = 'DATE_ACQUIRED',
                  all.x = F, all.y = T) %>% 
  filter(!is.na(avg_chlor_rfu)) %>% 
  filter(avg_chlor_rfu > 1000) %>%
  # filter(avg_chlor_rfu < 1000) %>%
  mutate(log_chl = log(avg_chlor_rfu),
         sqrt_chl = sqrt(avg_chlor_rfu),
         .after = avg_chlor_rfu)

corrplot(cor(chl_LS[,-1],use='pairwise'),type='lower', tl.cex = 0.9, tl.col = "black", diag = F)

cors <- cor(chl_LS[,-1],use='pairwise')[,1:3]
corrplot(t(cors), type = "upper", diag = F, tl.col = "black", cl.pos = "b", cl.ratio = 1)

# try ratio regressions ---------------------------------------------------
avg_cors = tibble("cor_value" = abs(cors[-c(1:3),1]), "name" = names(cors[-c(1:3),1])) %>% 
  filter(cor_value > 0.15)

chl_LS = chl_LS %>% 
  select(-sqrt_chl, -log_chl) %>% 
  pivot_longer(Aerosol:OC4) %>% 
  filter(name %in% avg_cors$name) %>% 
  pivot_wider(names_from = name, values_from = value)

chl_lm = lm(avg_chlor_rfu ~ ., data = chl_LS)
summary(chl_lm)
sub_lm = step(chl_lm)
summary(sub_lm)

corrplot(cor(chl_LS[,-1], use = 'pairwise'), type = "upper", diag = F)
