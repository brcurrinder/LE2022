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

ggplot(na.omit(buoy_long), aes(sampledate,avg))+
  geom_point()+
  facet_wrap(~variable,scales='free')


buoy_long %>% 
  filter(variable == "chlor_rfu") %>% 
  na.omit() %>% 
  mutate(indic = case_when(year(sampledate) < 2008 ~ "Set 1",
                           year(sampledate) > 2018 ~ "Set 3",
                           T ~ "Set 2")) %>% 
ggplot(aes(sampledate,avg))+
  geom_point()+
 # facet_wrap(~indic,scales='free')+
  labs(y = "avg_chlor_rfu")

#recreate wide table
buoy = buoy_long %>% 
  pivot_wider(names_from = variable, values_from = avg, names_prefix = "avg_") %>% 
#upate chl and phyco data post 2018
  mutate(avg_chlor_rfu = case_when(year(sampledate) > 2018 ~ 1000*avg_chlor_rfu,
                               T ~ avg_chlor_rfu),
         avg_phyco_rfu = case_when(year(sampledate) > 2018 ~ 1000*avg_phyco_rfu,
                               T ~ avg_phyco_rfu))


ggplot(buoy)+
  geom_point(aes(x = sampledate, y = avg_chlor_rfu))
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

## calculate algorithms

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

# select relevant measurements from the buoy
buoy_sub = buoy %>% 
  select(sampledate, avg_air_temp, avg_chlor_rfu, avg_phyco_rfu, 
         avg_do_wtemp, avg_pco2_ppm, avg_fdom, avg_turbidity)


# try predicting chl-a -------------------------------------------------------------------

# join chla data with algorithm and band values
chl_LS <- merge(x = buoy_sub[,c(1,3)], y = ls8_algs,
                  by.x = 'sampledate', by.y = 'DATE_ACQUIRED',
                  all.x = F, all.y = T) %>% 
  filter(!is.na(avg_chlor_rfu)) %>% 
  filter(avg_chlor_rfu > 1000) %>%
  # filter(avg_chlor_rfu < 1000) %>%
  mutate(log_chl = log(avg_chlor_rfu),
         sqrt_chl = sqrt(avg_chlor_rfu),
         #scale_chl = scale(avg_chlor_rfu),
         .after = avg_chlor_rfu)

# plot correlations
corrplot(cor(chl_LS[,-1],use='pairwise'),type='lower', tl.cex = 0.9, tl.col = "black", diag = F)

cors <- cor(chl_LS[,-1],use='pairwise')[,1:3]


#view subset
corrplot(t(cors), type = "upper", diag = F, tl.col = "black", cl.pos = "b", cl.ratio = 1)

## it seems not to matter what transformation of the chlorophyll values we use

# try regressions on most relevant ratios ---------------------------------------------------
avg_cors = cors %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column() %>% 
  select(rowname, cor_value = avg_chlor_rfu) %>% 
  mutate(cor_value = abs(cor_value)) %>% 
  filter(cor_value > 0.15)

chl_LS2 = chl_LS %>% 
  select(-sqrt_chl) %>% 
  pivot_longer(Aerosol:OC4) %>% 
  filter(name %in% avg_cors$rowname) %>% 
  pivot_wider(names_from = name, values_from = value)

#chl_lm1 = step(lm(log_chl ~.-avg_chlor_rfu -sqrt_chl -sampledate, data = chl_LS))

chl_lm = lm(log_chl ~ .-avg_chlor_rfu -sampledate, data = chl_LS2)
summary(chl_lm)
sub_lm = step(chl_lm)
summary(sub_lm)
plot(sub_lm)
#SABI, NDVI, G.BGR, N_R, R_B.R_N, chla_nasa
corrplot(cor(chl_LS2[,-1], use = 'pairwise'), type = "upper", method = "ellipse", 
         number.cex = 0.8, diag = F)

#var part
library(vegan)

#set row names
chl_LS2 = as.data.frame(chl_LS2)
rownames(chl_LS2) = chl_LS2$sampledate
labels = chl_LS2$sampledate

chl.dat = as.data.frame(chl_LS2[,3])
rownames(chl.dat) = labels
names(chl.dat) = "log_chl"
ls.dat = chl_LS2[,4:12]

edaPlot = rda(chl.dat ~., data = ls.dat)
plot(edaPlot, choices=c(1,2), type="text")

out = varpart(chl.dat, ~ SABI, ~ G.BGR, ~ N_R, ~ chla_nasa, data = ls.dat)
out$part$indfract
#out$part$indfract$Adj.R.square = out$part$indfract$Adj.R.square * 100 # plot as percentages
plot(out, bg = c("hotpink", "skyblue", "limegreen", "orange2"), digits = 0, Xnames = c("SABI", "NDVI", "G_R.N", "chla_nasa"))


chl_lm2 = lm(avg_chlor_rfu ~ SABI+G_R.N, data = chl_LS2)
summary(chl_lm2)
plot(chl_lm2, 2)
