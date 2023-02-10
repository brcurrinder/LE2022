### Make plots for metab saved csv output 

library(tidyverse)
library(lubridate)
library(ggpmisc)

#manually set working directory to metab folder 

##read in data 
data <- read_csv("C:/Users/dwh18/OneDrive/Desktop/MetabolismMendotaLE2022/MetabolismMendotaLE2022_5janAdded/Output/SimResultsMatrix.csv")
seasons <- read_csv("C:/Users/dwh18/OneDrive/Desktop/MetabolismMendotaLE2022/seasons_from_Robins_paper.csv")
head(seasons)


#seaons
head(seasons)

# seasons_summary_spring <- seasons %>% 
#   select(Year, Spring) %>% 
#   rename(value = Spring ) %>% 
#   summarise(min = min(value),
#             # quantile25 = quantile(value, probs = c(.25), na.rm = T),
#             mean = mean(value, na.rm = T),
#             median = median(value, na.rm = T),
#             # quantile75 = quantile(value, probs = c(.75), na.rm = T),
#             max = max(value),
#             sd = sd(value, na.rm = T)
#   ) 



#### buoy #### 
buoy <- read_csv("C:/Users/dwh18/OneDrive/Desktop/MetabolismMendotaLE2022/mendota_buoy/mendota_buoy/daily_buoy.csv")

head(buoy)

buoy %>% 
  ggplot(aes(x = sampledate, y = avg_chlor_rfu))+
  geom_point()+
  ggtitle("Chlorophyll (RFU)")+
  theme_classic()
summary(buoy$avg_chlor_rfu)

buoy %>% 
  ggplot(aes(x = sampledate, y = avg_phyco_rfu))+
  geom_point()+
  ggtitle("Phycocyanin (RFU)")+
  theme_classic()
summary(buoy$avg_phyco_rfu)


buoy %>% 
  ggplot(aes(x = sampledate, y = avg_fdom))+
  geom_point()+
  ggtitle("fDOM (RFU)")+
  theme_classic()
summary(buoy$avg_fdom)


buoy %>% 
  ggplot(aes(x = sampledate, y = avg_turbidity))+
  geom_point()+
  ggtitle("Turbidity (FNU)")+
  theme_classic()
summary(buoy$avg_turbidity)


buoy %>% 
  ggplot(aes(x = sampledate, y = avg_do_raw))+
  geom_point()+
  ggtitle("DO (mg/L)")+
  theme_classic()
summary(buoy$avg_do_raw)


#buoy, chla, phytos to metab 
head(data)
head(buoy)
joined <- left_join(data, buoy, by = c("SimDate" = "sampledate"))

# Data set title: North Temperate Lakes LTER: Chlorophyll - Madison Lakes Area 1995 - current.
inUrl5  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/38/28/66796c3bc77617e7cc95c4b09d4995c5" 
infile5 <- tempfile()
try(download.file(inUrl5,infile5,method="curl"))
chlorophyll = read_csv(infile5)
# Filter for Mendota 
chlorophyll = chlorophyll %>% filter(lakeid == 'ME') 

# Data set title: North Temperate Lakes LTER: Phytoplankton - Madison Lakes Area 1995 - current.
inUrl6  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/88/31/f2de15b2fff6ae962a04c150c0a1c510" 
infile6 <- tempfile()
try(download.file(inUrl6,infile6,method="curl"))
phytoplankton = read_csv(infile6)
# Filter for Mendota 
phytoplankton = phytoplankton %>% filter(lakeid == 'ME')

phytos_chla <- left_join(chlorophyll, phytoplankton, by = "sampledate")

joined <- left_join(joined, phytos_chla, by = c("SimDate" = "sampledate"))


joined %>% 
  # filter(EpiNPP_mgC_L >= 0.1) %>% 
  # filter(avg_chlor_rfu >= 50) %>% 
  ggplot(aes(x = avg_chlor_rfu, y = EpiNPP_mgC_L))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq()+
  # geom_smooth()+
  theme_classic()

joined %>% 
  filter(avg_phyco_rfu > 0) %>% 
  ggplot(aes(x = avg_phyco_rfu, y = EpiNPP_mgC_L))+
  stat_poly_line()+
  stat_poly_eq()+
  geom_point()+
  theme_classic()

joined %>% 
  ggplot(aes(x = cells_per_ml, y = EpiNPP_mgC_L))+
  stat_poly_line()+
  stat_poly_eq()+
  geom_point()+
  theme_classic()

joined %>% 
  ggplot(aes(x = mono_chl_spec, y = EpiNPP_mgC_L))+
  stat_poly_line()+
  stat_poly_eq()+
  geom_point()+
  ggtitle("Chl-a via spec")+
  theme_classic()

joined %>% 
  ggplot(aes(x = correct_chl_fluor, y = EpiNPP_mgC_L))+
  stat_poly_line()+
  stat_poly_eq()+
  geom_point()+
  ggtitle("Chl-a via flourometer")+
  theme_classic()

joined %>% 
  ggplot(aes(x = phaeo_spec, y = EpiNPP_mgC_L))+
  stat_poly_line()+
  stat_poly_eq()+
  geom_point()+
  ggtitle("Phaeo via spec")+
  theme_classic()

joined %>% 
  ggplot(aes(x = phaeo_fluor, y = EpiNPP_mgC_L))+
  stat_poly_line()+
  stat_poly_eq()+
  geom_point()+
  ggtitle("phaeo via flourometer")+
  theme_classic()

joined %>% 
  select(SimDate, EpiNPP_mgC_L, avg_chlor_rfu, avg_phyco_rfu, 
         mono_chl_spec, correct_chl_fluor, phaeo_spec, phaeo_fluor,
         cells_per_ml, biovolume_conc, biomass_conc) %>% 
  pivot_longer(cols = -c(1:2)) %>% 
  ggplot(aes(x = value, y = EpiNPP_mgC_L))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(small.p = T)+
  facet_wrap(~name, scales = "free_x")+
  theme_classic()


#look at metab data ----
data %>% 
  select(SimDate, EpiPOC_mgC_L, EpiDOC_mgC_L) %>% 
  pivot_longer(cols = c(2:3)) %>% 
  ggplot(aes(x = SimDate, y = value))+
    geom_point()+
    facet_wrap(~name, ncol = 1)+
    theme_classic()


data %>% 
  select(SimDate, EpiDO_mgO2_L, EpiDO_sat) %>% 
  pivot_longer(cols = c(2:3)) %>% 
  ggplot(aes(x = SimDate, y = value))+
  geom_point()+
  facet_wrap(~name, ncol = 1)+
  theme_classic()
  


data %>% 
  select(SimDate, EpiNPP_mgC_L, EpiR_mgC_L) %>% 
  pivot_longer(cols = c(2:3)) %>% 
  ggplot(aes(x = SimDate, y = value))+
  geom_point()+
  facet_wrap(~name, ncol = 1)+
  theme_classic()

data %>% 
  select(SimDate, EpiI_w_m2, Secchi_m) %>% 
  pivot_longer(cols = c(2:3)) %>% 
  ggplot(aes(x = SimDate, y = value))+
  geom_point()+
  facet_wrap(~name, ncol = 1)+
  theme_classic()


nppseasons <- data %>% 
  select(SimDate, EpiNPP_mgC_L) %>% 
  # mutate(Year = year(SimDate),
  #        monthday = paste(month(SimDate), day(SimDate), sep = "-")) %>% 
  ggplot(aes(x = SimDate, y = EpiNPP_mgC_L))+
  geom_point()+
  # facet_wrap(~Year)+
  theme_classic()

nppseasons


wintermonths <- c(12,1,2)
springmonths <- c(3,4,5)
summermonths <- c(6,7,8)
fallmonths <- c(9,10,11)

data_seasons <- data %>% 
  mutate(Month = month(SimDate)) %>% 
  mutate(Season = ifelse(Month %in% wintermonths, "Winter", NA),
         Season = ifelse(Month %in% springmonths, "Spring", Season),
         Season = ifelse(Month %in% summermonths, "Summer", Season),
         Season = ifelse(Month %in% fallmonths, "Fall", Season))

data_seasons %>% 
  select(SimDate, EpiNPP_mgC_L, Season) %>% 
  filter(SimDate >= ymd("2015-01-01")) %>%
  mutate(Year = year(SimDate)) %>% 
  mutate(FakeYeardate = ymd(paste("2024", month(SimDate), day(SimDate), sep = "-"))) %>% 
  # mutate(monthday = paste(month(SimDate), day(SimDate), sep = "-")) %>%
  ggplot(aes(x = FakeYeardate, y = EpiNPP_mgC_L, col = Season))+
  geom_point()+
  facet_wrap(~Year, ncol = 2)+
  scale_x_date(date_labels="%b")+
  xlab(element_blank())+
  theme_classic()

head(seasons)
head(data)
data$Year <- year(data$SimDate)
data_robins_seasons <- left_join(data, seasons, by = "Year")


####Metab DWH ran 12jan22
#create csv after running model in MetabolismWrapper.R
source('WriteResultsMatrix.R') #read in function needed 
getwd()  #need to make sure working directory is in main metab folder ie Metab..._5JanAdded, where Output is folder, thats where the function will write to
WriteResultsMatrix = function(OutputFileListName) #Output File List name needs to be the csv name generated from Wrapper 
#This gives csv named SimResultsMatrix.csv in Output folder, I renamed to date ran and my initals to put on google drive and read in below 


metab_2020_correctmet <- read_csv("C:/Users/dwh18/OneDrive/Desktop/MetabolismMendotaLE2022/MetabolismMendotaLE2022_5janAdded/Output/SimResultsMatrix.csv")
plot(metab_2020_correctmet$SimDate, metab_2020_correctmet$EpiNPP_mgC_L)

metab_2020_correctmet %>% 
  select(SimDate, EpiNPP_mgC_L) %>% 
  # mutate(Year = year(SimDate),
  #        monthday = paste(month(SimDate), day(SimDate), sep = "-")) %>% 
  ggplot(aes(x = SimDate, y = EpiNPP_mgC_L))+
  geom_point()+
  # facet_wrap(~Year)+
  theme_classic()

#pull buoy and metab for HS ----
#2016 aug 22 and 31
dates2016 <- c(ymd("2016-08-22"), ymd("2016-08-31"))
head(buoy)
buoyHS16 <- buoy %>% 
  filter(sampledate %in% dates2016)

head(data)
metabHS16 <- data %>% 
  filter(SimDate %in% dates2016) %>% 
  rename("sampledate" = "SimDate")


joinHS16 <- full_join(buoyHS16, metabHS16, by = c("sampledate"))
getwd()
# write.csv(joinHS16, "./Aug2016_buoy_and_metab_data.csv", row.names = F) 





