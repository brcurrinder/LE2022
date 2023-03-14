## Summarizing HS data ###
# GLEON LEXP 22 ###
#start date: 3/8/2023 ###

#load packages
library(tidyverse)

#load file-DESIS for now
setwd()
desis100<-read.csv('Data/DESIS/DESIS_all_clip2buoy_destriped.csv')

#Data summarize -mean, median, std dev
desis_stats<-desis100 %>%
  select(-c(pixel_ct,lon,lat)) %>% #drop pixel count, Long, Lat
  group_by(date) %>%
  summarise_all(list(avg=mean,med=median,
                     sd=sd)) #overall mean & median are relatively similar
#export desis summary file
write.csv(desis_stats,'Data/DESIS/DESIS_stats_use.csv')