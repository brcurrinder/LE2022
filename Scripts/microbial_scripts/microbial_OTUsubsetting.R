# Selecting top OTUs in the 16S microbial data in limony dataset
# April 2023
# AGS

pacman::p_load(tidyverse, devtools, lubridate, limony)

# Load the data_long data frame generated in the microbial_cleaning.R script
df <- data_long
# df is a data frame in the long format
# with each OTU and its abundance for a certain date
# and columns for all of its taxonomic levels

max <- df %>%
  group_by(seqID) %>%
  summarise(max = max(abund)) %>%
  ungroup() %>% arrange(desc(max))
max_subset <- max %>% slice_max(max, n = 50)

sum <- df %>%
  group_by(seqID) %>%
  summarise(sum = sum(abund)) %>%
  ungroup() %>% arrange(desc(sum))
sum_subset <- sum %>% slice_max(sum, n = 50)

mean <- df %>%
  group_by(seqID) %>%
  summarise(mean = mean(abund)) %>%
  ungroup() %>% arrange(desc(mean))
mean_subset <- mean %>% slice_max(mean, n = 50)

pres <- df %>%
  group_by(seqID) %>%
  summarise(pres = sum(abund > 0.0)) %>%
  ungroup() %>% arrange(desc(pres))
pres_subset <- pres %>% slice_max(pres, n = 50)