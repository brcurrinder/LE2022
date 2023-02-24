# use this script after clipping all DESIS destriped data to the buoy buffer
# takes all the individual files as inputs and returns one CSV

library(tidyverse)

flist = dir("Data/DESIS/L2W_clipped_to_buoy", full.names = T)

dataset = tibble()
for (data in flist){
    tempory <-read_csv(data)
    dataset <-bind_rows(dataset, tempory)
    rm(tempory)
  }

t = dataset %>% 
  left_join(dataset %>%
              group_by(date) %>% 
              summarise(pixel_ct = n())
  ) %>% 
  relocate(Rrs_909, .before = Rrs_912) %>% 
  relocate(pixel_ct, lon, lat, .after = date) %>% 
  select(-c(id, y, x, geometry, l2_flags, SPM_Nechad2010_645, ORIG_FID, BUFF_DIST))

write_csv(t, "Data/DESIS/DESIS_all_clip2buoy_destriped.csv")
