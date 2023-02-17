
#setup
library(tidyverse)
library(ggforce)

source_dir = "Data/DESIS/L2W_raw_from_SMCE/"
flist = dir(source_dir, pattern = ".csv", full.names = T)


# draw figures ------------------------------------------------------------

#for each member of flist
for(j in 1:length(flist)){

  # read in image name and file contents  
  img_name = str_remove(flist[j], source_dir) %>% str_remove(".csv")
  
  surface_data = read_csv(flist[j]) %>% 
    rename(id = `...1`)

  # Check for stripes in the blue wavelengths (400-500)
  ## prepare plotting data
  plot_data = surface_data %>% 
    select(id:Rrs_450) %>% 
    pivot_longer(Rrs_401:Rrs_450, names_to = "band", values_to = "Rrs") 
 
  ## draw and save plots
  for(i in 1:5){
    ggplot(plot_data)+
      geom_raster(aes(x,-y, fill=Rrs))+
      facet_wrap_paginate(~band, nrow = 2, ncol = 2, page = i)+
      labs(x = "", y = "", title = img_name)+
      theme(panel.background = element_blank(), panel.grid = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank())+
      scale_fill_continuous(name = "Rrs")
    
    ggsave(paste0(source_dir, "figures/", img_name,
                  "_", i, ".png"), width = 12, height = 9, units = "in" )
    print(i)
  }
}

# Write down striped wavelengths for each file
# [1] DESIS_HSI_002_2020_06_02_20_11_33: 401
# [2] maybe ice covered
# [3] 401, 404, 409?
# [4] 401?
# [5] 401
# [6] 401, 404, 409
# [7] OK
# [8] 401
# [9] - maybe need to drop everything below 437?


# remove striped wavelengths ----------------------------------------------

#For each image:
# - read in
# - remove striped wavelengths
# - filter flagged pixels (https://oceancolor.gsfc.nasa.gov/atbd/ocl2flags/)
# - save to /destriped file

target_dir = "Data/DESIS/L2W_destriped/"

#2020_06_02
# read in image name and file contents  
img_name = str_remove(flist[1], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[1]) %>% rename(id = `...1`) 

clean_surface = surface_data %>% 
  filter(l2_flags == 0) %>% 
  mutate(Rrs_401 = as.numeric(NA))

write_csv(clean_surface, paste0(target_dir, img_name, "_clean.csv"))

## 2020_12_25

# read in image name and file contents  
img_name = str_remove(flist[2], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[2]) %>% rename(id = `...1`) 

clean_surface = surface_data %>% 
  filter(l2_flags == 0) %>% 
  mutate(across(Rrs_401:Rrs_450, replace, values = as.numeric(NA)))

write_csv(clean_surface, paste0(target_dir, img_name, "_clean.csv"))

## 2021_02_12
# read in image name and file contents  
img_name = str_remove(flist[3], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[3]) %>% rename(id = `...1`) 

clean_surface = surface_data %>% 
  filter(l2_flags != 21,
         !is.na(Rrs_401)) 
# NO AVAIL PIXELS

## 2021_06_15 - missing columns
img_name = str_remove(flist[4], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[4]) %>% rename(id = `...1`)

clean_surface = surface_data %>% 
  filter(l2_flags == 0) %>% 
  mutate(across(Rrs_401, replace, values = as.numeric(NA)))

ggplot(clean_surface)+
  #geom_raster(aes(x=x, y=-y, fill = as_factor(l2_flags)))+
  geom_raster(aes(x=x, y=-y, fill = Rrs_445))

write_csv(clean_surface, paste0(target_dir, img_name, "_clean.csv"))

## 2021_08_15
# read in image name and file contents  
img_name = str_remove(flist[5], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[5]) %>% rename(id = `...1`)

clean_surface = surface_data %>% 
  filter(l2_flags == 0) %>% 
  mutate(across(c(Rrs_401, Rrs_448, Rrs_450), replace, values = as.numeric(NA)))

write_csv(clean_surface, paste0(target_dir, img_name, "_clean.csv"))

## 2022_06_23
# read in image name and file contents  
img_name = str_remove(flist[6], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[6]) %>% rename(id = `...1`)

clean_surface = surface_data %>% 
  filter(l2_flags == 0) %>% 
  mutate(across(c(Rrs_401, Rrs_448, Rrs_404), replace, values = as.numeric(NA)))

write_csv(clean_surface, paste0(target_dir, img_name, "_clean.csv"))

## 2022_06_02
# read in image name and file contents  
img_name = str_remove(flist[7], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[7]) %>% rename(id = `...1`)

clean_surface = surface_data %>% 
  filter(l2_flags == 0) %>% 
  mutate(across(Rrs_401, replace, values = as.numeric(NA)))

write_csv(clean_surface, paste0(target_dir, img_name, "_clean.csv"))

## 2022_10_10
# read in image name and file contents  
img_name = str_remove(flist[8], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[8]) %>% rename(id = `...1`)

clean_surface = surface_data %>% 
  filter(l2_flags == 0) %>% 
  mutate(across(Rrs_401, replace, values = as.numeric(NA)))

write_csv(clean_surface, paste0(target_dir, img_name, "_clean.csv"))

## 2022_08_31
# read in image name and file contents  
img_name = str_remove(flist[9], source_dir) %>% str_remove(".csv")
surface_data = read_csv(flist[9]) %>% rename(id = `...1`)

clean_surface = surface_data %>% 
  filter(l2_flags == 0) %>% 
  mutate(across(Rrs_401:Rrs_458, replace, values = as.numeric(NA)))

ggplot(clean_surface)+
  #geom_raster(aes(x=x, y=-y, fill = as_factor(l2_flags)))+
  geom_raster(aes(x=x, y=-y, fill = Rrs_665))

write_csv(clean_surface, paste0(target_dir, img_name, "_clean.csv"))

#   -----------------------------------------------------------------------
##  test w smaller subset

# index = sample(nrow(surface_data), nrow(surface_data)/10)
# plot_data %>% 
#   filter(band == "rhos_401") %>% 
#   ggplot()+
#   geom_raster(aes(x,#=lon, 
#                   -y, #=lat, 
#                   fill = rhos))

# surface_data %>% 
#   mutate(l2_flags= as_factor(l2_flags)) %>% 
#   filter(l2_flags==0,
#          # Rrs_865 > 0.003,
#          Rrs_865 < 0.02) %>% 
#   ggplot()+
#   #geom_raster(aes(x,-y, fill=l2_flags))+
#   geom_raster(aes(x,-y, fill=Rrs_858))+
#   labs(x = "", y = "", title = img_name)+
#   theme(panel.background = element_blank(), panel.grid = element_blank(),
#         axis.text = element_blank(), axis.ticks = element_blank())

