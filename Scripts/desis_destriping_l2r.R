
#setup
library(tidyverse)
library(ggforce)

source_dir = "Data/DESIS/L2R_raw_from_SMCE/"
flist = dir(source_dir, pattern = ".csv", full.names = T)

#for each member of flist
for(j in 1:length(flist)){

  # read in image name and file contents  
  img_name = str_remove(flist[j], source_dir) %>% str_remove(".csv")
  
  surface_data = read_csv(flist[j]) %>% 
    rename(id = `...1`) %>% 
    select(-contains("rhot"))

  # Check for stripes in the blue wavelengths (400-500)
  ## prepare plotting data
  plot_data = surface_data %>% 
    select(id:rhos_502) %>% 
    pivot_longer(rhos_401:rhos_502, names_to = "band", values_to = "rhos") 
 
  ## draw and save plots
  for(i in 1:10){
    ggplot(plot_data)+
      geom_raster(aes(x,-y, fill=rhos))+
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


#   -----------------------------------------------------------------------
##  test w smaller subset

# index = sample(nrow(surface_data), nrow(surface_data)/10)
# plot_data %>% 
#   filter(band == "rhos_401") %>% 
#   ggplot()+
#   geom_raster(aes(x,#=lon, 
#                   -y, #=lat, 
#                   fill = rhos))

