library(tidyverse)
library(sf)

# set up files
source_dir = "Data/DESIS/L2W_destriped/"
flist = dir(source_dir, pattern = ".csv", full.names = T)

# read shapefile
shape = sf::read_sf("Data/Buoy_100m_shp/Buoy_100m.shp")
#plot(shape) 

for(i in 1:length(flist)){
# read each SPDF
data = read_csv(flist[i]) %>% 
  st_as_sf(coords = c("lon", "lat"))
# set CRS
st_crs(data) = st_crs(shape) #EPSG:4326

#extract filename and date
img_name = (str_remove(flist[i], source_dir)) %>% str_remove("_Mendota_clean.csv")
t = (img_name %>% str_split("_"))[[1]][4:6]
tt = paste(t[1], t[2], t[3]) %>% 
  lubridate::as_date(format = "%Y %m %d")

# clip SPDF to shapefile
cutout = st_intersection(data, shape)
#plot(cutout)

# add sample date as column in results
result = cutout %>% 
  mutate(date = tt, .before = everything())
# write file
write_csv(result, paste0("Data/DESIS/L2W_clipped_to_buoy/", img_name, "_buoy.csv"))
}
