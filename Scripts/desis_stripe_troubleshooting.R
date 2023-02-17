library(tidyverse)

clipped_data = read_csv("Data/DESIS/L2R_raw_from_SMCE/DESIS_HSI_034_2022_08_31_13_45_35_Mendota.csv") %>% 
  rename(id = `...1`)

surface = clipped_data %>% 
  select(-contains("rhot"), -y, -x)

# toa = clipped_data %>% 
#   select(-contains("rhos"))

rm(clipped_data)

# ggplot(surface)+
#   geom_point(aes(x=lon, y=lat, color = rhos_401))


# Plot spectra

index = sample(nrow(surface), nrow(surface)/100)

sub = surface %>% 
  slice(index) %>% 
  select(id, rhos_401:rhos_1000) 

meanline = sub %>% 
  summarise(across(.cols = -id, .fns = mean, na.rm = T))%>% 
  pivot_longer(rhos_401:rhos_1000, names_prefix = "rhos_", 
               names_transform = as.numeric)

sub %>% 
  pivot_longer(rhos_401:rhos_1000, names_prefix = "rhos_", 
               names_transform = as.numeric) %>% 
  ggplot(aes(x = name, y = value))+
  geom_point(aes(color = id))+
  geom_line(data = meanline, lwd = 2)+
  guides(color = "none")
  lims(x=c(500,700))
  

# https://opg.optica.org/oe/fulltext.cfm?uri=oe-22-23-28058&id=303697&ibsearch=false