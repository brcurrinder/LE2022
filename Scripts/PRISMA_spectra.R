# PRISMA draw spectral curves

## Load libraries
library(tidyverse)

## Load 1 prisma scene
dat = read_csv("Data/PRISMA/lake/PRISMA_2021_06_05_17_02_29_lake.csv")
#dat = read_csv("Data/PRISMA/buoy/PRISMA_2021_06_05_17_02_29_buoy.csv")
data = dat %>% 
  select(-aot_550, -contains("rhot")) %>% 
  rename("id" = `...1`)
rm(dat)

## Make spectral curves
index = sample(1:nrow(data), nrow(data)/10)

spec_dat = data %>% 
  slice(index) %>% 
  mutate(group = rhos_807/rhos_505) %>% 
  select(-c(lon, lat), -geometry) %>%
  pivot_longer(cols = (-c(x, y, group, id)), names_prefix = "rhos_") %>% 
  mutate(name = as.numeric(name)) %>% 
  filter(name > 450, name < 850) %>% 
  na.omit()

ggplot(spec_dat, aes(x=x, y=y))+
  geom_tile(aes(fill = value))

ggplot(spec_dat, aes(x = name, y = value))+
  geom_line(aes(group = id, color = y))+
  labs(x="Wavelength", y="RS Reflectance", title = "2021-06-05 PRISMA around buoy")+
  geom_line(data = spec_dat %>% 
              group_by(name) %>% 
              summarize(value = mean(value)),
            lwd= 2)+
  lims(y = c(0.05, 0.1))+
  theme_bw()

spec_dat %>% 
  filter(y > 100, y < 200) %>% 
ggplot(aes(x = name, y = value))+
  geom_line(aes(group = id, color = y))+
  labs(x="Wavelength", y="RS Reflectance", title = "2021-06-05 PRISMA around buoy")+
  geom_line(data = spec_dat %>% 
              group_by(name) %>% 
              summarize(value = mean(value)),
            lwd= 2)+
  theme_bw()

ggplot(spec_dat)+geom_histogram(aes(x=group))+
  scale_y_log10()

ggplot(spec_dat, aes(x = name, y = value))+
  geom_line(aes(group = id, color = x))+
  labs(x="wavelength", y="reflectance")+
  facet_grid(rows = vars(group > 1, group > 10),
             cols = vars(group < 6, y < 150))
