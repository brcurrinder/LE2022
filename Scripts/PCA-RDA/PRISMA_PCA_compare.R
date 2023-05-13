# PRISMA PCA comparisons

## Load libraries
library(sp)
library(raster)
library(RStoolbox)
library(tidyverse)

# Load list of scenes

flist = dir("Data/PRISMA/lake/")
#Initialize storage for variances and components
vars = tibble()
lds = tibble()

#Loop over each file to run PCA and save top components
#currently: saving everything that represents at least 1% of variance, plus one

for(i in 1:length(flist)){
# save image date to identify observations
  img_date = flist[i] %>% substr(8, 17)
# read data
  dat = read_csv(paste0("Data/PRISMA/lake/", flist[i])) %>% 
  select(-aot_550, -contains("rhot")) %>% #remove unneeded columns (TOA radiance)
  rename("id" = `...1`) %>% # clean up column names
  sf::st_as_sf(coords = c("x", "y")) #convert to SpatialPointsDataFrame

#convert to raster
  dat_ras = raster(dat, resolution = c(1,1)) #each point becomes one raster pixel
  #trim to reasonable spectral range. this is for 450:850nm
  colnames = names(dat)[which(names(dat)=="rhos_453"):which(names(dat)=="rhos_849")]
  #fill the raster with data from each band
  dat_ras = rasterize(dat, dat_ras, field = colnames)
  
# run pca
  rpc = rasterPCA(dat_ras, spca = T #, nSamples = ncell(dat_ras)/10 #can uncomment this to subsample and run faster
                  )
  summary(rpc$model) #prints loadings
  screeplot(rpc$model) #shows variance represented by each component
  ggRGB(rpc$map, 1,2,3, stretch = "lin") #maps the top 3 components
  #ggsave(paste0("CO_ignore/",img_date,".png")) #change filepath if you want to save these images

# get variances represented by each PC
  variances = as_tibble(rpc$model$sdev^2, rownames = "component") %>% #convert SD to variance
    mutate(proportion = value/sum(value), #proportion of total variance represented
           img_date = img_date) #add image information in a separate column
# decide how many components to retain
  variances = variances %>% 
    filter(proportion >= 0.01) #just keep those representing at least 1% of variance
    #slice_head(n = (nrow(variances[variances$proportion >= 0.01,])+1)) #or keep those plus one
  
# get loadings for each variance
  loadings = rpc$model$loadings[,1:nrow(variances)] %>% #only keep the ones that have been retained
    as_tibble(rownames = "band") %>% 
    janitor::clean_names() %>%
    pivot_longer(cols = starts_with("comp")) %>% #this will cause an error if only retaining one PC - just comment out this line and fix it manually for that case
    mutate(band = as.numeric(str_remove(band, "rhos_")), #format wavelengths for plotting
           img_date = img_date) #tag with image date
  
# plot loadings  
  ggplot(loadings, aes(x = band, y = value))+
    theme_bw()+
    geom_line(aes(group = name, color = name))+
    scale_color_discrete(name = "Component")+
    labs(x = "Wavelength", y = "Loading")

# save information from this image
  vars = bind_rows(vars, variances)
  lds = bind_rows(lds, loadings)
}

# save information to files
# write_csv(lds, "CO_ignore/prisma_PCA_loadings.csv")
# write_csv(vars, "CO_ignore/prisma_PCA_variances.csv")

#plot all PC loadings

lds %>% 
  mutate(grp = paste(name, img_date),
         #following line is to group the panels
         cat = case_when(name %in% c("comp_1", "comp_2") ~ "PC 1, 2",
                         name %in% c("comp_4", "comp_5") ~ "PC 4, 5",
                         name %in% c("comp_6", "comp_7") ~ "PC 6, 7",
                         T ~ "PC 3")) %>% 
  ggplot(aes(x = band, y = value))+
  geom_line(aes(group = grp, color= img_date, lty = name))+
  scale_linetype_manual(values = c(1,2,1,2,1,2,1), guide = "none")+
  facet_wrap(~ cat)+
  scale_color_discrete(name = "Image Date")+
  labs(x = "Wavelength", y = "Loading", title = "PRISMA: all PCs over 1% variance")+
  theme_bw()
  
## scratch
vars = read_csv("CO_ignore/prisma_PCA_variances.csv") %>% 
  mutate(component = component %>% str_replace("Comp.", "comp_")) %>% 
  select(-value)

data = lds %>% 
  left_join(vars, by = c("name" = "component", "img_date"))

data %>% 
  mutate(grp = paste(name, img_date),
         #following line is to group the panels by shape
         cat = case_when(name %in% c("comp_1", "comp_2") ~ "PC 1, 2",
                         name %in% c("comp_4", "comp_5") ~ "PC 4, 5",
                         name %in% c("comp_6", "comp_7") ~ "PC 6, 7",
                         T ~ "PC 3"))%>% 
  ggplot(aes(x = band, y = value))+
  geom_line(aes(group = grp, color= img_date, alpha = sqrt(sqrt(proportion))))+
  facet_wrap(~ cat)+
  scale_color_discrete(name = "Image Date")+
  scale_alpha(name = "Prop. Var \nExplained", range = c(0, 1))+
  labs(x = "Wavelength", y = "Loading", title = "PRISMA: all PCs at/over 1% variance")+
  theme_bw()
