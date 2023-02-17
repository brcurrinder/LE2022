library(vegan)
#from LS_buoy, source ar_merge

chl_dat = ar_merge %>% 
  select(-c(avg_pco2_ppm:avg_spec_cond, avg_par, avg_par_below)) %>% 
  na.omit() %>% 
  left_join(chl_LS2) %>% 
  na.omit()

rownames(chl_dat) = chl_dat$sampledate

datelabels = chl_dat$sampledate

buoy.dat = chl_dat[,c(47,2:10)]
bands.dat = chl_dat[,11:18]
ratios.dat = chl_dat[,19:46]
algs.dat = chl_dat[,48:56]


edaPlot = rda(buoy.dat ~ ., data = algs.dat)
#blue.swir, aero.thermal, thermal, OC4, NDVI
plot(edaPlot, choices = c(1,2), type = "text")

# RDA of full model

out <- varpart(buoy.dat, ~Blue.SWIR1.r, ~Aerosol.Thermal.r, ~Thermal, ~NDVI, data = chl_dat)
#Let's see if our calculations match those that `varpart` calculated
out$part$fract #Fractions (ab, bc, abc)
out$part$indfract

plot(out, bg = c("hotpink","skyblue", "limegreen", "orange"), Xnames = c("Blue.SWIR1.r", "Aerosol.Thermal.r", "OC4", "NDVI"))

chl_dat %>% 
  filter(OC4 < 10^6) %>% 
ggplot() +
  geom_point(aes(x = OC4, y = avg_chlor_rfu))
