#from LS_buoy, source ar_merge

chl_dat = ar_merge %>% 
  select(-c(avg_pco2_ppm:avg_spec_cond, avg_par, avg_par_below)) %>% 
  na.omit()

rownames(chl_dat) = chl_dat$sampledate

datelabels = chl_dat$sampledate

buoy.dat = chl_dat[,2:10]
ls.dat = chl_dat[,11:46]
date.dat = chl_dat[,1]

library(vegan)
edaPlot = rda(buoy.dat ~ ., data = ls.dat)

plot(edaPlot, choices = c(1,2), type = "text")

# RDA of full model, gives us the fractions of [a+b+c] 
rda.all <- rda(buoy.dat ~ Thermal + Aerosol.NIR.r, data=ls.dat) 
# fractions [a+b] 
rda.MAT <- rda(buoy.dat ~ Thermal, data = ls.dat) 
#fractions [b+c] 
rda.MCMT <- rda(buoy.dat ~ Aerosol.NIR.r, data = ls.dat)
plot(rda.all, choices=c(1,2), type="text") #check the plot

#fractions [a+b+c] 
RsquareAdj(rda.all)

abc <- RsquareAdj(rda.all)$adj.r.squared # Extract the adjusted r-squared for the full model
#fraction [a+b] 
ab <- RsquareAdj(rda.MAT)$adj.r.squared
#fraction [b+c] 
bc <- RsquareAdj(rda.MCMT)$adj.r.squared
#individual fractions 
b <- ab + bc - abc 
a <- ab - b 
c <- bc - b
out <- varpart(buoy.dat, ~ Thermal, ~ Aerosol.NIR.r, data=ls.dat)
#Let's see if our calculations match those that `varpart` calculated
out$part$fract #Fractions (ab, bc, abc)
out$part$indfract

plot(out, bg = c("hotpink","skyblue"), Xnames = c("Thermal", "Aerosol.NIR.r"))
