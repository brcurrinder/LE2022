# Scrit to test various NPP functions
# Iz = I0 e - kz, e.g., eutrophic = 1200 * exp(-1*1*4) = 22; oligotrophic = 1200 * exp(-1*0.25*4) = 442

TPconc = c(0.01,0.060) # g/m3 (mg/L)
IP = c(0.02,0.1) # 0.1 for eutrophic, 0.02 fo oligotrophic

# Calibration points:
# Oligotrophic: 
# Eutrophic: Peak NPP should be about 0.4 gC/m3 for TP 50 ug/L at I=50

CalPoints_I = c(100, 50) # Maximum epi light oligotrophic, eutrophic
CalPoints_NPP = c(0.04, 0.4) # Oligotrophic NPP, eutrophic NPP

I = 0:100 # w/m2/s

NPP_1 = data.frame(NPP=rep(NA,length(I)))

# Cycle through the light levels and calculate NPP
# Assume 20C (Arrhenius = 1)
for (i in 1:length(I)){
  # First model
  NPP_1$NPP[i] = TPconc[2] * IP[2] * I[i]
}

plot(I,NPP_1$NPP,type='l',xlab="I (mmol/m2/s)",ylab="NPP (gC/m3/d)")
points(CalPoints_I,CalPoints_NPP,col='red')
grid()
