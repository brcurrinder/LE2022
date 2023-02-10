

# Function for P flux between sediment and water column

# DiffFlux (g/m3/d) = (C1 / (ZSed(m) * ZHyp(m))) * (PSed(g/m3) - PHyp(g/m3))

nDays = 365 * 100
Area = 19805000 #10000
C1 = 1*10^-7 #5*10^-7 # diffusivity in m2/d
ZSed = 0.1 # m
ZHyp = 4 # m
PSedInit = 120 # gP/m2 of sediment area
PHypInit = 0.1 # gP/m3 in hypolimnion, e.g., 0.06 gP/m3 = 0.06 mgP/L = 60 ugP/L
ExportFrac = 0.25 # Fraction of hypolimnetic export each year (due to, e.g., WRT)

PSed = data.frame(PSed = rep(NA,nDays))
PSed[1] = PSedInit * Area # 140 (g/m2) * Area(m) = gP
PHyp = data.frame(PHyp = rep(NA,nDays))
PHyp[1] = PHypInit * Area * ZHyp # mg/m3 * m3 = gP
DiffFlux = data.frame(Flux=rep(NA,nDays))
DiffFlux[1] = 0

for (i in 2:nDays){
  # all fluxes calculated in gP/d
  DiffFlux$Flux[i] = (C1 / (ZSed * ZHyp)) * (PSed$PSed[i-1]/(ZSed * Area) - PHyp$PHyp[i-1]/(ZHyp * Area)) # g/m3/d
  Flux_g = DiffFlux$Flux[i] * ZHyp * Area # Flux in g/d
  PSed$PSed[i] = PSed$PSed[i-1] - Flux_g # in g
  PHyp$PHyp[i] = PHyp$PHyp[i-1] + Flux_g # in g
  PHyp$PHyp[i] = PHyp$PHyp[i-1] + Flux_g - ExportFrac/365 * PHyp$PHyp[i-1]
}

plot(2:nDays/365,DiffFlux$Flux[2:nDays],type='l',xlab='Year',ylab='Flux (mgP/L/d)')
plot(2:nDays/365,DiffFlux$Flux[2:nDays]/ZHyp*1000,type='l',xlab='Year',ylab='Flux (mgP/m2/d)')
plot(1:nDays/365,PHyp$PHyp/(ZHyp*Area),type='l',xlab='Year',ylab='Hypo P (mgP/L)')
plot(1:nDays/365,PSed$PSed/Area,type='l',xlab='Year',ylab='Sed P (g/m2)')


