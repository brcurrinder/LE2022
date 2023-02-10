PlotMetabResultsHypoBudget = function(MR,myLake,SimSetup,LakeHypsometry,TargetYear){
  
  # Lake area
  LakeArea = max(LakeHypsometry$area)
  
  # Extract data frames from metabolism model
  StartDates = SimSetup$StartDate
  EndDates  = SimSetup$EndDate
  StatesEpiC = MR[[1]]
  StatesHypoC = MR[[2]]
  StatesEpiGas = MR[[3]]
  StatesHypoGas = MR[[4]]
  RatesEpiC = MR[[5]]
  RatesHypoC = MR[[6]]
  RatesEpiGas = MR[[7]]
  RatesHypoGas = MR[[8]]
  StatesEpiP = MR[[9]]
  StatesHypoP = MR[[10]]
  StatesSedP = MR[[11]]
  RatesEpiP = MR[[12]]
  RatesHypoP = MR[[13]]
  RatesSedP = MR[[14]]
  StatesSedC = MR[[15]]
  RatesSedC = MR[[16]]
  LakeVolumes = MR[[17]]
  
  LakeAreas = data.frame(TotalArea = max(LakeHypsometry))
  
  
  
  # get simulation date range
  # SimDateRange = c(min(as.Date(ModelRun[[1]]$StartDate)),max(as.Date(ModelRun[[1]]$EndDate)))
  SimDateRange = c(as.Date(StartDates),as.Date(EndDates))
  # generate sequence of dates for, e.g., plotting
  SimDates = seq(as.Date(SimDateRange[1]),as.Date(SimDateRange[2]),by='day')
  # remove first value, as output begins on day 2
  # If more than one simulation is strung together, a day gets dropped for each sim
  nTotDates = length(SimDates)-1
  SimDates = SimDates[2:nTotDates]
  
  
  ########################################
  # Plot results
  YearFrac = round(1:nDays/365,digits=2)
  iStratified = which(LakeVolumes$HypoVol>0)
  # Determine onset/offset of stratification
  CurrentlyStratified = FALSE
  S1=NULL
  S2=NULL
  for (iS in 1:length(LakeVolumes$HypoVol)){
    if (LakeVolumes$HypoVol[iS]>0 & !CurrentlyStratified){
      CurrentlyStratified = TRUE
      S1 = c(S1,iS)
    }
    if (LakeVolumes$HypoVol[iS]==0 & CurrentlyStratified){
      CurrentlyStratified = FALSE
      S2 = c(S2,iS)
    }
  }
  # Normalize to fraction of year
  S1 = S1/365
  S2 = S2/365
  
  # Calc strat ablines
  Years = unique(floor(YearFrac))
  
  # Get subrange based on target year
  TYBegin = paste(TargetYear,'-01-01',sep="")
  TYEnd = paste(TargetYear,'-12-31',sep="")
  iThisYear = which(as.POSIXct(SimDates)>=as.POSIXct(TYBegin) & as.POSIXct(SimDates)<=as.POSIXct(TYEnd))

  # Get the budget parts in C units
  RespirationSed = (RatesSedC$POCLR + RatesSedC$POCRR)/ LakeArea
  SedO2Demand = (RespirationSed /12 *1000 /365)/myLake$SedAreaProp # mmol O2,C /m2SA/day
  RespirationHypoDOC = (RatesHypoC$DOCRR + RatesHypoC$DOCLR)/LakeArea
  RespirationHypoPOC = (RatesHypoC$POCRR + RatesHypoC$POCLR)/LakeArea
  RespirationTotal = RespirationSed+RespirationHypoDOC+RespirationHypoPOC
  
  # Mass balance equation from Metabolism.R
  # StatesHypoGas$DO[i] = StatesHypoGas$DO[i-1] + RatesEpiGas$DOMix[i] - 
  #          RatesHypoC$DOCLR[i]*O2toCMR - RatesHypoC$DOCRR[i]*O2toCMR - 
  #          RatesHypoC$POCLR[i]*O2toCMR - RatesHypoC$POCRR[i]*O2toCMR - 
  #          RatesHypoGas$DOSed[i]
  
  DOMix = RatesEpiGas$DOMix/LakeArea
  RDOWQ = (RatesHypoC$DOCLR + RatesHypoC$DOCRR + 
    RatesHypoC$POCLR + RatesHypoC$POCRR) * 32/12 / LakeArea
  DORSed = RatesHypoGas$DOSed / LakeArea
  HypoDO = StatesHypoGas$DO/LakeVolumes$HypoVol
  
  par(mfrow=c(2,1),mai = c(0.15,0.75, 0.05, 0.05),oma = c(2,2,0.2,0.2), cex = 0.9)
  
  # Plot DO
  myYlim = c(-0.1,1)
  plot(as.POSIXct(SimDates[iThisYear]),HypoDO[iThisYear]/10,ylim=myYlim,type='l',
       ylab = 'O2 and fluxes(g/m3,gO2/m2/d)')
  lines(as.POSIXct(SimDates[iThisYear]),DOMix[iThisYear]/20,col='red')
  lines(as.POSIXct(SimDates[iThisYear]),RDOWQ[iThisYear],col='green')
  lines(as.POSIXct(SimDates[iThisYear]),DORSed[iThisYear],col='blue')
  legend("topleft",legend=c('O2/10','Mix/20','WQ','Sed'),lty=c(1,1,1,1,1),
         col=c('black','red','green','blue'))
  grid()
  
  # Plot Carbon
  plot(as.POSIXct(SimDates[iThisYear]),RespirationTotal[iThisYear],ylim=myYlim,type='l',
       ylab = 'C Fluxes (gC/m2/d)')
  lines(as.POSIXct(SimDates[iThisYear]),RespirationSed[iThisYear],col='blue')
  lines(as.POSIXct(SimDates[iThisYear]),RespirationHypoPOC[iThisYear],col='green')
  lines(as.POSIXct(SimDates[iThisYear]),RespirationHypoDOC[iThisYear],col='brown')
  legend("topleft",legend=c('RTot','RSed','RPOC','RDOC'),lty=c(1,1,1,1,1),col=c('black','blue','green','brown'))
  grid()
  
  # Plot cumulative entrainment
  plot(as.POSIXct(SimDates[iThisYear]),cumsum(DOMix[iThisYear]),ylim=c(0,120), type='l',
                                              ylab = 'Cumulative entrainment (mgDO/m2)')
  grid()
  
  myCumSumEntrO2 = data.frame(Date=SimDates[iThisYear],EntrCumSumO2=cumsum(DOMix[iThisYear]))
  myRTotO2 = data.frame(Date=SimDates[iThisYear],RTotO2=RespirationTotal[iThisYear]*32/12)
  
  print(myCumSumEntrO2)
  print(myRTotO2)
 
}

StratifiedShade <- function(S1,S2,myYLim){
  for (iMix in 1:length(S1)){
    #polygon(c(1,2,2,1),c(99.2,99.2,100,100),col=rgb(0, 0, 0,0.2))
    polygon(c(S1[iMix],S2[iMix],S2[iMix],S1[iMix]),c(min(myYLim),min(myYLim),max(myYLim),max(myYLim)),
            col=rgb(0, 0, 0,0.1),border=NA)
    grid()
  }
}