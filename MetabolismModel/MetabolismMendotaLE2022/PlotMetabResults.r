PlotMetabResults = function(MR,myLake,LakeHypsometry){
  
  # Lake area
  LakeArea = max(LakeHypsometry$area)
  
  # Extract data frames from metabolism model
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
  
  # S1 = (Years) + myLake$StratBegin/365
  # S2 = (Years) + myLake$StratEnd/365
  # S1 = (Years) + min(iStratified)/365
  # S2 = (Years) + max(iStratified)/365

  ######################
  # Plot OC rates
  # Calculate a bunch of stuff
  RSettling = (RatesSedC$POCRSettling)*365 / LakeArea
  LSettling = (RatesSedC$POCLSettling)*365 / LakeArea
  Settling = RSettling + LSettling
  PropRSettling = round(mean(RSettling/Settling,na.rm=TRUE),digits=2)
  PropLSettling = round(mean(LSettling/Settling,na.rm=TRUE),digits=2)
  RespirationSed = (RatesSedC$POCLR + RatesSedC$POCRR)*365 / LakeArea
  Burial = (RatesSedC$POCLBurial + RatesSedC$POCRBurial)*365 / LakeArea
  Net = Settling - RespirationSed - Burial
  SedimentCarbon = StatesSedC$POCR/LakeAreas$TotalArea+StatesSedC$POCL/LakeAreas$TotalArea # gC/m2LA
  SedO2Demand = (RespirationSed /12 *1000 /365)/myLake$SedAreaProp # mmol O2,C /m2SA/day
  RespirationHypoDOC = (RatesHypoC$DOCRR + RatesHypoC$DOCLR)*365/LakeArea
  RespirationHypoPOC = (RatesHypoC$POCRR + RatesHypoC$POCLR)*365/LakeArea
  
  # Hypo
  par(mfrow=c(3,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myYLim = c(min(c(RespirationHypoDOC,RespirationHypoPOC,RespirationSed),na.rm=TRUE),
             max(c(RespirationHypoDOC,RespirationHypoPOC,RespirationSed),na.rm=TRUE))
  plot(YearFrac,RespirationSed,ylim=myYLim,col='blue',type='l',ylab="Hypo R (gC/m2LA/y)",main=myLake$LakeName)
  lines(YearFrac,RespirationHypoDOC,col='red')
  lines(YearFrac,RespirationHypoPOC,col='green')
  StratifiedShade(S1,S2,myYLim)
  legend('topleft',legend=c('Rsed','Rdoc','Rpoc'),lty=c(1,1,1),col=c('blue','red','green'))
  
  # Sed
  myYLim = c(min(c(Settling,RespirationSed,Burial,Net),na.rm=TRUE),
             max(c(Settling,RespirationSed,Burial,Net),na.rm=TRUE))
  plot(YearFrac,Settling,ylim=myYLim,col='green',type='l',ylab='OCTot-Sed (gC/m2LA/y)',main="")
  lines(YearFrac,RespirationSed,col='blue')
  lines(YearFrac,Burial,col='brown')
  lines(YearFrac,Net,col='grey',lty=2,lwd=2)
  abline(h=0,col='black',lty=1)
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('Settling','Resp','Burial','Net'),'topleft',lty=c(1,1,1,2),col=c('green','blue','brown','grey'),cex=0.7)
  myLegTxt = paste('Proportion settling POCR: ',PropRSettling,', POCL: ',PropLSettling)
  legend(legend=myLegTxt,'topright')
  # Plot sediment oxygen demand
  par(mai = c(0.4,0.8, 0, 0.1),cex = 0.9)
  myYLim = c(0,max(SedO2Demand,na.rm=TRUE))
  plot(YearFrac,SedO2Demand,ylim=myYLim,type='l',ylab='SedO2demand (mmolO2/m2SA/day)')
  StratifiedShade(S1,S2,myYLim)

  # Plot OC and P retention
  # OCExport = RatesEpiC$DOCLOutflow + RatesEpiC$DOCROutflow +
  #               RatesEpiC$POCLOutflow + RatesEpiC$POCROutflow 
  # Allocthony = RatesEpiC$DOCLInflow + RatesEpiC$DOCRInflow +
  #                 RatesEpiC$POCLInflow + RatesEpiC$POCRInflow
  # Autochthony = RatesEpiC$GPP
  
  # OCRetention = 1 - (OCExport / (Allocthony+Autochthony))
  # PRetention = 1 - (RatesEpiP$TPOutflow / RatesEpiP$TPInflow) # Daily values
  # myYLim = c(min(OCRetention,PRetention,na.rm=TRUE),max(OCRetention,PRetention,na.rm=TRUE))
  # if (is.finite(myYLim[1]) & is.finite(myYLim[2])){
  #   plot(YearFrac,PRetention,ylim=myYLim,type='l',col='blue',ylab='Retention (Export/Load')
  #   lines(YearFrac,OCRetention,col='brown')
  #   StratifiedShade(S1,S2,myYLim)
  #   legend('topleft',c('P','OC_tot'),lty=c(1,1),col=c('blue','brown'))
  # }

  ######################
  #DOCTot
  # Calc y lims
  YLower = min((StatesEpiC$DOCR+StatesEpiC$POCR+StatesEpiC$DOCL+StatesEpiC$POCL)/LakeVolumes$EpiVol) - 2
  YUpper = max((StatesEpiC$DOCR+StatesEpiC$POCR+StatesEpiC$DOCL+StatesEpiC$POCL)/LakeVolumes$EpiVol) + 2
  myYLim = c(YLower,YUpper)
  par(mfrow=c(4,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  # Epi
  plot(YearFrac,(StatesEpiC$DOCR)/LakeVolumes$EpiVol,ylim=myYLim,type='l',ylab='DOCTot-epi (gC/m3)',main=myLake$LakeName)
  lines(YearFrac,(StatesEpiC$DOCR+StatesEpiC$POCR)/LakeVolumes$EpiVol,col='grey')
  lines(YearFrac,(StatesEpiC$DOCR+StatesEpiC$POCR+StatesEpiC$DOCL)/LakeVolumes$EpiVol,col='brown')
  lines(YearFrac,(StatesEpiC$DOCR+StatesEpiC$POCR+StatesEpiC$DOCL+StatesEpiC$POCL)/LakeVolumes$EpiVol,col='green')
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('DOCR','POCR','DOCL','POCL'),'topleft',lty=c(1,1,1,1),col=c('black','grey','brown','green'),cex=0.7)
  # Hypo
  YLower = 0
  YUpper = max((StatesHypoC$DOCR[iStratified]+StatesHypoC$POCR[iStratified]+StatesHypoC$DOCL[iStratified]+StatesHypoC$POCL[iStratified]) /
                 LakeVolumes$HypoVol[iStratified]) + 2
  myYLim = c(YLower,YUpper)
  plot(YearFrac,(StatesHypoC$DOCR)/LakeVolumes$HypoVol,ylim=myYLim,type='l',ylab='DOCTot-Hypo (gC/m3)')
  lines(YearFrac,(StatesHypoC$DOCR+StatesHypoC$POCR)/LakeVolumes$HypoVol,col='grey')
  lines(YearFrac,(StatesHypoC$DOCR+StatesHypoC$POCR+StatesHypoC$DOCL)/LakeVolumes$HypoVol,col='brown')
  lines(YearFrac,(StatesHypoC$DOCR+StatesHypoC$POCR+StatesHypoC$DOCL+StatesHypoC$POCL)/LakeVolumes$HypoVol,col='green')
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('DOCR','POCR','DOCL','POCL'),'topleft',lty=c(1,1,1,1),col=c('black','grey','brown','green'),cex=0.7)
  # Sed
  myYLim = c(min(StatesSedC$POCR/(LakeAreas$TotalArea*myLake$SedAreaProp)),
          max(max(StatesSedC$POCR/(LakeAreas$TotalArea*myLake$SedAreaProp)+StatesSedC$POCL/(LakeAreas$TotalArea*myLake$SedAreaProp))))
  plot(YearFrac,(StatesSedC$POCR)/(LakeAreas$TotalArea*myLake$SedAreaProp),ylim=myYLim,col='brown',type='l',ylab='DOCTot-Sed (gC/m2SA)')
  lines(YearFrac,(StatesSedC$POCR+StatesSedC$POCL)/(LakeAreas$TotalArea*myLake$SedAreaProp),col='green')
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('POCR','POCL'),'topleft',lty=c(1,1),col=c('brown','green'),cex=0.7)
  # Secchi
  myYLim = c(min(StatesEpiC$Secchi),max(StatesEpiC$Secchi))
  plot(YearFrac,(StatesEpiC$Secchi),ylim=myYLim,type='l',ylab='Secchi (m)')
  StratifiedShade(S1,S2,myYLim)
  
  ######################
  #DO
  # Calc y lims
  YLower = 0
  YUpper = 15
  myYLim = c(YLower,YUpper)
  
  par(mfrow=c(3,1), mai = c(0.3,0.8, 0.3, 0.1),cex = 0.9)
  plot(YearFrac,(StatesEpiGas$DO)/LakeVolumes$EpiVol,ylim=myYLim,type='l',xlab='',ylab='DO-epi (gO2/m3)',main=myLake$LakeName)
  lines(YearFrac,(StatesEpiGas$DOsat),lty=2,col='grey')
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('DO','DOsat'),'topleft',lty=c(1,2),col=c('black','grey','brown','green'),cex=0.7)
  plot(YearFrac,(StatesHypoGas$DO)/LakeVolumes$HypoVol,ylim=myYLim,type='l',xlab='',ylab='DO-Hypo (gO2/m3)')
  lines(YearFrac,(StatesHypoGas$DOsat),lty=2,col='grey')
  StratifiedShade(S1,S2,myYLim)
  myYLim = range(RatesEpiGas$DOFatm/LakeAreas$TotalArea,na.rm=TRUE)
  plot(YearFrac,c(0,RatesEpiGas$DOFatm[2:length(RatesEpiGas$DOFatm)]/LakeAreas$TotalArea),type='l',
                  ylim=myYLim,xlab='',ylab='Cum DO Flux (gO2/m2LA;+=in)')
  abline(h=0,lty=2)
  StratifiedShade(S1,S2,myYLim)
  

  # Phosphorus
  # Calc y lims
  myYLim = c(0,max(RatesEpiP$TPSettling/LakeAreas$TotalArea*365,
                   RatesEpiP$TPInflow/LakeAreas$TotalArea*365,
                   #RatesEpiP$TPRebind/LakeAreas$TotalArea*365,
                   RatesEpiP$TPOutflow/LakeAreas$TotalArea*365,na.rm=TRUE))
  
  ######################
  # Plot the rates
  par(mfrow=c(3,1), mai = c(0.3,0.8, 0.3, 0.1),cex = 0.9)
  # Epi
  plot(YearFrac,(RatesEpiP$TPSettling)/LakeAreas$TotalArea*365,ylim=myYLim,type='l',xlab='',ylab='P-epi rates (gP/m2LA/y)',main=myLake$LakeName)
  lines(YearFrac,RatesEpiP$TPRebind/LakeAreas$TotalArea*365,lty=2,col='black')
  lines(YearFrac,RatesEpiP$TPRecycling/LakeAreas$TotalArea*365,lty=1,col='blue')
  lines(YearFrac,RatesEpiP$TPRelease/LakeAreas$TotalArea*365,lty=2,col='blue')
  lines(YearFrac,(RatesEpiP$TPInflow)/LakeAreas$TotalArea*365,col='green')
  lines(YearFrac,(RatesEpiP$TPOutflow)/LakeAreas$TotalArea*365,col='red')
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('Settling','Rebind','Recycling','Release','Inflow','Outflow'),'topleft',lty=c(1,2,1,2,1,1),col=c('black','black','blue','blue','green','red'),cex=0.7)
  # Hypo
  plot(YearFrac,(RatesHypoP$TPRecycling)/LakeAreas$TotalArea*365,ylim=myYLim,type='l',col='blue',xlab='',ylab='P-hypo rates (gP/m2LA/y)',main="")
  lines(YearFrac,RatesHypoP$TPRelease/LakeAreas$TotalArea*365,lty=2,col='blue')
  lines(YearFrac,(RatesHypoP$TPSettling)/LakeAreas$TotalArea*365,col='black')
  lines(YearFrac,RatesHypoP$TPRebind/LakeAreas$TotalArea*365,lty=2,col='black')
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('Recycling','Release','Settling','Rebinding'),'topleft',lty=c(1,2,1,2),col=c('blue','blue','black','black'),cex=0.7)
  # Sed
  plot(YearFrac,(RatesSedP$TPSettling)/LakeAreas$TotalArea*365,ylim=myYLim,type='l',xlab='',ylab='P-sed rates (gP/m2LA/y)',main="")
  lines(YearFrac,(RatesSedP$TPRebind)/LakeAreas$TotalArea*365,lty=2,col='black')
  lines(YearFrac,(RatesSedP$TPRecycling)/LakeAreas$TotalArea*365,col='blue')
  lines(YearFrac,RatesSedP$TPRelease/LakeAreas$TotalArea*365,lty=2,col='blue')
  
  lines(YearFrac,(RatesSedP$TPBurial)/LakeAreas$TotalArea*365,col='red')
  lines(YearFrac,(RatesSedP$TPBBurial)/LakeAreas$TotalArea*365,lty=2,col='red')
  legend(legend=c('Settling','Rebinding','Recycling','Release','Burial','BurialBound'),'topleft',
         lty=c(1,2,1,2,1,2),col=c('black','black','blue','blue','red','red'),cex=0.7)
  StratifiedShade(S1,S2,myYLim)
  
  ######################
  # Plot the states
  par(mfrow=c(3,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  # Epi
  plot(YearFrac,(StatesEpiP$TP)/LakeVolumes$EpiVol,type='l',xlab='',ylab='P-epi (gP/m3)',main=myLake$LakeName)
  StratifiedShade(S1,S2,myYLim)
  # Hypo
  plot(YearFrac,(StatesHypoP$TP)/LakeVolumes$HypoVol,type='l',xlab='',ylab='P-hypo (gP/m3)')
  StratifiedShade(S1,S2,myYLim)

  # Sed
  myRange = range((StatesSedP$TP+StatesSedP$TPB)/(LakeAreas$TotalArea*myLake$SedAreaProp))
  myDiff = diff(myRange)
  myYLim = c(myRange[1]-0.5*myDiff,myRange[2]+0.5*myDiff)
  plot(YearFrac,(StatesSedP$TP+StatesSedP$TPB)/(LakeAreas$TotalArea*myLake$SedAreaProp),ylim=myYLim,type='l',xlab='',ylab='P-sed (gP/m2)')
  grid(col='grey')
  StratifiedShade(S1,S2,myYLim)
  
}

StratifiedShade <- function(S1,S2,myYLim){
  for (iMix in 1:length(S1)){
    #polygon(c(1,2,2,1),c(99.2,99.2,100,100),col=rgb(0, 0, 0,0.2))
    polygon(c(S1[iMix],S2[iMix],S2[iMix],S1[iMix]),c(min(myYLim),min(myYLim),max(myYLim),max(myYLim)),
            col=rgb(0, 0, 0,0.1),border=NA)
    grid()
  }
}