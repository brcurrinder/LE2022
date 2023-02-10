PlotMetabResults = function(MetabResults,myLake,LakeVolumes){
  
  # Extract data frames from metabolism model
  StatesEpiC = MetabResults[[1]]
  StatesHypoC = MetabResults[[2]]
  StatesEpiGas = MetabResults[[3]]
  StatesHypoGas = MetabResults[[4]]
  RatesEpiC = MetabResults[[5]]
  RatesHypoC = MetabResults[[6]]
  RatesEpiGas = MetabResults[[7]]
  RatesHypoGas = MetabResults[[8]]
  StatesEpiP = MetabResults[[9]]
  StatesHypoP = MetabResults[[10]]
  StatesSedP = MetabResults[[11]]
  RatesEpiP = MetabResults[[12]]
  RatesHypoP = MetabResults[[13]]
  RatesSedP = MetabResults[[14]]
  StatesSedC = MetabResults[[15]]
  RatesSedC = MetabResults[[16]]
  
  ########################################
  # Plot results
  YearFrac = round(1:nDays/365,digits=2)
  iStratified = which(LakeVolumes$HypoVol>0)
  
  # Calc strat ablines
  Years = unique(floor(YearFrac))
  
  S1 = (Years) + myLake$StratBegin/365
  S2 = (Years) + myLake$StratEnd/365

  # Plot rates
  # Labile epi
  myYLim = c(min(RatesEpiC$GPP/LakeVolumes$EpiVol,na.rm=TRUE),max(RatesEpiC$GPP/LakeVolumes$EpiVol,na.rm=TRUE))
  par(mfrow=c(3,1), mai = c(0.1,0.8, 0.3, 0.1),cex = 0.9)
  plot(YearFrac,RatesEpiC$GPP/LakeVolumes$EpiVol,ylim=myYLim,type='l',ylab='GPP (gC/m3/d)',main=myLake$LakeName)
  StratifiedShade(S1,S2,myYLim)
  plot(YearFrac,RatesEpiC$DOCLR/LakeVolumes$EpiVol,ylim=myYLim,,type='l',ylab='DOCLR (gC/m3/d)')
  StratifiedShade(S1,S2,myYLim)
  plot(YearFrac,RatesEpiC$DOCRInflow/LakeVolumes$EpiVol,ylim=myYLim,,type='l',ylab='Inflow (gC/m3/d)')
  StratifiedShade(S1,S2,myYLim)
  
  # Recalcitrant epi
  myYLim = c(min(RatesEpiC$DOCLR/LakeVolumes$EpiVol,na.rm=TRUE),max(RatesEpiC$DOCLR/LakeVolumes$EpiVol,na.rm=TRUE))
  plot(YearFrac,RatesEpiC$DOCLR/LakeVolumes$EpiVol,type='l',ylab='DOCLR (gC/m3/d)',main=myLake$LakeName)
  StratifiedShade(S1,S2,myYLim)
  plot(YearFrac,RatesEpiC$DOCRInflow/LakeVolumes$EpiVol,type='l',ylab='DOCR Inflow (gC/m3/d)')
  StratifiedShade(S1,S2,myYLim)
  plot(YearFrac,RatesEpiC$DOCROutflow/LakeVolumes$EpiVol,type='l',ylab='DOCR Outflow (gC/m3/d)')
  StratifiedShade(S1,S2,myYLim)
  
  # Plot states
  
  # POCL
  # Calc y lims
  YLower = 0
  YUpper = max(StatesEpiC$POCL/LakeVolumes$EpiVol) + 1
  myYLim = c(YLower,YUpper)
  par(mfrow=c(3,1), mai = c(0.1,0.8, 0.3, 0.1),cex = 0.9)
  plot(YearFrac,StatesEpiC$POCL/LakeVolumes$EpiVol,ylim=myYLim,type='l',ylab='POCL-epi (gC/m3)',main=myLake$LakeName)
  StratifiedShade(S1,S2,myYLim)
  plot(YearFrac,StatesHypoC$POCL/LakeVolumes$HypoVol,ylim=myYLim,type='l',ylab='POCL-hyp (gC/m3)')
  StratifiedShade(S1,S2,myYLim)
  
  # POCR
  # Calc y lims
  YLower = 0
  YUpper = max(StatesEpiC$POCR/LakeVolumes$EpiVol) + 1
  myYLim = c(YLower,YUpper)
  par(mfrow=c(3,1), mai = c(0.1,0.8, 0.3, 0.1),cex = 0.9)
  plot(YearFrac,StatesEpiC$POCR/LakeVolumes$EpiVol,ylim=myYLim,type='l',ylab='POCR-epi (gC/m3)',main=myLake$LakeName)
  StratifiedShade(S1,S2,myYLim)
  plot(YearFrac,StatesHypoC$POCR/LakeVolumes$HypoVol,ylim=myYLim,type='l',ylab='POCR-hyp (gC/m3)')
  StratifiedShade(S1,S2,myYLim)
  
  # DOCL
  # Calc y lims
  YLower = 0
  YUpper = max(StatesEpiC$DOCL/LakeVolumes$EpiVol) + 1
  myYLim = c(YLower,YUpper)
  par(mfrow=c(3,1), mai = c(0.1,0.8, 0.3, 0.1),cex = 0.9)
  plot(YearFrac,StatesEpiC$DOCL/LakeVolumes$EpiVol,ylim=myYLim,type='l',ylab='DOCL-epi (gC/m3)',main=myLake$LakeName)
  StratifiedShade(S1,S2,myYLim)
  plot(YearFrac,StatesHypoC$DOCL/LakeVolumes$HypoVol,ylim=myYLim,type='l',ylab='DOCL-hyp (gC/m3)')
  StratifiedShade(S1,S2,myYLim)
  
  #DOCR
  # Calc y lims
  YLower = min(StatesEpiC$DOCR/LakeVolumes$EpiVol) - 1
  YUpper = max(StatesEpiC$DOCR/LakeVolumes$EpiVol) + 1
  myYLim = c(YLower,YUpper)
  par(mfrow=c(3,1), mai = c(0.1,0.8, 0.3, 0.1),cex = 0.9)
  plot(YearFrac,StatesEpiC$DOCR/LakeVolumes$EpiVol,ylim=myYLim,type='l',ylab='DOCR-epi (gC/m3)',main=myLake$LakeName)
  StratifiedShade(S1,S2,myYLim)
  plot(YearFrac,StatesHypoC$DOCR/LakeVolumes$HypoVol,ylim=myYLim,type='l',ylab='DOCR-hyp (gC/m3)')
  StratifiedShade(S1,S2,myYLim)
  
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
  legend(legend=c('POCR','POCL'),'topleft',lty=c(1,1),col=c('brown','green'),cex=0.7)
  # Sed
  myYLim = c(min(StatesSedC$POCR/LakeAreas$TotalArea),
          max(max(StatesSedC$POCR/LakeAreas$TotalArea+StatesSedC$POCL/LakeAreas$TotalArea)))
  plot(YearFrac,(StatesSedC$POCR)/LakeAreas$TotalArea,ylim=myYLim,col='brown',type='l',ylab='DOCTot-Sed (gC/m2LA)')
  lines(YearFrac,(StatesSedC$POCR+StatesSedC$POCL)/LakeAreas$TotalArea,col='green')
  StratifiedShade(S1,S2,myYLim)
  # Secchi
  myYLim = c(min(StatesEpiC$Secchi),max(StatesEpiC$Secchi))
  plot(YearFrac,(StatesEpiC$Secchi),ylim=myYLim,type='l',ylab='Secchi (m)')
  StratifiedShade(S1,S2,myYLim)
  
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
  plot(YearFrac,(StatesHypoGas$DO)/LakeVolumes$HypoVol,ylim=myYLim,type='l',xlab='',ylab='DO-Hypo (gO2/m3)')
  lines(YearFrac,(StatesHypoGas$DOsat),lty=2,col='grey')
  StratifiedShade(S1,S2,myYLim)
  

  # Phosphorus
  # Calc y lims
  myYLim = c(0,max(RatesEpiP$TPSettling/LakeAreas$TotalArea*365,
                   RatesEpiP$TPInflow/LakeAreas$TotalArea*365,
                   RatesEpiP$TPOutflow/LakeAreas$TotalArea*365,na.rm=TRUE))
  # Plot the rates
  par(mfrow=c(3,1), mai = c(0.3,0.8, 0.3, 0.1),cex = 0.9)
  # Epi
  plot(YearFrac,(RatesEpiP$TPSettling)/LakeAreas$TotalArea*365,ylim=myYLim,type='l',xlab='',ylab='P-epi rates (gP/m2/y)',main=myLake$LakeName)
  lines(YearFrac,(RatesEpiP$TPInflow)/LakeAreas$TotalArea*365,col='green')
  lines(YearFrac,(RatesEpiP$TPOutflow)/LakeAreas$TotalArea*365,col='red')
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('Settling','Inflow','Outflow'),'topleft',lty=c(1,1,1),col=c('black','green','red'),cex=0.7)
  # Hypo
  plot(YearFrac,(RatesHypoP$TPRecycling)/LakeAreas$TotalArea*365,type='l',xlab='',ylab='P-hypo rates (gP/m2/y)',main="")
  lines(YearFrac,(RatesHypoP$TPSettling)/LakeAreas$TotalArea*365,col='green')
  StratifiedShade(S1,S2,myYLim)
  legend(legend=c('Recycling','Settling'),'topleft',lty=c(1,1),col=c('black','green'),cex=0.7)
  # Sed
  myYLim = c(0,max(RatesSedP$TPSettling/LakeAreas$TotalArea*365,
                   RatesSedP$TPRecycling/LakeAreas$TotalArea*365,
                   RatesSedP$TPBurial/LakeAreas$TotalArea*365,na.rm=TRUE))
  plot(YearFrac,(RatesSedP$TPSettling)/LakeAreas$TotalArea*365,ylim=myYLim,type='l',xlab='',ylab='P-sed rates (gP/m2/y)',main="")
  lines(YearFrac,(RatesSedP$TPRecycling)/LakeAreas$TotalArea*365,col='green')
  lines(YearFrac,(RatesSedP$TPBurial)/LakeAreas$TotalArea*365,col='red')
  legend(legend=c('Settling','Recycling','Burial'),'topleft',lty=c(1,1),col=c('black','green','red'),cex=0.7)
  StratifiedShade(S1,S2,myYLim)
  
  # Plot the states
  par(mfrow=c(3,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  # Epi
  plot(YearFrac,(StatesEpiP$TP)/LakeVolumes$EpiVol,type='l',xlab='',ylab='P-epi (gP/m3)',main=myLake$LakeName)
  StratifiedShade(S1,S2,myYLim)
  # Hypo
  plot(YearFrac,(StatesHypoP$TP)/LakeVolumes$HypoVol,type='l',xlab='',ylab='P-hypo (gP/m3)')
  StratifiedShade(S1,S2,myYLim)

  # Sed
  myRange = range((StatesSedP$TP)/LakeAreas$TotalArea)
  myDiff = diff(myRange)
  myYLim = c(myRange[1]-0.5*myDiff,myRange[2]+0.5*myDiff)
  plot(YearFrac,(StatesSedP$TP)/LakeAreas$TotalArea,ylim=myYLim,type='l',xlab='',ylab='P-sed (gP/m2)')
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