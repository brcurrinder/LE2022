PlotMetabResultsAllSimsAnnualSummaries = function(OutputFileNames){

  # ModelRun[[1]] = SimSetup
  # ModelRun[[2]] = myLake
  # ModelRun[[3]] = myParameters
  # ModelRun[[4]] = LakeVolumes
  # ModelRun[[5]] = LakeTemps
  # ModelRun[[6]] = LakeAreas
  # ModelRun[[7]] = MetabResults[[i]]
  # ModelRun[[8]] = MetabSummaries[[i]]

  # MetabSummary[[1]] = MassBalanceC
  # MetabSummary[[2]] = MassBalanceP
  # MetabSummary[[3]] = CarbonSummary
  # MetabSummary[[4]] = PhosphorusSummary
  # MetabSummary[[5]] = DOSummary
  # MetabSummary[[6]] = SecchiSummary
  nFiles = dim(OutputFileNames)[1]
  
  # myLake = ModelRun[[2]]
  
  StatesEpiC    = NULL
  StatesHypoC   = NULL
  StatesEpiGas  = NULL
  StatesHypoGas = NULL
  RatesEpiC     = NULL
  RatesHypoC    = NULL
  RatesEpiGas   = NULL
  RatesHypoGas  = NULL
  StatesEpiP    = NULL
  StatesHypoP   = NULL
  StatesSedP    = NULL
  RatesEpiP     = NULL
  RatesHypoP    = NULL
  RatesSedP     = NULL
  StatesSedC    = NULL
  RatesSedC     = NULL
  LakeVolumes   = NULL
  nDays = 0

  for (i in 1:nFiles){
    load(OutputFileNames$Names[i])
    
    #Load time series results
    MR = ModelRun[[7]]
    Hypso = ModelRun[[4]]
    myLake = ModelRun[[2]]
    
    # Lake area
    LakeArea = max(Hypso$area)
    nDays = nDays + ModelRun[[1]]$nDays[1]
    
    # Extract data frames from metabolism model
    StatesEpiC    = rbind(StatesEpiC,MR[[1]])
    StatesHypoC   = rbind(StatesHypoC,MR[[2]])
    StatesEpiGas  = rbind(StatesEpiGas,MR[[3]])
    StatesHypoGas = rbind(StatesHypoGas,MR[[4]])
    RatesEpiC     = rbind(RatesEpiC,MR[[5]])
    RatesHypoC    = rbind(RatesHypoC,MR[[6]])
    RatesEpiGas   = rbind(RatesEpiGas,MR[[7]])
    RatesHypoGas  = rbind(RatesHypoGas,MR[[8]])
    StatesEpiP    = rbind(StatesEpiP,MR[[9]])
    StatesHypoP   = rbind(StatesHypoP,MR[[10]])
    StatesSedP    = rbind(StatesSedP,MR[[11]])
    RatesEpiP     = rbind(RatesEpiP,MR[[12]])
    RatesHypoP    = rbind(RatesHypoP,MR[[13]])
    RatesSedP     = rbind(RatesSedP,MR[[14]])
    StatesSedC    = rbind(StatesSedC,MR[[15]])
    RatesSedC     = rbind(RatesSedC,MR[[16]])
    LakeVolumes   = rbind(LakeVolumes,MR[[17]])
    
    LakeAreas = data.frame(TotalArea = max(Hypso))
    
  }
  ########################################
  # Plot results
  YearFrac = round(1:nDays/365,digits=3)
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
  # Indeces for years can be determined as, e.g. the first year, which(floor(YearFrac)==0)

  ######################
  # Plot OC rates
  # Calculate a bunch of stuff
  RSettling = NULL
  LSettling = NULL
  Settling = NULL
  PropRSettling = NULL
  PropLSettling = NULL
  RespirationSed = NULL
  Burial = NULL
  Net = NULL
  SedimentCarbon = NULL
  SedO2Demand = NULL
  RespirationHypoDOC = NULL
  RespirationHypoPOC = NULL
  
  for (i in 1:length(Years)){
    iThisYear = which(floor(YearFrac)==Years[i])
    RSettling[i] = mean(RatesSedC$POCRSettling[iThisYear],na.rm=TRUE)*365 / LakeArea
    
    LSettling[i] = mean(RatesSedC$POCLSettling[iThisYear],na.rm=TRUE)*365 / LakeArea
    Settling[i] = RSettling[i] + LSettling[i]
    PropRSettling[i] = round(mean(RSettling[i]/Settling[i],na.rm=TRUE),digits=2)
    PropLSettling[i] = round(mean(LSettling[i]/Settling[i],na.rm=TRUE),digits=2)
    RespirationSed[i] = mean(RatesSedC$POCLR[iThisYear] + RatesSedC$POCRR[iThisYear],na.rm=TRUE)*365 / LakeArea
    Burial[i] = mean(RatesSedC$POCLBurial[iThisYear] + RatesSedC$POCRBurial[iThisYear],na.rm=TRUE)*365 / LakeArea
    Net[i] = Settling[i] - RespirationSed[i] - Burial[i]
    SedimentCarbon[i] = mean(StatesSedC$POCR[iThisYear]/LakeAreas$TotalArea+StatesSedC$POCL[iThisYear]/LakeAreas$TotalArea,na.rm=TRUE) # gC/m2LA
    SedO2Demand[i] = (RespirationSed[i] /12 *1000 /365)/myLake$SedAreaProp # mmol O2,C /m2SA/day
    RespirationHypoDOC[i] = mean(RatesHypoC$DOCRR[iThisYear] + RatesHypoC$DOCLR[iThisYear],na.rm=TRUE)*365/LakeArea
    RespirationHypoPOC[i] = mean(RatesHypoC$POCRR[iThisYear] + RatesHypoC$POCLR[iThisYear],na.rm=TRUE)*365/LakeArea
  }

  # Hypo
  par(mfrow=c(3,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myYLim = c(min(c(RespirationHypoDOC,RespirationHypoPOC,RespirationSed),na.rm=TRUE),
             max(c(RespirationHypoDOC,RespirationHypoPOC,RespirationSed),na.rm=TRUE))
  plot(Years,RespirationSed,ylim=myYLim,col='blue',type='l',ylab="Hypo R (gC/m2LA/y)",main=myLake$LakeName)
  grid(col='grey')
  lines(Years,RespirationHypoDOC,col='red')
  lines(Years,RespirationHypoPOC,col='green')
  #StratifiedShade(S1,S2,myYLim)
  legend('topleft',legend=c('Rsed','Rdoc','Rpoc'),lty=c(1,1,1),col=c('blue','red','green'))
  
  # Sed
  myYLim = c(min(c(Settling,RespirationSed,Burial,Net),na.rm=TRUE),
             max(c(Settling,RespirationSed,Burial,Net),na.rm=TRUE))
  plot(Years,Settling,ylim=myYLim,col='green',type='l',ylab='OCTot-Sed (gC/m2LA/y)',main="")
  grid(col='grey')
  lines(Years,RespirationSed,col='blue')
  lines(Years,Burial,col='brown')
  lines(Years,Net,col='grey',lty=2,lwd=2)
  abline(h=0,col='black',lty=1)
  # StratifiedShade(S1,S2,myYLim)
  legend(legend=c('Settling','Resp','Burial','Net'),'topleft',lty=c(1,1,1,2),col=c('green','blue','brown','grey'),cex=0.7)
  #myLegTxt = paste('Proportion settling POCR: ',PropRSettling,', POCL: ',PropLSettling)
  #legend(legend=myLegTxt,'topright')
  # Plot sediment oxygen demand
  par(mai = c(0.4,0.8, 0, 0.1),cex = 0.9)
  myYLim = c(0,max(SedO2Demand,na.rm=TRUE))
  plot(Years,SedO2Demand,ylim=myYLim,type='l',ylab='SedO2demand (mmolO2/m2SA/day)')
  grid(col='grey')
  #StratifiedShade(S1,S2,myYLim)

  ######################
  EDOCR = NULL
  EDOCL = NULL
  EPOCR = NULL
  EPOCL = NULL
  HDOCR = NULL
  HDOCL = NULL
  HPOCR = NULL
  HPOCL = NULL
  SPOCR = NULL
  SPOCL = NULL
  ESecchiMean = NULL
  ESecchiMin = NULL
  ESecchiMax = NULL
  
  for (i in 1:length(Years)){
    iThisYear = which(floor(YearFrac)==Years[i])
      EDOCR[i] = mean(StatesEpiC$DOCR[iThisYear]/LakeVolumes$EpiVol[iThisYear],na.rm=TRUE)
      EDOCL[i] = mean(StatesEpiC$DOCL[iThisYear]/LakeVolumes$EpiVol[iThisYear],na.rm=TRUE)
      EPOCR[i] = mean(StatesEpiC$POCR[iThisYear]/LakeVolumes$EpiVol[iThisYear],na.rm=TRUE)
      EPOCL[i] = mean(StatesEpiC$POCL[iThisYear]/LakeVolumes$EpiVol[iThisYear],na.rm=TRUE)
      
      # This works for hypolimnion, because hypo values that are zero (when lake is mixed)
      # are divided by hypolimnion volume, which is zero. The result is NA during mixed periods
      HDOCR[i] = mean(StatesHypoC$DOCR[iThisYear]/LakeVolumes$HypoVol[iThisYear],na.rm=TRUE)
      HDOCL[i] = mean(StatesHypoC$DOCL[iThisYear]/LakeVolumes$HypoVol[iThisYear],na.rm=TRUE)
      HPOCR[i] = mean(StatesHypoC$POCR[iThisYear]/LakeVolumes$HypoVol[iThisYear],na.rm=TRUE)
      HPOCL[i] = mean(StatesHypoC$POCL[iThisYear]/LakeVolumes$HypoVol[iThisYear],na.rm=TRUE)

      SPOCR[i] = mean(StatesSedC$POCR[iThisYear]/(LakeAreas$TotalArea*myLake$SedAreaProp),na.rm=TRUE)
      SPOCL[i] = mean(StatesSedC$POCL[iThisYear]/(LakeAreas$TotalArea*myLake$SedAreaProp),na.rm=TRUE)
      
      ESecchiMean[i] = mean(StatesEpiC$Secchi[iThisYear],na.rm=TRUE)
      ESecchiMin[i] = min(StatesEpiC$Secchi[iThisYear],na.rm=TRUE)
      ESecchiMax[i] = max(StatesEpiC$Secchi[iThisYear],na.rm=TRUE)
      
    }
  
  #DOCTot
  # Calc y lims
  YLower = min(EDOCR+EPOCR+EDOCL+EPOCL) - 2
  YUpper = max(EDOCR+EPOCR+EDOCL+EPOCL) + 2
  myYLim = c(YLower,YUpper)
  par(mfrow=c(4,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  # Epi
  plot(Years,(EDOCR),ylim=myYLim,type='l',ylab='DOCTot-epi (gC/m3)',main=myLake$LakeName)
  grid(col='grey')
  lines(Years,(EDOCR+EPOCR),col='grey')
  lines(Years,(EDOCR+EPOCR+EDOCL),col='brown')
  lines(Years,(EDOCR+EPOCR+EDOCL+EPOCL),col='green')
  # StratifiedShade(S1,S2,myYLim)
  legend(legend=c('DOCR','POCR','DOCL','POCL'),'topleft',lty=c(1,1,1,1),col=c('black','grey','brown','green'),cex=0.7)
  # Hypo
  YLower = 0
  YUpper = max(HDOCR+HPOCR+HDOCL+HPOCL,na.rm=TRUE) + 2
  myYLim = c(YLower,YUpper)
  plot(Years,(HDOCR),ylim=myYLim,type='l',ylab='DOCTot-Hypo (gC/m3)')
  grid(col='grey')
  lines(Years,(HDOCR+HPOCR),col='grey')
  lines(Years,(HDOCR+HPOCR+HDOCL),col='brown')
  lines(Years,(HDOCR+HPOCR+HDOCL+HPOCL),col='green')
  # StratifiedShade(S1,S2,myYLim)
  legend(legend=c('DOCR','POCR','DOCL','POCL'),'topleft',lty=c(1,1,1,1),col=c('black','grey','brown','green'),cex=0.7)
  # Sed
  myYLim = c(min(SPOCR),max(max(SPOCR+SPOCL)))
  plot(Years,(SPOCR),ylim=myYLim,col='brown',type='l',ylab='DOCTot-Sed (gC/m2SA)')
  grid(col='grey')
  lines(Years,(SPOCR+SPOCL),col='green')
  # StratifiedShade(S1,S2,myYLim)
  legend(legend=c('POCR','POCL'),'topleft',lty=c(1,1),col=c('brown','green'),cex=0.7)
  # Secchi
  myYLim = c(min(ESecchiMin),max(ESecchiMax))
  plot(Years,(ESecchiMean),ylim=myYLim,type='l',ylab='Secchi (m)')
  grid(col='grey')
  lines(Years,ESecchiMin,lty=2)
  lines(Years,ESecchiMax,lty=2)
  legend(legend=c('Mean','Range'),'topleft',lty=c(1,2),col=c('black','black'),cex=0.7)
  # StratifiedShade(S1,S2,myYLim)
  
  ######################
  #DO
  EDO = NULL
  EDOsat = NULL
  EDOFatm = NULL
  HDO = NULL
  HDOmin = NULL
  HDOmax = NULL
  HDOsat = NULL
  for (i in 1:length(Years)){
    iThisYear = which(floor(YearFrac)==Years[i])
    EDO[i] = mean(StatesEpiGas$DO[iThisYear]/LakeVolumes$EpiVol[iThisYear],na.rm=TRUE)
    EDOsat[i] = mean(StatesEpiGas$DOsat[iThisYear],na.rm=TRUE)
    EDOFatm[i] = mean(RatesEpiGas$DOFatm[iThisYear]/LakeAreas$TotalArea,na.rm=TRUE)
    
    HDO[i] = mean(StatesHypoGas$DO[iThisYear]/LakeVolumes$HypoVol[iThisYear],na.rm=TRUE)
    HDOmin[i] = min(StatesHypoGas$DO[iThisYear]/LakeVolumes$HypoVol[iThisYear],na.rm=TRUE)
    HDOmax[i] = max(StatesHypoGas$DO[iThisYear]/LakeVolumes$HypoVol[iThisYear],na.rm=TRUE)
    
    HDOsat[i] = mean(StatesHypoGas$DOsat[iThisYear],na.rm=TRUE)
    
  }
  
  # Calc y lims
  YLower = 0
  YUpper = 15
  myYLim = c(YLower,YUpper)
  
  par(mfrow=c(3,1), mai = c(0.3,0.8, 0.3, 0.1),cex = 0.9)
  plot(Years,(EDO),ylim=myYLim,type='l',xlab='',ylab='DO-epi (gO2/m3)',main=myLake$LakeName)
  grid(col='grey')
  lines(Years,(EDOsat),lty=2,col='grey')
  # StratifiedShade(S1,S2,myYLim)
  legend(legend=c('DO','DOsat'),'topleft',lty=c(1,2),col=c('black','grey','brown','green'),cex=0.7)
  plot(Years,(HDO),ylim=myYLim,type='l',xlab='',ylab='DO-Hypo (gO2/m3)')
  grid(col='grey')
  lines(Years,HDOmin,lty=2)
  lines(Years,(HDOmax),lty=2)
  legend(legend=c('DO','DOmin,max'),'topleft',lty=c(1,2),col=c('black','black','brown','green'),cex=0.7)
  #StratifiedShade(S1,S2,myYLim)
  myYLim = range(EDOFatm)
  plot(Years,EDOFatm,type='l',
       ylim=myYLim,xlab='',ylab='DO Flux (gO2/m2LA/d)')
  abline(h=0,lty=2)
  grid(col='grey')
  #StratifiedShade(S1,S2,myYLim)

  # Phosphorus
  
  ETPSettling = NULL
  ETPRebind = NULL
  ETPRecycling = NULL
  ETPRelease = NULL
  ETPInflow = NULL
  ETPOutflow = NULL
  ETPMix = NULL
  
  HTPRecycling = NULL
  HTPRelease = NULL
  HTPSettling = NULL
  HTPRebind = NULL
  
  STPRecycling = NULL
  STPRelease = NULL
  STPSettling = NULL
  STPRebind = NULL
  STPBurial = NULL
  STPBBurial = NULL
  
  for (i in 1:length(Years)){
    iThisYear = which(floor(YearFrac)==Years[i])
    # Note the mass balance for the epilimnion is:
    # ETPInflow*0.7+ ETPRelease + ETPRecycling = ETPSettling + ETPRebind + ETPOutflow + ETPMix + dTPWC; (ETPMix is negative)
    ETPSettling[i] = mean(RatesEpiP$TPSettling[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    ETPRebind[i] = mean(RatesEpiP$TPRebind[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    ETPRecycling[i] = mean(RatesEpiP$TPRecycling[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    ETPRelease[i] = mean(RatesEpiP$TPRelease[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    ETPInflow[i] = mean(RatesEpiP$TPInflow[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    ETPOutflow[i] = mean(RatesEpiP$TPOutflow[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    ETPMix[i] = mean(RatesEpiP$TPMix[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    
    HTPRecycling[i] = mean(RatesHypoP$TPRecycling[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    HTPRelease[i] = mean(RatesHypoP$TPRelease[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    HTPSettling[i] = mean(RatesHypoP$TPSettling[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    HTPRebind[i] = mean(RatesHypoP$TPRebind[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    
    STPRecycling[i] = mean(RatesSedP$TPRecycling[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    STPRelease[i] = mean(RatesSedP$TPRelease[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    STPSettling[i] = mean(RatesSedP$TPSettling[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    STPRebind[i] = mean(RatesSedP$TPRebind[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    STPBurial[i] = mean(RatesSedP$TPBurial[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    STPBBurial[i] = mean(RatesSedP$TPBBurial[iThisYear]/LakeAreas$TotalArea*365,na.rm=TRUE)
    
  }
  
  # Calc y lims
  myYLim = c(0,max(ETPSettling, ETPInflow, ETPOutflow))
  
  ######################
  # Plot the rates
  par(mfrow=c(3,1), mai = c(0.3,0.8, 0.3, 0.1),cex = 0.9)
  # Epi
  plot(Years,(ETPSettling),ylim=myYLim,type='l',xlab='',ylab='P-epi rates (gP/m2LA/y)',main=myLake$LakeName)
  grid(col='grey')
  lines(Years,ETPRebind,lty=2,col='black')
  lines(Years,ETPRecycling,lty=1,col='blue')
  lines(Years,ETPRelease,lty=2,col='blue')
  lines(Years,(ETPInflow),col='green')
  lines(Years,(ETPOutflow),col='red')
  # StratifiedShade(S1,S2,myYLim)
  legend(legend=c('Settling','Rebind','Recycling','Release','Inflow','Outflow'),'topleft',lty=c(1,2,1,2,1,1),col=c('black','black','blue','blue','green','red'),cex=0.7)
  # Hypo
  plot(Years,(HTPRecycling),ylim=myYLim,type='l',col='blue',xlab='',ylab='P-hypo rates (gP/m2LA/y)',main="")
  grid(col='grey')
  lines(Years,HTPRelease,lty=2,col='blue')
  lines(Years,(HTPSettling),col='black')
  lines(Years,HTPRebind,lty=2,col='black')
  # StratifiedShade(S1,S2,myYLim)
  legend(legend=c('Recycling','Release','Settling','Rebinding'),'topleft',lty=c(1,2,1,2),col=c('blue','blue','black','black'),cex=0.7)
  # Sed
  plot(Years,(STPSettling),ylim=myYLim,type='l',xlab='',ylab='P-sed rates (gP/m2LA/y)',main="")
  grid(col='grey')
  lines(Years,(STPRebind),lty=2,col='black')
  lines(Years,(STPRecycling),col='blue')
  lines(Years,STPRelease,lty=2,col='blue')
  
  lines(Years,(STPBurial),col='red')
  lines(Years,(STPBBurial),lty=2,col='red')
  legend(legend=c('Settling','Rebinding','Recycling','Release','Burial','BurialBound'),'topleft',
         lty=c(1,2,1,2,1,2),col=c('black','black','blue','blue','red','red'),cex=0.7)
  # StratifiedShade(S1,S2,myYLim)
  
  ######################
  # Plot the states
  ETP = NULL
  ETPmin = NULL
  ETPmax = NULL
  HTP = NULL
  STP = NULL
  STPB = NULL
  
  for (i in 1:length(Years)){
    iThisYear = which(floor(YearFrac)==Years[i])
    ETP[i] = mean(StatesEpiP$TP[iThisYear]/LakeVolumes$EpiVol[iThisYear],na.rm=TRUE)
    ETPmin[i] = min(StatesEpiP$TP[iThisYear]/LakeVolumes$EpiVol[iThisYear],na.rm=TRUE)
    ETPmax[i] = max(StatesEpiP$TP[iThisYear]/LakeVolumes$EpiVol[iThisYear],na.rm=TRUE)
    HTP[i] = mean(StatesHypoP$TP[iThisYear]/LakeVolumes$HypoVol[iThisYear],na.rm=TRUE)
    STP[i] = mean(StatesSedP$TP[iThisYear]/(LakeAreas$TotalArea*myLake$SedAreaProp),na.rm=TRUE)  
    STPB[i] = mean(StatesSedP$TPB[iThisYear]/(LakeAreas$TotalArea*myLake$SedAreaProp),na.rm=TRUE) 
  }
  
  myYLim = c(min(ETPmin),max(ETPmax))
  par(mfrow=c(3,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  # Epi
  plot(Years,(ETP),type='l',ylim = myYLim,xlab='',ylab='P-epi (gP/m3)',main=myLake$LakeName)
  grid(col='grey')
  lines(Years,ETPmin,lty = 2)
  lines(Years,ETPmax, lty = 2)
  legend(legend=c('Mean','Range'),'topleft',lty=c(1,2),col=c('black','black'),cex=0.7)
  # StratifiedShade(S1,S2,myYLim)
  # Hypo
  plot(Years,(HTP),type='l',xlab='',ylab='P-hypo (gP/m3)')
  grid(col='grey')
  # StratifiedShade(S1,S2,myYLim)
  
  # Sed
  myRange = range((STP+STPB))
  myDiff = diff(myRange)
  myYLim = c(myRange[1]-0.5*myDiff,myRange[2]+0.5*myDiff)
  plot(Years,(STP+STPB),ylim=myYLim,type='l',xlab='',ylab='P-sed (gP/m2SA)')
  grid(col='grey')
  # StratifiedShade(S1,S2,myYLim)
  
}

StratifiedShade <- function(S1,S2,myYLim){
  for (iMix in 1:length(S1)){
    #polygon(c(1,2,2,1),c(99.2,99.2,100,100),col=rgb(0, 0, 0,0.2))
    polygon(c(S1[iMix],S2[iMix],S2[iMix],S1[iMix]),c(min(myYLim),min(myYLim),max(myYLim),max(myYLim)),
            col=rgb(0, 0, 0,0.1),border=NA)
    grid()
  }
}
