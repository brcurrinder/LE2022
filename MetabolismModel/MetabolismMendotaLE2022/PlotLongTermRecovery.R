PlotLongTermRecovery = function(OutputFileNames,TurnDay=NA,PrintLegend){

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

  # User defined water quality thresholds
  TrophicThresholds = GetTrophicThresholds()
  nCol = dim(TrophicThresholds)[2]
  # Order thresholds by eutrophic-Meso, Meso-Oligo
  ChlThresholds         = TrophicThresholds[which(TrophicThresholds$Variable=='Chl'),c(3,2)]
  EpiPThresholds        = TrophicThresholds[which(TrophicThresholds$Variable=='TP'),c(3,2)]
  SecchiThresholds      = 1.7 / TrophicThresholds[which(TrophicThresholds$Variable=='Secchi'),c(3,2)]
  nAnoxiaDaysThresholds = TrophicThresholds[which(TrophicThresholds$Variable=='AnoxicDays'),c(3,2)]
  EpiPOCThresholds      = TrophicThresholds[which(TrophicThresholds$Variable=='EpiPOC'),c(3,2)]
  AnoxicThreshold = 0.2

  nFiles = dim(OutputFileNames)[1]

  PStratified = data.frame(Epimin=rep(NA,nFiles),Epimean=rep(NA,nFiles),Epimax=rep(NA,nFiles),
                           Hypomin=rep(NA,nFiles),Hypomean=rep(NA,nFiles),Hypomax=rep(NA,nFiles),
                           Sedmin=rep(NA,nFiles),Sedmean=rep(NA,nFiles),Sedmax=rep(NA,nFiles),SedAvailableP=rep(NA,nFiles),
                           PRetention=rep(NA,nFiles))
  CStratified = data.frame(EpiDOCmean=rep(NA,nFiles),EpiPOCmean=rep(NA,nFiles),
                           HypoDOCmean=rep(NA,nFiles),HypoPOCmean=rep(NA,nFiles),
                           SedPOCLmean=rep(NA,nFiles),SedPOCRmean=rep(NA,nFiles),
                           SedAvailableOC=rep(NA,nFiles),EpiOC2PRatioStratified=rep(NA,nFiles),SedOC2PRatio=rep(NA,nFiles),
                           OCRetention=rep(NA,nFiles))
  CRates = data.frame(OCSedimentR=rep(NA,nFiles))
  PRates = data.frame(PRecyclingRelease=rep(NA,nFiles),PRecycling=rep(NA,nFiles),PRelease=rep(NA,nFiles),
                      PSettling=rep(NA,nFiles),PRebind=rep(NA,nFiles),PBurial=rep(NA,nFiles),dSedPCum=rep(NA,nFiles),PLoads=rep(NA,nFiles))
  Misc = data.frame(Secchimin=rep(NA,nFiles),Secchimean=rep(NA,nFiles),Secchimax=rep(NA,nFiles),HypoDOmin=rep(NA,nFiles),HypoDOminOnset=rep(NA,nFiles),HypoDOAnoxicDays=rep(NA,nFiles))

  RunInfo = data.frame(nYears=rep(NA,nFiles))

  for (i in 1:nFiles){
    load(OutputFileNames$Names[i])
    
    #Load time series results
    
    # Everything from here down is summary info for this simulation
    SimSetup = ModelRun[[1]]
    myLake = ModelRun[[2]]
    myParameters = ModelRun[[3]]
    MS = ModelRun[[8]]
    MassBalanceC = MS[[1]]
    MassBalanceP = MS[[2]]
    CarbonSummary = MS[[3]]
    PhosphorusSummary = MS[[4]]
    DOSummary = MS[[5]]
    SecchiSummary = MS[[6]]

    RunInfo$nYears[i] = SimSetup$nDays[1]/365

    # Load up the Secchi and Misc
    Misc$Secchimin[i] = SecchiSummary$Secchimin
    Misc$Secchimean[i] = SecchiSummary$Secchimean
    Misc$Secchimax[i] = SecchiSummary$Secchimax
    
    # Calculation that should have been done in the summarize code
    LakeVolumes   = ModelRun[[7]][[17]]
    iStratified = which(LakeVolumes$HypoVol>0)
    StatesHypoGas = ModelRun[[7]][[4]]
    Misc$HypoDOAnoxicDays[i] = length(which(StatesHypoGas$DO[iStratified]/LakeVolumes$HypoVol[iStratified] < AnoxicThreshold)) / RunInfo$nYears[i]

    Misc$HypoDOminOnset[i] = DOSummary$HypoDOminOnset
    Misc$HypoDOmin[i] = DOSummary$HypoDOmin

    # Load up the phosphorus
    PStratified$Epimin[i] = PhosphorusSummary$EpiPStratifiedmin
    PStratified$Epimean[i] = PhosphorusSummary$EpiPStratifiedmean
    PStratified$Epimax[i] = PhosphorusSummary$EpiPStratifiedmax

    PStratified$Hypomin[i] = PhosphorusSummary$HypoPStratifiedmin
    PStratified$Hypomean[i] = PhosphorusSummary$HypoPStratifiedmean
    PStratified$Hypomax[i] = PhosphorusSummary$HypoPStratifiedmax

    PStratified$Sedmin[i] = PhosphorusSummary$SedPmin
    PStratified$Sedmean[i] = PhosphorusSummary$SedPmean
    PStratified$Sedmax[i] = PhosphorusSummary$SedPmax
    # Calculate the sediment available P in mgP/gSed
    gSedAreal = myLake$SedBulkDensity * myLake$ActiveSedDepth
    PStratified$SedAvailableP[i] = PhosphorusSummary$SedPmean / gSedAreal * 1000 # mgP/gSed
    PStratified$PRetention[i] = MassBalanceP$PRetention

    PRates$PSettlingRebind[i] = MassBalanceP$PSettlingRebind
    PRates$PSettling[i] = MassBalanceP$PSettling
    PRates$PRebind[i] = MassBalanceP$PRebind
    PRates$PRecyclingRelease[i] = MassBalanceP$PRecyclingRelease
    PRates$PRecycling[i] = MassBalanceP$PRecycling
    PRates$PRelease[i] = MassBalanceP$PRelease
    PRates$PBurial[i] = MassBalanceP$PBurial+MassBalanceP$PBBurial
    PRates$dSedPCum[i] = (PRates$PSettlingRebind[i] - PRates$PRecyclingRelease[i] - PRates$PBurial[i]) * RunInfo$nYears[i]
    PRates$PLoads[i] = MassBalanceP$PLoad

    # Load up the Organic carbon
    CRates$OCSedimentR[i] = MassBalanceC$OCSedimentR

    CStratified$EpiDOCmean[i] = CarbonSummary$EpiDOCStratifiedmean
    CStratified$EpiPOCmean[i] = CarbonSummary$EpiPOCStratifiedmean
    CStratified$HypoDOCmean[i] = CarbonSummary$HypoDOCStratifiedmean
    CStratified$HypoPOCmean[i] = CarbonSummary$HypoPOCStratifiedmean
    CStratified$SedPOCLmean[i] = CarbonSummary$SedPOCLmean
    CStratified$SedPOCRmean[i] = CarbonSummary$SedPOCRmean
    CStratified$OCRetention[i] = MassBalanceC$OCRetention
    # Calculate the sediment available OC in mgOC/gSed
    gSedAreal = myLake$SedBulkDensity * myLake$ActiveSedDepth
    CStratified$SedAvailableOC[i] = CarbonSummary$SedPOCmean / gSedAreal * 1000 # mgC/gSed

    # Calculate molar rations of C:P
    EpiCarbonStratified = (CarbonSummary$EpiPOCStratifiedmean + CarbonSummary$EpiDOCStratifiedmean) / 12 # molar units
    EpiPhosphStratified = PhosphorusSummary$EpiPStratifiedmean / 31 # molar units
    CStratified$EpiOC2PRatioStratified[i] = EpiCarbonStratified/EpiPhosphStratified
    SedCarbon = (CarbonSummary$SedPOCRmean + CarbonSummary$SedPOCLmean) / 12 # molar units
    SedPhosph =  PhosphorusSummary$SedPmean / 31 # molar units
    CStratified$SedOC2PRatio[i] = SedCarbon/SedPhosph

  }

  myTime = cumsum(RunInfo$nYears)-10
  # Turn time around if necessary to show hysteresis
  if (!is.na(TurnDay)){
    newTime = myTime
    StartT = TurnDay
    StopT = length(myTime)
    newTime[StartT:StopT] = (TurnDay-1) - (newTime[StartT:StopT] - (TurnDay-1))
    myTime = newTime
  }
  
  # Calculate water quality transitions
  # Secchi
  LEC = 1.7/Misc$Secchimin
  LECThresh =    c(which(LEC <= SecchiThresholds$MesotrophicUpper)[1],
                   which(LEC <= SecchiThresholds$OligotrophicUpper)[1]) * mean(RunInfo$nYears)-10
  AnoxiaThresh = c(which(Misc$HypoDOAnoxicDays <= nAnoxiaDaysThresholds$MesotrophicUpper)[1],
                   which(Misc$HypoDOAnoxicDays <= nAnoxiaDaysThresholds$OligotrophicUpper)[1]) * mean(RunInfo$nYears)-10
  EpiPOCThresh = c(which(CStratified$EpiPOCmean <= EpiPOCThresholds$MesotrophicUpper)[1],
                   which(CStratified$EpiPOCmean <= EpiPOCThresholds$OligotrophicUpper)[1]) * mean(RunInfo$nYears)-10
  EpiPThresh =   c(which(PStratified$Epimean <= EpiPThresholds$MesotrophicUpper)[1],
                   which(PStratified$Epimean <= EpiPThresholds$OligotrophicUpper)[1]) * mean(RunInfo$nYears)-10

  ###############################
  # First recovery plot
  par(mfrow=c(1,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myCols = c('black','red','red','green','green','grey','grey')
  myLtys = c(1,1,2,1,2,1,2)
  myLwds = c(2,2,2,2,2,2,2,2)
  myLabels = c('Secchi','AnoxicDays','SedO2Demand','EpiPOC','HypoPOC','EpiTP','SedTP')
  myYLim = c(0,1)
  iVar = 0 # Variable counter
  # myTime[StartT:StopT],SedO2Demand[StartT:StopT]
  # Phosphorus
  
  # Secchi
  myVar = 1.7/ Misc$Secchimin
  iVar = iVar+1
  # zi = (xi – min(x)) / (max(x) – min(x))
  myVarScaled = (myVar - min(myVar)) / (max(myVar) - min(myVar))
  plot(myTime,myVarScaled,type='l',ylim=myYLim,col=myCols[iVar],lwd=myLwds[iVar],
       ylab='Normalized across recovery range')
  # abline(h=0,lty=2)
  # abline(h=1,lty=2)
  
  # Hypolimnetic DO
  myVar = Misc$HypoDOAnoxicDays
  iVar = iVar+1
  myVarScaled = (myVar - min(myVar)) / (max(myVar) - min(myVar))
  lines(myTime,myVarScaled,col=myCols[iVar],lty=myLtys[iVar],lwd=myLwds[iVar])

  # SedO2Demand
  SedO2Demand = CRates$OCSedimentR /12 *1000 /365 # mmol O2,C /m2/day
  myVar = SedO2Demand
  iVar = iVar+1
  myVarScaled = (myVar - min(myVar)) / (max(myVar) - min(myVar))
  lines(myTime,myVarScaled,col=myCols[iVar],lty=myLtys[iVar],lwd=myLwds[iVar])

  # Epi OC
  # EpiTOC = CStratified$EpiDOCmean+CStratified$EpiPOCmean
  myVar = CStratified$EpiPOCmean
  iVar = iVar+1
  myVarScaled = (myVar - min(myVar)) / (max(myVar) - min(myVar))
  lines(myTime,myVarScaled,col=myCols[iVar],lty=myLtys[iVar],lwd=myLwds[iVar])
  
  # Hypo POC
  myVar = CStratified$HypoPOCmean
  iVar = iVar+1
  myVarScaled = (myVar - min(myVar)) / (max(myVar) - min(myVar))
  lines(myTime,myVarScaled,col=myCols[iVar],lty=myLtys[iVar],lwd=myLwds[iVar])

  myVar = PStratified$Epimean
  iVar = iVar+1
  myVarScaled = (myVar - min(myVar)) / (max(myVar) - min(myVar))
  lines(myTime,myVarScaled,col=myCols[iVar],lty=myLtys[iVar],lwd=myLwds[iVar])
  
  myVar = PStratified$Sedmean
  iVar = iVar+1
  myVarScaled = (myVar - min(myVar)) / (max(myVar) - min(myVar))
  lines(myTime,myVarScaled,col=myCols[iVar],lty=myLtys[iVar],lwd=myLwds[iVar])

  # grid()
  # Plot the thresholds
  abline(v=0,col='black',lty=2)
  abline(v=EpiPThresh,col='brown',lty=2)
  abline(v=LECThresh,col='blue',lty=2)
  abline(v=AnoxiaThresh,col='grey',lty=2)
  #abline(v=c(EpiPThresh[1],AnoxiaThresh[1],LECThresh[1]),col='grey',lty=3,lwd=3)
  #abline(v=c(EpiPThresh[2],AnoxiaThresh[2],LECThresh[2]),col='blue',lty=3,lwd=3)
  # points(c(EpiPThresh[1],AnoxiaThresh[1],LECThresh[1]),rep(0,3),col='green',pch=22)
  # points(c(EpiPThresh[2],AnoxiaThresh[2],LECThresh[2]),rep(0,3),col='blue',pch=22)
  myLabels = c(myLabels,"EpiTPThresh","LECThresh","AnoxiaThresh")
  myCols = c(myCols,'brown','blue','grey')
  myLtys = c(myLtys,2,2,2)
  
  if (PrintLegend){
    legend('topright',myLabels,lty=myLtys,col=myCols)
  }else{
    # Don't print the legend
  }
  
}
