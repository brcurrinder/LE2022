PlotLongTermResults = function(OutputFileNames,TurnDay=NA){

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
  Misc = data.frame(Secchimin=rep(NA,nFiles),Secchimean=rep(NA,nFiles),Secchimax=rep(NA,nFiles),HypoDOmin=rep(NA,nFiles),HypoDOminOnset=rep(NA,nFiles))

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

  myTime = cumsum(RunInfo$nYears)
  # Turn time around if necessary to show hysteresis
  if (!is.na(TurnDay)){
    newTime = myTime
    iTurn = which(newTime==TurnDay)+1
    StartT = iTurn-1
    StopT = length(myTime)
    #newTime[StartT:StopT] = (TurnDay-1) - (newTime[StartT:StopT] - (TurnDay-1))
    newTime[iTurn:length(newTime)] = (TurnDay) + (TurnDay - newTime[iTurn:length(newTime)])
    myTime = newTime
  }

  ###############################
  # Plot cools stuff, like Secchi
  par(mfrow=c(3,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myYLim = c(0,max(Misc$Secchimax))
  plot(myTime,Misc$Secchimax,type='l',ylim = myYLim,lty=2,col='blue',
       xlab='Years',ylab='Secchi (m)',main=paste(myLake$LakeName,sep=""))
  lines(myTime,Misc$Secchimean,col='blue')
  lines(myTime,Misc$Secchimin,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],Misc$Secchimax[StartT:StopT],lty=2,col='red')
    lines(myTime[StartT:StopT],Misc$Secchimean[StartT:StopT],lwd=2,col='red')
    lines(myTime[StartT:StopT],Misc$Secchimin[StartT:StopT],lty=2,col='red')
  }

  grid()
  # Hypo DO min onset
  par(mai = c(0.4,0.8, 0, 0.1),cex = 0.9)
  #myYLim = c(min(Misc$HypoDOmin),max(Misc$HypoDOmin))
  myYLim = c(0,max(Misc$HypoDOmin))
  plot(myTime,Misc$HypoDOmin,type='l',col='blue',ylim=myYLim,ylab='Minimum Hypo DO (g/m3)')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],Misc$HypoDOmin[StartT:StopT],lty=1,lwd=2,col='red')
  }
  grid()
  # Onset day for min hypo
  # myYLim = c(min(Misc$HypoDOminOnset),max(Misc$HypoDOminOnset))
  myYLim = c(0,max(Misc$HypoDOminOnset))
  plot(myTime,Misc$HypoDOminOnset,type='l',col='blue',ylim=myYLim,ylab='HypoDOminOnset (nDays after strat)')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],Misc$HypoDOminOnset[StartT:StopT],lty=1,lwd=2,col='red')
  }
  grid()

  ###############################
  # Plot More cool stuff
  SedO2Demand = CRates$OCSedimentR /12 *1000 /365 # mmol O2,C /m2/day
  par(mfrow=c(3,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myYLim = c(0,max(SedO2Demand))
  plot(myTime,SedO2Demand,type='l',ylim = myYLim,col='blue',
       xlab='Years',ylab='Sed O2 demand (mmol O2,C/m2LA/day)',main=paste(myLake$LakeName,sep=""))
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],SedO2Demand[StartT:StopT],lty=1,lwd=2,col='red')
  }
  grid()

  par(mai = c(0.4,0.8, 0, 0.1),cex = 0.9)
  myYLim = c(min(CStratified$SedAvailableOC,CStratified$SedOC2PRatio),max(CStratified$SedAvailableOC,CStratified$SedOC2PRatio))
  plot(myTime,CStratified$SedAvailableOC,type='l',ylim = myYLim,col='blue',
       xlab='Years',ylab='OC:P molar ratio',main="")
  grid()
  lines(myTime,CStratified$SedOC2PRatio,col='black')
  legend('topleft',c('Mixed surface','Sediment'),lty=c(1,1),col=c('blue','black'))

  # Plot OC and P retention
  myYLim = c(min(CStratified$OCRetention,PStratified$PRetention,na.rm=TRUE),max(CStratified$OCRetention,PStratified$PRetention,na.rm=TRUE))
  if (is.finite(myYLim[1]) & is.finite(myYLim[2])){
    plot(myTime,PStratified$PRetention,ylim=myYLim,type='l',col='blue',ylab='Retention (Export/Load')
    lines(myTime,CStratified$OCRetention,col='brown')
    grid()
    legend('topleft',c('P','OC_tot'),lty=c(1,1),col=c('blue','brown'))
  }

  # Epi OC all by itself
  par(mfrow=c(1,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myYLim = c(min(CStratified$EpiDOCmean),max(CStratified$EpiDOCmean+CStratified$EpiPOCmean))
  #myYLim[1] = 0
  plot(myTime,CStratified$EpiDOCmean+CStratified$EpiPOCmean,type='l',lwd=3,ylim = myYLim,col='green',
       xlab='Years',ylab='Epi TOC (g/m3)',main="")
  grid()
  #lines(myTime,CStratified$EpiDOCmean+CStratified$EpiPOCmean,col='green')
  if (!is.na(TurnDay)){
    #lines(myTime[StartT:StopT],CStratified$EpiDOCmean[StartT:StopT],lty=1,lwd=2,col='brown')
    lines(myTime[StartT:StopT],CStratified$EpiDOCmean[StartT:StopT]+CStratified$EpiPOCmean[StartT:StopT],lty=1,lwd=4,col='brown')
  }
  #legend('topleft',c('DOC','POC'),lty=c(1,1),col=c('brown','green'))
  
  
  ###############################
  # Plot OC Stratified
  par(mfrow=c(4,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myYLim = c(min(CStratified$EpiDOCmean),max(CStratified$EpiDOCmean+CStratified$EpiPOCmean))
  plot(myTime,CStratified$EpiDOCmean,type='l',ylim = myYLim,col='brown',
       xlab='Years',ylab='Epi OC (g/m3)',main=paste(myLake$LakeName,', Stratified',sep=""))
  grid()
  lines(myTime,CStratified$EpiDOCmean+CStratified$EpiPOCmean,col='green')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],CStratified$EpiDOCmean[StartT:StopT],lty=1,lwd=2,col='brown')
    lines(myTime[StartT:StopT],CStratified$EpiDOCmean[StartT:StopT]+CStratified$EpiPOCmean[StartT:StopT],lty=1,lwd=2,col='green')
  }
  legend('topleft',c('DOC','POC'),lty=c(1,1),col=c('brown','green'))

  myYLim = c(min(CStratified$HypoDOCmean),max(CStratified$HypoDOCmean+CStratified$HypoPOCmean))
  par(mai = c(0.4,0.8, 0, 0.1),cex = 0.9)
  plot(myTime,CStratified$HypoDOCmean,type='l',ylim = myYLim,col='brown',
       xlab='Years',ylab='Hypo OC (g/m3)',main="")
  grid()
  lines(myTime,CStratified$HypoDOCmean+CStratified$HypoPOCmean,col='green')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],CStratified$HypoDOCmean[StartT:StopT],lty=1,lwd=2,col='brown')
    lines(myTime[StartT:StopT],CStratified$HypoDOCmean[StartT:StopT]+CStratified$HypoPOCmean[StartT:StopT],lty=1,lwd=2,col='green')
  }

  myYLim = c(min(CStratified$SedPOCRmean),max(CStratified$SedPOCRmean+CStratified$SedPOCLmean))
  plot(myTime,CStratified$SedPOCRmean,type='l',ylim = myYLim,col='brown',
       xlab='Years',ylab='Sediment OC (g/m2SA)',main="")
  grid()
  lines(myTime,CStratified$SedPOCRmean+CStratified$SedPOCLmean,col='green')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],CStratified$SedPOCRmean[StartT:StopT],lty=1,lwd=2,col='brown')
    lines(myTime[StartT:StopT],CStratified$SedPOCRmean[StartT:StopT]+CStratified$SedPOCLmean[StartT:StopT],lty=1,lwd=2,col='green')
  }
  legend('topleft',c('POCR','POCL'),lty=c(1,1),col=c('brown','green'))

  myYLim = c(0,max(CStratified$SedAvailableOC,0))
  plot(myTime,CStratified$SedAvailableOC,ylim=myYLim,type='l',col='black',
       xlab='Years',ylab='Available OC (mgC/gSed)',main="")
  grid()
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],CStratified$SedAvailableOC[StartT:StopT],lty=1,lwd=2,col='red')
  }
  legend("topleft",legend=c('Oligo','Meso','Eutro'),lty=c(2,2,2),col=c('blue','grey','green'))

  ###############################
  # Plot P rates
  par(mfrow=c(2,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myYLim = c(min(0,min(PRates),min(cumsum(PRates$dSedPCum))),max(PRates,max(cumsum(PRates$dSedPCum))))
  myYLim = c(0,3)
  plot(myTime,PRates$PSettlingRebind,type='l',ylim = myYLim,col='green',
       xlab='Years',ylab='Sediment P rates (g/m2LA/y)',main=paste(myLake$LakeName,sep=""))
  grid()
  lines(myTime,PRates$PRecyclingRelease,col='red')
  lines(myTime,PRates$PBurial,col='blue')
  lines(myTime,PRates$dSedPCum,col='black')
  lines(myTime,cumsum(PRates$dSedPCum),col='black',lty=2,lwd=2)
  legend('topleft',c('Settling','Recycling','Burial(incl bound)','dSed','cumsum(dSed)'),lty=c(1,1,1,1,2),col=c('green','red','blue','black','black'))

  myYLim = c(min(PRates),max(PRates))
  myYLim = c(0,3)
  plot(myTime,PRates$PSettling,type='l',ylim = myYLim,col='black',
       xlab='Years',ylab='Sediment P rates (g/m2LA/y)',main=paste(myLake$LakeName,sep=""))
  grid()
  lines(myTime,PRates$PRebind,lty=2,col='black')
  lines(myTime,PRates$PRecycling,col='blue')
  lines(myTime,PRates$PRelease,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],PRates$PSettling[StartT:StopT],lty=1,lwd=2,col='black')
    lines(myTime[StartT:StopT],PRates$PRebind[StartT:StopT],lty=2,lwd=2,col='black')
    lines(myTime[StartT:StopT],PRates$PRecycling[StartT:StopT],lty=1,lwd=2,col='blue')
    lines(myTime[StartT:StopT],PRates$PRelease[StartT:StopT],lty=2,lwd=2,col='blue')
  }
  legend('topleft',c('Settling','Rebind','Recycling','Release'),lty=c(1,2,1,2),col=c('black','black','blue','blue'))

  ###############################
  # Plot P Stratified
  par(mfrow=c(4,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  myYLim = c(min(PStratified$Epimin),max(PStratified$Epimax))
  plot(myTime,PStratified$Epimax,type='l',ylim = myYLim,lty=2,col='blue',
       xlab='Years',ylab='Epi P (g/m3)',main=paste(myLake$LakeName,', Stratified',sep=""))
  grid()
  lines(myTime,PStratified$Epimean,col='blue')
  lines(myTime,PStratified$Epimin,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],PStratified$Epimax[StartT:StopT],lty=2,lwd=1,col='red')
    lines(myTime[StartT:StopT],PStratified$Epimean[StartT:StopT],lty=1,lwd=2,col='red')
    lines(myTime[StartT:StopT],PStratified$Epimin[StartT:StopT],lty=2,lwd=1,col='red')
  }
  legend('topleft',c('Mean','Min, Max'),lty=c(1,2),col=c('blue','blue'))

  myYLim = c(min(PStratified$Hypomin),max(PStratified$Hypomax))
  par(mai = c(0.4,0.8, 0, 0.1),cex = 0.9)
  plot(myTime,PStratified$Hypomax,type='l',ylim = myYLim,lty=2,col='blue',
       xlab='Years',ylab='Hypo P (g/m3)',main="")
  grid()
  lines(myTime,PStratified$Hypomean,col='blue')
  lines(myTime,PStratified$Hypomin,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],PStratified$Hypomax[StartT:StopT],lty=2,lwd=1,col='red')
    lines(myTime[StartT:StopT],PStratified$Hypomean[StartT:StopT],lty=1,lwd=2,col='red')
    lines(myTime[StartT:StopT],PStratified$Hypomin[StartT:StopT],lty=2,lwd=1,col='red')
  }

  myYLim = c(min(PStratified$Sedmin),max(PStratified$Sedmax))
  plot(myTime,PStratified$Sedmax,type='l',ylim = myYLim,lty=2,col='blue',
       xlab='Years',ylab='Sed P (g/m2SA)',main="")
  grid()
  lines(myTime,PStratified$Sedmean,col='black')
  lines(myTime,PStratified$Sedmin,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],PStratified$Sedmean[StartT:StopT],lty=1,lwd=2,col='red')
  }

  myYLim = c(0,max(PStratified$SedAvailableP,4))
  plot(myTime,PStratified$SedAvailableP,ylim=myYLim,type='l',col='black',
       xlab='Years',ylab='Available P (mgP/gSed (SA))',main="")
  grid()
  abline(h=1,lty=2,lwd=2,col='blue')
  abline(h=2,lty=2,lwd=2,col='grey')
  abline(h=3,lty=3,lwd=2,col='green')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],PStratified$SedAvailableP[StartT:StopT],lty=1,lwd=2,col='red')
  }
  legend("topleft",legend=c('Oligo','Meso','Eutro'),lty=c(2,2,2),col=c('blue','grey','green'))
  
  #####################################
  # Plots to be printed
  ###############################
  # Plot P Stratified
  
  par(mfrow=c(3,1),mai = c(0.05,0.05, 0.08, 0.05),oma = c(2,2,0.2,0.2), cex = 0.7)
  FileName = paste('./Figures/TPLongTerm.png',sep="")

  myYLim = c(min(PStratified$Epimin),max(PStratified$Epimax))
  plot(myTime,PStratified$Epimax,type='l',ylim = myYLim,lty=2,col='blue',
       xlab='Years',ylab='Epi P (g/m3)',main="",xaxt='n')
  grid()
  lines(myTime,PStratified$Epimean,col='blue',lwd=2)
  lines(myTime,PStratified$Epimin,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],PStratified$Epimax[StartT:StopT],lty=2,lwd=1,col='red')
    lines(myTime[StartT:StopT],PStratified$Epimean[StartT:StopT],lty=1,lwd=2,col='red')
    lines(myTime[StartT:StopT],PStratified$Epimin[StartT:StopT],lty=2,lwd=1,col='red')
  }
  #legend('topleft',c('Mean','Min, Max'),lty=c(1,2),col=c('blue','blue'))
  
  myYLim = c(min(PStratified$Hypomin),max(PStratified$Hypomax))
  plot(myTime,PStratified$Hypomax,type='l',ylim = myYLim,lty=2,col='blue',
       xlab='Years',ylab='Hypo P (g/m3)',main="",xaxt='n')
  grid()
  lines(myTime,PStratified$Hypomean,col='blue',lwd=2)
  lines(myTime,PStratified$Hypomin,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],PStratified$Hypomax[StartT:StopT],lty=2,lwd=1,col='red')
    lines(myTime[StartT:StopT],PStratified$Hypomean[StartT:StopT],lty=1,lwd=2,col='red')
    lines(myTime[StartT:StopT],PStratified$Hypomin[StartT:StopT],lty=2,lwd=1,col='red')
  }
  
  myYLim = c(min(PStratified$Sedmin),max(PStratified$Sedmax))
  plot(myTime,PStratified$Sedmax,type='l',ylim = myYLim,lty=2,col='blue',
       xlab='Years',ylab='Sed P (g/m2SA)',main="")
  grid()
  lines(myTime,PStratified$Sedmean,col='blue',lwd=2)
  lines(myTime,PStratified$Sedmin,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],PStratified$Sedmean[StartT:StopT],lty=1,lwd=2,col='red')
  }
  dev.copy(png,FileName,width=2,height=4,units="in",res=300)
  dev.off()
  ###############################
  # Plot OC Stratified
  
  par(mfrow=c(3,1),mai = c(0.05,0.05, 0.08, 0.05),oma = c(2,2,0.2,0.2), cex = 0.7)
  FileName = paste('./Figures/OCLongTerm.png',sep="")
  
  myYLim = c(min(CStratified$EpiDOCmean),max(CStratified$EpiDOCmean+CStratified$EpiPOCmean))
  plot(myTime,CStratified$EpiDOCmean,type='l',ylim = myYLim,col='brown',
       xlab='Years',ylab='Epi OC (g/m3)',main="",lwd=2,xaxt='n')
  grid()
  lines(myTime,CStratified$EpiDOCmean+CStratified$EpiPOCmean,col='green',lwd=2)
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],CStratified$EpiDOCmean[StartT:StopT],lty=1,lwd=2,col='brown')
    lines(myTime[StartT:StopT],CStratified$EpiDOCmean[StartT:StopT]+CStratified$EpiPOCmean[StartT:StopT],lty=1,lwd=2,col='green')
  }
  #legend('topleft',c('DOC','POC'),lty=c(1,1),col=c('brown','green'))
  
  myYLim = c(min(CStratified$HypoDOCmean),max(CStratified$HypoDOCmean+CStratified$HypoPOCmean))
  plot(myTime,CStratified$HypoDOCmean,type='l',ylim = myYLim,col='brown',
       xlab='Years',ylab='Hypo OC (g/m3)',main="",lwd=2,xaxt='n')
  grid()
  lines(myTime,CStratified$HypoDOCmean+CStratified$HypoPOCmean,col='green',lwd=2)
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],CStratified$HypoDOCmean[StartT:StopT],lty=1,lwd=2,col='brown')
    lines(myTime[StartT:StopT],CStratified$HypoDOCmean[StartT:StopT]+CStratified$HypoPOCmean[StartT:StopT],lty=1,lwd=2,col='green')
  }
  
  myYLim = c(min(CStratified$SedPOCRmean),max(CStratified$SedPOCRmean+CStratified$SedPOCLmean))
  plot(myTime,CStratified$SedPOCRmean,type='l',ylim = myYLim,col='brown',
       xlab='Years',ylab='Sediment OC (g/m2SA)',main="",lwd=2)
  grid()
  lines(myTime,CStratified$SedPOCRmean+CStratified$SedPOCLmean,col='green',lwd=2)
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],CStratified$SedPOCRmean[StartT:StopT],lty=1,lwd=2,col='brown')
    lines(myTime[StartT:StopT],CStratified$SedPOCRmean[StartT:StopT]+CStratified$SedPOCLmean[StartT:StopT],lty=1,lwd=2,col='green')
  }
  #legend('topleft',c('POCR','POCL'),lty=c(1,1),col=c('brown','green'))

  dev.copy(png,FileName,width=2,height=4,units="in",res=300)
  dev.off()
  
  ###############################
  # Plot DO and Secchi
  
  par(mfrow=c(3,1),mai = c(0.05,0.05, 0.08, 0.05),oma = c(2,2,0.2,0.2), cex = 0.7)
  FileName = paste('./Figures/DOandSecchiLongTerm.png',sep="")
  
  myYLim = c(0,max(Misc$Secchimax))
  plot(myTime,Misc$Secchimax,type='l',ylim = myYLim,lty=2,col='blue',
       xlab='Years',ylab='Secchi (m)',main="",xaxt='n')
  lines(myTime,Misc$Secchimean,col='blue',lwd=2)
  lines(myTime,Misc$Secchimin,lty=2,col='blue')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],Misc$Secchimax[StartT:StopT],lty=2,col='red')
    lines(myTime[StartT:StopT],Misc$Secchimean[StartT:StopT],lwd=2,col='red')
    lines(myTime[StartT:StopT],Misc$Secchimin[StartT:StopT],lty=2,col='red')
  }
  
  grid()
  # Hypo DO min onset
  #myYLim = c(min(Misc$HypoDOmin),max(Misc$HypoDOmin))
  myYLim = c(0,max(Misc$HypoDOmin))
  plot(myTime,Misc$HypoDOmin,type='l',col='blue',ylim=myYLim,ylab='Minimum Hypo DO (g/m3)',lwd=2,xaxt='n')
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],Misc$HypoDOmin[StartT:StopT],lty=1,lwd=2,col='red')
  }
  grid()
  # Onset day for min hypo
  # myYLim = c(min(Misc$HypoDOminOnset),max(Misc$HypoDOminOnset))
  myYLim = c(0,max(Misc$HypoDOminOnset))
  plot(myTime,Misc$HypoDOminOnset,type='l',col='blue',ylim=myYLim,ylab='HypoDOminOnset (nDays after strat)',lwd=2)
  if (!is.na(TurnDay)){
    lines(myTime[StartT:StopT],Misc$HypoDOminOnset[StartT:StopT],lty=1,lwd=2,col='red')
  }
  grid()
  
  
  dev.copy(png,FileName,width=2,height=4,units="in",res=300)
  dev.off()
}
