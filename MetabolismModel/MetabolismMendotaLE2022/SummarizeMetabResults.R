SummarizeMetabResults = function(MetabResults,myLake,LakeHypsometry,nYears2Save,Print){
  
  AnoxicThreshold = 1 # Below this value, we consider the system anoxic
  # Setup indeces for date ranges
  nRecs2Save = nYears2Save*365
  ibOutput = nDays-nRecs2Save+1 # Beginning index
  ieOutput = nDays # Ending index
  FirstDay=1 # This gets changed if there's only on year because of NAs
  
  # Calculate lake area and volume
  LakeArea = max(LakeHypsometry$area)
  
  # Extract data frames from metabolism model and shorten to the desired analysis period
  StatesEpiC = MetabResults[[1]][ibOutput:ieOutput,]
  StatesHypoC = MetabResults[[2]][ibOutput:ieOutput,]
  StatesEpiGas = MetabResults[[3]][ibOutput:ieOutput,]
  StatesHypoGas = MetabResults[[4]][ibOutput:ieOutput,]
  RatesEpiC = MetabResults[[5]][ibOutput:ieOutput,]
  RatesHypoC = MetabResults[[6]][ibOutput:ieOutput,]
  RatesEpiGas = MetabResults[[7]][ibOutput:ieOutput,]
  RatesHypoGas = MetabResults[[8]][ibOutput:ieOutput,]
  
  # Change starting index for the year if the first rate value = NA
  # This happens for the very first record of the whole simulation
  if (is.na(RatesEpiC$DOCRR[1])){
    FirstDay=2
  }
  
  StatesEpiP = MetabResults[[9]][ibOutput:ieOutput,]
  StatesHypoP = MetabResults[[10]][ibOutput:ieOutput,]
  StatesSedP = MetabResults[[11]][ibOutput:ieOutput,] 
  RatesEpiP = MetabResults[[12]][ibOutput:ieOutput,] 
  RatesHypoP = MetabResults[[13]][ibOutput:ieOutput,] 
  RatesSedP = MetabResults[[14]][ibOutput:ieOutput,]
  
  StatesSedC = MetabResults[[15]][ibOutput:ieOutput,]
  RatesSedC = MetabResults[[16]][ibOutput:ieOutput,]
  
  LakeVolumes = MetabResults[[17]][ibOutput:ieOutput,]
  
  # Setup return data frame
  MassBalanceC = data.frame(Allocthony=NA)
  MassBalanceP = data.frame(PLoad=NA)
  

  ##################################
  # Calculate the phosphorus budget
  MassBalanceP$PLoad = sum(RatesEpiP$TPInflow,na.rm=TRUE) / LakeArea /nRecs2Save*365 # gP/m2LA/y
  MassBalanceP$PExport = sum(RatesEpiP$TPOutflow,na.rm=TRUE) / LakeArea /nRecs2Save*365 # gP/m2LA/y
  MassBalanceP$PBurial = sum(RatesSedP$TPBurial,na.rm=TRUE) / LakeArea /nRecs2Save*365
  MassBalanceP$PBBurial = sum(RatesSedP$TPBBurial,na.rm=TRUE) / LakeArea /nRecs2Save*365
  # Settling and recycling were summed and in main model and saved in sediment
  MassBalanceP$PSettlingRebind = sum((RatesSedP$TPSettling+RatesSedP$TPRebind),na.rm=TRUE) / LakeArea /nRecs2Save*365
  MassBalanceP$PSettling = sum((RatesSedP$TPSettling),na.rm=TRUE) / LakeArea /nRecs2Save*365
  MassBalanceP$PRebind = sum((RatesSedP$TPRebind),na.rm=TRUE) / LakeArea /nRecs2Save*365
  MassBalanceP$PRecyclingRelease = sum((RatesSedP$TPRecycling+RatesSedP$TPRelease),na.rm=TRUE) / LakeArea /nRecs2Save*365
  MassBalanceP$PRecycling = sum((RatesSedP$TPRecycling),na.rm=TRUE) / LakeArea /nRecs2Save*365
  MassBalanceP$PRelease = sum((RatesSedP$TPRelease),na.rm=TRUE) / LakeArea /nRecs2Save*365
  
  
  # Change in water column values
  MassBalanceP$dPwc = ((StatesEpiP$TP[nRecs2Save] + StatesHypoP$TP[nRecs2Save]) - 
    (StatesEpiP$TP[FirstDay] + StatesHypoP$TP[FirstDay])) / LakeArea
  MassBalanceP$dPsed = ((StatesSedP$TP[nRecs2Save]+ StatesSedP$TPB[nRecs2Save]) - (StatesSedP$TP[FirstDay]+StatesSedP$TPB[FirstDay])) / LakeArea
  MassBalanceP$PBalance = MassBalanceP$PLoad - MassBalanceP$PExport - MassBalanceP$PBurial - MassBalanceP$PBBurial
  MassBalanceP$PRetention = 1 - (MassBalanceP$PExport / MassBalanceP$PLoad)
  ##################################
  # Calculate phosphorus characteristics
  iStratified = which(LakeVolumes$HypoVol>0)
  iMixed = which(LakeVolumes$HypoVol==0)
  
  EpiPMixedmin = min(StatesEpiP$TP[iMixed]/LakeVolumes$EpiVol[iMixed])
  EpiPMixedmean = mean(StatesEpiP$TP[iMixed]/LakeVolumes$EpiVol[iMixed])
  EpiPMixedmax = max(StatesEpiP$TP[iMixed]/LakeVolumes$EpiVol[iMixed])

  EpiPStratifiedmin = min(StatesEpiP$TP[iStratified]/LakeVolumes$EpiVol[iStratified])
  EpiPStratifiedmean = mean(StatesEpiP$TP[iStratified]/LakeVolumes$EpiVol[iStratified])
  EpiPStratifiedmax = max(StatesEpiP$TP[iStratified]/LakeVolumes$EpiVol[iStratified])

  HypoPStratifiedmin = min(StatesHypoP$TP[iStratified]/LakeVolumes$HypoVol[iStratified])
  HypoPStratifiedmean = mean(StatesHypoP$TP[iStratified]/LakeVolumes$HypoVol[iStratified])
  HypoPStratifiedmax = max(StatesHypoP$TP[iStratified]/LakeVolumes$HypoVol[iStratified])
  
  SedPmin = min((StatesSedP$TP+StatesSedP$TPB)/(LakeArea*myLake$SedAreaProp))
  SedPmean = mean((StatesSedP$TP+StatesSedP$TPB)/(LakeArea*myLake$SedAreaProp))
  SedPmax = max((StatesSedP$TP+StatesSedP$TPB)/(LakeArea*myLake$SedAreaProp))
  
  PhosphorusSummary = data.frame(EpiPMixedmin=EpiPMixedmin,EpiPMixedmean=EpiPMixedmean,EpiPMixedmax=EpiPMixedmax,
                                 EpiPStratifiedmin=EpiPStratifiedmin,EpiPStratifiedmean=EpiPStratifiedmean,EpiPStratifiedmax=EpiPStratifiedmax,
                                 HypoPStratifiedmin=HypoPStratifiedmin,HypoPStratifiedmean=HypoPStratifiedmean,HypoPStratifiedmax=HypoPStratifiedmax,
                                 SedPmin=SedPmin,SedPmean=SedPmean,SedPmax=SedPmax)
  
  #mySoilDepEquil = SedimentDepositionEquilibrium(myLake,MassBalanceP$PRetention)
  
  ##################################
  # Calculate the carbon budget
  MassBalanceC$Allocthony = (sum(RatesEpiC$DOCLInflow,na.rm=TRUE) + sum(RatesEpiC$DOCRInflow,na.rm=TRUE) +
                  sum(RatesEpiC$POCLInflow,na.rm=TRUE) + sum(RatesEpiC$POCRInflow,na.rm=TRUE)) / LakeArea /nRecs2Save*365 # gOC/m2LA/y
  MassBalanceC$Autochthony = (sum(RatesEpiC$GPP,na.rm=TRUE)) / LakeArea /nRecs2Save*365 # gOC/m2LA/y
  MassBalanceC$OCExport = (sum(RatesEpiC$DOCLOutflow,na.rm=TRUE) + sum(RatesEpiC$DOCROutflow,na.rm=TRUE) +
                sum(RatesEpiC$POCLOutflow,na.rm=TRUE) + sum(RatesEpiC$POCROutflow,na.rm=TRUE)) / LakeArea /nRecs2Save*365 # gOC/m2LA/y
  MassBalanceC$OCR = (sum(RatesEpiC$DOCLR,na.rm=TRUE) + sum(RatesEpiC$DOCRR,na.rm=TRUE) + 
           sum(RatesEpiC$POCLR,na.rm=TRUE) + sum(RatesEpiC$POCRR,na.rm=TRUE) +
           sum(RatesHypoC$DOCLR,na.rm=TRUE) + sum(RatesHypoC$DOCRR,na.rm=TRUE) + 
           sum(RatesHypoC$POCLR,na.rm=TRUE) + sum(RatesHypoC$POCRR,na.rm=TRUE)) / LakeArea /nRecs2Save*365 # gOC/m2LA/y
  MassBalanceC$OCSedimentation = (sum(RatesHypoC$POCLSettling,na.rm=TRUE) + sum(RatesHypoC$POCRSettling,na.rm=TRUE)) / LakeArea /nRecs2Save*365
  # Add sedimentation from the epilimnion when the lake is not stratified
  iNotStratified = which(LakeVolumes$HypoVol==0)
  SumAllEpiSeds = RatesEpiC$POCLSettling + RatesEpiC$POCRSettling
  EpiSedimentation = sum(SumAllEpiSeds[iNotStratified],na.rm=TRUE) / LakeArea
  MassBalanceC$OCSedimentation = MassBalanceC$OCSedimentation + EpiSedimentation
  MassBalanceC$OCSedimentR = (sum(RatesSedC$POCLR,na.rm=TRUE) + sum(RatesSedC$POCRR,na.rm=TRUE)) / LakeArea
  MassBalanceC$OCBurial = (sum(RatesSedC$POCLBurial,na.rm=TRUE) + sum(RatesSedC$POCRBurial,na.rm=TRUE)) / LakeArea             
  # MassBalanceC$SedimentR = (sum(RatesEpiGas$DOSed) + sum(RatesHypoGas$DOSed)) / LakeArea
  
  # Calculate the carbon balance
  MassBalanceC$OCBalance = MassBalanceC$Allocthony + MassBalanceC$Autochthony - MassBalanceC$OCExport -
    MassBalanceC$OCR - MassBalanceC$OCSedimentR - MassBalanceC$OCBurial
  
  # To calculate change in water column pool of C
  MassBalanceC$EndC = (StatesEpiC$DOCL[nRecs2Save] + StatesEpiC$DOCR[nRecs2Save] + StatesEpiC$POCL[nRecs2Save] + StatesEpiC$POCR[nRecs2Save] + 
            StatesHypoC$DOCL[nRecs2Save] + StatesHypoC$DOCR[nRecs2Save] + StatesHypoC$POCL[nRecs2Save] + StatesHypoC$POCR[nRecs2Save]) / LakeArea
  MassBalanceC$BeginC = (StatesEpiC$DOCL[FirstDay] + StatesEpiC$DOCR[FirstDay] + StatesEpiC$POCL[FirstDay] + StatesEpiC$POCR[FirstDay] + 
              StatesHypoC$DOCL[FirstDay] + StatesHypoC$DOCR[FirstDay] + StatesHypoC$POCL[FirstDay] + StatesHypoC$POCR[FirstDay]) / LakeArea
  # To calculate change in sediment pool of OC
  MassBalanceC$EndSedC = (StatesSedC$POCL[nRecs2Save] + StatesSedC$POCR[nRecs2Save]) / LakeArea
  MassBalanceC$BeginSedC = (StatesSedC$POCL[FirstDay] + StatesSedC$POCR[FirstDay]) / LakeArea
  
  MassBalanceC$dOCwc = MassBalanceC$EndC - MassBalanceC$BeginC
  MassBalanceC$dOCsed = MassBalanceC$EndSedC  - MassBalanceC$BeginSedC
  
  MassBalanceC$Burial = MassBalanceC$OCBurial
  MassBalanceC$OCRetention = 1 - (MassBalanceC$OCExport / (MassBalanceC$Allocthony+MassBalanceC$Autochthony))
  
  ##################################
  # Calculate Carbon characteristics
  iStratified = which(LakeVolumes$HypoVol>0)
  
  # Only use days when the lake is stratified
  EpiDOCMixedmin = min(StatesEpiC$DOCR[iMixed]/LakeVolumes$EpiVol[iMixed]+StatesEpiC$DOCL[iMixed]/LakeVolumes$EpiVol[iMixed])
  EpiDOCMixedmean = mean(StatesEpiC$DOCR[iMixed]/LakeVolumes$EpiVol[iMixed]+StatesEpiC$DOCL[iMixed]/LakeVolumes$EpiVol[iMixed])
  EpiDOCMixedmax = max(StatesEpiC$DOCR[iMixed]/LakeVolumes$EpiVol[iMixed]+StatesEpiC$DOCL[iMixed]/LakeVolumes$EpiVol[iMixed])
  
  EpiPOCMixedmin = min(StatesEpiC$POCR[iMixed]/LakeVolumes$EpiVol[iMixed]+StatesEpiC$POCL[iMixed]/LakeVolumes$EpiVol[iMixed])
  EpiPOCMixedmean = mean(StatesEpiC$POCR[iMixed]/LakeVolumes$EpiVol[iMixed]+StatesEpiC$POCL[iMixed]/LakeVolumes$EpiVol[iMixed])
  EpiPOCMixedmax = max(StatesEpiC$POCR[iMixed]/LakeVolumes$EpiVol[iMixed]+StatesEpiC$POCL[iMixed]/LakeVolumes$EpiVol[iMixed])
  
  EpiDOCStratifiedmin = min(StatesEpiC$DOCR[iStratified]/LakeVolumes$EpiVol[iStratified]+StatesEpiC$DOCL[iStratified]/LakeVolumes$EpiVol[iStratified])
  EpiDOCStratifiedmean = mean(StatesEpiC$DOCR[iStratified]/LakeVolumes$EpiVol[iStratified]+StatesEpiC$DOCL[iStratified]/LakeVolumes$EpiVol[iStratified])
  EpiDOCStratifiedmax = max(StatesEpiC$DOCR[iStratified]/LakeVolumes$EpiVol[iStratified]+StatesEpiC$DOCL[iStratified]/LakeVolumes$EpiVol[iStratified])
  
  EpiPOCStratifiedmin = min(StatesEpiC$POCR[iStratified]/LakeVolumes$EpiVol[iStratified]+StatesEpiC$POCL[iStratified]/LakeVolumes$EpiVol[iStratified])
  EpiPOCStratifiedmean = mean(StatesEpiC$POCR[iStratified]/LakeVolumes$EpiVol[iStratified]+StatesEpiC$POCL[iStratified]/LakeVolumes$EpiVol[iStratified])
  EpiPOCStratifiedmax = max(StatesEpiC$POCR[iStratified]/LakeVolumes$EpiVol[iStratified]+StatesEpiC$POCL[iStratified]/LakeVolumes$EpiVol[iStratified])
  
  HypoDOCStratifiedmin = min(StatesHypoC$DOCR[iStratified]/LakeVolumes$HypoVol[iStratified]+StatesHypoC$DOCL[iStratified]/LakeVolumes$HypoVol[iStratified])
  HypoDOCStratifiedmean = mean(StatesHypoC$DOCR[iStratified]/LakeVolumes$HypoVol[iStratified]+StatesHypoC$DOCL[iStratified]/LakeVolumes$HypoVol[iStratified])
  HypoDOCStratifiedmax = max(StatesHypoC$DOCR[iStratified]/LakeVolumes$HypoVol[iStratified]+StatesHypoC$DOCL[iStratified]/LakeVolumes$HypoVol[iStratified])
  
  HypoPOCStratifiedmin = min(StatesHypoC$POCR[iStratified]/LakeVolumes$HypoVol[iStratified]+StatesHypoC$POCL[iStratified]/LakeVolumes$HypoVol[iStratified])
  HypoPOCStratifiedmean = mean(StatesHypoC$POCR[iStratified]/LakeVolumes$HypoVol[iStratified]+StatesHypoC$POCL[iStratified]/LakeVolumes$HypoVol[iStratified])
  HypoPOCStratifiedmax = max(StatesHypoC$POCR[iStratified]/LakeVolumes$HypoVol[iStratified]+StatesHypoC$POCL[iStratified]/LakeVolumes$HypoVol[iStratified])
  
  SedPOCRmin = min(StatesSedC$POCR/(LakeArea*myLake$SedAreaProp))
  SedPOCRmean = mean(StatesSedC$POCR/(LakeArea*myLake$SedAreaProp))
  SedPOCRmax = max(StatesSedC$POCR/(LakeArea*myLake$SedAreaProp))
  
  SedPOCLmin = min(StatesSedC$POCL/(LakeArea*myLake$SedAreaProp))
  SedPOCLmean = mean(StatesSedC$POCL/(LakeArea*myLake$SedAreaProp))
  SedPOCLmax = max(StatesSedC$POCL/(LakeArea*myLake$SedAreaProp))
  
  SedPOCmin = min(StatesSedC$POCR/(LakeArea*myLake$SedAreaProp) + StatesSedC$POCL/(LakeArea*myLake$SedAreaProp))
  SedPOCmean = mean(StatesSedC$POCR/(LakeArea*myLake$SedAreaProp) + StatesSedC$POCL/(LakeArea*myLake$SedAreaProp))
  SedPOCmax = max(StatesSedC$POCR/(LakeArea*myLake$SedAreaProp) + StatesSedC$POCL/(LakeArea*myLake$SedAreaProp))
  
  CarbonSummary = data.frame(EpiDOCMixedmin=EpiDOCMixedmin,EpiDOCMixedmean=EpiDOCMixedmean,EpiDOCMixedmax=EpiDOCMixedmax,
                             EpiPOCMixedmin=EpiPOCMixedmin,EpiPOCMixedmean=EpiPOCMixedmean,EpiPOCMixedmax=EpiPOCMixedmax,
                             EpiDOCStratifiedmin=EpiDOCStratifiedmin,EpiDOCStratifiedmean=EpiDOCStratifiedmean,EpiDOCStratifiedmax=EpiDOCStratifiedmax,
                             EpiPOCStratifiedmin=EpiPOCStratifiedmin,EpiPOCStratifiedmean=EpiPOCStratifiedmean,EpiPOCStratifiedmax=EpiPOCStratifiedmax,
                             HypoDOCStratifiedmin=HypoDOCStratifiedmin,HypoDOCStratifiedmean=HypoDOCStratifiedmean,HypoDOCStratifiedmax=HypoDOCStratifiedmax,
                             HypoPOCStratifiedmin=HypoPOCStratifiedmin,HypoPOCStratifiedmean=HypoPOCStratifiedmean,HypoPOCStratifiedmax=HypoPOCStratifiedmax,
                             SedPOCRmin=SedPOCRmin,SedPOCRmean=SedPOCRmean,SedPOCRmax=SedPOCRmax,
                             SedPOCLmin=SedPOCLmin,SedPOCLmean=SedPOCLmean,SedPOCLmax=SedPOCLmax,
                             SedPOCmin=SedPOCmin,SedPOCmean=SedPOCmean,SedPOCmax=SedPOCmax)
  
  ##################################
  # Calculate Anoxia characteristics
  
  # Only use days when the lake is stratified
  iStratified = which(LakeVolumes$HypoVol>0)
  EpiDOmean = mean(StatesEpiGas$DO[iStratified]/mean(LakeVolumes$EpiVol[iStratified]))
  EpiDOmeanSatPercent = mean(StatesEpiGas$DO[iStratified]/mean(LakeVolumes$EpiVol[iStratified]) /
                               StatesEpiGas$DOsat[iStratified]) * 100
  
  EpiDOmax = max(StatesEpiGas$DO[iStratified])
  # Get the index of the first occurence of onset AFTER stratification begins
  EpiDOmaxOnset = min(which(StatesEpiGas$DO[iStratified]==EpiDOmax))
  # Adjust to concentration
  EV = LakeVolumes$EpiVol[iStratified[EpiDOmaxOnset]]
  EpiDOmax = EpiDOmax/EV
  # Use that index to get the actual day of year
  EpiDOmaxSatPercent = (StatesEpiGas$DO[iStratified[EpiDOmaxOnset]]/EV) /
    StatesEpiGas$DOsat[iStratified[EpiDOmaxOnset]] * 100
  
  HypoDOmin = min(StatesHypoGas$DO[iStratified])
  HypoDOAnoxicDays = length(which(StatesHypoGas$DO[iStratified] < AnoxicThreshold))
  # Get the index of the first occurence of onset AFTER stratification begins
  HypoDOminOnset = min(which(StatesHypoGas$DO[iStratified]==HypoDOmin))
  # Adjust to concentration
  HV = LakeVolumes$HypoVol[iStratified[HypoDOminOnset]]
  HypoDOmin = HypoDOmin/HV
  # Use that index to get the actual day of year
  HypoDOminSatPercent = (StatesHypoGas$DO[iStratified[HypoDOminOnset]]/HV) /
    StatesHypoGas$DOsat[iStratified[HypoDOminOnset]] * 100
  
  DOSummary = data.frame(EpiDOmean=EpiDOmean,EpiDOmeanSatPercent,
                         EpiDOmax=EpiDOmax,EpiDOmaxOnset=EpiDOmaxOnset,EpiDOmaxSatPercent=EpiDOmaxSatPercent,
                         HypoDOmin=HypoDOmin,HypoDOminOnset,HypoDOminSatPercent,HypoDOAnoxicDays)
  
  ##################################
  # Calculate Clarity characteristics
  
  # Only use days when the lake is stratified
  Secchimax = max(StatesEpiC$Secchi)
  Secchimin = min(StatesEpiC$Secchi)
  Secchimean = mean(StatesEpiC$Secchi)

  SecchiSummary = data.frame(Secchimean=Secchimean,Secchimin=Secchimin,Secchimax=Secchimax)
  
  ##################################
  # Print summary
  SigDig=3
  if (Print){
    print('================================================')
    print('OC budget (gC/m2LA/y)')
    print(paste('    Allochthony: ',round(MassBalanceC$Allocthony),sep=""))
    print(paste('    Autochthony: ',round(MassBalanceC$Autochthony),sep=""))
    print(paste('         Export: ',round(MassBalanceC$OCExport),sep=""))
    print(paste(' WC Respiration: ',round(MassBalanceC$OCR),sep=""))
    print(paste('Sed Respiration: ',round(MassBalanceC$OCSedimentR),sep=""))
    print(paste('         Burial: ',round(MassBalanceC$OCBurial),sep=""))
    print('---------------------')
    print(paste('Balance (dOC = All+Aut-Rwc-Exp-Rsed-Burial) '))
    print(paste(round(MassBalanceC$dOCwc+MassBalanceC$dOCsed),' = ',round(MassBalanceC$OCBalance),sep=""))
    print(paste('dOC in water col: ',round(MassBalanceC$dOCwc),sep=""))
    print(paste('dOC in sediments: ',round(MassBalanceC$dOCsed),sep=""))
    print(paste('Net OC burial (Burial+dSed): ',round(MassBalanceC$OCBurial+MassBalanceC$dOCsed)))
    print(paste('LakeOC retention: ',round(MassBalanceC$OCRetention,digits=SigDig),sep=""))
    MassBalanceC$OCRetention
  }
  
  if (Print){
    print('================================================')
    print('P budget (gP/m2LA/y)')
    print(paste('           Load: ',round(MassBalanceP$PLoad,digits=SigDig),sep=""))
    print(paste('         Export: ',round(MassBalanceP$PExport,digits=SigDig),sep=""))
    print(paste('         Burial: ',round((MassBalanceP$PBurial+MassBalanceP$PBBurial),digits=SigDig),sep=""))
    print('---------------------')
    print(paste('Balance (dP = Load-Export-Burial) '))
    print(paste(round(MassBalanceP$dPwc+MassBalanceP$dPsed,digits=SigDig),' = ',round(MassBalanceP$PBalance,digits=SigDig),sep=""))
    print(paste('dP in water col: ',round(MassBalanceP$dPwc,digits=SigDig),sep=""))
    print(paste('dP in sediments: ',round(MassBalanceP$dPsed,digits=SigDig),sep=""))
    print(paste('Net P burial (Burial+dSed): ',round(MassBalanceP$PBurial+MassBalanceP$PBBurial+MassBalanceP$dPsed,digits=SigDig)))
    print(paste('       Settling: ',round(MassBalanceP$PSettlingRebind,digits=SigDig),sep=""))
    print(paste('      Recycling: ',round(MassBalanceP$PRecyclingRelease,digits=SigDig),sep=""))
    print(paste('LakeP retention: ',round(MassBalanceP$PRetention,digits=SigDig),sep=""))
    #print(paste('Equilibrium soil deposition (mm/y): ',round(mySoilDepEquil,digits = SigDig),sep=""))
  }

  ##################################
  # Prepare the return structure
  MetabSummary = list()
  MetabSummary[[1]] = MassBalanceC
  MetabSummary[[2]] = MassBalanceP
  MetabSummary[[3]] = CarbonSummary
  MetabSummary[[4]] = PhosphorusSummary
  MetabSummary[[5]] = DOSummary
  MetabSummary[[6]] = SecchiSummary
  MetabSummary[[7]] = StatesEpiC
  MetabSummary[[8]] = StatesEpiP
  MetabSummary[[9]] = StatesSedC
  
  return(MetabSummary)
}
