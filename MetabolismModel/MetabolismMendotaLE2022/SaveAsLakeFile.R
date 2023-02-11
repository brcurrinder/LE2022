SaveAsLakeFile = function(InputFile,OutputFile,LakeHypsometry){
  
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
  # MetabSummary[[7]] = StatesEpiC
  # MetabSummary[[8]] = StatesEpiP
  
  load(InputFile)
  
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
  
  StatesEpiC = MS[[7]]
  StatesEpiP = MS[[8]]
  StatesSedC = MS[[9]]
  
  Output = myLake
  # Update states
  nRecs = length(StatesEpiC$DOCR)
  LakeVolume = max(LakeHypsometry$area) * myLake$Zmean
  Output$DOCRinit = round(StatesEpiC$DOCR[nRecs] / LakeVolume,digits=3)
  Output$DOCLinit = round(StatesEpiC$DOCL[nRecs] / LakeVolume,digits=3)
  Output$POCRinit = round(StatesEpiC$POCR[nRecs] / LakeVolume,digits=3)
  Output$POCLinit = round(StatesEpiC$POCL[nRecs] / LakeVolume,digits=3)
  Output$SedPOCRinitProp = round(StatesSedC$POCR[nRecs] / (StatesSedC$POCR[nRecs] + StatesSedC$POCL[nRecs]),digits=3)
  Output$SedPOCLinitProp = 1-Output$SedPOCRinitProp
  
  
  gSedAreal = myLake$SedBulkDensity * myLake$ActiveSedDepth
  Output$SedAvailOC = round(CarbonSummary$SedPOCmean / gSedAreal * 1000,digits=3) # mgC/gSed
  Output$Pinit = round(StatesEpiP$TP[nRecs]/LakeVolume,digits=3)
  Output$SedAvailPinit = round(PhosphorusSummary$SedPmean / gSedAreal * 1000,digits=3) # mgP/gSed

  OutputTranspose <- as.data.frame(as.matrix(t(Output)))
  
  write.csv(OutputTranspose,OutputFile,quote=FALSE)
  
  return(Output)
  
}