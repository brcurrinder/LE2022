ComparePredictionsObservations = function(OutputFileListName){
  # ComparePredictionsObservations
  # For now, use just the last file in the file list that comes from OutputFileListName
  # Inputs
  #    CSV file containing the file list of the simulation output
  #    Date range over which to make the comparison
  
  # Load Observational data
  # myObs = read.csv("./ObservationalData/Mendota_observations_v1.csv")
  # myObsSecchi = read.csv("./ObservationalData/Mendota_secchi.csv")
  myObs = read.csv("./ObservationalData/ME_observations_updated.csv")
  myObsSecchi = read.csv("./ObservationalData/Mendota_secchi.csv")
  
  # Load the simulation results
  # OutputFileListName = "./Output/Mendota_1980_Present_1_1.csv"
  myFileList = read.csv(OutputFileListName)
  nFiles = dim(myFileList)[1]
  
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
  LakeTemps     = NULL
  nDays = 0
  StartDates    = NULL
  EndDates      = NULL
  
  # Use only the last simulation
  for (i in 1:nFiles){
    #load(OutputFileNames$Names[i])
    load(myFileList$Names[i])
    
    #Load time series results
    MR = ModelRun[[7]]
    Hypso = ModelRun[[4]]
    
    # Lake area
    LakeArea = max(Hypso$area)
    nDays = nDays + ModelRun[[1]]$nDays[1]

    # Extract data frames from metabolism model
    StartDates    = rbind(StartDates,ModelRun[[1]]$StartDate)
    EndDates      = rbind(EndDates,ModelRun[[1]]$EndDate)
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
    LakeTemps     = rbind(LakeTemps,MR[[19]])
    
    LakeAreas = data.frame(TotalArea = max(LakeHypsometry))
    
  }

  # get simulation date range
  # SimDateRange = c(min(as.Date(ModelRun[[1]]$StartDate)),max(as.Date(ModelRun[[1]]$EndDate)))
  SimDateRange = c(min(as.Date(StartDates)),max(as.Date(EndDates)))
  # generate sequence of dates for, e.g., plotting
  SimDates = seq(as.Date(SimDateRange[1]),as.Date(SimDateRange[2]),by='day')
  # remove first value, as output begins on day 2
  # If more than one simulation is strung together, a day gets dropped for each sim
  nTotDates = length(SimDates)-(nFiles-1)
  SimDates = SimDates[2:nTotDates]
  
  #################################################
  # For each observation and simulation combination
  # par(mfrow=c(4,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  par(mfrow=c(2,1),mai = c(0.05,0.05, 0.05, 0.05),oma = c(2,2,0.2,0.2), cex = 0.5)
  FileName = paste('./Figures/DOCalibration.png',sep="")
  
  #################
  # Epilimnetic DO
  myLayer = 'epi'
  SimVector = StatesEpiGas$DO/LakeVolumes$EpiVol
  OtherVector = StatesEpiGas$DOsat
  ObsVector = myObs$o2
  ObsDates = myObs$sampledate
  ObsLayerVector = myObs$layer
  stratFlagVector = myObs$stratFlag
  myYLim = NULL # if null, will calculate from the data
  myYLabel = 'Epi DO (mg/L)'
  # myLegendTxt = c('DO Model','DO obs','DO sat')
  myLegendTxt = NA
  PlotXLab = FALSE
  myLegendLines = c(1,NA,2)
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','grey')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)

  #################
  # Hypo DO
  myLayer = 'hypo'
  SimVector = StatesHypoGas$DO/LakeVolumes$HypoVol
  OtherVector = NULL
  ObsVector = myObs$o2
  ObsDates = myObs$sampledate
  ObsLayerVector = myObs$layer
  stratFlagVector = myObs$stratFlag
  myYLim = NULL # if null, will calculate from the data
  myYLabel = 'Hypo DO (mg/L)'
  myLegendTxt = c('DO Model','DO obs',NULL)
  myLegendTxt = NA
  PlotXLab = TRUE
  myLegendLines = c(1,NA,NA)
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','black')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)
  dev.copy(png,FileName,width=2,height=2.5,units="in",res=300)
  dev.off()
  
  #################
  # Epilimnetic TP
  par(mfrow=c(2,1),mai = c(0.05,0.05, 0.05, 0.05),oma = c(2,2,0.2,0.2), cex = 0.5)
  FileName = paste('./Figures/TPCalibration.png',sep="")
  myLayer = 'epi'
  SimVector = StatesEpiP$TP/LakeVolumes$EpiVol*1000
  OtherVector = NULL
  ObsVector = myObs$tp
  ObsDates = myObs$sampledate
  ObsLayerVector = myObs$layer
  stratFlagVector = myObs$stratFlag
  myYLim = NULL # if null, will calculate from the data
  myYLim = c(0,200)
  myYLabel = 'Epi TP (ug/L)'
  myLegendTxt = NA
  PlotXLab = FALSE
  #myLegendTxt = c('TP Model','TP obs',NULL)
  myLegendLines = c(1,NA,NA)
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','black')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)
  #################
  # Hypo TP
  myLayer = 'hypo'
  SimVector = StatesHypoP$TP/LakeVolumes$HypoVol*1000
  OtherVector = NULL
  ObsVector = myObs$tp
  ObsDates = myObs$sampledate
  ObsLayerVector = myObs$layer
  stratFlagVector = myObs$stratFlag
  myYLim = NULL # if null, will calculate from the data
  myYLabel = 'Hypo TP (ug/L)'
  myLegendTxt = NA
  #myLegendTxt = c('TP Model','TP obs',NULL)
  PlotXLab = TRUE
  myLegendLines = c(1,NA,NA)
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','black')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)
  dev.copy(png,FileName,width=2,height=2.5,units="in",res=300)
  dev.off()
  
  #################
  # Epi Temp 
  par(mfrow=c(3,1),mai = c(0.05,0.05, 0.05, 0.05),oma = c(2,2,0.2,0.2), cex = 0.5)
  FileName = paste('./Figures/TempCalibration.png',sep="")
  
  myLayer = 'epi'
  SimVector = LakeTemps$EpiT
  OtherVector = NULL
  ObsVector = myObs$temp
  ObsDates = myObs$sampledate
  ObsLayerVector = myObs$layer
  stratFlagVector = myObs$stratFlag
  myYLim = NULL # if null, will calculate from the data
  myYLabel = 'Epi T (C)'
  myLegendTxt = NA
  #myLegendTxt = c('Temp Model','Temp obs',NULL)
  PlotXLab = FALSE
  myLegendLines = c(1,NA,NA)
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','black')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)
  #################
  # Hypo Temp 
  myLayer = 'hypo'
  SimVector = LakeTemps$HypoT
  OtherVector = NULL
  ObsVector = myObs$temp
  ObsDates = myObs$sampledate
  ObsLayerVector = myObs$layer
  stratFlagVector = myObs$stratFlag
  myYLim = NULL # if null, will calculate from the data
  myYLabel = 'Hypo T (C)'
  myLegendTxt = NA
  #myLegendTxt = c('Temp Model','Temp obs',NULL)
  myLegendLines = c(1,NA,NA)
  PlotXLab = FALSE
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','black')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)
  
  #################
  # Secchi depth 
  myLayer = 'all'
  SimVector = StatesEpiC$Secchi
  OtherVector = NULL
  ObsVector = myObsSecchi$secnview
  ObsDates = myObsSecchi$sampledate
  ObsLayerVector = rep('all',length(ObsVector))
  myYLim = NULL # if null, will calculate from the data
  myYLabel = 'Secchi depth (m)'
  myLegendTxt = NA
  #myLegendTxt = c('Secchi Model','Secchi obs',NULL)
  PlotXLab = TRUE
  myLegendLines = c(1,NA,NA)
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','black')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)
  
  dev.copy(png,FileName,width=2,height=2.5,units="in",res=300)
  dev.off()
  
  par(mfrow=c(2,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  #################
  # Epi DOC 
  myLayer = 'epi'
  SimVector = (StatesEpiC$DOCR + StatesEpiC$DOCL)/LakeVolumes$EpiVol
  OtherVector = NULL
  ObsVector = myObs$doc
  ObsDates = myObs$sampledate
  ObsLayerVector = myObs$layer
  myYLim = NULL # if null, will calculate from the data
  myYLabel = 'Epi DOC (g/m3)'
  myLegendTxt = NA
  #myLegendTxt = c('DOC Model','DOC obs',NULL)
  myLegendLines = c(1,NA,NA)
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','black')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)
  #################
  # Hypo DOC 
  myLayer = 'hypo'
  SimVector = (StatesHypoC$DOCR + StatesHypoC$DOCL)/LakeVolumes$HypoVol
  OtherVector = NULL
  ObsVector = myObs$doc
  ObsDates = myObs$sampledate
  ObsLayerVector = myObs$layer
  myYLim = NULL # if null, will calculate from the data
  myYLabel = 'Hypo DOC (g/m3)'
  myLegendTxt = NA
  #myLegendTxt = c('DOC Model','DOC obs',NULL)
  myLegendLines = c(1,NA,NA)
  myPCH = c(NA,1,NA)
  myColSeq = c('black','red','black')
  PlotComparison(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,myLegendTxt,
                 myLegendLines,myPCH,myColSeq,PlotXLab)

}

PlotComparison = function(myLayer,SimDateRange,SimDates,SimVector,OtherVector,ObsVector,stratFlagVector,ObsDates,ObsLayerVector,myYLim,myYLabel,
                          myLegendTxt,myLegendLines,myPCH,myColSeq,PlotXLab){
  # Reduce Observational data range to the simulation date range, just for easier plotting
  # if hypo, get just hypo, but if epi, get epi while stratified and all while not
  if (myLayer=='hypo' | myLayer=='all'){
    iSimDates = which(as.Date(ObsDates) >= as.Date(SimDateRange[1]) & 
                        as.Date(ObsDates) <= as.Date(SimDateRange[2]) &
                        ObsLayerVector == myLayer)
  }else{
    iSimDates = which(as.Date(ObsDates) >= as.Date(SimDateRange[1]) & 
                        as.Date(ObsDates) <= as.Date(SimDateRange[2]) &
                        (ObsLayerVector == myLayer | (ObsLayerVector == 'all' & stratFlagVector=='FALSE')))
  }
  ObsVectorSubset = ObsVector[iSimDates]
  ObsDatesSubset = ObsDates[iSimDates]
  # Get the intersection of dates between simulation and observations
  IntersectDates = as.Date(lubridate::intersect(as.Date(SimDates),as.Date(ObsDatesSubset)),origin="1970-01-01")
  # Get these records for observation data
  iIntersectObs = which(as.Date(ObsDatesSubset)==as.Date(IntersectDates))
  # Get those records for simulation data
  iIntersectSim = which(as.Date(SimDates) %in% as.Date(IntersectDates))
  
  # Calculate the residuals
  Residuals = ObsVectorSubset[iIntersectObs] - SimVector[iIntersectSim]
  
  # Plot it
  if (is.null(myYLim)){
    myYLim = range(ObsVectorSubset,SimVector,na.rm=TRUE)
  }
  if(PlotXLab){
    plot(as.POSIXct(SimDates),SimVector,ylim=myYLim,type='l',ylab=myYLabel,col=myColSeq[1])
    
  }else{
    plot(as.POSIXct(SimDates),SimVector,ylim=myYLim,type='l',ylab=myYLabel,col=myColSeq[1],xaxt='n')
  }
  
  if (!is.null(OtherVector)){
    lines(as.POSIXct(SimDates),OtherVector,lty=2,col=myColSeq[3])
  }
  points(as.POSIXct(ObsDatesSubset),ObsVectorSubset,col=myColSeq[2])
  if (!is.na(myLegendTxt)){
    legend('topleft',myLegendTxt,lty=myLegendLines,pch=myPCH,col=myColSeq)
  }
  PlotResids = FALSE
  if (PlotResids){
    # Plot residuals
    plot(as.POSIXct(IntersectDates),Residuals)
    abline(h=0,lty=2)
  }
}
