WriteResultsMatrix = function(OutputFileListName){
  # Writes the results of a simulation to disk
  # Inputs
  #    CSV file containing the file list of the simulation output
  
  # Date range over which to make the comparison
  DateRange2Write = c("2010-01-01","2020-12-31")
  # Where to write the output\
  OutputFileName = "./Output/SimResultsMatrix.csv"
  PlotOutput = TRUE # Generates plots of output variables

  # Load the simulation results
  # e.g., OutputFileListName = "./Output/Mendota_1980_Present_1_1.csv"
  myFileList = read.csv(OutputFileListName)
  nFiles = dim(myFileList)[1] 
  # There will be only 1 file if only 1 lake is defined and number of repeats of the sim is 1

  # Setup easy-to-read structures to hold the data  
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
  
  # Cycle through the files and transfer simulation to above variables
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
    # The end of Metabolism.R has a key to the following list
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
  
  # Note that most variables are in total mass for the lake
  # Easier to understand if they are divided by lake area or by respective volumes

  # get simulation date range
  # SimDateRange = c(min(as.Date(ModelRun[[1]]$StartDate)),max(as.Date(ModelRun[[1]]$EndDate)))
  SimDateRange = c(min(as.Date(StartDates)),max(as.Date(EndDates)))
  # generate sequence of dates for, e.g., plotting
  SimDates = seq(as.Date(SimDateRange[1]),as.Date(SimDateRange[2]),by='day')
  # remove first value, as output begins on day 2
  # If more than one simulation is strung together, a day gets dropped for each sim
  nTotDates = length(SimDates)-(nFiles-1)
  SimDates = SimDates[2:nTotDates]
  
  # Get indeces of the desired date range to write
  iWr = which(as.Date(SimDates) >= as.Date(DateRange2Write[1]) 
        & as.Date(SimDates) <= as.Date(DateRange2Write[2]) )
  
  # Create the output data frame
  Output = data.frame(SimDate = SimDates[iWr],
                      EpiPOC_mgC_L  = (StatesEpiC$POCL[iWr]+StatesEpiC$POCR[iWr])/LakeVolumes$EpiVol[iWr],
                      EpiDOC_mgC_L  = (StatesEpiC$DOCL[iWr]+StatesEpiC$DOCR[iWr])/LakeVolumes$EpiVol[iWr],
                      EpiDO_mgO2_L  = (StatesEpiGas$DO[iWr])/LakeVolumes$EpiVol[iWr],
                      EpiDO_sat     = (StatesEpiGas$DOsat[iWr]),
                      EpiNPP_mgC_L  = (RatesEpiC$GPP[iWr])/LakeVolumes$EpiVol[iWr],
                      EpiR_mgC_L    = (RatesEpiC$DOCLR[iWr]+RatesEpiC$DOCRR[iWr])/LakeVolumes$EpiVol[iWr],
                      EpiT_C        = LakeTemps$EpiT[iWr],
                      EpiI_w_m2     = LakeTemps$MeanEpiIrradiance[iWr],
                      Secchi_m      = StatesEpiC$Secchi[iWr])
  
  # Write output to disc
  write_csv(Output,OutputFileName)
  
  if (PlotOutput){
    #################################################
    # par(mfrow=c(4,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
    par(mfrow=c(3,1),mai = c(0.55,0.75, 0.05, 0.05),oma = c(2,2,0.2,0.2), cex = 0.9)
    
    nCols = dim(Output)[2]
    for (i in 2:nCols){
      myYLab = colnames(Output)[i]
      plot(Output[,1],Output[,i],ylab=myYLab, type = 'l')
    }
  } # End plotting
  
  # End function
}

  
