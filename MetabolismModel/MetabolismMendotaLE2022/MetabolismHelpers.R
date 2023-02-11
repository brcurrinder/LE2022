# Help functions for metabolism

GetTrophicThresholds <- function(FileName=NULL){
  # if Whichparameters = NA, then return them all
  # Assume parameter names are the first column
  FDCol = 3 # First Data Column
  # Open the lake file
  if (is.null(FileName)){
    FileName = './Config/trophicThresholds.csv'
  }
  # Read into a structure
  Thresholds = read.csv(FileName,header = TRUE, stringsAsFactors = FALSE)
  return(Thresholds)
}

GetParameters <- function(WhichParams){
  # if Whichparameters = NA, then return them all
  # Assume parameter names are the first column
  FDCol = 4 # First Data Column
  # Open the lake file
  ParamFile = './Config/Parameters.csv'
  # Read into a structure
  AllParamsets = read.csv(ParamFile,header = TRUE, stringsAsFactors = FALSE)
  nVars = length(AllParamsets$Parameter)
  nParams = dim(AllParamsets)[2]-3 # First three cols are not data
  # Build the data frame
  dfText = "myParamset = data.frame("
  for (i in 1:nVars){
    if (i<nVars){
      dfText = paste(dfText,AllParamsets$Parameter[i],'=rep(NA,',nParams,'),',sep="")
    }else{
      dfText = paste(dfText,AllParamsets$Parameter[i],'=rep(NA,',nParams,'))',sep="")
    }
  }
  # Create the data frame
  eval(parse(text=dfText))
  
  # Load the data into the DF
  for (j in 4:dim(AllParamsets)[2]){ # Across the number of Params
    for (i in 1:length(AllParamsets$Parameter)){ # and the variables
      #print(AllParams$Parameter[i])
      if (i==1){
        eval(parse(text = paste("myParamset$",AllParamsets$Parameter[i],"[j-3]=",'"',AllParamsets[i,j],'"', sep="")))
      }else{
        eval(parse(text = paste("myParamset$",AllParamsets$Parameter[i],"[j-3]=",AllParamsets[i,j], sep="")))
      }
    }
  }
  if (is.na(WhichParams)){
    return(myParamset)
  }else{
    return(myParamset[WhichParams,])
  }
}

GetLakeConfig <- function(WhichLake){
  # if WhichLake = NA, then return them all
  # Assume parameter names are the first column
  FDCol = 5 # First Data Column
  # Open the lake file
  LakeFile = './Config/Lakes.csv'
  # Read into a structure
  AllLakes = read.csv(LakeFile,header = TRUE, stringsAsFactors = FALSE)
  nVars = length(AllLakes$Variable)
  nLakes = dim(AllLakes)[2]-(FDCol-1) # First three cols are not data
  # Build the data frame
  dfText = "myLake = data.frame("
  for (i in 1:nVars){
    if (i<nVars){
      dfText = paste(dfText,AllLakes$Variable[i],'=rep(NA,',nLakes,'),',sep="")
    }else{
      dfText = paste(dfText,AllLakes$Variable[i],'=rep(NA,',nLakes,'))',sep="")
    }
  }
  # Create the data frame
  eval(parse(text=dfText))
  
  # Load the data into the DF
  for (j in FDCol:dim(AllLakes)[2]){ # Across the number of lakes
    for (i in 1:length(AllLakes$Variable)){ # and the variables
      #print(AllLakes$Variable[i])
      if (i==1){
        eval(parse(text = paste("myLake$",AllLakes$Variable[i],"[j-4]=",'"',AllLakes[i,j],'"', sep="")))
      }else{
        eval(parse(text = paste("myLake$",AllLakes$Variable[i],"[j-4]=",AllLakes[i,j], sep="")))
      }
    }
  }
  if (is.na(WhichLake)){
    return(myLake)
  }else{
    return(myLake[WhichLake,])
  }
}

GetRandomLakes = function(BaseLakes,Lakes2Sample,nSamples){
  # BaseLakes is same structure as GetLakes returns
  # Lakes2Sample identifies which of the base lakes to sample
  # nSamples is the number of samples per lake
  nBaseLakes = dim(BaseLakes)[2]
  NewLakes = BaseLakes[0,]
  rownames(NewLakes)=c()
  for (i in 1:length(Lakes2Sample)){
    iThisLake = Lakes2Sample[i]
    ThisLakes = BaseLakes[iThisLake,]
    for (j in 1:nSamples){
      # Replicate the origina lake
      ThisLakes[j+1,] = ThisLakes[1,]
      # Change the name of the present lake
      ThisLakes$LakeName[j+1] = paste(ThisLakes$LakeName[j+1],"_rand_",j,sep="")
      # Sample lake characteristics
      
    } # End samples
    rownames(ThisLakes) = c()
    NewLakes = rbind(NewLakes,ThisLakes)
  } # End lakes
  
  return(NewLakes)
}

GetLakeThermal <- function(nDays,Area,Zmax,Zmean,Ztherm,StratBegin,StratEnd,Irradiance,MinTemp,MaxEpiTemp,MaxHypoTemp){
  # Setup lake volume time series
  # Loop through time, assuming 1 is first day of the year
  # MinTemp = 4
  # MaxEpiTemp = 20
  # MaxHypoTemp = 12
  source('./LakeShapeAreaVolume.R')
  
  # Use irradiance pattern to control lake temperature pattern
  IMix = min(c(Irradiance$SopenMean[StratBegin],Irradiance$SopenMean[StratEnd]))
  IMax = max(Irradiance$SopenMean)
  
  # Get lake shape, based on Carpenter 1983 lake shape algorithm
  myLayers = LakeShapeAreaVolume(Zmax,Zmean,Ztherm,FALSE)
  
  LakeVolumes = data.frame(DOY=rep(NA,nDays),EpiVol=rep(NA,nDays),HypoVol=rep(NA,nDays),TotalVol=rep(NA,nDays))
  LakeAreas = data.frame(DOY=rep(NA,nDays),EpiArea=rep(NA,nDays),HypoArea=rep(NA,nDays),TotalArea=rep(NA,nDays))
  LakeTemps = data.frame(DOY=rep(NA,nDays),EpiT=rep(NA,nDays),HypoT=rep(NA,nDays))
  DOY = 0
  for (i in 1:nDays){
    DOY = DOY + 1
    # day of year between 1 and 365
    if (DOY==366){
      DOY = 1
    }
    IToday = Irradiance$SopenMean[DOY]
    IScaler = (IToday - IMix) / (IMax - IMix)
    if (IScaler < 0){
      IScaler = 0
    }
    if (DOY < StratBegin | DOY > StratEnd){
      LakeVolumes$EpiVol[i] = Area * Zmean
      LakeAreas$EpiArea[i] = Area
      LakeVolumes$HypoVol[i] = 0
      LakeAreas$HypoArea[i] = 0
      LakeTemps$EpiT[i] = MinTemp 
      LakeTemps$HypoT[i] = MinTemp 
    }else{
      LakeVolumes$EpiVol[i] = Area * Zmean * myLayers$VolumeAboveZ
      LakeAreas$EpiArea[i] = Area * myLayers$AreaAboveZ
      LakeVolumes$HypoVol[i] = Area * Zmean * (1-myLayers$VolumeAboveZ)
      LakeAreas$HypoArea[i] = Area * (1-myLayers$AreaAboveZ)
      LakeTemps$EpiT[i] = MinTemp + (MaxEpiTemp - MinTemp) * IScaler
      LakeTemps$HypoT[i] = MinTemp + (MaxHypoTemp - MinTemp) * IScaler
    }
    LakeVolumes$TotalVol[i] = LakeVolumes$EpiVol[i] + LakeVolumes$HypoVol[i]
    LakeAreas$TotalArea[i] = LakeAreas$EpiArea[i] + LakeAreas$HypoArea[i]
    LakeVolumes$DOY[i] = DOY
    LakeAreas$DOY[i] = DOY
    LakeTemps$DOY[i] = DOY
  }
  # plot(1:nDays,LakeVolumes$EpiVol,type='l',ylim=c(10e5,3e7))
  # lines(1:nDays,LakeVolumes$HypoVol,type='l',col='red')
  # lines(1:nDays,LakeVolumes$TotalVol,type='l',col='green')
  LakeVolAreaTemp = list()
  LakeVolAreaTemp[[1]] = LakeVolumes
  LakeVolAreaTemp[[2]] = LakeAreas
  LakeVolAreaTemp[[3]] = LakeTemps
  return(LakeVolAreaTemp)
}

GetIrradiance <- function(myLat,myLon,myElevation,IceOn,IceOff){
  # install.packages('solrad',repos="https://cloud.r-project.org")
  library(solrad)

  # IDOY <- seq(100,101,0.25) #third value controls how many points per day are sampled.
  IDOY <- seq(1,366,0.05) #third value controls how many points per day are sampled.
  Sdirect <- DirectRadiation(IDOY, Lat = myLat, Lon=myLon, SLon=myLon, DS=0, Elevation = myElevation, Slope = 5, Aspect = 10) #surface radiation
  Sopen <- OpenRadiation(IDOY, Lat = myLat, Lon=myLon, SLon=myLon, DS=0, Elevation = myElevation) #open air radiation
  # plot(IDOY, Sopen,,ylim=c(0,1000),type='l')
  # lines(IDOY, Sdirect,col='red')
  # Get daily averages
  Irradiance = data.frame(DOY=seq(1,365),SopenMean=rep(NA,365),SopenMax=rep(NA,365))
  for (i in 1:365){
    iToday = which(floor(IDOY)==i)
    Irradiance$SopenMean[i] = mean(Sopen[iToday])
    Irradiance$SopenMax[i] = max(Sopen[iToday])
  }
  # plot(Irradiance$DOY,Irradiance$SopenMax,ylim=c(0,1000),type='l')
  # lines(Irradiance$DOY,Irradiance$SopenMean,col='red')
  # find and replace any NaNs
  iBad = which(is.nan(Irradiance$SopenMean))
  if (length(iBad)>0){
    iGood = which(!is.nan(Irradiance$SopenMean))
    xout = 1:365
    myApprox = approx(Irradiance$DOY[iGood],Irradiance$SopenMean[iGood],xout)
    Irradiance$SopenMean = myApprox$y
    iGood = which(!is.nan(Irradiance$SopenMax))
    myApprox = approx(iGood,Irradiance$SopenMax[iGood],xout)
    Irradiance$SopenMax = myApprox$y
  }
  # Reduce light for ice covered days
  if (!is.na(IceOff) & !is.na(IceOn)){
    iNoLight = c(1:IceOff,IceOn:365)
    Irradiance$SopenMean[iNoLight] = 0
    Irradiance$SopenMax[iNoLight] = 0
  }
  return(Irradiance)
}

WriteResults2Disk = function(Output,nYears2Save,OutputFile,Append){
  # Write results to disk
  nYears2Save
  nRecs2Save = nYears2Save*365
  OutputLength = dim(Output)[1]
  StartOutput = OutputLength-nRecs2Save+1
  write.csv(Output[StartOutput:OutputLength,],file=OutputFile)
}

VollenweiderPLoads <- function(myLake){
  # Calculate P load according to Vollenwieder
  # TPlake = Load / Zmean(1/RT + Sigma); Load = TPLake * (Zmean(1/RT + Sigma))
  #   e.g., Lake Mendota Load = 0.14 g/m3 (12.8 m (0.25 y-1 + 0.7 y-1)) = 1.22 g/m2LA/y
  PLoadVoll = myLake$Pinit * (myLake$Zmean*(1/(myLake$RT/365) + myParameters$PRetention))
  return(PLoadVoll)
}

HansonOCLoads <- function(myLake,myParameters,LakeVolumes,LakeTemps){
  # Estimate OC load according to equilibrium
  # Add a factor for POC load at the end
  # dOC/dt = Load - Outflow - R
  # Load (g/m2/y) = Outflow + R
  # Load = [DOC] * Qout + [DOC]/Vol * Cresp
  # Load = [DOC] * (Qout + Cresp/Vol) / LakeArea
  
  # Assuming that input/output happens through the epilimnion,
  #  calculate the weighted average for the epilimnion
  iStratDays = which(LakeVolumes$HypoVol>0)
  iMixedDays = which(LakeVolumes$HypoVol==0)
  StratDays = myLake$StratEnd - myLake$StratBegin + 1
  MixedDays = 365-StratDays
  StratDaysWeight = StratDays/365
  MixedDaysWeight = MixedDays/365
  DOC = myLake$DOCRinit
  LakeArea = myLake$Area
  LakeVolume = LakeArea * myLake$Zmean
  RT = myLake$RT/365
  Cresp = myParameters$RDOCR*365
  
  # Total load will be weighted average of the two
  # Mixed days loads
  Outflow = 1/(RT) * LakeVolume
  TAdjEpi = myParameters$RTheta^(mean(LakeTemps$EpiT[iMixedDays])-myParameters$Tbase)
  LoadMixed = DOC * (Outflow + Cresp*TAdjEpi*LakeVolume) / LakeArea
  
  # Stratified days load
  EpiVol = mean(LakeVolumes$EpiVol[iStratDays])
  # RTstratified = RT*(StratDays/365)*EpiVol/LakeVolume
  # RTstratified = RT
  # Outflow = 1/(RTstratified) * LakeVolume
  TAdjEpi = myParameters$RTheta^(mean(LakeTemps$EpiT[iStratDays])-myParameters$Tbase)
  LoadStratified = DOC * (Outflow + Cresp*TAdjEpi*EpiVol) / LakeArea
  
  Load = LoadMixed*MixedDaysWeight + LoadStratified*StratDaysWeight
  Load = Load / myLake$DOCRinflow
  return(Load)
}

LakePEquilibria <- function(myLake,PRetention){
  # This function calculates equilibrium values for a number of states/fluxes relevant to lake P cycling
  # 
  # Calculate what sediment deposition should be (mm/y), assuming load, Zsed, retention, P initial condition
  # SoilDepEquil (mm/y) = PLoad (gP/m2LA/y)/PropSed * Sigma * Zsed (m) * 1/Pinit (gP/m2) * 1000mm/m
  
  nCalcs = 10
  i=1
  
  OutT = data.frame(WRT=rep(NA,nCalcs),Zmean=rep(NA,nCalcs),SedProp=rep(NA,nCalcs),Zsed=rep(NA,nCalcs),SedDeposition=rep(NA,nCalcs),
                    PLoad=rep(NA,nCalcs),Pwc=rep(NA,nCalcs),PSed=rep(NA,nCalcs),PSedDens=rep(NA,nCalcs),BulkDens=rep(NA,nCalcs),PRetention=rep(NA,nCalcs))
  
  OutT$WRT[1] = myLake$RT/365
  OutT$Zmean[1] = myLake$Zmean
  OutT$SedProp[1] = myLake$SedAreaProp
  OutT$Zsed[1] = myLake$ActiveSedDepth
  OutT$SedDeposition[1] = myLake$SoilDeposition
  OutT$PLoad[1] = myLake$PLoad
  OutT$Pwc[1] = myLake$Pinit
  OutT$PSedDens[1] = myLake$SedAvailPinit
  OutT$BulkDens[1] = myLake$SedBulkDensity
  OutT$PRetention[1] = PRetention

  # Vollenwider retention
  i=i+1
  #PRetention = (OutT$PLoad[1] / OutT$Zmean[1]) / OutT$Pwc[1] - 1/OutT$WRT[1]
  #PRetentionEq = (PLoad/Pwc - 1) / WRT
  OutT$PLoad[i] = 1
  OutT$Zmean[i] = 1
  OutT$Pwc[i] = 1
  OutT$WRT[i] = 1
  
  # Vollenweider PwcEq
  i=i+1
  OutT$Pwc[i] = OutT$PLoad[1] / (OutT$Zmean[1]*(1/OutT$WRT[1]+OutT$PRetention[1])) # g/m3
  #myEquilstxt$PwcEqtxt = 'PwcEq = f(PLoad, Zmean, WRT, PRetention)'
  OutT$PRetention[i] = 1
  OutT$PLoad[i] = 1
  OutT$Zmean[i] = 1
  OutT$WRT[i] = 1
  
  # Vollenweider load
  i=i+1
  OutT$PLoad[i] = OutT$Pwc[1] * (OutT$Zmean[1]*(1/OutT$WRT[1]+OutT$PRetention[1])) # g/m2
  #myEquilstxt$PLoadEqtxt = 'PLoadEq = f(Pwc, Zmean, WRT, PRetention)'
  OutT$Pwc[i] = 1
  OutT$PRetention[i] = 1
  OutT$Zmean[i] = 1
  OutT$WRT[i] = 1
  
  # Sed deposition equilibrium estimated from steady state diff eq
  i=i+1
  PSedinit = OutT$PSedDens[1]/1000 * OutT$BulkDens[1] * OutT$Zsed[1] # gP/m2 SA
  print(PSedinit)
  OutT$SedDeposition[i] = OutT$PLoad[1]/OutT$SedProp[1] * OutT$PRetention[1] * OutT$Zsed[1] * 1/PSedinit * 1000 # mm/y
  #myEquilstxt$SedDepositionEqtxt = 'SedDepositionEq = f(Psedinit=f(PSedDens,BulkDens,Zsed),PLoad,SedProp,PRetention)'
  OutT$PSedDens[i] = 1
  OutT$BulkDens[i] = 1
  OutT$Zsed[i] = 1
  OutT$PLoad[i] = 1
  OutT$SedProp[i] = 1
  OutT$PRetention[i] = 1
  
  #PSed (g/m2) = (PLoad/SedProp*PRetention) / (SedDeposition/Zsed)
  i=i+1
  
  #SedimentPConcAreal = SedAvailPinit * (1/1000) * SedBulkDensity * ActiveSedDepth # gP/m2
  
  OutT$PSed[i] = (OutT$PLoad[1]/OutT$SedProp[1]*OutT$PRetention[1]) / ((OutT$SedDeposition[1]/1000)/OutT$Zsed[1])
  print(OutT$PSed[i])
  #myEquilstxt$PsedEqtxt = 'PSedEq = f(PLoad,SedProp,PRetention,SedDeposition,Zsed)'
  OutT$PLoad[i] = 1
  OutT$SedProp[i] = 1
  OutT$PRetention[i] = 1
  OutT$SedDeposition[i] = 1
  OutT$Zsed[i] = 1
  
  i=i+1
  OutT$PSedDens[i] = OutT$PSed[i-1] * 1000 / (OutT$BulkDens[1] * OutT$Zsed[1])
  #myEquilstxt$PSedDensEqtxt = 'PSedDensEq = f(PSedEq,BulkDens,Zsed)'
  OutT$PSed[i] = 1
  OutT$BulkDens[i] = 1
  OutT$Zsed[i] = 1

  #reactable(round(OutT,digits=3))
  return(OutT)
  
}

LakeOCEquilibria <- function(myLake,myParameters,WCAllocthProp=NA){
  CTwc = 0.6
  CTsed = 0.3
  # WCAllocthProp is the proportion of observed watercolumn OC that is allocthonous
  # This matters a lot for estimates of sedimentation, sed OC, and burial
  if (is.na(WCAllocthProp)){
    WCAllocthProp = 0.8
  }
  # This function calculates equilibrium values for a number of states/fluxes relevant to lake OC cycling
  # 
  # Water column, loads
  # dOC/dt = OCLoad + NEP – Outflow – R – Sedimentation
  # OCLoad = Outflow + R + Sedimentation – NEP
  # Outflow = 1/RT * OCwc
  # R = RDOC20*CTwc * OCwc
  # Sedimentation = 0.5 * RDOC20*CTwc * OCwc
  # NEP = 0
  # 
  # OCLoad = (1/RT + RDOC20*CTwc + 0.5*RDOC20*CTwc)*Zmean*OCwc
  # where RDOC20 = 0.365 y-1; CTwc = Arrhenius water column temperature adjustment 
  # OCwc = OCLoad/[(1/RT + RDOC20*CTwc  + 0.5*RDOC20*CTwc )*Zmean]
  # Sedimentation = OCLoad + NEP – Outflow – R
  # Sedimentation = OCLoad + NEP – 1/RT*OCwc – RDOC20*CTwc *OCwc
  # 
  # Sediments
  # dOCsed/dt = Sedimentation – R – Burial 
  # Burial = OCsed*Deposition/Zsed
  # R = OCsed * RDOC20sed*CTsed; where RDOC20sed = 0.00365 y-1
  # Sedimentation = OCLoad + NEP – 1/RT*OCwc – RDOC20*CTsed *OCwc
  # 0 = OCLoad + NEP – 1/RT*OCwc – RDOC20*CTsed *OCwc – OCsed*Deposition/Zsed
  # OCsed = [OCLoad + NEP – 1/RT*OCwc – RDOC20*CT *OCwc] / (Deposition/Zsed)
  # 
  nCalcs = 6
  i=1
  
  OutT = data.frame(WRT=rep(NA,nCalcs),Zmean=rep(NA,nCalcs),SedProp=rep(NA,nCalcs),Zsed=rep(NA,nCalcs),SedDeposition=rep(NA,nCalcs),
                    OCLoad=rep(NA,nCalcs),OCwc=rep(NA,nCalcs),OCSedimentationSA=rep(NA,nCalcs),OCsedSA=rep(NA,nCalcs),OCSedDens=rep(NA,nCalcs),BulkDens=rep(NA,nCalcs))
  
  OutT$WRT[1] = myLake$RT/365
  OutT$Zmean[1] = myLake$Zmean
  OutT$SedProp[1] = myLake$SedAreaProp
  OutT$Zsed[1] = myLake$ActiveSedDepth
  OutT$SedDeposition[1] = myLake$SoilDeposition
  OutT$OCLoad[1] = myLake$OCLoad
  OutT$OCwc[1] = myLake$DOCRinit
  OutT$OCSedDens[1] = myLake$SedAvailOC
  OutT$BulkDens[1] = myLake$SedBulkDensity
  RDOC20 = myParameters$RDOCR*365
  RDOC20sed = myParameters$RPOCRSed*365
  DOCRinflow = myLake$DOCRinflow
  OutT$OCsedSA[1] = OutT$OCSedDens[1]/1000 * OutT$BulkDens[1] * OutT$Zsed[1]
  
  ######################################
  # Calculate the organic carbon budget, assuming equilibrium in OCwc and OCsed
  # Assume we know load and outflow
  # dOCwc/dt = Allochthony + netAutochthony - Sedimentation - R - Outflow
  # Allochthony - Outflow = Sedimentation + R - Autochthony; because outflow, R, and Sed are based on
  #   water column OC (which includes allochthonous OC), Autochthony is subsumed
  # R (g/m2LA/y) = OCwc * RDOC20 * CTwc * Zmean
  eqRwc = OutT$OCwc[1]*RDOC20*CTwc*OutT$Zmean[1] * WCAllocthProp
  # Outflow = OCwc * 1/WRT * Zmean
  eqOutflow = OutT$OCwc[1]/OutT$WRT[1]*OutT$Zmean[1] * WCAllocthProp
  # Sedimentation = Load - Outflow - R
  eqSedimentation = OutT$OCLoad[1] - eqOutflow - eqRwc
  dOCwc = OutT$OCLoad[1] - (eqRwc+eqOutflow+eqSedimentation)

  ######################################
  # Calcuate equilibrium values in various ways
  # OCLoad = (1/RT + RDOC20*CTwc + 0.5*RDOC20*CTwc)*Zmean*OCwc
  # Version 2 uses POC input as the sedimentation value
  i=i+1
  # Calculate load based on observed water column OC
  # Note that it's different from assumed load
  #OutT$OCLoad[i] = (1/OutT$WRT[1] + RDOC20*CTwc + 0.5*RDOC20*CTwc)*OutT$Zmean[1]*OutT$OCwc[1]
  eqOCLoad = ((1/OutT$WRT[1] + RDOC20*CTwc)*OutT$Zmean[1]*OutT$OCwc[1])
  OutT$OCLoad[i] = eqOCLoad
  OutT$WRT[i] = 1
  OutT$Zmean[i] = 1
  OutT$OCwc[i] = 1
  
  # Water column OC
  # OCwc = OCLoad/[(1/RT + RDOC20*CTwc  + 0.5*RDOC20*CTwc )*Zmean]
  i=i+1
  eqOCwc = OutT$OCLoad[1] / ((1/OutT$WRT[1] + RDOC20*CTwc + 0.5*RDOC20*CTwc)*OutT$Zmean[1])
  OutT$OCwc[i] = eqOCwc
  OutT$OCLoad[i] = 1
  OutT$WRT[i] = 1
  OutT$Zmean[i] = 1
  
  # Sedimentation of OC
  # In units of g/m2/y in sediment area
  # Sedimentation(1) = OCLoad + NEP – 1/RT*OCwc – RDOC20*CTwc *OCwc
  i=i+1
  eqOCSedimentationSA = (OutT$OCLoad[1] - OutT$OCwc[1]/OutT$WRT[1]*OutT$Zmean[1]* WCAllocthProp - 
                           OutT$OCwc[1]*RDOC20*CTwc*OutT$Zmean[1]* WCAllocthProp)/OutT$SedProp[1]
  OutT$OCSedimentationSA[i] = eqOCSedimentationSA
  OutT$Zmean[i] = 1
  OutT$OCLoad[i] = 1
  OutT$OCwc[i] = 1

  # Sed deposition equilibrium estimated from steady state diff eq
  # For sediments, estimate the proportion that's POCL
  maxPOCLsedProp = 0.8
  kPOCLsedProp = 0.08
  POCLsedProp = (maxPOCLsedProp*myLake$Pinit)/(kPOCLsedProp+myLake$Pinit)
  POCRsedProp = 1

  # OCsed in areal units
  i=i+1
  # dOCsed/dt = Sedimentation - Rsed - Burial
  # Rsed + Burial = Sedimentation; OCsed = Sedimentation / (Rlossrate + Buriallossrate)
  eqOCsedSA = ((eqSedimentation/OutT$SedProp[1]) / (RDOC20sed*CTsed + (OutT$SedDeposition[1]/1000)/OutT$Zsed[1]))
  OutT$OCsedSA[i] = eqOCsedSA
  OutT$OCLoad[i] = 1
  OutT$SedProp[i] = 1
  OutT$OCwc[i] = 1
  OutT$WRT[i] = 1
  OutT$SedDeposition[i] = 1
  OutT$Zsed[i] = 1

  # OC Sed density
  i=i+1
  eqOCSedDensSA = eqOCsedSA * 1000 / (OutT$BulkDens[1] * OutT$Zsed[1])
  OutT$OCSedDens[i] = eqOCSedDensSA
  OutT$OCsedSA[i] = 1
  OutT$BulkDens[i] = 1
  OutT$Zsed[i] = 1
  
  # Burial, based on previous OCsed and deposition
  eqBurialSA = eqOCsedSA * (OutT$SedDeposition[1]/1000)/OutT$Zsed[1]  
  
  # Rsed, based on previous OCsed
  eqRsedSA = eqOCsedSA * RDOC20sed*CTsed 
  
  dOCsed = eqSedimentation - eqRsedSA*OutT$SedProp[1] - eqBurialSA*OutT$SedProp[1]
  
  #################################
  # Print the budget
  print('===================================================')
  print('Estimates of lake OC from equilibrium modeling.')
  print('These underestimate autochthonous contributions.')
  print('===================================================')
  print('Water column OC budget (gC/m2LA/y)')
  print('Load = R + Sedimentation* + Outflow')
  print(paste(round(OutT$OCLoad[1],digits=3),' = ',round(eqRwc,digits=3),' + ',round(eqSedimentation,digits=3),' + ',round(eqOutflow,digits=3),sep=""))
  print(paste('Water column budget balance: ',round(dOCwc,digits=3),sep=""))
  print('*Sedimentation calculated by difference; should add net autochthony')
  print('------------------------')
  print(paste('Sediment values assume POCL proportion: ',round(POCLsedProp,digits=3),sep=""))
  print('Sediment OC budget (gC/m2LA/y)')
  print('dOCsed/dt = Sedimentation - R - Burial')
  print('Sedimentation = Rsed + Burial')
  print(paste(round(eqSedimentation,digits=3),' = ',round(eqRsedSA*OutT$SedProp[1],digits=3)," + ",round(eqBurialSA*OutT$SedProp[1],digits=3),sep=""))
  print(paste('Sediment budget balance: ',round(dOCsed,digits=3),sep=""))
  print('In sediment area units (gC/m2SA/y)')
  print(paste(round(eqSedimentation/OutT$SedProp[1],digits=3),' = ',round(eqRsedSA,digits=3)," + ",round(eqBurialSA,digits=3),sep=""))
  print('---------------------------------------------------')
  
  return(OutT)
  
}

MeanEpiLight <- function(z1,z2,a1,a2,ir,mu){
  deps = seq(z1, z2 , 0.1)
  irr =  ir * exp(-(mu)*deps)
  area = approx(c(z1,z2), c(a1,a2), deps)$y
  intarea = pracma::trapz((deps),(area))
  
  return((area %*% irr)/ (sum(area)))
}

DensityByYear <- function(SimulationDates,Variable,Label){
  PlotIt = TRUE
  uYears = unique(format(SimulationDates,'%Y'))
  myDensities = list()
  # par(mfrow=c(length(uYears),1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  par(mfrow = c(1,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  for (i in 1:length(uYears)){
    thisYearDensities = list()
    thisYear = uYears[i]
    ithisYear = which(format(SimulationDates,'%Y')==thisYear)
    thisDensity = density(na.omit(Variable[ithisYear]))
    if (PlotIt){
      plot(density(na.omit(Variable[ithisYear])),xlim=c(-4,4),main="",xlab="")
    }
    thisYearDensities[[1]] = thisYear
    thisYearDensities[[2]] = thisDensity
    myDensities[[i]] = thisYearDensities
  }
  
  for (i in 1:length(uYears)){
    thisYearDensities = list()
    thisYear = uYears[i]
    ithisYear = which(format(SimulationDates,'%Y')==thisYear)
    thisDensity = density(diff(na.omit(Variable[ithisYear])))
    if (PlotIt){
      plot(density(diff(na.omit(Variable[ithisYear]))),xlim=c(-1,1),main="",xlab="")
    }
    thisYearDensities[[1]] = thisYear
    thisYearDensities[[2]] = thisDensity
    myDensities[[i]] = thisYearDensities
  }
  return(myDensities)
}



