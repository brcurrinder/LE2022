Metabolism = function(nDays,nYears2Save,myLake,myParameters,Irradiance,
                      WindSpeed,Discharge,OCLoad,TPLoad,LakeHypsometry,LoadFromState,StateFile,Verbose){
  # Metabolism again
  # Assume all boundary data are of the same length and cover the same date range
  # This allows the code to run through data as indices rather than checking all the dates
  
  # source('MetabolismHelpers.R')
  # Load lake metabolizer for dissolved gas functions
  # library(LakeMetabolizer)
  O2toCMR = 32/12
  AnoxiaThreshold = 1.0 # mgO2/L, used for sed P release, rebinding
  MaxEpiDepthProp = 0.7 # Proportion of total depth that can be epilimnion during stratification
  #  controls spikes in hypo concentrations during late in stratification
  KO2IceConditions = 0.6 # KO2 when the lake is ice covered
  AlbedoIce = 0.9 # albedo under ice on conditions
  ModelClearwaterPhase = TRUE # Whether or not to shutdown GPP for part of spring
  
  StartIndex = 2 # Starts the model on day 2 so we have previous day's states;
  # however, this gets switched to 1 automatically if this is not the first simulation
  
  ####################
  # General
  LakeAreaInit = max(LakeHypsometry$area)
  ActiveSedDepth = myLake$ActiveSedDepth
  Albedo = myLake$Albedo
  
  LakeVolumes = data.frame(EpiVol=rep(NA,nDays), HypoVol=rep(NA,nDays), TotalVol=rep(sum(LakeHypsometry$volume),nDays),ThermoclineDepth=rep(NA,nDays))
  LakeVolumes$EpiVol[1] = LakeVolumes$TotalVol[1] # Assumes we are starting in mixed conditions
  LakeVolumes$HypoVol[1] = LakeVolumes$TotalVol[1] - LakeVolumes$EpiVol[1]
  
  LakeAreas = data.frame(EpiArea=rep(NA,nDays), HypoArea=rep(NA,nDays), TotalArea=rep(max(LakeHypsometry$area),nDays))
  LakeAreas$EpiArea[1] = LakeAreas$TotalArea[1] # Assumes we are starting in mixed conditions
  LakeAreas$HypoArea[1] = LakeAreas$TotalArea[1] - LakeAreas$EpiArea[1]
  
  LakeTemps = data.frame(EpiT=rep(NA,nDays), HypoT=rep(NA,nDays), IceCovered=rep(NA,nDays), 
                         IceThickness=rep(FALSE,nDays),MeanEpiIrradiance=rep(FALSE,nDays))
  LakeTemps$EpiT[1] = 4
  LakeTemps$HypoT[1] = 4
  LakeTemps$IceCovered[1] = TRUE
  LakeTemps$IceThickness[1] = 0.5
  LakeTemps$MeanEpiIrradiance[1] = 0
  ####################
  # Internal rates
  # Primary production
  CP = myParameters$CP
  GPP2POC = myParameters$GPP2POC
  GPPTheta = myParameters$GPPTheta
  # Michaelis Menten scaler for DO
  # plot(DO,DOmaxProp*DO/(kDO+DO),type='l')
  DOmaxProp = myParameters$DOmaxProp
  kDO = myParameters$kDO
  mmDOmin = myParameters$mmDOmin
  # #kDO = 0.1 # g/m3 half saturation coefficient for Michaelis Menten equation
  # # Respiration
  RDOCR = myParameters$RDOCR
  RDOCL = myParameters$RDOCL
  RPOCR = myParameters$RPOCR
  RPOCL = myParameters$RPOCL
  RPOCRSed = myParameters$RPOCRSed
  RPOCLSed = myParameters$RPOCLSed
  Tbase = myParameters$Tbase
  RTheta = myParameters$RTheta
  SedAvailOC = myLake$SedAvailOC
  SedAvailOCBurial = myLake$SedAvailOCBurial

  # P cycling
  SedAreaProp = myLake$SedAreaProp
  SedAvailPBurial = myLake$SedAvailPBurial
  # PRetention = myParameters$PRetention # proportion of P load retained /y
  SedBulkDensity = myLake$SedBulkDensity # g/m3 standard buld for clay/silt
  # SedWaterProp = myLake$SedWaterProp # Proportion of sediment volume that's water
  SedAvailPinit = myLake$SedAvailPinit # mgP/g sediment dry weight
  SedBoundP = myLake$SedBoundP # converted to 
  SoilDeposition = myLake$SoilDeposition # mm/y
  SedPOCRinitProp = myLake$SedPOCRinitProp
  SedPOCLinitProp = myLake$SedPOCLinitProp
  #Sedimentation Rate ((12.5)*(1/1000)*(1/.3) # pats dissertation)
  TPcSed = myParameters$TPcSed # 1/d, 0.0137 for Mendota, Hanson et al. 2020
  TPb = myParameters$TPb # g dry weight m-2 d-1; second constant in Nurnberg Recycling Eq
  TPa = myParameters$TPa #pars [3] #original is -4.3, first constant in Nurnberg Recycling Eq
  TPRecyclingTheta = myParameters$TPRecyclingTheta # Theta for sedimentation
  TPSedimentationTheta = myParameters$TPSedimentationTheta # Theta for sedimentation
  TPSedimentationBaseT = myParameters$TPSedimentationBaseT # C
  TPcRecycle = myParameters$TPcRecycle # Recycling coefficient, model 1 1/d
  TPkRecycle = myParameters$TPkRecycle # Recycling coefficient, model 2 mm/d
  TPRecyclingBaseT = myParameters$TPRecyclingBaseT # C
  TPcRelease = myParameters$TPcRelease # P release under anoxia
  TPcRebind = myParameters$TPcRebind # P rebinding under anoxia
  
  # Settling
  POCLSettling = myParameters$POCLSettling
  POCRSettling = myParameters$POCRSettling

  # Clarity
  LECwater = myParameters$LECwater
  LECPOCL = myParameters$LECPOCL
  LECPOCR = myParameters$LECPOCR
  LECDOCR = myParameters$LECDOCR
  LECDOCL = myParameters$LECDOCL
  LEC2Secchi = myParameters$LEC2Secchi
  # LEC = LECwater + LECPOCL*0.2 + LECPOCR*0.1 + LECDOCR*3 + LECDOCL*0.3
  # Secchi = LEC2Secchi/LEC
  
  #############################
  # Initialize states and rates data frames
  
  # States
  StatesEpiC  = data.frame(SimDay=rep(NA,nDays),DOCR = rep(NA,nDays),DOCL=rep(NA,nDays),
                           POCL=rep(NA,nDays),POCR=rep(NA,nDays),Secchi=rep(NA,nDays))
  StatesHypoC = data.frame(SimDay=rep(NA,nDays),DOCR = rep(NA,nDays),DOCL=rep(NA,nDays),
                           POCL=rep(NA,nDays),POCR=rep(NA,nDays))
  StatesSedC  = data.frame(SimDay=rep(NA,nDays),DOCR = rep(NA,nDays),DOCL=rep(NA,nDays),
                           POCL=rep(NA,nDays),POCR=rep(NA,nDays))
  StatesEpiGas = data.frame(SimDay=rep(NA,nDays),DO=rep(NA,nDays),DOsat=rep(NA,nDays),CO2=rep(NA,nDays),CH4=rep(NA,nDays))
  StatesHypoGas = data.frame(SimDay=rep(NA,nDays),DO=rep(NA,nDays),DOsat=rep(NA,nDays),CO2=rep(NA,nDays),CH4=rep(NA,nDays))
  
  StatesEpiP  = data.frame(SimDay=rep(NA,nDays),TP = rep(NA,nDays))
  StatesHypoP  = data.frame(SimDay=rep(NA,nDays),TP = rep(NA,nDays))
  StatesSedP  = data.frame(SimDay=rep(NA,nDays),TP = rep(NA,nDays),TPB = rep(NA,nDays))
  
  # Set initial values
  StatesEpiC$DOCR[1] = myLake$DOCRinit * LakeAreaInit * myLake$Zmean
  StatesEpiC$DOCL[1] = myLake$DOCLinit * LakeAreaInit * myLake$Zmean
  StatesEpiC$POCL[1] = myLake$POCLinit * LakeAreaInit * myLake$Zmean
  StatesEpiC$POCR[1] = myLake$POCRinit * LakeAreaInit * myLake$Zmean
  StatesEpiP$TP[1]   = myLake$Pinit * LakeAreaInit * myLake$Zmean
  
  ##################
  # General sediment
  SedimentArea = SedAreaProp * LakeAreaInit
  SedimentVolume = SedimentArea * ActiveSedDepth #m3
  
  ##################
  # Sediment OC
  # Assume 90% of initial sediment is POCR, the rest is POCL
  # Example from phosphorus: SedAvailPinit * (1/1000) * BulkDensity  * SedDepth
  CsedEmpiricalAreal = SedAvailOC/1000 * SedBulkDensity  * ActiveSedDepth
  # Adjust sediment OC according to initial P concentration
  StatesSedC$POCR[1] = CsedEmpiricalAreal * SedimentArea * SedPOCRinitProp
  StatesSedC$POCL[1] = CsedEmpiricalAreal * SedimentArea * SedPOCLinitProp

  ##################
  # Phosphorus
  SedimentPConcAreal = SedAvailPinit * (1/1000) * SedBulkDensity * ActiveSedDepth # gP/m2
  # First fill up bound P if appropriate
  StatesSedP$TPB[1] =  min(SedimentPConcAreal*SedimentArea, SedBoundP*(1/1000)*SedBulkDensity * ActiveSedDepth * SedimentArea) # gP boundP
  StatesSedP$TP[1] = SedimentPConcAreal*SedimentArea - StatesSedP$TPB[1] # g
  
  # Initialize a variable for maximum bound P for this simulation
  SedTPBoundMaxgrams = SedBoundP * (1/1000) * SedBulkDensity * ActiveSedDepth * SedimentArea

  ##################
  # LEC & Secchi
  LEC = LECwater + LECPOCL*StatesEpiC$POCL[1]/LakeVolumes$EpiVol[1] + LECPOCR*StatesEpiC$POCR[1]/LakeVolumes$EpiVol[1] + 
    LECDOCR*StatesEpiC$DOCR[1]/LakeVolumes$EpiVol[1] + LECDOCL*StatesEpiC$DOCL[1]/LakeVolumes$EpiVol[1]
  StatesEpiC$Secchi[1] = LEC2Secchi/LEC
  
  StatesHypoC$DOCR[1] = 0 # myLake$DOCRinit * LakeAreaInit * myLake$Zmean
  StatesHypoC$DOCL[1] = 0 # myLake$DOCLinit * LakeAreaInit * myLake$Zmean
  StatesHypoC$POCL[1] = 0 # myLake$POCinit * LakeAreaInit * myLake$Zmean
  StatesHypoC$POCR[1] = 0 # myLake$POCinit * LakeAreaInit * myLake$Zmean

  # Fill in saturation values for gases and set gas initial conditions
  TforEpiDO = data.frame(datetime=1,wtr=LakeTemps$EpiT[1])
  StatesEpiGas$DOsat[1] = o2.at.sat(TforEpiDO,altitude = myLake$Elevation)$do.sat # Concentration units of g/m3
  StatesEpiGas$DO[1] = StatesEpiGas$DOsat[1] * LakeAreaInit * myLake$Zmean
  StatesHypoGas$DO[1] = StatesHypoGas$DOsat[1] * LakeAreaInit * myLake$Zmean
  
  # Overwrite initial states with end of previous state if LoadFromState
  if (LoadFromState){ # Only if TRUE
    StartIndex = 2 # switched from day 2, because we have states and rates from end of previous sim
    load(StateFile)
    nRecs = dim(ModelRun[[7]][[1]])[1]
    StatesEpiC$DOCR[1] = ModelRun[[7]][[1]]$DOCR[nRecs]
    StatesEpiC$DOCL[1] = ModelRun[[7]][[1]]$DOCL[nRecs]
    StatesEpiC$POCL[1] = ModelRun[[7]][[1]]$POCL[nRecs]
    StatesEpiC$POCR[1] = ModelRun[[7]][[1]]$POCR[nRecs]
    StatesEpiP$TP[1]   = ModelRun[[7]][[9]]$TP[nRecs]
    
    StatesSedC$POCR[1] = ModelRun[[7]][[15]]$POCR[nRecs]
    StatesSedC$POCL[1] = ModelRun[[7]][[15]]$POCL[nRecs]
    
    StatesSedP$TP[1] = ModelRun[[7]][[11]]$TP[nRecs]
    StatesSedP$TPB[1] = ModelRun[[7]][[11]]$TPB[nRecs]
    
    StatesEpiC$Secchi[1] = ModelRun[[7]][[1]]$Secchi[nRecs]
    
    StatesHypoC$DOCR[1] = ModelRun[[7]][[2]]$DOCR[nRecs]
    StatesHypoC$DOCL[1] = ModelRun[[7]][[2]]$DOCL[nRecs]
    StatesHypoC$POCL[1] = ModelRun[[7]][[2]]$POCL[nRecs]
    StatesHypoC$POCR[1] = ModelRun[[7]][[2]]$POCR[nRecs]
    
    StatesEpiGas$DO[1] = ModelRun[[7]][[3]]$DO[nRecs]
    StatesHypoGas$DO[1] = ModelRun[[7]][[4]]$DO[nRecs]
    rm(ModelRun)
  }
  
  ##################
  # Setup rates
  RatesEpiC = data.frame(SimDay=rep(NA,nDays),
                         DOCLMix=rep(NA,nDays),GPP=rep(NA,nDays),DOCLR=rep(NA,nDays),DOCLInflow=rep(NA,nDays),DOCLOutflow=rep(NA,nDays),
                         DOCRMix=rep(NA,nDays),DOCRR=rep(NA,nDays),DOCRInflow=rep(NA,nDays),DOCROutflow=rep(NA,nDays),
                         POCLMix=rep(NA,nDays),POCLR=rep(NA,nDays),POCLInflow=rep(NA,nDays),POCLOutflow=rep(NA,nDays),POCLSettling=rep(NA,nDays),
                         POCRMix=rep(NA,nDays),POCRR=rep(NA,nDays),POCRInflow=rep(NA,nDays),POCROutflow=rep(NA,nDays),POCRSettling=rep(NA,nDays))
  RatesHypoC = data.frame(SimDay=rep(NA,nDays),
                          DOCLMix=rep(NA,nDays),GPP=rep(NA,nDays),DOCLR=rep(NA,nDays),DOCLInflow=rep(NA,nDays),DOCLOutflow=rep(NA,nDays),
                          DOCRMix=rep(NA,nDays),DOCRR=rep(NA,nDays),DOCRInflow=rep(NA,nDays),DOCROutflow=rep(NA,nDays),
                          POCLMix=rep(NA,nDays),POCLR=rep(NA,nDays),POCLInflow=rep(NA,nDays),POCLOutflow=rep(NA,nDays),POCLSettling=rep(NA,nDays),
                          POCRMix=rep(NA,nDays),POCRR=rep(NA,nDays),POCRInflow=rep(NA,nDays),POCROutflow=rep(NA,nDays),POCRSettling=rep(NA,nDays))
  RatesSedC = data.frame(SimDay=rep(NA,nDays),
                          DOCLMix=rep(NA,nDays),GPP=rep(NA,nDays),DOCLR=rep(NA,nDays),DOCLInflow=rep(NA,nDays),DOCLOutflow=rep(NA,nDays),
                          DOCRMix=rep(NA,nDays),DOCRR=rep(NA,nDays),DOCRInflow=rep(NA,nDays),DOCROutflow=rep(NA,nDays),
                          POCLMix=rep(NA,nDays),POCLR=rep(NA,nDays),POCLInflow=rep(NA,nDays),POCLOutflow=rep(NA,nDays),POCLSettling=rep(NA,nDays),
                          POCRMix=rep(NA,nDays),POCRR=rep(NA,nDays),POCRInflow=rep(NA,nDays),POCROutflow=rep(NA,nDays),POCRSettling=rep(NA,nDays),
                         POCLBurial=rep(NA,nDays),POCRBurial=rep(NA,nDays))
  RatesEpiP = data.frame(SimDay=rep(NA,nDays),
                         TPInflow=rep(NA,nDays),TPOutflow=rep(NA,nDays),
                         TPMix=rep(NA,nDays),TPSettling=rep(NA,nDays),TPRecycling=rep(NA,nDays),
                         TPRebind=rep(NA,nDays),TPRelease=rep(NA,nDays))
  RatesHypoP = data.frame(SimDay=rep(NA,nDays),
                         TPInflow=rep(NA,nDays),TPOutflow=rep(NA,nDays),
                         TPMix=rep(NA,nDays),TPSettling=rep(NA,nDays),TPRecycling=rep(NA,nDays),
                         TPRebind=rep(NA,nDays),TPRelease=rep(NA,nDays))
  RatesSedP = data.frame(SimDay=rep(NA,nDays),
                         TPInflow=rep(NA,nDays),TPOutflow=rep(NA,nDays),
                         TPMix=rep(NA,nDays),TPSettling=rep(NA,nDays),
                         TPcFsedP=rep(NA,nDays),TPRecycling=rep(NA,nDays),TPBurial=rep(NA,nDays),
                         TPRebind=rep(NA,nDays),TPRelease=rep(NA,nDays),TPBBurial=rep(NA,nDays))

  RatesEpiGas  = data.frame(SimDay=rep(NA,nDays),DOK=rep(NA,nDays),DOFatm=rep(NA,nDays),DOMix=rep(NA,nDays),DOSed=rep(NA,nDays))                      
  RatesHypoGas = data.frame(SimDay=rep(NA,nDays),DOK=rep(NA,nDays),DOFatm=rep(NA,nDays),DOMix=rep(NA,nDays),DOSed=rep(NA,nDays))
  
  ##################
  # Get gas exchange values (KO2) and save to rates vector
  # WindForGasEx = data.frame(datetime=LakeTemps$DOY,wnd=rep(myLake$WindSpeed,nDays))
  WindForGasEx = data.frame(datetime=WindSpeed$myDays,wnd=rep(WindSpeed$u10))
  K600 = k.cole(WindForGasEx)
  # KO2 = k600.2.kGAS(data.frame(datetime=K600$datetime,k600=K600$k600,wtr=LakeTemps$EpiT),gas="O2")
  # RatesEpiGas$DOK = KO2$k.gas
  
  ###################
  # Initialize lake physics
  u = u_ini
  # Start time is an index of the interpolated met data, with index of 1 being the first day
  # increments are by seconds, 
  # startTime = as.numeric((Irradiance$myDays[2] - Irradiance$myDays[1]) * 86400) # current - index 1 in the driver data in seconds  )
  startTime = 1 # current - index 1 in the driver data in seconds  )
  # End time in this implementation is 1 day after start time
  endTime = startTime + hydrodynamic_timestep - 1 # 
  ice = TRUE
  Hi = LakeTemps$IceThickness[1]
  iceT = 0
  supercooled = 0
  kd_light = 0.15
  # matrix_range = max(1, (startTime/dt)):((endTime/dt)+1)
  matrix_range_start = max(1, round(startTime/dt) + 1)
  matrix_range_end = round(endTime/dt)
  # Run the following only once to get the met data
  meteo = get_interp_drivers(meteo_all=meteo_all, total_runtime=nDays, 
                             hydrodynamic_timestep=hydrodynamic_timestep, dt=dt)
  # daily_meteo$Shoftwave_Radiation_Downwelling_wattPerMeterSquared = daily_meteo$Shoftwave_Radiation_Downwelling_wattPerMeterSquared * 0.9
  
  # Run physics model for a day
  res <-  run_thermalmodel(u = u, 
                           startTime = startTime, 
                           endTime =  endTime, 
                           ice = ice, 
                           Hi = Hi, 
                           iceT = iceT,
                           supercooled = supercooled,
                           kd_light = kd_light,
                           zmax = zmax,
                           nx = nx,
                           dt = dt,
                           dx = dx,
                           area = hyps_all[[1]], # area
                           depth = hyps_all[[2]], # depth
                           volume = hyps_all[[3]], # volume
                           # daily_meteo = meteo_all[[1]],
                           daily_meteo = meteo[,matrix_range_start:matrix_range_end],
                           # secview = meteo_all[[2]],
                           Cd = 0.0008,
                           densThresh = 1e-3)
  
  ###############
  # Setup a progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = nDays, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  ########################################################
  ########################################################
  # Run the model
  # Loop through time, assuming 1 is first day of the year
  SimRunTime = rep(NA,nDays)
  DayOfStratification = 0
  for (i in StartIndex:nDays){
    SimDay = i
    # if (SimDay == 5475){
    #   TPcRecycle = 0.00075
    # }
    if (Verbose){
      setTxtProgressBar(pb, i)
    }

    StatesEpiC$SimDay[i] = i
    StatesHypoC$SimDay[i] = i
    StatesSedC$SimDay[i] = i
    StatesEpiGas$SimDay[i] = i
    StatesHypoGas$SimDay[i] = i
    StatesEpiP$SimDay[i] = i
    StatesHypoP$SimDay[i] = i
    StatesSedP$SimDay[i] = i
    
    ###########################################
    # Robert's lake physics
    total_runtime <- 1 # number of days to run the model
    # Set the inputs for physics to outputs from last time step
    u = res$temp[, ncol(res$temp)]
    startTime = as.numeric((Irradiance$myDays[i] - Irradiance$myDays[1]) * 86400) # current - index 1 in the driver data in seconds  )
    # End time in this implementation is 1 day after start time
    endTime = startTime + hydrodynamic_timestep - 1 #
    ice = res$iceflag
    Hi = res$icethickness 
    iceT = res$icemovAvg
    supercooled = res$supercooled
    kd_light = 1.7/StatesEpiC$Secchi[i-1]
    #matrix_range = max(1, (startTime/dt)):((endTime/dt)+1)
    matrix_range_start = max(1, round(startTime/dt) + 1)
    matrix_range_end = round(endTime/dt)
    res <-  run_thermalmodel(u = u, 
                             startTime = startTime, # 86400  25833600
                             endTime =  endTime, # 172799  25919999
                             ice = ice, 
                             Hi = Hi, 
                             iceT = iceT,
                             supercooled = supercooled,
                             kd_light = kd_light,
                             zmax = zmax,
                             nx = nx,
                             dt = dt,
                             dx = dx,
                             area = hyps_all[[1]], # area
                             depth = hyps_all[[2]], # depth
                             volume = hyps_all[[3]], # volume
                             # daily_meteo = meteo_all[[1]],
                             daily_meteo = meteo[,matrix_range_start:matrix_range_end],
                             # secview = meteo_all[[2]],
                             Cd = 0.0008, #0.0008,
                             Hgeo = 0.1,
                             densThresh = 1e-3)
    
    # Store values used for the P-model, these are daily averages by stratum and for whole lake
    LakePhysicsAverages <- res$average %>% mutate(datetime = as.POSIXct(startingDate) + time,
                                                  Date = as.Date(datetime, format = "%m/%d/%Y")) %>% summarise_all(mean)
    
    # meanT = mean(c(LakePhysicsAverages$epi,LakePhysicsAverages$hypo),na.rm=TRUE)
    # Save ice values
    LakeTemps$IceCovered[i] = res$iceflag
    LakeTemps$IceThickness[i] = res$icethickness
    if (LakeTemps$IceCovered[i]){
      Albedo = AlbedoIce
      # print(Hi)
      # TPcSed = myParameters$TPcSed * 0.5
      # TPcRebind = myParameters$TPcRebind * 0.5
      mmDOepi = max((DOmaxProp*StatesEpiGas$DO[i-1]/LakeVolumes$EpiVol[i-1])/(kDO+StatesEpiGas$DO[i-1]/LakeVolumes$EpiVol[i-1]),
                   mmDOmin) * 1
    }else{
      Albedo = myLake$Albedo
      TPcSed = myParameters$TPcSed
      TPcRebind = myParameters$TPcRebind
      mmDOepi = max((DOmaxProp*StatesEpiGas$DO[i-1]/LakeVolumes$EpiVol[i-1])/(kDO+StatesEpiGas$DO[i-1]/LakeVolumes$EpiVol[i-1]),mmDOmin)
    }

    # Set flag for stratification
    # Set lake strata volumes and areas
    # if(LakePhysicsAverages$stratFlag==1 | LakeTemps$IceCovered[i]){
    if(LakePhysicsAverages$stratFlag==1){
      isStratified = TRUE
      
      DayOfStratification = DayOfStratification + 1
      # ThermoclineDepth = round(LakePhysicsAverages$thermoclineDep)
      ThermoclineDepth = LakePhysicsAverages$thermoclineDep
      # Do not allow thermocline depth to exceed maximum proportion to control for hypo spikes
      ThermoclineDepth = min(ThermoclineDepth, MaxEpiDepthProp*(max(LakeHypsometry$depth)))
      # MeanEpiThickness = ThermoclineDepth # meters, this is used for particle settling
      MeanEpiThickness = sum(LakeHypsometry$volume[1:round(ThermoclineDepth)]) / max(LakeHypsometry$area) # meters, this is used for particle settling
      # MeanHypoThickness = max(LakeHypsometry$depth) - MeanEpiThickness
      MeanHypoThickness = sum(LakeHypsometry$volume[round(ThermoclineDepth):length(LakeHypsometry$volume)]) / LakeHypsometry$area[round(ThermoclineDepth)] # meters, this is used for particle settling
      #print(ThermoclineDepth)
      LakeVolumes$ThermoclineDepth[i] = ThermoclineDepth
      # print(paste('Thermocline depth: ', LakePhysicsAverages$thermoclineDep))
      # Lake volumes
      # LakeVolumes$EpiVol[i] = sum(LakeHypsometry$volume[1:ThermoclineDepth])
      LakeVolumes$EpiVol[i] = sum(LakeHypsometry$volume[1:floor(ThermoclineDepth)]) + 
        (ThermoclineDepth - floor(ThermoclineDepth)) * LakeHypsometry$volume[floor(ThermoclineDepth)+1]
      LakeVolumes$HypoVol[i] = LakeVolumes$TotalVol[i] - LakeVolumes$EpiVol[i]
      iThermoclineDepth = which(LakeHypsometry$depth >= ThermoclineDepth)[1]
      # lake areas
      LakeAreas$HypoArea[i] = LakeHypsometry$area[iThermoclineDepth]
      LakeAreas$EpiArea[i] = LakeAreas$TotalArea[i] - LakeAreas$HypoArea[i]
      LakeTemps$EpiT[i] = LakePhysicsAverages$epi
      LakeTemps$HypoT[i] = LakePhysicsAverages$hyp
      # Light at depth
      LakeTemps$MeanEpiIrradiance[i] = (1-Albedo) * MeanEpiLight(z1=0,z2=ThermoclineDepth,a1=max(LakeHypsometry$area),
             a2=LakeAreas$HypoArea[i],ir=Irradiance$I[i],mu=kd_light)
    }else{
      isStratified = FALSE
      ThermoclineDepth = 0
      DayOfStratification = 0
      LakeVolumes$EpiVol[i] = LakeVolumes$TotalVol[i]
      LakeVolumes$HypoVol[i] = 0
      LakeVolumes$ThermoclineDepth[i] = 0
      LakeAreas$EpiArea[i] = LakeAreas$TotalArea[i]
      LakeAreas$HypoArea[i] = 0
      #############################
      # Because of the way Robert's lake physics work, need to calculate my own average water temp
      LakeTemps$EpiT[i] = mean(res$temp,na.rm = TRUE) # Would be better to volume weight this
      #############################
      #LakeTemps$EpiT[i] = LakePhysicsAverages$epi
      LakeTemps$HypoT[i] = NA
      LakeTemps$MeanEpiIrradiance[i] = (1-Albedo) * MeanEpiLight(z1=0,z2=max(LakeHypsometry$depth),a1=max(LakeHypsometry$area),
             a2=max(LakeHypsometry$area),ir=Irradiance$I[i],mu=kd_light)
      # max(LakeHypsometry$depth) # meters, this is used for particle settling
      MeanEpiThickness = sum(LakeHypsometry$volume) / max(LakeHypsometry$area)
      MeanHypoThickness = 0
    }
    
    # print(paste(MeanEpiThickness,MeanHypoThickness,MeanEpiThickness+MeanHypoThickness))
    # Update gas saturation values
    # Fill in saturation values for gases and set gas initial conditions
    TforEpiDO = data.frame(datetime=1,wtr=LakeTemps$EpiT[i])
    StatesEpiGas$DOsat[i] = o2.at.sat(TforEpiDO,altitude = myLake$Elevation)$do.sat # Concentration units of g/m3
    
    KO2 = k600.2.kGAS(data.frame(datetime=K600$datetime[i],k600=K600$k600[i],wtr=max(LakeTemps$EpiT[i],4)),gas="O2")
    if (LakeTemps$EpiT[i] < 4 & LakeTemps$IceCovered[i]){
      # We have ice
      RatesEpiGas$DOK[i] = KO2IceConditions
    }else{
      RatesEpiGas$DOK[i] = KO2$k.gas
    }
      
    ###########################################
    # This first section it doesn't matter whether the lake is stratified
    
    ###############
    # hydrologic loads and outflow
    OutflowVol = Discharge$Flow[i]
    # Partition loads by predefined speciation
    RatesEpiC$DOCLInflow[i] = OCLoad[i] * myLake$DOCLinflow
    RatesEpiC$DOCLOutflow[i] = OutflowVol * StatesEpiC$DOCL[i-1]/LakeVolumes$EpiVol[i-1] # m3/d * g/m3
    RatesEpiC$DOCRInflow[i] = OCLoad[i] * myLake$DOCRinflow
    RatesEpiC$DOCROutflow[i] = OutflowVol * StatesEpiC$DOCR[i-1]/LakeVolumes$EpiVol[i-1]
    RatesEpiC$POCLInflow[i] = OCLoad[i] * myLake$POCLinflow
    RatesEpiC$POCLOutflow[i] = OutflowVol * StatesEpiC$POCL[i-1]/LakeVolumes$EpiVol[i-1]
    RatesEpiC$POCRInflow[i] = OCLoad[i] * myLake$POCRinflow
    RatesEpiC$POCROutflow[i] = OutflowVol * StatesEpiC$POCR[i-1]/LakeVolumes$EpiVol[i-1]
    RatesEpiP$TPInflow[i] = TPLoad[i]
    # if (SimDay >= 5475){
    #   RatesEpiP$TPInflow[i] = TPLoad[i]*0.20
    # }
    RatesEpiP$TPOutflow[i] = OutflowVol * StatesEpiP$TP[i-1]/LakeVolumes$EpiVol[i-1]
    
    ###############
    # Epi rates
    # Assume GPP includes both water column and benthic
    # RatesEpiC$GPP[i] = StatesEpiP$TP[i-1]/LakeVolumes$EpiVol[i-1] * CP * GPPTheta^(LakeTemps$EpiT[i-1]-Tbase) * 
      # Irradiance$I[i] * StatesEpiC$Secchi[i-1]*2/myLake$Zmean * LakeVolumes$EpiVol[i-1] # g C
    RatesEpiC$GPP[i] = (StatesEpiP$TP[i-1]/LakeVolumes$EpiVol[i-1]) * CP * LakeTemps$MeanEpiIrradiance[i] *
       GPPTheta^(LakeTemps$EpiT[i-1]-Tbase) * LakeVolumes$EpiVol[i-1] # g C
    
    # Kill NPP during clearwater phase, if that's desired
    if (ModelClearwaterPhase){
      DayOfYear = SimDay %% 365
      if (DayOfYear >120 & DayOfYear < 165){
        RatesEpiC$GPP[i] = RatesEpiC$GPP[i] * 0.001
      }
    }
    
    # TPConcentration = (StatesEpiP$TP[i-1]/LakeVolumes$EpiVol[i-1])
    # RatesEpiC$GPP[i] = TPConcentration * (TPConcentration*2-0.01) * LakeTemps$MeanEpiIrradiance[i] *
    #   GPPTheta^(LakeTemps$EpiT[i-1]-Tbase) * LakeVolumes$EpiVol[i-1] # g C
    
    # Respiration
    RatesEpiC$DOCLR[i] = StatesEpiC$DOCL[i-1] * RDOCL * RTheta^(LakeTemps$EpiT[i]-Tbase) # g C
    RatesEpiC$DOCRR[i] = StatesEpiC$DOCR[i-1] * RDOCR * RTheta^(LakeTemps$EpiT[i]-Tbase) # g C
    RatesEpiC$POCLR[i] = StatesEpiC$POCL[i-1] * RPOCL * RTheta^(LakeTemps$EpiT[i]-Tbase) # g C
    RatesEpiC$POCRR[i] = StatesEpiC$POCR[i-1] * RPOCR * RTheta^(LakeTemps$EpiT[i]-Tbase) # g C
    # Settling, based on velocity and depth
    RatesEpiC$POCLSettling[i] = min(StatesEpiC$POCL[i-1],StatesEpiC$POCL[i-1] * POCLSettling / MeanEpiThickness) #* (1-RTheta^(LakeTemps$EpiT[i]-Tbase))
    #print(POCLSettling * (1-RTheta^(LakeTemps$EpiT[i]-Tbase)))
    RatesEpiC$POCRSettling[i] = min(StatesEpiC$POCR[i-1],StatesEpiC$POCR[i-1] * POCRSettling / MeanEpiThickness)
    RatesEpiP$TPSettling[i] = min(StatesEpiP$TP[i-1],StatesEpiP$TP[i-1] * TPcSed / MeanEpiThickness * 
      TPSedimentationTheta^(LakeTemps$EpiT[i]-TPSedimentationBaseT))

    # Epi gas exchange
    RatesEpiGas$DOFatm[i] = RatesEpiGas$DOK[i] * (StatesEpiGas$DOsat[i-1] - StatesEpiGas$DO[i-1]/LakeVolumes$EpiVol[i-1]) *
      LakeAreaInit # g O2
    
    ###############
    # Sediments
    # Deal with sediment area and proportion of that in the epi versus hypo
    PropSedInHypo = min(1,LakeAreas$HypoArea[i-1]/SedimentArea) # in case hypo area > sed area (extremely unusual)
    PropSedInEpi = 1-PropSedInHypo
    
    # Soil burial
    # SedBuried is used for the parameter version of sediment burial but not for the empirical version
    SedBuried =  myLake$SedBulkDensity * LakeAreaInit * SedAreaProp * SoilDeposition * (1/1000) *  (1/365)  # g/d
    # Calculate the coefficient that adjusts rates based on low O2, and take the larger of that value or assumed minimul
    # mmDOepi = max((DOmaxProp*StatesEpiGas$DO[i-1]/LakeVolumes$EpiVol[i-1])/(kDO+StatesEpiGas$DO[i-1]/LakeVolumes$EpiVol[i-1]),mmDOmin)

    ###############
    # Calculate sediment P recycling rates
    # Applied separately to Epi and Hypo
    # Two models: (1) Nurmberg equation; (2) First order decay
    # Model (1)
    # Next find the recycling rate for Nurnberg equation
    # Current mass P per mass sed, used for both Epi and Hypo
    #TPsed = StatesSedP$TP[i-1]*(1/SedBulkDensity)*(1/SedimentVolume)*1000 #mg of P/ g of soil
    #                = gP * (m^3Soil/1.33e6g) * 1/m^3 * 1000mg/g
    # This ends up being nearly equivalent to SedAvailPinit, which was used for the
    # start condition of the sediments. 
    # Nurmburg equation, do not allow negative
    #RatesSedP$TPcFsedP[i] = max((TPa + TPb * (TPsed)) / 1000,0) # g/m^2*d
    
    # Model (2)
    RatesEpiP$TPRecycling[i] = PropSedInEpi * (StatesSedP$TP[i-1] * TPcRecycle *
        TPRecyclingTheta^(LakeTemps$EpiT[i]-TPRecyclingBaseT)) #recycling in g/day

    # Epilimnetic recycling
    # RatesEpiP$TPRecycling[i] = 0
    # Model (3) is transfer coefficient
    # RatesEpiP$TPRecycling[i] = PropSedInEpi * (TPkRecycle/1000 * (StatesSedP$TP[i-1]/(LakeAreas$EpiArea[i-1]*ActiveSedDepth) - StatesEpiP$TP[i-1]/LakeVolumes$EpiVol[i-1]) * 
    # LakeAreas$EpiArea[i-1] * TPRecyclingTheta^(LakeTemps$EpiT[DOY]-TPRecyclingBaseT) ) #recycling in g/day
    
    ###############
    # Calculate sediment P release and rebinding rates
    if (StatesEpiGas$DO[i-1]/LakeVolumes$EpiVol[i-1] <= AnoxiaThreshold ){ # Anoxic
      # Release
      RatesEpiP$TPRelease[i] = PropSedInEpi * StatesSedP$TPB[i-1] * TPcRelease *
        TPRecyclingTheta^(LakeTemps$EpiT[i]-TPRecyclingBaseT) #recycling in g/day
      RatesEpiP$TPRebind[i] = 0
    }else{ # Oxic
      # Rebind if room available in bound P pool
      # Calculate static seds available for binding
      AvailSedBoundEpi = max(0,PropSedInEpi * (SedTPBoundMaxgrams - StatesSedP$TPB[i-1]))  # in grams
      RebindEpi = PropSedInEpi * (StatesEpiP$TP[i-1] * TPcRebind) *
        TPSedimentationTheta^(LakeTemps$EpiT[i]-TPRecyclingBaseT) # Potential rebinding
      RatesEpiP$TPRebind[i] = min(RebindEpi,AvailSedBoundEpi)
      RatesEpiP$TPRelease[i] = 0
    }
    
    # Respiration of POC
    sedPOCLRepi = PropSedInEpi * StatesSedC$POCL[i-1] * RPOCLSed * RTheta^(LakeTemps$EpiT[i]-Tbase) * mmDOepi # g C
    sedPOCRRepi = PropSedInEpi * StatesSedC$POCR[i-1] * RPOCRSed * RTheta^(LakeTemps$EpiT[i]-Tbase) * mmDOepi # g C

    # Sediment gas exchange, already takes into account epi area proportion
    RatesEpiGas$DOSed[i]  = (sedPOCLRepi + sedPOCRRepi) * O2toCMR # g O2

    # Burial of POC
    # Calculate according to parameter
    POCParamBurial = SedAvailOCBurial/1000 * SedBuried
    # Calculate according to available POC in sediments
    POCLEmpBurial = SoilDeposition * (1/365) * (1/1000) * StatesSedC$POCL[i-1]/ActiveSedDepth  #Burial term in g of P/d
    POCREmpBurial = SoilDeposition * (1/365) * (1/1000) * StatesSedC$POCR[i-1]/ActiveSedDepth
    POCEmpBurialTotal = POCLEmpBurial + POCREmpBurial
    # Choose the lesser of the two
    # POCBurialActual = min(POCParamBurial,POCEmpBurialTotal)
    POCBurialActual = POCEmpBurialTotal
    # Determine relative proportions of POC species in sediments
    POCRproportion = StatesSedC$POCR[i-1] / (StatesSedC$POCR[i-1]+StatesSedC$POCL[i-1])
    # Apply burial to the two POC species
    # The following calculation was necessary when 2 different burial types were possible
    # RatesSedC$POCRBurial[i] = POCBurialActual * POCRproportion # Burial in g of C/d for recalcitrant
    # RatesSedC$POCLBurial[i] = POCBurialActual * (1-POCRproportion) # Burial in g of C/d for recalcitrant
    # The following assumes empirical burial
    RatesSedC$POCRBurial[i] = POCREmpBurial # Burial in g of C/d for recalcitrant
    RatesSedC$POCLBurial[i] = POCLEmpBurial # Burial in g of C/d for recalcitrant

    # Calculate permanent burial of the sediments
    # Do the same for phosphorus
    PParamBurial = SedAvailPBurial/1000 * SedBuried
    # Now calculate directly from data
    PEmpBurial = SoilDeposition * (1/365) * (1/1000) * StatesSedP$TP[i-1]/ActiveSedDepth #Burial term in g of P/d
    PEmpBurialBound = SoilDeposition * (1/365) * (1/1000) * StatesSedP$TPB[i-1]/ActiveSedDepth #Burial term in g of P/d
    # Total empirical burial
    PEmpBurialTotal = PEmpBurial+PEmpBurialBound
    PEmpBurialBoundProp = PEmpBurialBound / PEmpBurialTotal
    # Use the lesser of the two
    # print('_______________')
    # print(PEmpBurialTotal/myLake$Area)
    # print(PParamBurial/myLake$Area)
    # PBurialActual = min(PParamBurial,PEmpBurialTotal)
    PBurialActual = PEmpBurialTotal
    #PBurialActual = PEmpBurialTotal
    # Divid into relative proportion of bound and not bound
    # The following calculation was necessary when 2 different burial types were possible
    # RatesSedP$TPBurial[i]  = PBurialActual * (1-PEmpBurialBoundProp)
    # RatesSedP$TPBBurial[i] = PBurialActual * (PEmpBurialBoundProp) 
    # The following assignments assumes the empirical burial
    RatesSedP$TPBurial[i]  = PEmpBurial
    RatesSedP$TPBBurial[i] = PEmpBurialBound 

    ###############
    # Hypo rates
    #if (DOY >myLake$StratBegin & DOY <=myLake$StratEnd){ # Stratified, note ">" for begin, which allows indexing to work
    if (isStratified & DayOfStratification > 1){ # Stratified, note ">" for begin, which allows indexing to work
      # Calculate Michaelis Menten scaler (0-1)
      mmDOhypo = max((DOmaxProp*StatesHypoGas$DO[i-1]/LakeVolumes$HypoVol[i-1])/(kDO+StatesHypoGas$DO[i-1]/LakeVolumes$HypoVol[i-1]),mmDOmin)
      RatesHypoC$DOCLR[i] = StatesHypoC$DOCL[i-1] * RDOCL * RTheta^(LakeTemps$HypoT[i]-Tbase) * mmDOhypo # g C
      RatesHypoC$DOCRR[i] = StatesHypoC$DOCR[i-1] * RDOCR * RTheta^(LakeTemps$HypoT[i]-Tbase) * mmDOhypo # g C
      RatesHypoC$POCLR[i] = StatesHypoC$POCL[i-1] * RPOCL * RTheta^(LakeTemps$HypoT[i]-Tbase) * mmDOhypo # g C
      RatesHypoC$POCRR[i] = StatesHypoC$POCR[i-1] * RPOCR * RTheta^(LakeTemps$HypoT[i]-Tbase) * mmDOhypo # g C
      RatesHypoC$POCLSettling[i] = min(StatesHypoC$POCL[i-1],StatesHypoC$POCL[i-1] * POCLSettling / MeanHypoThickness) # g C
      RatesHypoC$POCRSettling[i] = min(StatesHypoC$POCR[i-1],StatesHypoC$POCR[i-1] * POCRSettling / MeanHypoThickness) # g C
      # There are two possible approaches to recycling
      #(1) Areal
      # RatesHypoP$TPRecycling[i] = LakeAreas$HypoArea[i-1] * RatesSedP$TPcFsedP[i] *
      #   TPRecyclingTheta^(LakeTemps$HypoT[DOY]-TPRecyclingBaseT) #recycling in g/day
      #(2) Volumetric
      RatesHypoP$TPRecycling[i] = PropSedInHypo * (StatesSedP$TP[i-1] * TPcRecycle * #mmDOhypo *
         TPRecyclingTheta^(LakeTemps$HypoT[i]-TPRecyclingBaseT)) * mmDOhypo #recycling in g/day
      #(3) Piston velocity areal, positive is out of the sediments
      # RatesHypoP$TPRecycling[i] = PropSedInHypo * (TPkRecycle/1000 * (StatesSedP$TP[i-1]/(LakeAreas$HypoArea[i-1]*ActiveSedDepth) - StatesHypoP$TP[i-1]/LakeVolumes$HypoVol[i-1]) * 
        #LakeAreas$HypoArea[i-1] * TPRecyclingTheta^(LakeTemps$HypoT[DOY]-TPRecyclingBaseT) ) #recycling in g/day
      
      ###############
      # Calculate sediment P release and rebinding rates
      if (StatesHypoGas$DO[i-1]/LakeVolumes$HypoVol[i-1] <= AnoxiaThreshold ){ # Anoxic
        # Release
        RatesHypoP$TPRelease[i] = PropSedInHypo * StatesSedP$TPB[i-1] * TPcRelease *
          TPRecyclingTheta^(LakeTemps$HypoT[i]-TPRecyclingBaseT) #recycling in g/day
        RatesHypoP$TPRebind[i] = 0
      }else{ # Oxic
        # Rebind if room available in bound P pool
        # Calculate static seds available for binding
        AvailSedBoundHypo = max(0,PropSedInHypo * (SedTPBoundMaxgrams - StatesSedP$TPB[i-1]))  # in grams
        RebindHypo = PropSedInHypo * (StatesHypoP$TP[i-1] * TPcRebind) *
          TPSedimentationTheta^(LakeTemps$HypoT[i]-TPRecyclingBaseT) # Potential rebinding
        RatesHypoP$TPRebind[i] = min(RebindHypo,AvailSedBoundHypo)
        RatesHypoP$TPRelease[i] = 0
      }
    
      # Sediment respiration in the hypo
      sedPOCLRhypo = PropSedInHypo * StatesSedC$POCL[i-1] * RPOCLSed * RTheta^(LakeTemps$HypoT[i]-Tbase) * mmDOhypo # g C
      sedPOCRRhypo = PropSedInHypo * StatesSedC$POCR[i-1] * RPOCRSed * RTheta^(LakeTemps$HypoT[i]-Tbase) * mmDOhypo # g C
      
      # Carbon respiration already accounts for T and low DO above
      RatesHypoGas$DOSed[i] = (sedPOCLRhypo + sedPOCRRhypo) * O2toCMR # g O2
      RatesEpiGas$DOSed[i]  = (sedPOCLRepi + sedPOCRRepi) * O2toCMR # g O2
      
    }else{ # Mixed
      # Calculate Michaelis Menten scaler (0-1)
      sedPOCLRhypo = 0
      sedPOCRRhypo = 0
      RatesHypoC$DOCLR[i] = 0
      RatesHypoC$DOCRR[i] = 0
      RatesHypoC$POCLR[i] = 0
      RatesHypoC$POCRR[i] = 0
      RatesHypoC$POCLSettling[i] = 0
      RatesHypoC$POCRSettling[i] = 0
      RatesHypoP$TPSettling[i] = 0
      RatesSedP$TPcFsedP[i] = 0
      RatesHypoP$TPRecycling[i] = 0
      RatesHypoP$TPRebind[i] = 0
      RatesHypoP$TPRelease[i] = 0
      RatesHypoGas$DOSed[i] = 0
    }
    
    # Add the Epi and Hypo sedimentation
    RatesSedC$POCLR[i] = sedPOCLRepi + sedPOCLRhypo # g C 
    RatesSedC$POCRR[i] = sedPOCRRepi + sedPOCRRhypo # g C 

    ###############
    # Mixing
    # Check for change in stratification
    # Adjust rates as necessary
    dEpiVol = LakeVolumes$EpiVol[i] - LakeVolumes$EpiVol[i-1]
    MixProp = dEpiVol/LakeVolumes$EpiVol[i-1] # Proport lost due to entrainment
    RatesEpiC$DOCLMix[i] = 0 # If neither of the two following cases are true
    RatesEpiC$DOCRMix[i] = 0 # If neither of the two following cases are true
    RatesEpiC$POCLMix[i] = 0 # If neither of the two following cases are true
    RatesEpiC$POCRMix[i] = 0 # If neither of the two following cases are true
    RatesEpiGas$DOMix[i] = 0 #
    RatesEpiP$TPMix [i] = 0
    if (MixProp!=0){ # If there is mixing
      if (MixProp<0){ # Epi loses water and mass
        # Lose these masses from states from the epi
        RatesEpiC$DOCLMix[i] = StatesEpiC$DOCL[i-1] * (abs(MixProp))
        RatesEpiC$DOCRMix[i] = StatesEpiC$DOCR[i-1] * (abs(MixProp))
        RatesEpiC$POCLMix[i] = StatesEpiC$POCL[i-1] * (abs(MixProp))
        RatesEpiC$POCRMix[i] = StatesEpiC$POCR[i-1] * (abs(MixProp)) 
        RatesEpiP$TPMix[i]   = StatesEpiP$TP[i-1]   * (abs(MixProp))
        RatesEpiGas$DOMix[i] = StatesEpiGas$DO[i-1] * (abs(MixProp))
      }else{ # Epi gains water, hypo loses water
        # Calculate in terms of hypolimnion - easier to understand
          dHypoVol = LakeVolumes$HypoVol[i] - LakeVolumes$HypoVol[i-1]
          MixProp = dHypoVol/LakeVolumes$HypoVol[i-1] # Proport lost due to entrainment
          # Lose this from the Hypo
          RatesEpiC$DOCLMix[i] = -1*(StatesHypoC$DOCL[i-1] * (abs(MixProp)))
          RatesEpiC$DOCRMix[i] = -1*(StatesHypoC$DOCR[i-1] * (abs(MixProp)))
          RatesEpiC$POCLMix[i] = -1*(StatesHypoC$POCL[i-1] * (abs(MixProp)))
          RatesEpiC$POCRMix[i] = -1*(StatesHypoC$POCR[i-1] * (abs(MixProp)))
          RatesEpiP$TPMix[i]   = -1*(StatesHypoP$TP[i-1]   * (abs(MixProp)))
          RatesEpiGas$DOMix[i] = -1*(StatesHypoGas$DO[i-1] * (abs(MixProp)))
      }
    }

    ##################
    # Mass balance
    
    # Epilimnion, always exists, whether mixed or stratified
    StatesEpiC$DOCL[i] = StatesEpiC$DOCL[i-1] + RatesEpiC$DOCLInflow[i] + RatesEpiC$GPP[i]*(1-GPP2POC) - 
      RatesEpiC$DOCLOutflow[i] - RatesEpiC$DOCLMix[i]  - RatesEpiC$DOCLR[i]
    
    StatesEpiC$DOCR[i] = StatesEpiC$DOCR[i-1] + RatesEpiC$DOCRInflow[i] - 
      RatesEpiC$DOCROutflow[i] - RatesEpiC$DOCRMix[i] - RatesEpiC$DOCRR[i]
    
    StatesEpiC$POCL[i] = StatesEpiC$POCL[i-1] + RatesEpiC$POCLInflow[i] + RatesEpiC$GPP[i]*GPP2POC - 
      RatesEpiC$POCLR[i] - RatesEpiC$POCLOutflow[i] - RatesEpiC$POCLMix[i]  - RatesEpiC$POCLSettling[i]
    
    StatesEpiC$POCR[i] = StatesEpiC$POCR[i-1] + RatesEpiC$POCRInflow[i] - 
      RatesEpiC$POCRR[i] - RatesEpiC$POCROutflow[i] - RatesEpiC$POCRMix[i]  - RatesEpiC$POCRSettling[i]
    
    StatesEpiP$TP[i] = StatesEpiP$TP[i-1] + RatesEpiP$TPInflow[i]*(1-myLake$PLoadAsSed) + 
      RatesEpiP$TPRecycling[i] + RatesEpiP$TPRelease[i] -
      RatesEpiP$TPOutflow[i] - RatesEpiP$TPMix[i] - RatesEpiP$TPSettling[i] - RatesEpiP$TPRebind[i]

    StatesEpiGas$DO[i] = StatesEpiGas$DO[i-1] + RatesEpiC$GPP[i]*O2toCMR + RatesEpiGas$DOFatm[i] - RatesEpiGas$DOMix[i] -
      RatesEpiC$DOCLR[i]*O2toCMR - RatesEpiC$DOCRR[i]*O2toCMR - RatesEpiC$POCLR[i]*O2toCMR - RatesEpiC$POCRR[i]*O2toCMR - 
      RatesEpiGas$DOSed[i]
    
    # Hypolimnion, only when stratified, indexing is tricky
    # if (DOY >=myLake$StratBegin & DOY <=myLake$StratEnd){ # Stratified, note ">" for begin, which allows indexing to work
    if (isStratified){ # Stratified, note ">" for begin, which allows indexing to work
      # if (DOY==myLake$StratBegin){ # On the first day of stratification, add the epi DOC
      if (DayOfStratification == 1){ # On the first day of stratification, add the epi DOC
        StatesHypoC$DOCL[i] = RatesEpiC$DOCLMix[i]
        StatesHypoC$DOCR[i] = RatesEpiC$DOCRMix[i]
        StatesHypoC$POCL[i] = RatesEpiC$POCLMix[i]
        StatesHypoC$POCR[i] = RatesEpiC$POCRMix[i]
        StatesHypoP$TP[i] = RatesEpiP$TPMix[i]
        StatesHypoGas$DO[i] = RatesEpiGas$DOMix[i]
        # To not lose continuity in sediments
        StatesSedC$POCR[i] = StatesSedC$POCR[i-1] 
        StatesSedC$POCL[i] = StatesSedC$POCL[i-1] 
      }else{ # Do regular mass balance
        StatesHypoC$DOCL[i] = StatesHypoC$DOCL[i-1] + RatesEpiC$DOCLMix[i] - RatesHypoC$DOCLR[i]
        StatesHypoC$DOCR[i] = StatesHypoC$DOCR[i-1] + RatesEpiC$DOCRMix[i] - RatesHypoC$DOCRR[i]
        StatesHypoC$POCL[i] = StatesHypoC$POCL[i-1] + RatesEpiC$POCLMix[i] + RatesEpiC$POCLSettling[i] - RatesHypoC$POCLR[i] - RatesHypoC$POCLSettling[i]
        StatesHypoC$POCR[i] = StatesHypoC$POCR[i-1] + RatesEpiC$POCRMix[i] + RatesEpiC$POCRSettling[i] - RatesHypoC$POCRR[i] - RatesHypoC$POCRSettling[i] 
        StatesHypoP$TP[i] = StatesHypoP$TP[i-1] + RatesEpiP$TPMix[i] + RatesHypoP$TPRecycling[i] +
          + RatesHypoP$TPRelease[i] - RatesHypoP$TPRebind[i]
        
        StatesHypoGas$DO[i] = StatesHypoGas$DO[i-1] + RatesEpiGas$DOMix[i] - 
          RatesHypoC$DOCLR[i]*O2toCMR - RatesHypoC$DOCRR[i]*O2toCMR - RatesHypoC$POCLR[i]*O2toCMR - RatesHypoC$POCRR[i]*O2toCMR - RatesHypoGas$DOSed[i]
        
        # Sediment OC receives only hypo
        RatesSedC$POCRSettling[i] = RatesHypoC$POCRSettling[i]
        RatesSedC$POCLSettling[i] = RatesHypoC$POCLSettling[i]
        StatesSedC$POCR[i] = StatesSedC$POCR[i-1] + RatesSedC$POCRSettling[i] - RatesSedC$POCRR[i] - RatesSedC$POCRBurial[i]
        StatesSedC$POCL[i] = StatesSedC$POCL[i-1] + RatesSedC$POCLSettling[i] - RatesSedC$POCLR[i] - RatesSedC$POCLBurial[i]
      }
    }else{ # Mixed
      StatesHypoC$DOCL[i] = 0
      StatesHypoC$DOCR[i] = 0
      StatesHypoC$POCL[i] = 0
      StatesHypoC$POCR[i] = 0
      StatesHypoP$TP[i] = 0
      StatesHypoGas$DO[i] = 0
      # Sediment OC receives only epi settling
      RatesSedC$POCRSettling[i] = RatesEpiC$POCRSettling[i]
      RatesSedC$POCLSettling[i] = RatesEpiC$POCLSettling[i]
      StatesSedC$POCR[i] = StatesSedC$POCR[i-1] + RatesSedC$POCRSettling[i] - RatesSedC$POCRR[i] - RatesSedC$POCRBurial[i]
      StatesSedC$POCL[i] = StatesSedC$POCL[i-1] + RatesSedC$POCLSettling[i] - RatesSedC$POCLR[i] - RatesSedC$POCLBurial[i]
    }

    # Notes about this mass balance: We assume epiPsettling goes straight to sediments; no hypo P settling
    # Sum inflow settling and epi settling, epi and hypo recycling for convenience
    RatesSedP$TPSettling[i] = RatesEpiP$TPInflow[i]*(myLake$PLoadAsSed) + RatesEpiP$TPSettling[i]
    RatesSedP$TPRecycling[i] = RatesHypoP$TPRecycling[i] + RatesEpiP$TPRecycling[i]
    
    RatesSedP$TPRelease[i] = + RatesEpiP$TPRelease[i] + + RatesHypoP$TPRelease[i]
    RatesSedP$TPRebind[i]  = + RatesEpiP$TPRebind[i] + + RatesHypoP$TPRebind[i]

    StatesSedP$TP[i] = StatesSedP$TP[i-1] + RatesEpiP$TPInflow[i]*(myLake$PLoadAsSed) + RatesEpiP$TPSettling[i] - 
      RatesHypoP$TPRecycling[i] - RatesEpiP$TPRecycling[i] - RatesSedP$TPBurial[i]
    StatesSedP$TPB[i] = StatesSedP$TPB[i-1] + RatesSedP$TPRebind[i] - RatesSedP$TPRelease[i] - RatesSedP$TPBBurial[i]
    
    # Do not allow states to go negative
    StatesHypoGas$DO[i] = max(StatesHypoGas$DO[i],0)
    
    # Secchi
    LEC = LECwater + LECPOCL*StatesEpiC$POCL[i]/LakeVolumes$EpiVol[i] + LECPOCR*StatesEpiC$POCR[i]/LakeVolumes$EpiVol[i] + 
      LECDOCR*StatesEpiC$DOCR[i]/LakeVolumes$EpiVol[i] + LECDOCL*StatesEpiC$DOCL[i]/LakeVolumes$EpiVol[i]
    StatesEpiC$Secchi[i] = LEC2Secchi/LEC
  }
  
  ####################
  # Return results
  
  MetabolismResults = list()
  MetabolismResults[[1]] = StatesEpiC
  MetabolismResults[[2]] = StatesHypoC
  MetabolismResults[[3]] = StatesEpiGas
  MetabolismResults[[4]] = StatesHypoGas
  MetabolismResults[[5]] = RatesEpiC
  MetabolismResults[[6]] = RatesHypoC
  MetabolismResults[[7]] = RatesEpiGas
  MetabolismResults[[8]] = RatesHypoGas
  MetabolismResults[[9]] = StatesEpiP
  MetabolismResults[[10]] = StatesHypoP
  MetabolismResults[[11]] = StatesSedP
  MetabolismResults[[12]] = RatesEpiP
  MetabolismResults[[13]] = RatesHypoP
  MetabolismResults[[14]] = RatesSedP
  MetabolismResults[[15]] = StatesSedC
  MetabolismResults[[16]] = RatesSedC
  MetabolismResults[[17]] = LakeVolumes
  MetabolismResults[[18]] = LakeAreas
  MetabolismResults[[19]] = LakeTemps
  
  close(pb) # Closes the progress bar

  return(MetabolismResults)
}
