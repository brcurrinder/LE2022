# Wrapper file for the metabolism model
# Paul C Hanson, October 2020

#install.packages('neonstore',repos="https://CRAN.R-project.org")
#remotes::install_github("aemon-j/gotmtools")
source('MetabolismHelpers.R')
source('Metabolism.R')
source('PlotMetabResultsLastSim.R')
source('PlotLongTermResults.R')
source('PlotMetabResultsAllSims.R')
source('PlotMetabResultsAllSimsAnnualSummaries.R')
source('PlotMetabResultsHypoBudget.R')
source('SummarizeMetabResults.R')
source('ExploreLakePhosphorus.R')
source('SaveAsLakeFile.R')
# source('./LakePhysics/1D_HeatMixing/1D_HeatMixing_functions_metabolism.R')
source('./LakePhysics/1D_HeatMixing/1D_LakeModel_functions.R')
source('ComparePredictionsObservations.R')
source('PlotLongTermRecovery.R')
# Load lake metabolizer for dissolved gas functions
library(LakeMetabolizer)
library(tidyverse)

########################################
# Simulation Setup
# Can run multiple lakes in series to reflect
# changing lake condition, e.g., nutrient loads
# Simulations 
# Each lake, it's associated parameter set, met data, start and end dates

# the number of iterations through loops depends first on length(LakeIDs)
# LakeIDs = c(3,2,3) #,3) # Can be a vector, lakes loaded from text file
LakeIDs = c(2,4)
LakeIDs = 2
ParameterIDs = c(2,2) # Can be a vector, params loaded from text file
nReruns = c(1,20) # Can be a vector, If re-running the same inputs but continuing states
# ComparePredictionsObservations(OutputFileListName)

# Input data
# StartDates = rep("2000-01-01 UTC",3) # StartDates = rep("1995-01-01 UTC",3) # StartDates = rep("2009-01-01 UTC",3)
# EndDates   = rep("2009-12-31 UTC",3) # EndDates   = rep("2015-12-31 UTC",3) # EndDates   = rep("2013-12-31 UTC",3)
# StartDates = rep("2000-01-01 UTC",3) # StartDates = rep("1995-01-01 UTC",3) # StartDates = rep("2009-01-01 UTC",3)
# EndDates   = rep("2007-12-31 UTC",3) # EndDates   = rep("2015-12-31 UTC",3) # EndDates   = rep("2013-12-31 UTC",3)
# StartDates = c("1995-01-01 UTC","2010-01-01") 
# EndDates   = c("2015-12-31 UTC","2015-12-31")
StartDates = c("1995-01-01 UTC","2014-01-01") 
EndDates   = c("2020-12-31 UTC","2015-12-31")


InitTemps  = rep(3,3) # Can be a vector, Initial lake temperature on, e.g., Jan 1
HypsoFiles = rep('bc/LakeEnsemblR_bathymetry_standard.csv',3)
# MeteoFiles = rep('bc/LakeEnsemblR_meteo_standard.csv',3)
# ___________________________________________________________________________________________
# Note: the script "CleanMetFile.R" was used to create this cleaned met file
# MeteoFiles = c('bc/NLDAS2_Mendota_1979_2016_cell_5_GLMReady_cut_timezonecorr_CLEANED.csv',
#                'bc/NLDAS2_Mendota_2010_2016_cell_5_GLMReady_cut_timezonecorr_CLEANED.csv',
#                'bc/NLDAS2_Mendota_2010_2016_cell_5_GLMReady_cut_timezonecorr_CLEANED.csv')
#MeteoFiles = rep('bc/NLDAS2_Mendota_1979_2016_forKludg_CLEANED.csv',3) #before 12 Jan 2022
MeteoFiles = rep('bc/NLDAS2_Mendota_1995_2020_forKludg_DWH_12jan22.csv',3) 

# ___________________________________________________________________________________________
# HydroFiles = rep('DriverData/hydro_inputs.csv',3)
HydroFiles = rep('DriverData/hydro_inputs_1979_2020_CopyMethod.csv',3)
SecchiFiles =rep('bc/light.csv',3) # placeholder needed for Robert's function

# Setup for some of the output
TurnDay = 55 # Can be a vector, this is for graphing purposes and plots days following TurnDay back toward the origin
PrintSummary = TRUE # Whether to print the summary to screen
PlotResults = FALSE # Whether to plot the time series of results
PlotLTResults = TRUE # Whether to plot the long term results
SaveResults = TRUE # whether to save results to disk
nYears2Save = 1 # The final n years to be saved
OutputDir = './Output'

# Quick check on setup of nLakes, params, and input data
nLakes = length(LakeIDs)
if (length(ParameterIDs)<nLakes | length(nReruns)<nLakes | length(InitTemps)<nLakes |
    length(HypsoFiles)<nLakes | length(MeteoFiles)<nLakes | length(SecchiFiles)<nLakes){
  print('Length mismatch in inputs and number of lakes')
  stop()
}

########################################
# Some setup for the physical model
# For now, these should not change
nx = 25 # number of layers we will have
dt = 3600*2 # 24 hours times 60 min/hour times 60 seconds/min, * 2 changes time step to 2 hours
hydrodynamic_timestep = 12 * dt # for the lake physics model; hard code as 86400 as long as physics is for 1 day at a time
WindFactor = 0.8 # Energy transfer adjustment
SWFactor = 1 # SW adjustment
########################################
# Load the lakes
# First load all pre-defined lake configurations
AllLakes = GetLakeConfig(NA)

########################################
# Run the models

# Setup lists for output
MetabResults = list() # For multiple lakes
MetabSummaries = list() # For multiple lakes
myLakes = list()
OutputFileNames = data.frame(Names=rep(NA,sum(sum(nReruns[1:length(LakeIDs)]))))
FileNameCounter = 0
TotalptmSim = proc.time()

# For each lake in the list of LakeIDs, run the model
for (i in 1:length(LakeIDs)){
  Totalptm = proc.time()
  LakeID = LakeIDs[i]
  StartDate = StartDates[i]
  startingDate = StartDate # global variable for the physics model
  EndDate = EndDates[i]
  nDays = as.numeric(difftime(as.Date(EndDate), as.Date(StartDate), units = 'days'))
  SimulationDates = as.Date(StartDate):as.Date(EndDate)
  SimulationDates = as.Date(SimulationDates,origin = '1970-01-01')
  
  ########################################
  # Load parameters
  myParameters = GetParameters(ParameterIDs[i])
  # Save simulation setup informaion
  SimSetup = data.frame(nDays=nDays,nReruns=nReruns[i],LakeID=LakeID,ParameterIDs=myParameters,
                        nYears2Save=nYears2Save,StartDate=StartDate,EndDate=EndDate)
  ########################################
  # Get Lake Configuration
  myLake = AllLakes[LakeIDs[i],]
  myLakes[[i]] = myLake
  print('__________________________________________________')
  print('<><><><><><><><><><><><><><><><><><><><><><><><><>')
  print(paste('Running lake ',i,' of ',length(LakeIDs),', ',myLake$LakeName),sep="")
  print('------------------------------------------')

  ## Some temporary hard coded stuff ro Lake Mendota
  zmax = myLake$Zmax # maximum lake depth
  dx = zmax/nx # spatial step
  
  # Load hypsography from file
  hyps_all <- get_hypsography(hypsofile = toString(HypsoFiles[i]),
                              dx = dx, nx = nx)
  # Place in a structure
  LakeHypsometry = data.frame(depth = 1:length(hyps_all[[1]]), area = hyps_all[[1]], volume = hyps_all[[3]])
  # Define initial temperature profile
  # u_ini <- initial_profile(initfile = 'bc/obs.txt', nx = nx, dx = dx,
                           # depth = hyps_all[[3]],
                           # processed_meteo = meteo_all[[1]])
  u_ini = rep(InitTemps[i],as.numeric(dim(LakeHypsometry)[1]))
  
  ########################################
  # Load atmospheric boundary conditions and hydrology
  LoadMetData = TRUE
  # but do not reload if it's already loaded
  if (i>1){
    if (StartDate==StartDates[i-1] & EndDate==EndDates[i-1] &
            MeteoFiles[i]==MeteoFiles[i-1]){
      LoadMetData = FALSE
    }
  }
  # LoadMetData = FALSE
  
  if (LoadMetData){ # and hydrology data
    #########################
    # Meteorlogy
    print('Loading atmospheric driving data...')
    meteo_all <- provide_meteorology(meteofile = toString(MeteoFiles[i]),
                                     secchifile = toString(SecchiFiles[i]), 
                                     windfactor = WindFactor)
    # Subset meteo_all to the simulation date range for this lake
    # meteo_all is a list of 2, each with it's unique date field
    # Get dates from the first list and truncate to day (remove time)
    d = as.POSIXct(meteo_all[[1]]$datetime,format = '%Y-%m-%d')
    d2 = format(d,format='%Y-%m-%d')
    iSimDates = which(d2>=as.POSIXct(StartDate,format = '%Y-%m-%d') & d2<=as.POSIXct(EndDate,format = '%Y-%m-%d'))
    # Now truncate the first list
    newList = meteo_all[[1]][iSimDates,]
    meteo_all[[1]] = meteo_all[[1]][iSimDates,]
    # Repeat for the second list
    d = as.POSIXct(meteo_all[[2]]$sampledate,format = '%Y-%m-%d')
    d2 = format(d,format='%Y-%m-%d')
    iSimDates = which(d2>=as.POSIXct(StartDate,format = '%Y-%m-%d') & d2<=as.POSIXct(EndDate,format = '%Y-%m-%d'))
    meteo_all[[2]] = meteo_all[[2]][iSimDates,]
    
    # Summarize to daily means
    myIrradiance = data.frame(myDate=meteo_all[[1]]$datetime,I=meteo_all[[1]]$Shortwave_Radiation_Downwelling_wattPerMeterSquared)
    myIrradiance = aggregate(x=myIrradiance[c("I")],FUN=mean,
                             by = list(myDays = as.Date(myIrradiance$myDate,format="%Y-%m-%d")))
    # Adjust irradiance by factor (as in wind)
    myIrradiance$I = myIrradiance$I * SWFactor
    
    myWindSpeed = data.frame(myDate=meteo_all[[1]]$datetime,u10=meteo_all[[1]]$Ten_Meter_Elevation_Wind_Speed_meterPerSecond)
    myWindSpeed = aggregate(x=myWindSpeed[c("u10")],FUN=mean,
                            by = list(myDays = as.Date(myWindSpeed$myDate,format="%Y-%m-%d")))
    #########################
    # Hydrology
    print('Loading hydrology data...')
    FlowFactor = 0.6
    Hydro_all <- read.csv(HydroFiles[i],)
    Hydro_all$total_inflow_volume = Hydro_all$total_inflow_volume*FlowFactor
    # Subset meteo_all to the simulation date range for this lake
    # Get dates from the first list and truncate to day (remove time)
    d = as.POSIXct(Hydro_all$time,format = '%Y-%m-%d')
    d2 = format(d,format='%Y-%m-%d')
    iSimDates = which(d2>=as.POSIXct(StartDate,format = '%Y-%m-%d') & d2<=as.POSIXct(EndDate,format = '%Y-%m-%d'))
    # Now truncate the first list
    Hydro_all = Hydro_all[iSimDates,]
    Discharge = data.frame(Date = Hydro_all$time,Flow = Hydro_all$total_inflow_volume)
  }else{
    print("Re-using driver data from previous lake")
  }
  
  # Setup the loads, which can change if we change lake file
  OCLoad = (myLake$OCLoad/365 * Discharge$Flow) / mean(Discharge$Flow) * myLake$Area
  TPLoad = (myLake$PLoad/365 * Discharge$Flow) / mean(Discharge$Flow) * myLake$Area
  
  ########################################
  # Run the model
  for (j in 1:nReruns[i]){
    # Following is a hack
    # if (i==2 & j==1){
    #   # Adjusting flow for second lake
    #   print('Adjusting flow as hack')
    #   Discharge$Flow = Discharge$Flow / 0.6
    # }
    Iterationptm = proc.time()
    print('================================================')
    print(paste('Lake ',i,' of ',length(LakeIDs),', iteration ',j,' of ',nReruns[i],sep=""))
    # Create outputfile name
    FileNameCounter = FileNameCounter+1
    OutputFileNames$Names[FileNameCounter] = paste(OutputDir,'/',trimws(myLake$LakeName),'_',i,'_',j,'.RData',sep="")
    LoadFromState=FALSE # Assume this is the case unless otherwise set
    StateFile = ""
    if (FileNameCounter>1){
      LoadFromState=TRUE
      StateFile = OutputFileNames$Names[FileNameCounter-1]
    }
    MetabResults[[1]] = Metabolism(nDays,nYears2Save,myLake,myParameters,myIrradiance,myWindSpeed,
                                   Discharge,OCLoad,TPLoad,LakeHypsometry,LoadFromState,StateFile,Verbose=TRUE)
    
    TestPlots = FALSE
    if (TestPlots){
      par(mfrow=c(3,1), mai = c(0.7,0.8, 0.2, 0.1),cex = 0.9)
      # GPP
      plot(MetabResults[[1]][[3]]$SimDay,MetabResults[[1]][[5]]$GPP/MetabResults[[1]][[17]]$EpiVol,type='l',
           xlab = 'Simulation day', ylab = 'NPP (gC/m3/d)')
      # Plot irradiance for NPP
      plot(1:nDays,(MetabResults[[1]][[19]]$MeanEpiIrradiance),type='l',xlab='Sim day',ylab='mean epi I (mmol/m2/s)')
      plot(MetabResults[[1]][[19]]$MeanEpiIrradiance,MetabResults[[1]][[5]]$GPP/MetabResults[[1]][[17]]$EpiVol,
           xlab='I (mmol/m2/s)',ylab='NPP (gC/m3/d)')
      
      # Dissolved oxygen
      par(mfrow=c(3,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
      plot(MetabResults[[1]][[3]]$SimDay,MetabResults[[1]][[3]]$DO/MetabResults[[1]][[17]]$EpiVol,type='l',ylim=c(0,15),
           xlab = 'Simulation day', ylab = 'DO (mg/L)')
      lines(MetabResults[[1]][[4]]$SimDay,MetabResults[[1]][[4]]$DO/MetabResults[[1]][[17]]$HypoVol,type='l',col='red')
      # Phosphorus
      plot(MetabResults[[1]][[9]]$SimDay,MetabResults[[1]][[9]]$TP/MetabResults[[1]][[17]]$EpiVol,type='l',#ylim=c(0.01,0.02),
           xlab = 'Simulation day', ylab = 'TP (ug/L)')
      lines(MetabResults[[1]][[10]]$SimDay,MetabResults[[1]][[10]]$TP/MetabResults[[1]][[17]]$HypoVol,type='l',col='red')
      # Water temperatures
      plot(MetabResults[[1]][[10]]$SimDay,MetabResults[[1]][[19]]$EpiT,type='l',xlab='Sim day',ylab='T')
      lines(MetabResults[[1]][[10]]$SimDay,MetabResults[[1]][[19]]$HypoT,type='l',col='red')
      lines(MetabResults[[1]][[10]]$SimDay,MetabResults[[1]][[19]]$IceThickness*10,type='l',col='blue')
      par(mfrow=c(2,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
      # Hypolimnetic volume as proportion of max
      plot(MetabResults[[1]][[10]]$SimDay,MetabResults[[1]][[17]]$HypoVol/max(MetabResults[[1]][[17]]$HypoVol),type='l',
           xlab='Sib day',ylab='Hypo volume / max(hypo volume)')
      # Ratio of hypolimnetic sediment area to hypolimnetic volume
      plot(MetabResults[[1]][[10]]$SimDay,MetabResults[[1]][[18]]$HypoArea/MetabResults[[1]][[17]]$HypoVol,ylim=c(0,1),type='l',
           xlab = 'Sim day',ylab='Hypo sed area: hypo vol')
    }
    ########################################
    # Summarize results
    MetabSummaries[[1]] = SummarizeMetabResults(MetabResults[[1]],myLake,LakeHypsometry,nYears2Save,PrintSummary)
    # MetabSummaries[[1]] = SummarizeMetabResultsFullSim(MetabResults[[1]],myLake,LakeHypsometry,nYears2Save,PrintSummary)
    #MetabSummaries[[Lake]][[1]]$States*,[[2]]$Rates*,[[3]]$Massbalance*

    ########################################
    # Plot results
    if (PlotResults){
      PlotMetabResultsLastSim(MetabResults[[1]],myLake,LakeHypsometry)
    }
    ########################################
    # Write results to disk
    if (SaveResults){
      ModelRun = list()
      ModelRun[[1]] = SimSetup
      ModelRun[[2]] = myLake
      ModelRun[[3]] = myParameters
      ModelRun[[4]] = LakeHypsometry
      ModelRun[[5]] = NULL
      ModelRun[[6]] = NULL
      ModelRun[[7]] = MetabResults[[1]]
      ModelRun[[8]] = MetabSummaries[[1]]
      save(ModelRun,file=OutputFileNames$Names[FileNameCounter])
    }
    RunTime = proc.time()-Iterationptm
    print('------------------------------------------')
    print(paste('Iteration run time: ',round(RunTime[3]/60,digits=3),' min',sep=""))
  } # End of reruns
  RunTime = proc.time()-Totalptm
  print('------------------------------------------')
  print(paste('Lake run time: ',round(RunTime[3]/60,digits=3),' min',sep=""))
}
RunTime = proc.time()-TotalptmSim
print('------------------------------------------')
print(paste('Total run time: ',round(RunTime[3]/60,digits=3),' min',sep=""))

# Plot results across iterations
if (PlotLTResults){
  # OutputFileNames = read.csv('./Output_1000y_eutrophic_equilibrium/Mendota_1980_Present_2_100.csv')
  # OutputFileNames = read.csv('./Output/Mendota_No_P_Load_2_60.csv')
  # OutputFileNames = read.csv('./Output/Mendota_zero_Pload_2_20.csv')
  # PlotLongTermResults(OutputFileNames,TurnDay=NA)
  # PlotLongTermRecovery(OutputFileNames,TurnDay=NA,TRUE)
  PlotMetabResultsLastSim(MetabResults[[1]],myLake,LakeHypsometry)
  # PlotMetabResultsHypoBudget(MetabResults[[1]],myLake,SimSetup,LakeHypsometry,2002)
  # PlotMetabResultsAllSims(OutputFileNames)
  # PlotMetabResultsAllSimsAnnualSummaries(OutputFileNames)
  # To test any one simulation output:
    # OutputFileListName = "./Output_1000y_eutrophic_equilibrium/Mendota_1980_Present_2_100.csv"
  # ComparePredictionsObservations(OutputFileListName)
  # PlotAnnualCycles(OutputFileListName,StartDate,EndDate,"eutrophic")
  # PlotAnnualCycles(OutputFileListName,"2005-01-01 UTC","2006-12-31 UTC","eutrophic")
}

# Save final state values, which can be handy for long burn-in runs
nFiles = dim(OutputFileNames)[1]
SaveAsLakeFile(OutputFileNames$Names[nFiles],'./Output/0_NewState.csv',LakeHypsometry)

# Save the file list for possible manual reloading later
save(ModelRun,file=OutputFileNames$Names[FileNameCounter])
OutputFileListName = paste(OutputDir,'/',trimws(myLake$LakeName),'_',i,'_',j,'.csv',sep="")
write.csv(OutputFileNames,OutputFileListName,quote=FALSE,row.names=FALSE)


GenerateCoolDistributions = FALSE
if (GenerateCoolDistributions){
  myDO = MetabResults[[1]][[3]]$DO/MetabResults[[1]][[17]]$EpiVol
  myDOsat = MetabResults[[1]][[3]]$DOsat
  Variable = myDO - myDOsat
  i2009 = which(as.POSIXct(SimulationDates)>=2008)
  DensityByYear(SimulationDates[i2009],Variable[i2009],Label)
}

#plot(density(Variable),xlim=c(-4,4),main="",ylab='Modeled density (dev sat) all years')

plot(MetabResults[[1]][[11]]$SimDay, MetabResults[[1]][[11]]$TP/(myLake$Area*myLake$SedAreaProp),type='l',
     ylab='Sed TP unbound (gP/m2SA)')
plot(MetabResults[[1]][[11]]$SimDay, MetabResults[[1]][[11]]$TPB/(myLake$Area*myLake$SedAreaProp),
     type='l',ylab='Sed TP bound (gP/m2Sa)')

# Print carbon budget based on mean values

OC_autoch_total = (sum(MetabResults[[1]][[5]]$GPP,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365

OC_alloch_total = OC_export_total = (sum(MetabResults[[1]][[5]]$DOCLInflow,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365 +
  (sum(MetabResults[[1]][[5]]$DOCRInflow,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365 +
  (sum(MetabResults[[1]][[5]]$POCLInflow,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365 +
  (sum(MetabResults[[1]][[5]]$POCLInflow,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365

OC_export_total = (sum(MetabResults[[1]][[5]]$DOCLOutflow,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365 +
  (sum(MetabResults[[1]][[5]]$DOCROutflow,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365 +
  (sum(MetabResults[[1]][[5]]$POCLOutflow,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365 +
  (sum(MetabResults[[1]][[5]]$POCLOutflow,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365

# Approximate sedimentation + burial, because we do not have net difference in water column
OC_netsedimentation = OC_autoch_total+OC_alloch_total - OC_export_total

OC_sedimentR = (sum(MetabResults[[1]][[5]]$GPP,na.rm=TRUE)/myLake$Area)/length((MetabResults[[1]][[5]]$GPP/myLake$Area)) * 365

