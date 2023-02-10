PlotAnnualCycles = function(OutputFileListName,StartDate=NA,EndDate=NA,TrophicState){
  # Code for converting R output to a huge text file of the simulation
  
  # Set print parameters depending on screen or png
  PrintIt = TRUE
  if (PrintIt){
    myLWD = 0.5
    myLTY2 = 1
    myLTYabline = 3
  }else{
    myLWD = 1
    myLTY2 = 1
    myLTYabline = 2
  }
  # Set dates for visualization
  if (is.na(StartDate)){
    StartDate = "2000-01-01 UTC" #Earliest start date
  }
  if (is.na(EndDate)){
    EndDate   = "2005-12-31 UTC" #Latest end date
  }
  
  # ComparePredictionsObservations
  # For now, use just the last file in the file list that comes from OutputFileListName
  # Inputs
  #    CSV file containing the file list of the simulation output
  #    Date range over which to make the comparison
  
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
  
  # Use only the last simulation
  for (i in nFiles:nFiles){
    #load(OutputFileNames$Names[i])
    load(myFileList$Names[i])
    
    #Load time series results
    MR = ModelRun[[7]]
    Hypso = ModelRun[[4]]
    
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
    #LakeAreas     = rbind(LakeVolumes,MR[[18]])
    LakeTemps     = rbind(LakeTemps,MR[[19]])
    
    LakeAreas = data.frame(TotalArea = max(LakeHypsometry))
    
  }
  
  # get simulation date range
  SimDateRange = c(min(as.Date(ModelRun[[1]]$StartDate)),max(as.Date(ModelRun[[1]]$EndDate)))
  # generate sequence of dates for, e.g., plotting
  SimDates = seq(as.Date(SimDateRange[1]),as.Date(SimDateRange[2]),by='day')
  # remove first value, as output begins on day 2
  SimDates = SimDates[2:length(SimDates)]
  
  # Subset to date range
  iPlotDates = which(as.POSIXct(SimDates)>=as.POSIXct(StartDate,format = '%Y-%m-%d') & 
                       as.POSIXct(SimDates)<=as.POSIXct(EndDate,format = '%Y-%m-%d'))
  
  # Get stratified periods
  CurrentlyStratified = FALSE
  S1=NULL
  S2=NULL
  myHypoVol = LakeVolumes$HypoVol[iPlotDates]
  for (iS in 1:length(myHypoVol)){
    if (myHypoVol[iS]>0 & !CurrentlyStratified){
      CurrentlyStratified = TRUE
      S1 = c(S1,iS)
    }
    if (myHypoVol[iS]==0 & CurrentlyStratified){
      CurrentlyStratified = FALSE
      S2 = c(S2,iS)
    }
  }
  myDates = SimDates[iPlotDates]

  #########################
  # Calculate Ice
  # use epi Y limits for plotting ice
  myYLim = c(0,30)

  myIce = LakeTemps$IceThickness
  myIce[!LakeTemps$IceCovered] = NA
  myIce[LakeTemps$IceCovered] = myYLim[2]
  # Calculate ice durations
  uYears = unique(format(as.POSIXct(myDates),format="%Y"))
  IceEachYear = data.frame(uYears=uYears,IceCoveredDays=rep(NA,length(uYears)),
                           IceOnDate=rep(as.Date('1990-01-01', format = "%Y-%m-%d"),length(uYears)),
                           IceOffDate=rep(as.Date('1990-01-01', format = "%Y-%m-%d"),length(uYears)),
                           JulianDayOn=rep(NA,length(uYears)),
                           JulianDayOff=rep(NA,length(uYears)))
  myIceCovered = FALSE
  IceCoveredDays = 0
  iYears = 1
  for (iC in 1:length(LakeTemps$IceCovered)){
    # Convert to julian day and use that instead of 180 next
    if (LakeTemps$IceCovered[iC] & iC > 180){
      # skip the first year of ice cover because it's partial season, thus the 180
      myIceCovered = TRUE
      IceCoveredDays = IceCoveredDays+1
      if (iYears<length(uYears) & IceCoveredDays==1){
        IceEachYear$IceOnDate[iYears+1] = as.Date(SimDates[iC], format = "%Y-%m-%d")
        tmp <- as.Date(IceEachYear$IceOnDate[iYears+1], format = "%d%b%y")
        IceEachYear$JulianDayOn[iYears+1] = format(tmp, "%j")
      }
    }else{
      if (myIceCovered & iYears < length(uYears)){ # the first day of ice off!
        # Check to make sure it's not another short ice period from same year
        ThisYear = format(as.Date(SimDates[iC], format = "%Y-%m-%d"),'%Y')
        LastSavedYear = format(IceEachYear$IceOffDate[iYears],'%Y')
        if (ThisYear==LastSavedYear){
          # Just another ice period for this year
        }else{
          IceEachYear$IceCoveredDays[iYears+1] = IceCoveredDays
          IceEachYear$IceOffDate[iYears+1] = as.Date(SimDates[iC], format = "%Y-%m-%d")
          tmp <- as.Date(IceEachYear$IceOffDate[iYears+1], format = "%d%b%y")
          IceEachYear$JulianDayOff[iYears+1] = format(tmp, "%j")
          myIceCovered = FALSE
          IceCoveredDays = 0
          iYears = iYears + 1
        }
      }
    }
  }

  nDaysAvgIceDuration = mean(IceEachYear$IceCoveredDays,na.rm=TRUE)
  sdIceDuration = sd(IceEachYear$IceCoveredDays,na.rm=TRUE)
  meanJulianOn = mean(as.numeric(IceEachYear$JulianDayOn),na.rm=TRUE)
  meanJulianOff = mean(as.numeric(IceEachYear$JulianDayOff),na.rm=TRUE)
  # print ice results to screen
  print('ICCCCCCCCCCEEEEEEE')
  print(paste('Mean ice duration: ',nDaysAvgIceDuration,sep=""))
  print(paste('Std  ice duration: ',round(sdIceDuration),sep=""))
  print(paste('Mean ice on   DOY: ',round(meanJulianOn),sep=""))
  print(paste('Mean ice off  DOY: ',round(meanJulianOff),sep=""))
  print('Ice cover by year...')
  print(IceEachYear)
  
  # # Plots, LTER ice data loaded here
  ObsIce = read_csv('./DriverData/NTLIceDataForComparison.csv')
  # par(mfrow=c(2,1),lend=2,mai = c(0.05,0.75, 0.4, 0.05),oma = c(2,1,0.2,0.2), cex = 0.8)
  # plot(IceEachYear$uYears[1:length(uYears)-1],IceEachYear$IceCoveredDays[2:length(uYears)],type='b',pch=19,ylim=c(0,140),
  #      main="",ylab='Ice covered days',col='blue')
  # lines(ObsIce$Year,ObsIce$IceDays,lty=2, type="b", col='blue',pch=1)
  # legend(legend=c('Simulated','Observed'),'bottomleft',
  #        lty=c(2,2),pch=c(19,1,1,1),col=c('blue','blue','blue','red','blue'),cex=0.7)
  # plot(IceEachYear$uYears[1:length(uYears)-1],IceEachYear$JulianDayOff[2:length(uYears)],type='b',pch=19,ylim=c(0,140),
  #      main="",ylab='Ice off day of year',col='red')
  # lines(ObsIce$Year,ObsIce$IceOffDOY,lty=2, type="b", col='red',pch=1)
  # legend(legend=c('Simulated','Observed'),'bottomleft',
  #        lty=c(2,2),pch=c(19,1,1,1),col=c('red','red','blue','red','blue'),cex=0.7)
  
  plot(IceEachYear$uYears,IceEachYear$IceCoveredDays,type='b',pch=19,ylim=c(0,140),
       main="",ylab='Ice covered days',col='blue')
  lines(ObsIce$Year,ObsIce$IceDays,lty=2, type="b", col='blue',pch=1)
  legend(legend=c('Simulated','Observed'),'bottomleft',
         lty=c(2,2),pch=c(19,1,1,1),col=c('blue','blue','blue','red','blue'),cex=0.7)
  plot(IceEachYear$uYears,IceEachYear$JulianDayOff,type='b',pch=19,ylim=c(0,140),
       main="",ylab='Ice off day of year',col='red')
  lines(ObsIce$Year,ObsIce$IceOffDOY,lty=2, type="b", col='red',pch=1)
  legend(legend=c('Simulated','Observed'),'bottomleft',
         lty=c(2,2),pch=c(19,1,1,1),col=c('red','red','blue','red','blue'),cex=0.7)
  
  
  # Make the plots
  ################################
  # Temperature, DO, and stratification
  # par(mfrow=c(3,1),mai = c(0.05,0.05, 0.05, 0.05),oma = c(2,2,0.2,0.2), cex = 0.5)
  par(mfrow=c(3,1),lend=2,mai = c(0.05,0.35, 0.08, 0.05),oma = c(2,1,0.2,0.2), cex = 0.4)
  FileName = paste('./Figures/AnnualCycleTemp.png',sep="")
  
  # epilimnion
  plot(SimDates[iPlotDates],LakeTemps$EpiT[iPlotDates],type='l', ylim=myYLim,
       main="",xlab="",ylab="T,DO,DOsat",xaxt='n')
  lines(SimDates[iPlotDates],myIce[iPlotDates],type='l',lwd=4)
  lines(SimDates[iPlotDates],(StatesEpiGas$DO[iPlotDates])/LakeVolumes$EpiVol[iPlotDates],
        col='blue')
  lines(SimDates[iPlotDates],(StatesEpiGas$DOsat[iPlotDates]),lty=3,col='blue')
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  # hypolimnion
  myYLim = c(0,20)
  plot(SimDates[iPlotDates],LakeTemps$HypoT[iPlotDates],type='l', ylim=myYLim,
       main="",xlab="",ylab="hypo",xaxt='n')
  lines(SimDates[iPlotDates],(StatesHypoGas$DO[iPlotDates])/LakeVolumes$HypoVol[iPlotDates],
        col='blue')
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  # thermocline depth
  myYLim = c(0,25)
  plot(SimDates[iPlotDates],LakeVolumes$ThermoclineDepth[iPlotDates],type='l', ylim=myYLim,
       main="",xlab="",ylab="Z depth")
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  
  if (PrintIt){
    dev.copy(png,FileName,width=2,height=2.5,units="in",res=300)
    dev.off()
  }
  ################################
  # Phosphorus states
  par(mfrow=c(3,1),mai = c(0.05,0.35, 0.08, 0.05),oma = c(2,1,0.2,0.2), cex = 0.4)
  FileName = paste('./Figures/AnnualCycleTPStates.png',sep="")
  
  # Epi
  myYLim = c(0,200)
  if (toupper(TrophicState)!='EUTROPHIC'){
    myYLim = c(0,20)
  }
  plot(SimDates[iPlotDates],StatesEpiP$TP[iPlotDates]/LakeVolumes$EpiVol[iPlotDates]*1000,type='l', ylim=myYLim,
       main="",xlab="",ylab="TP (ug/L)",xaxt='n')
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  # Phosphorus hypolimnion
  myYLim = c(0,2000)
  if (toupper(TrophicState)!='EUTROPHIC'){
    myYLim = c(0,20)
  }
  plot(SimDates[iPlotDates],StatesHypoP$TP[iPlotDates]/LakeVolumes$HypoVol[iPlotDates]*1000,type='l', ylim=myYLim,
       main="",xlab="",ylab="hypo",xaxt='n')
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  # Sediments
  myYLim = c(10,50)
  if (toupper(TrophicState)!='EUTROPHIC'){
    myYLim = c(0,40)
  }
  plot(SimDates[iPlotDates],StatesSedP$TP[iPlotDates]/LakeArea,type='l', ylim=myYLim,
       main="",xlab="",ylab="Sed (gP/m2LA)")
  lines(SimDates[iPlotDates],StatesSedP$TPB[iPlotDates]/LakeArea,lty=2, ylim=myYLim,
       main="",xlab="",ylab="")
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
 
  if (PrintIt){
    dev.copy(png,FileName,width=2,height=2.5,units="in",res=300)
    dev.off()
  }
  # Phosphorus rates
  par(mfrow=c(3,1),mai = c(0.05,0.35, 0.08, 0.05),oma = c(2,1,0.2,0.2), cex = 0.4)
  FileName = paste('./Figures/AnnualCycleTPRates.png',sep="")
  
  # Epilimnion
  myYLim = c(-20,20)
  if (toupper(TrophicState)!='EUTROPHIC'){
    myYLim = c(-1,1)
  }
  # Sources
  plot(SimDates[iPlotDates],RatesEpiP$TPInflow[iPlotDates]/LakeAreas$TotalArea*365,col='black',ylim=myYLim,
       type='l',xlab='',ylab='P-epi rates (gP/m2LA/y)',xaxt='n',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],RatesEpiP$TPRecycling[iPlotDates]/LakeAreas$TotalArea*365,col='red',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],RatesEpiP$TPRelease[iPlotDates]/LakeAreas$TotalArea*365,col='blue',lty=myLTY2,lwd=myLWD)
  # lines(SimDates[iPlotDates],RatesEpiP$TPMix[iPlotDates]/LakeAreas$TotalArea*365,col='blue',lty=myLTY2,lwd=myLWD)
  # Sinks
  lines(SimDates[iPlotDates],-RatesEpiP$TPOutflow[iPlotDates]/LakeAreas$TotalArea*365,col='black',lty=myLTY2,lwd=myLWD)
  lines( SimDates[iPlotDates],-RatesEpiP$TPSettling[iPlotDates]/LakeAreas$TotalArea*365,col='red',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-RatesEpiP$TPRebind[iPlotDates]/LakeAreas$TotalArea*365,col='blue',lty=myLTY2,lwd=myLWD)
  abline(h=0,col='grey',lty=myLTYabline,lwd=myLWD)
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  if (!PrintIt){
    legend(legend=c('Inflow','Recycling','Release'),'topleft',
           lty=c(1,1,1,1,1,1),col=c('black','red','blue','black','red','blue'),cex=1.2)
    legend(legend=c('Outflow','Settling','Rebinding'),'bottomleft',
           lty=c(1,1,1,1,1,1),col=c('black','red','blue','black','red','blue'),cex=1.2)
           
  }
  # Hypolimnion
  myYLim = c(-20,20)
  if (toupper(TrophicState)!='EUTROPHIC'){
    myYLim = c(-1,1)
  }
  # Sources
  plot(SimDates[iPlotDates],(RatesHypoP$TPSettling[iPlotDates])/LakeAreas$TotalArea*365,
       ylim=myYLim,type='l',xlab='',ylab='P-Hypo rates (gP/m2LA/y)',main="",xaxt='n',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],RatesHypoP$TPRecycling[iPlotDates]/LakeAreas$TotalArea*365,col='red',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],RatesHypoP$TPRelease[iPlotDates]/LakeAreas$TotalArea*365,col='blue',lty=myLTY2,lwd=myLWD)
  # Sinks
  lines(SimDates[iPlotDates],-RatesHypoP$TPRebind[iPlotDates]/LakeAreas$TotalArea*365,col='blue',lty=myLTY2,lwd=myLWD)
  abline(h=0,col='grey',lty=myLTYabline,lwd=myLWD)
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  if (!PrintIt){
    legend(legend=c('Recycling','Release'),'topleft',
           lty=c(1,1,1,1,1,1),col=c('red','blue','blue'),cex=1.2)
    legend(legend=c('Rebind'),'bottomleft',
           lty=c(1,1,1,1,1,1),col=c('blue'),cex=1.2)
  }
  
  # Sediments
  # Sources
  plot( SimDates[iPlotDates],RatesSedP$TPSettling[iPlotDates]/LakeAreas$TotalArea*365,ylim=myYLim,type='l',
       xlab='',ylab='P-sed rates (gP/m2LA/y)',main="",col='red',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],RatesSedP$TPRebind[iPlotDates]/LakeAreas$TotalArea*365,col='blue',lty=myLTY2,lwd=myLWD)
  # Sinks
  lines(SimDates[iPlotDates],-RatesSedP$TPRecycling[iPlotDates]/LakeAreas$TotalArea*365,col='red',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-RatesSedP$TPRelease[iPlotDates]/LakeAreas$TotalArea*365,col='blue',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-(RatesSedP$TPBurial[iPlotDates]+RatesSedP$TPBurial[iPlotDates])/LakeAreas$TotalArea*365,
                              col='black',lty=myLTY2,lwd=myLWD)
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  abline(h=0,col='grey',lty=myLTYabline,lwd=myLWD)
  if (!PrintIt){
    legend(legend=c('Settling','Rebind'),'topleft',
           lty=c(1,1,1,1,1,1),col=c('red','blue','red','blue','black'),cex=1.2)
    legend(legend=c('Recycling','Release','Burial'),'bottomleft',
           lty=c(1,1,1,1,1,1),col=c('red','blue','black'),cex=1.2)
  }
  
  if (PrintIt){
    dev.copy(png,FileName,width=2,height=2.5,units="in",res=300)
    dev.off()
  }
  ######################
  # Organic carbon
  # Calc y lims
  myYLim = c(0,10)
  if (toupper(TrophicState)!='EUTROPHIC'){
    myYLim = c(0,10)
  }
  # par(mfrow=c(4,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
  par(mfrow=c(3,1),mai = c(0.05,0.35, 0.08, 0.05),oma = c(2,1,0.2,0.2), cex = 0.4)
  FileName = paste('./Figures/AnnualCycleOCStates.png',sep="")
  # Epi
  plot(SimDates[iPlotDates],(StatesEpiC$DOCR[iPlotDates])/LakeVolumes$EpiVol[iPlotDates],ylim=myYLim,type='l',
       ylab='DOCTot-epi (gC/m3)',main="",xaxt='n')
  lines(SimDates[iPlotDates],StatesEpiC$Secchi[iPlotDates],col='black')
  lines(SimDates[iPlotDates],(StatesEpiC$DOCR[iPlotDates]+StatesEpiC$POCR[iPlotDates])/LakeVolumes$EpiVol[iPlotDates],col='grey')
  lines(SimDates[iPlotDates],(StatesEpiC$DOCR[iPlotDates]+StatesEpiC$POCR[iPlotDates]+StatesEpiC$DOCL[iPlotDates])/LakeVolumes$EpiVol[iPlotDates],col='brown')
  lines(SimDates[iPlotDates],(StatesEpiC$DOCR[iPlotDates]+StatesEpiC$POCR[iPlotDates]+StatesEpiC$DOCL[iPlotDates]+StatesEpiC$POCL[iPlotDates])/LakeVolumes$EpiVol[iPlotDates],col='green')
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  if (!PrintIt){
    legend(legend=c('DOCR','POCR','DOCL','POCL'),'topleft',lty=c(1,1,1,1),col=c('black','grey','brown','green'),cex=0.7)
  }
  # Hypo
  plot(SimDates[iPlotDates],(StatesHypoC$DOCR[iPlotDates])/LakeVolumes$HypoVol[iPlotDates],
       ylim=myYLim,type='l',ylab='DOCTot-Hypo (gC/m3)',xaxt='n')
  lines(SimDates[iPlotDates],(StatesHypoC$DOCR[iPlotDates]+StatesHypoC$POCR[iPlotDates])/LakeVolumes$HypoVol[iPlotDates],col='grey')
  lines(SimDates[iPlotDates],(StatesHypoC$DOCR[iPlotDates]+StatesHypoC$POCR[iPlotDates]+StatesHypoC$DOCL[iPlotDates])/LakeVolumes$HypoVol[iPlotDates],col='brown')
  lines(SimDates[iPlotDates],(StatesHypoC$DOCR[iPlotDates]+StatesHypoC$POCR[iPlotDates]+StatesHypoC$DOCL[iPlotDates]+StatesHypoC$POCL[iPlotDates])/LakeVolumes$HypoVol[iPlotDates],col='green')
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  if (!PrintIt){
    legend(legend=c('DOCR','POCR','DOCL','POCL'),'topleft',lty=c(1,1,1,1),col=c('black','grey','brown','green'),cex=0.7)
  }
  # Sed
  myYLim = c(min(StatesSedC$POCR/(LakeAreas$TotalArea*myLake$SedAreaProp)),
             max(max(StatesSedC$POCR/(LakeAreas$TotalArea*myLake$SedAreaProp)+StatesSedC$POCL/(LakeAreas$TotalArea*myLake$SedAreaProp))))
  plot(SimDates[iPlotDates],(StatesSedC$POCR[iPlotDates])/(LakeAreas$TotalArea*myLake$SedAreaProp),ylim=myYLim,col='brown',type='l',ylab='DOCTot-Sed (gC/m2SA)')
  lines(SimDates[iPlotDates],(StatesSedC$POCR[iPlotDates]+StatesSedC$POCL[iPlotDates])/(LakeAreas$TotalArea*myLake$SedAreaProp),col='green')
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  if (!PrintIt){
    legend(legend=c('POCR','POCL'),'topleft',lty=c(1,1),col=c('brown','green'),cex=0.7) 
  }
  if (PrintIt){
    dev.copy(png,FileName,width=2,height=2.5,units="in",res=300)
    dev.off()
  }
  # Organic carbon rates
  # par(mfrow=c(4,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)

  # Calculate a bunch of stuff
  RSettling = (RatesSedC$POCRSettling)*365 / LakeArea
  LSettling = (RatesSedC$POCLSettling)*365 / LakeArea
  Settling = RSettling + LSettling
  PropRSettling = round(mean(RSettling/Settling,na.rm=TRUE),digits=2)
  PropLSettling = round(mean(LSettling/Settling,na.rm=TRUE),digits=2)
  RespirationSed = (RatesSedC$POCLR + RatesSedC$POCRR)*365 / LakeArea
  Burial = (RatesSedC$POCLBurial + RatesSedC$POCRBurial)*365 / LakeArea
  Net = Settling - RespirationSed - Burial
  SedimentCarbon = StatesSedC$POCR/LakeAreas$TotalArea+StatesSedC$POCL/LakeAreas$TotalArea # gC/m2LA
  SedO2Demand = (RespirationSed /12 *1000 /365)/myLake$SedAreaProp # mmol O2,C /m2SA/day
  RespirationHypoDOC = (RatesHypoC$DOCRR + RatesHypoC$DOCLR)*365/LakeArea
  RespirationHypoPOC = (RatesHypoC$POCRR + RatesHypoC$POCLR)*365/LakeArea

  par(mfrow=c(3,1),mai = c(0.05,0.35, 0.08, 0.05),oma = c(2,1,0.2,0.2), cex = 0.4)
  FileName = paste('./Figures/AnnualCycleOCRates.png',sep="")
  myYLim = c(-500,1000)
  if (toupper(TrophicState)!='EUTROPHIC'){
    myYLim = c(-100,100)
  }
  # Epi
  # sources
  plot(SimDates[iPlotDates],RatesEpiC$GPP[iPlotDates]*365/LakeArea,ylim=myYLim,
       col='green',type='l',ylab="Epi (gC/m2LA/y)",xaxt='n',lwd=myLWD)
  abline(h=0,col='grey',lty=myLTYabline,lwd=myLWD)
  lines(SimDates[iPlotDates],(RatesEpiC$DOCLInflow[iPlotDates]+RatesEpiC$DOCRInflow[iPlotDates]+
        RatesEpiC$POCRInflow[iPlotDates]+RatesEpiC$POCLInflow[iPlotDates])*365/LakeArea,
        col='black',lwd=myLWD)
  # sinks
  lines(SimDates[iPlotDates],-(RatesEpiC$DOCLOutflow[iPlotDates]+RatesEpiC$DOCROutflow[iPlotDates]+
        RatesEpiC$POCROutflow[iPlotDates]+RatesEpiC$POCLOutflow[iPlotDates])*365/LakeArea,
        col='black',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-(RatesEpiC$DOCLR[iPlotDates]+RatesEpiC$DOCRR[iPlotDates])*365/LakeArea,col='red',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-(RatesEpiC$POCLR[iPlotDates]+RatesEpiC$POCRR[iPlotDates])*365/LakeArea,col='blue',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-RatesEpiC$POCLSettling[iPlotDates]*365/LakeArea,col='green',lty=myLTY2,lwd=myLWD)
  abline(h=0,col='grey',lty=myLTYabline,lwd=myLWD)
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  if (!PrintIt){
    legend('topleft',legend=c('NPP','Inflow'),
         lty=c(1,1,1,1,1,1),col=c('green','black','black','red','blue','green'))
    legend('bottomleft',legend=c('Outflow','RDOC','RPOC','Settling'),
           lty=c(1,1,1,1,1,1),col=c('black','red','blue','green'))
  }
  # Hypo
  myYLim=c(-200,200)
  if (toupper(TrophicState)!='EUTROPHIC'){
    myYLim = c(-50,50)
  }
  # Settling to the hypo only occurs during stratification
  EpiSettlingStrat = (RatesEpiC$POCLSettling[iPlotDates]+RatesEpiC$POCRSettling[iPlotDates])
  iEpiSettlingStrat = which(myHypoVol==0)
  EpiSettlingStrat[iEpiSettlingStrat] = NA

  plot(SimDates[iPlotDates],EpiSettlingStrat*365/LakeArea,
       ylim=myYLim,col='green',type='l',ylab="Hypo (gC/m2LA/y)",xaxt='n',lwd=myLWD)
  lines(SimDates[iPlotDates],-(RatesHypoC$DOCLR[iPlotDates]+RatesHypoC$DOCRR[iPlotDates])*365/LakeArea,
        col='red',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-(RatesHypoC$POCLR[iPlotDates]+RatesHypoC$POCRR[iPlotDates])*365/LakeArea,
        col='blue',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-RatesHypoC$POCLSettling[iPlotDates]*365/LakeArea,
        col='green',lty=myLTY2,lwd=myLWD)
  abline(h=0,col='grey',lty=myLTYabline,lwd=myLWD)
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  if (!PrintIt){
    legend('topleft',legend=c('Settling'),
           lty=c(1,1,1,1,1,1),col=c('green','red','blue','green'))
    legend('bottomleft',legend=c('RDOC','RPOC','Settling'),
           lty=c(1,1,1,1,1,1),col=c('red','blue','green'))
  }
  
  # Sed
  plot(SimDates[iPlotDates],Settling[iPlotDates],ylim=myYLim,
       col='green',type='l',ylab='OCTot-Sed (gC/m2LA/y)',main="",lwd=myLWD)
  lines(SimDates[iPlotDates],-RespirationSed[iPlotDates],col='blue',lty=myLTY2,lwd=myLWD)
  lines(SimDates[iPlotDates],-Burial[iPlotDates],col='brown',lty=myLTY2,lwd=myLWD)
  # lines(SimDates[iPlotDates],Net[iPlotDates],col='grey',lty=myLTY2,lwd=2)
  abline(h=0,col='grey',lty=myLTYabline,lwd=myLWD)
  StratifiedShade(myDates[S1],myDates[S2],myYLim)
  if (!PrintIt){
    legend(legend=c('Settling'),'topleft',lty=c(1,1,1,2),col=c('green','blue','brown'),cex=0.7)
    legend(legend=c('Resp','Burial'),'bottomleft',lty=c(1,1,1,2),col=c('blue','brown'),cex=0.7)
  }
  
  if (PrintIt){
    dev.copy(png,FileName,width=2,height=2.5,units="in",res=300)
    dev.off()
  }
  
}

StratifiedShade <- function(S1,S2,myYLim){
  for (iMix in 1:length(S1)){
    #polygon(c(1,2,2,1),c(99.2,99.2,100,100),col=rgb(0, 0, 0,0.2))
    polygon(c(S1[iMix],S2[iMix],S2[iMix],S1[iMix]),c(min(myYLim),min(myYLim),max(myYLim),max(myYLim)),
            col=rgb(0, 0, 0,0.1),border=NA)
    polygon(c(as.POSIXct(S1[iMix]),as.POSIXct(S2[iMix]),as.POSIXct(S2[iMix]),as.POSIXct(S1[iMix])),c(min(myYLim),min(myYLim),max(myYLim),max(myYLim)),
            col=rgb(0, 0, 0,0.1),border=NA)
  }
}
