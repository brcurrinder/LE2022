# Test various phosphorus assumptions for a lake
ExploreLakePhosphorus <- function(myLake,myParameters,DesiredSedimentPAreal){
  # This function calls SedimentEquilibrium to do the 4 sediment calculations.
  # There are four lake characteristics that are important and highly uncertain
  # and that likely covary by lake trophic state. They are:
  # Phosphorus load (PLoad), Depth of active sediments (ActiveSedDepth), 
  # Available P in sediments (SedAvailPinit), Soil depostion rate (SoilDeposition)
  myLakeTest = myLake
  
  nTests = 1000 # Number of samples
  nBest = 100 # Number of best results to print
  nParams2Test = 4 # Number of parameters to test
  Range = c(0.25,2.0) # Range of each parameter as proportion of default
  
  # Randomly sample parameter space
  Tests = data.frame(PLoad=rep(NA,nTests),ActiveSedDepth=rep(NA,nTests),SoilDeposition=rep(NA,nTests))
  Tests$PLoad = runif(nTests,myLake$PLoad*Range[1],myLake$PLoad*Range[2])
  Tests$ActiveSedDepth = runif(nTests,myLake$ActiveSedDepth*Range[1],myLake$ActiveSedDepth*Range[2])
  Tests$SedAvailPinit = runif(nTests,myLake$SedAvailPinit*Range[1],myLake$SedAvailPinit*Range[2])
  Tests$SoilDeposition= runif(nTests,myLake$SoilDeposition*Range[1],myLake$SoilDeposition*Range[2])
  
  TAdjEpi = myParameters$TPSedimentationTheta^(mean(LakeTemps$EpiT)-myParameters$TPSedimentationBaseT)
  TAdjHypo = myParameters$TPRecyclingTheta^(mean(LakeTemps$HypoT)-myParameters$TPRecyclingBaseT)
  
  # Cycle through tests
  PSedEstimates = data.frame(EmpSeds=rep(NA,nTests),EquilSeds=rep(NA,nTests),
                          EmpLoads=rep(NA,nTests),EquilLoads=rep(NA,nTests))
  for (i in 1:nTests){
    myLakeTest$PLoad = Tests$PLoad[i]
    myLakeTest$ActiveSedDepth = Tests$ActiveSedDepth[i]
    myLakeTest$SedAvailPinit = Tests$SedAvailPinit[i]
    myLakeTest$SoilDeposition = Tests$SoilDeposition[i]
    # Get estimate of areal concentration of sediment P
    # Four estimates returned: 1) Empirical seds, 2) Equil from seds, 
    #     3) Empirical from loads, 4) Equil from loads
    PSedEstimates[i,] = SedimentEquilibrium(myLakeTest,myParameters,TAdjEpi,TAdjHypo)
  }
  PSedEstimates$Means = rowMeans(PSedEstimates[,1:nParams2Test])
  # Calculate range for each row
  range(PSedEstimates[,1:nParams2Test])
  for (i in 1:nTests){
    PSedEstimates$Range[i] = abs(diff(range(PSedEstimates[i,1:nParams2Test])))
  }
  # For indexing to match parameter ranges before sorting
  PSedEstimates$OriginalOrder = 1:nTests
  # Sort so we can pull the best results by range
  SortedEstimates = PSedEstimates[order(PSedEstimates$Range),]
  MeanSedEstimate = mean(SortedEstimates$Means[1:nBest])
  # Get best results
  inBest = SortedEstimates$OriginalOrder[1:nBest]
  BestVals = Tests[inBest,]
  # Figure out how to change parameters to get best result
  NewParameters = SedimentDepthOrAvailP(myLake,MeanSedEstimate)
  print('=======================================================')
  print('The following output is the result of random sampling over parameterr space.')
  print('The output are from combinations of parameters that give the tightest cluster')
  print('for estimates of sediment P concentration.')
  print(paste('Means of top',nBest,'values from random sampling (each their own units)'))
  print(colMeans(BestVals))
  print('Mean of the sediment concentrations of P from above parameter estimates (gP/m2LA)')
  print(MeanSedEstimate)
  print('___________')
  print('However, you should run the model to equilibrium and note the final sediment')
  print('phosphorus concentration. The desired P value you passed to this function was:')
  print(paste('    ',DesiredSedimentPAreal,' gP/m2LA',sep=""))
  print('Given that desired sediment P, CHANGE ONLY 1 of these parameters to the new value')
  print(SedimentDepthOrAvailP(myLake,DesiredSedimentPAreal))
  print('=======================================================')
  
  # Plot results 
  par(mfrow=c(2,1), mai = c(0.8,0.8,0.5, 0.1),cex = 0.9)
  myYLim = range(BestVals)
  plot(SortedEstimates$Range[1:nBest],BestVals$PLoad,ylim=myYLim,type='l',col='black',
       xlab='Range among tests',ylab='Values of parameters (mixed units)',
       main=paste('Param values for',nBest,'best tests'))
  lines(SortedEstimates$Range[1:nBest],BestVals$ActiveSedDepth*10,col='green')
  lines(SortedEstimates$Range[1:nBest],BestVals$SoilDeposition,col='blue')
  lines(SortedEstimates$Range[1:nBest],BestVals$SedAvailPinit,col='red')
  LegendText = c("PLoad","SedDepth*10","SoilDeposition","SedAvailPinit")
  legend("topright",LegendText,lty=c(1,1,1,1),col=c('black','green','blue','red'),cex=0.7)
  # Same but rank order
  par(mai = c(0.8,0.8, 0.1, 0.1),cex = 0.9)
  plot(1:nBest,BestVals$PLoad,ylim=myYLim,type='l',col='black',
       xlab='Rank order',ylab='Values of parameters (mixed units)',main="")
  lines(1:nBest,BestVals$ActiveSedDepth*10,col='green')
  lines(1:nBest,BestVals$SoilDeposition,col='blue')
  lines(1:nBest,BestVals$SedAvailPinit,col='red')
  LegendText = c("PLoad","SedDepth*10","SoilDeposition","SedAvailPinit")
  legend("topright",LegendText,lty=c(1,1,1,1),col=c('black','green','blue','red'),cex=0.7)
  
  # Check for correlation among test results
  # plot(PSedEstimates[,1:nParams2Test],main="Correlation between test results")
  
  # Plot the test results as a function of the range in results
  par(mfrow=c(1,1), mai = c(1.0,0.8, 0.1, 0.1),cex = 0.9)
  myYLim = range(PSedEstimates[,1:nParams2Test])
  plot(PSedEstimates$Range,PSedEstimates[,1],ylim=myYLim,xlab='Range in results',ylab='Sed P (g/m2)')
  points(PSedEstimates$Range,PSedEstimates[,2],col='green')
  points(PSedEstimates$Range,PSedEstimates[,3],col='red')
  points(PSedEstimates$Range,PSedEstimates[,4],col='blue')
  points(PSedEstimates$Range,PSedEstimates$Means,pch=15,col='black')
  legend("topleft",c('EmpSed','EqSed','EmpLoad','EqLoad','Mean'),pch=c(1,1,1,1,15),col=c('black','green','red','blue','black'))
  
}
