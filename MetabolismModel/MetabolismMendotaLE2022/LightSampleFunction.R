# Sample light distribution function

library(pracma)
# source('MetabolismHelpers.R')
# source('Metabolism.R')
# source('PlotMetabResults.R')
# source('PlotLongTermResults.R')
# source('SummarizeMetabResults.R')
# source('ExploreLakePhosphorus.R')
# source('SaveAsLakeFile.R')
# source('./LakePhysics/1D_HeatMixing/1D_HeatMixing_functions.R')
# # Load lake metabolizer for dissolved gas functions
# library(LakeMetabolizer)
# library(tidyverse)
# 
# HypsoFiles = data.frame(hypsofile = 'bc/LakeEnsemblR_bathymetry_standard.csv')
# hyps_all <- get_hypsography(hypsofile = toString(HypsoFiles[i]),
#                             dx = dx, nx = nx)
# LakeHypsometry = data.frame(depth = 1:length(hyps_all[[1]]), area = hyps_all[[1]], volume = hyps_all[[3]])

# z1 = top of the layer you are calculating for (surface = 0)
#   z2 = bottom depth of layer
# a1 = upper area of the layer (surface)
# a2 = area of the bottom of the layer
# ir = incident light at top of layer
# mu = total LEC
ir = 1000
mu = seq(0,2,0.1)
myLight = rep(NA,length(mu))
myLight2 = rep(NA,length(mu))

z1 = 0
z2 = 5
a1 = 39850000
a2 = 0
#a2 = 30000000

for (i in 1:length(mu)){
  myLight[i] = light_calc(z1,z2,a1,a2,ir,mu[i])
  myLight2[i] = light_calc_simple(z1,z2,a1,a2,ir,mu[i])
}

#par(mfrow=c(2,1), mai = c(0.4,0.8, 0.3, 0.1),cex = 0.9)
plot(mu,myLight,type='l',xlab = 'mu',ylab = 'light',ylim = c(0,ir))
lines(mu,myLight2,type='l',xlab = 'mu',ylab = 'light', col='red')
legend(legend=c('Robert','Simple (I at Z)'),'topright',lty=c(1,1),col=c('black','red'),cex=0.7)

light_calc <- function(z1,z2,a1,a2,ir,mu){
    deps = seq(z1, z2 , 0.1)
    irr =  ir * exp(-(mu)*deps)
    area = approx(c(z1,z2), c(a1,a2), deps)$y
    intarea = pracma::trapz((deps),(area))
    
    # print(range(irr))
    # print(mean(irr))
    # 
    # print((area %*% irr)/ (sum(area)))
    
    # sum(pracma::cumtrapz(rev(area),rev(irr)))/a1
    
    # return(pracma::trapz(rev(area),rev(irr))/a1)
    return((area %*% irr)/ (sum(area)))
  }

light_calc_simple <- function(z1,z2,a1,a2,ir,mu){
  # light = log(ir) - mu* z2
  # light = exp(light)
  light = ir*exp(-mu*z2)
  return(light)
}

# PP_calc <- function(I,PBmax,alpha){
#   #PB(I)=PBmaxtanh(𝛼𝐁IPBmax)+RB
#   PB_I = I * IP
#   PB_I2 = Pmax(1-exp(-IP*I/Pmax))
# }
