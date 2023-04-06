# Need to run as administrator to use RSToolbox

# PCA Analysis of AVIRIS-NG scene

# Load libraries
library(sp)
library(raster)
library(RStoolbox)

# Load raster data
setwd('C:/Users/ssharp/Box/UCD EDL/GLEON LE2022/Data/AVIRIS/Lake_clipped')
filename <- 'ang20160831t201002_rfl_v1n2_img_clipped_lake'
AVIRIS <- stack(filename)

rPCA <- rasterPCA(AVIRIS,
          nSamples=NULL,
          nComp=length(AVIRIS@layers),
          spca=FALSE,
          maskCheck=TRUE)
