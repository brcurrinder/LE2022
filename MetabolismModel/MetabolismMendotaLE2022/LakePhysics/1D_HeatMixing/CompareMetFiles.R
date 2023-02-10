# Script to compare two meteolological driving files

File1 = 'bc/LakeEnsemblR_meteo_standard.csv'
File2 = "bc/NLDAS2_Mendota_1979_2016_cell_5_GLMReady_cut_timezonecorr_CLEANED.csv"

StartDate = "2009-01-01 UTC"
EndDate   = "2013-12-31 UTC"
nColsCompare = 8
par(mfrow=c(2,1))

Met1 = read_csv(File1)
Met2 = read_csv(File2)

# Subset files to sim dates
d = as.POSIXct(Met1$datetime,format = '%Y-%m-%d')
d2 = format(d,format='%Y-%m-%d')
iSimDates = which(d2>=as.POSIXct(StartDate,format = '%Y-%m-%d') & d2<=as.POSIXct(EndDate,format = '%Y-%m-%d'))
Met1 = Met1[iSimDates,]

d = as.POSIXct(Met2$datetime,format = '%Y-%m-%d')
d2 = format(d,format='%Y-%m-%d')
iSimDates = which(d2>=as.POSIXct(StartDate,format = '%Y-%m-%d') & d2<=as.POSIXct(EndDate,format = '%Y-%m-%d'))
Met2 = Met2[iSimDates,]

for (i in 2:nColsCompare){
  plot(Met1$datetime,Met1[[i]],type='l',ylab = colnames(Met1[i]))
  plot(Met2$datetime,Met2[[i]],type='l',ylab = colnames(Met1[i]))
}
