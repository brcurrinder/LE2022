ReadME for metabolism model

This folder contains the metabolism model Paul initially shared to the google drive
The main model is in the folder named 'MetabolismMendotaLE2022'
The other folder 'PaulNLDAS2' is a folder containing NDLAS2 data used to extend our meteoroloyg data through 2020
The script LineUpMeteo_DWH.R is the script DWH used to append the new NDLAS data to the prior version and is now set up to use in the metab model
Mendota_ancillary.R is a script that provides edi link and download insturctions for various LTER data sets
MetabDataExploration_DWH.R is a script DWH has used to explore inital metabolism model output data
seasons_from_Robins_paper.csv - is a csv containing the table defining seasons from Robin's microbial dataset 


TO RUN METAB MODEL
within R studio 
Set working directory to the root of the 'MetabolismMendotaLE2022' folder
Run the script, “MetabolsimWrapper.R”
Model is set to run about 20 years, which will take ~30 minutes, depending on your computer
To save model outputs; run the script 'WriteResultsMatrix.csv'
after running this script you can run these two lines of code to save the output as a csv 
1 - source('WriteResultsMatrix.R') #reads in function needed if you didn't run the script  
getwd()  #need to make sure working directory is in main metab folder, where Output is folder, thats where the function will write to
WriteResultsMatrix = function(OutputFileListName) #Output File List name needs to be the csv name generated from Wrapper; this should be generated after running the model 
#This gives csv named SimResultsMatrix.csv in Output folder, I renamed to date ran and my initals to put on google drive and read in below 
