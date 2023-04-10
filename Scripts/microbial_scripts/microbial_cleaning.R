# Cleaning of 16S microbial data in limony dataset
# November 2022, cleaned up March 2023
# AGS

pacman::p_load(tidyverse, devtools, lubridate, limony)

# Load the package and data
devtools::install_github("rrohwer/limony", build_vignettes = TRUE, force = TRUE)
library(limony)
data("limony")
#browseVignettes(package = "limony") # browse the vignette

# Remove chloroplasts and renormalize the data
flat.list <- flatten.to.single.level(my.list = limony)
index <- get.taxon.indexes(my.list = flat.list, taxa = "Chloroplast")
flat.no.chloro <- subset.by.taxa(my.list = flat.list, keep.index = -index, renormalize = TRUE)

# Remove non-standard samples that were pre-filtered
index <- get.sample.indexes(my.list = flat.no.chloro, remove.prefiltered = TRUE)
no.pf <- subset.by.sample(my.list = flat.no.chloro, keep.index = index)
no.chloro <- group.by.tax(my.list = no.pf)

# If you want to remove other non-standard, can use these options
index.only.standard.depth <- get.sample.indexes(my.list = no.chloro, remove.diff.depths = TRUE) # index to remove them
no.dd <- subset.by.sample(my.list = no.chloro, keep.index = index.only.standard.depth)
no.dd <- group.by.tax(my.list = no.dd)
#index.nonstandard.depth <- get.sample.indexes(my.list = no.chloro, only.diff.depths = TRUE) # index to keep only them, like you could just check that they are not outliers
# The depths are in the sample names, D0 is surface, D8 is 1-8m integrated, D15 is 1-15m, and D alone means the two replicates were different depths.
index.only.standard.locations <- get.sample.indexes(my.list = no.dd, remove.diff.loc = TRUE) # index to remove them
no.dl <- subset.by.sample(my.list = no.dd, keep.index = index.only.standard.locations)
no.dl <- group.by.tax(my.list = no.dl)
#index.nonstandard.locations <- get.sample.indexes(my.list = no.chloro, only.diff.loc = TRUE) # index to keep only them, like you could just check that they are not outliers
# the locations are also in the sample names, DC I'm not sure where it was, but was in mendota, the s ones are different transects in winter, not quite at the buoy location

# The other non-standard you can check is the DNA concentration, to remove low-yield extractions
# update the package and load a new data that I added to it
# data("key")
# ?key
# index.no.low.dna <- get.sample.indexes(my.list = no.chloro, my.key = key, low.yield.ng.uL = 3, remove.low.dna = TRUE) 
# index.low.dna <- get.sample.indexes(my.list = no.chloro, my.key = key, low.yield.ng.uL = 3, only.low.dna  = TRUE) 

# Get OTU table
otu.table <- no.dl$av$seqID
otu.taxonomy <- no.dl$names$seqID
sample.dates <- convert.sample.names.to.dates(sample.names = colnames(otu.table))

# Combine OTU table with OTU taxonomy table
otu.table <- as.data.frame(otu.table)
otu.taxomony <- as.data.frame(otu.taxonomy)
otu.table <- tibble::rownames_to_column(otu.table, "seqID")
data <- merge(otu.taxonomy, otu.table, by = "seqID")

# Make data long format
data_long <- gather(data, date, abund, "23Apr2001":"24Apr2019")

#--------------------------------------
# Option 1: summing abundances by phyla
# Sum abundances for each phyla on each date
data_summed <- data_long %>%
  group_by(date, phylum) %>%
  summarise(abund =  sum(abund, na.rm = TRUE)) %>%
  ungroup()
# This code is fine, but you can also get this directly from the limony list
# phyla.table <- no.dl$av$Phylum
# phyla.tax <- no.dl$names$Phylum
#--------------------------------------

#--------------------------------------
# Option 2: summing abundances by OTU
# Sum abundances for each OTU on each date
data_summed <- data_long %>%
  group_by(date, seqID) %>%
  summarise(abund =  sum(abund, na.rm = TRUE)) %>%
  ungroup()
#--------------------------------------

# Remove dates that have an extension after the date
#data_summed <- data_summed[!grepl("\\.", data_summed$date), ]
# Notes from Robin about this: Yes, this is true, but myself I would keep most of them and maybe just check that they are not outliers in your analysis later. I would remove prefiltered from cyano analysis though. And also, you will have already removed them above, with the package functions.

# Add day, month, year, and doy columns
data_summed <- data_summed %>%
  mutate(sampledate = dmy(date),
         day = day(sampledate),
         month = month(sampledate),
         year = year(sampledate),
         doy = yday(sampledate))
data_summed$year <- as.character(data_summed$year) # for plotting

# Create .csv file of phyla abundances across all samples, for further analysis in microbial_analysis.R and random_forest.R scripts.
# write.csv(data_summed, file = "Data/Microbial/phyla_relabund.csv")

# Create .csv file of OTU abundances across all samples, for further analysis in microbial_analysis.R and random_forest.R scripts.
# write.csv(data_summed, file = "Data/Microbial/OTU_relabund.csv")




# -------------------------------------------
# Below is some code that Robin sent me for playing around with season definitions. Don't worry about this for now, but revisit if we end up using the season boundaries instead of DOY.

# Just remember, when you subset samples from limony list, also subset the same rows from the key.
index.may <- get.sample.indexes(my.list = no.dl, start.YY.MM.DD = "00-5-1", end.YY.MM.DD = "19-5-31", dates.are.season.range = TRUE)
limony.may <- subset.by.sample(my.list = no.dl, keep.index = index.may)
key.may <- key[index.may, ]
all.equal(key.may$in.R.colnames, colnames(limony.may$av$Kingdom)) 

# Also you can subset by the season definitions in Robin's paper
index.clearwater <- get.sample.indexes(my.list = no.dl, my.key = key, selected.seasons = "clearwater")
limony.clear <- subset.by.sample(my.list = no.dl, keep.index = index.clearwater)
key.clear <- key[index.clearwater, ]

# And then you can say, plot aligned by end of clearwater instead of day of year. aligning by these seasons made a lot of patterns more clear
flat.list <- flatten.to.single.level(my.list = no.dl)
index <- get.taxon.indexes(my.list = flat.list, taxa = "Aphanizomenon")
aphan <- subset.by.taxa(my.list = flat.list, keep.index = index)
aphan <- group.by.tax(my.list = aphan)
aphan <- cbind(key, "Aphanizomenon" = aphan$av$`Genus/Clade`[1, ,drop=T])

ggplot(data = aphan, aes(x = yDay, y = Aphanizomenon))+
  geom_line(aes(group = Year, color = Invasion)) # sorta looks like aphanizom is appearing earlier post-invasion
ggplot(data = aphan, aes(x = Days.Since.Clearwater, y = Aphanizomenon))+
  geom_line(aes(group = Year, color = Invasion)) # more obvious when aligned by season

# *** also this object now, season starts from my paper ***
data("seasons")
seasons
ave.ice.off <- mean(yday(seasons$Spring) - 1) # "spring start" is the day after ice-off
ggplot(data = aphan, aes(x = yDay, y = Aphanizomenon))+
  geom_line(aes(group = Year, color = Invasion))+
  geom_vline(xintercept = ave.ice.off)