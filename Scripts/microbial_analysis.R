# Analysis of 16S microbial data in limony dataset
# November 2022, cleaned up January 2023
# AGS

pacman::p_load(tidyverse, devtools, lubridate, limony)

# Load the package and data
devtools::install_github("rrohwer/limony", build_vignettes = TRUE)
library(limony)
data("limony")
#browseVignettes(package = "limony") # browse the vignette

# Remove chloroplasts and renormalize the data
flat.list <- flatten.to.single.level(my.list = limony)
index <- get.taxon.indexes(my.list = flat.list, taxa = "Chloroplast")
flat.no.chloro <- subset.by.taxa(my.list = flat.list, keep.index = -index)

# Remove non-standard samples that were pre-filtered
index <- get.sample.indexes(my.list = flat.no.chloro, remove.prefiltered = TRUE)
no.pf <- subset.by.sample(my.list = flat.no.chloro, keep.index = index)
no.chloro <- group.by.tax(my.list = no.pf)

# Get OTU table
otu.table <- no.chloro$av$seqID
otu.taxonomy <- no.chloro$names$seqID
sample.dates <- convert.sample.names.to.dates(sample.names = colnames(otu.table))

# Combine OTU table with OTU taxonomy table
otu.table <- as.data.frame(otu.table)
otu.taxomony <- as.data.frame(otu.taxonomy)
otu.table <- tibble::rownames_to_column(otu.table, "seqID")
data <- merge(otu.taxonomy, otu.table, by = "seqID")

# Make data long format
data_long <- gather(data, date, abund, "23Apr2001":"24Apr2019")

# Sum abundances for each phyla on each date
data_summed <- data_long %>% 
  group_by(date, phylum) %>%
  summarise(abund =  sum(abund, na.rm = TRUE)) %>%
  ungroup()

# Remove dates that have an extension after the date -- check in with Robin about this -- I am pretty sure these are non-standard samples that should be removed
data_summed <- data_summed[!grepl("\\.", data_summed$date), ]

# Add day, month, year, and doy columns
data_summed <- data_summed %>%
  mutate(sampledate = dmy(date),
       day = day(sampledate),
       month = month(sampledate),
       year = year(sampledate),
       doy = yday(sampledate))
data_summed$year <- as.character(data_summed$year) # for plotting

#--------------------------------------------------------------
# Plot % cyanobacteria over time, make csv for time series fig
#-------------------------------------------------------------

# Subset to just look at Cyanobacteria
cyano <- subset(data_summed, phylum == "Cyanobacteria") 

# Plot with a different line for each year
cyano %>%
  ggplot(aes(x = doy, y = abund, col = year)) + 
  geom_line(stat = "identity") +
  theme_classic() +
  labs(x = "Day of year", y = "Cyanobacterial relative abundance (%)", title = "")

# Plot across all years
cyano %>%
  ggplot(aes(x = sampledate, y = abund)) + 
  geom_line(stat = "identity") +
  theme_classic() +
  labs(x = "Date", y = "Cyanobacterial relative abundance (%)", title = "")

# Write .csv for time series figure
rel_abund_cyano <- cyano %>%
  #select(-c(date, phylum, sampledate, day, month)) %>%
  select(-c(phylum)) %>%
  rename("Cyanobacteria_rel_abund_16S" = "abund")
write.csv(rel_abund_cyano, "Data/prelim_time_series_figure/Cyanobacteria_rel_abund_16S.csv", row.names = F, quote=T)

#----------------------------------------------------------------------
# For fun: looking at relative abundances of top phyla across all dates
#----------------------------------------------------------------------

# Subset to only keep phyla that have a % relative abundance >2 at some point
data_summed_subset <- subset(data_summed, abund > 2)

# Subset by year (if desired)
data_summed_subset <- subset(data_summed_subset, year > 2015)

# Stacked bar plot
data_summed_subset %>%
  ggplot(aes(x = doy, y = abund, fill = phylum)) +
  geom_bar(stat = "identity", position="stack") +
  labs(x = "Day of year", y = "Relative abundance (%)", title = "") +
  facet_wrap(vars(year), ncol = 1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        title=element_text(size=14),
        legend.title=element_text(size=14,face="bold"))


#----------------------------------------------------------------------
# For random forest analysis
#----------------------------------------------------------------------

library(randomForest)
library(pdp)

data_summed2 <- data_summed %>%
  select(-c(date, day, month, year, doy)) 

summarized <- data_summed2 %>%
  pivot_wider(names_from = phylum, values_from = abund)

metabolism <- read.csv("Data/Metabolism_Model_Output/SimResultsMatrix_MetabData_run12jan23.csv") # load metab data

metabolism <- metabolism %>% rename("sampledate" = "SimDate")

metabolism$sampledate <- as.Date(metabolism$sampledate)

rf_input <- merge(x = summarized, y = metabolism, by.x = "sampledate")

rf_input <- rf_input %>%
  select(1:73, 78) %>%
  select(-contains("-")) # columns with dashes were causing us trouble

# Create training and test sets
rf_training <- sample_frac(tbl=rf_input,size=0.7,replace=FALSE)
rf_test <- rf_input[!(rf_input$sampledate %in% rf_training$sampledate),]

# Remove the date from the training set
rf_training <- rf_training %>% select(-c(sampledate))

# Run the random forest
rf <- randomForest(data=rf_training, EpiNPP_mgC_L~., localImp=TRUE, err.rate=TRUE)

rf
plot(rf) # what it's supposed to look like

# Plot predicted vs. observed NPP
rf_test$predict <- predict(rf, rf_test)
rf_test %>%
  ggplot(aes(x= EpiNPP_mgC_L, y = predict)) +
  geom_point() +
  geom_abline()

# Variable importance plot
varImpPlot(rf)
# node purity - not many people use
# %IncMSE - % increase in MSE if that variable is not included

# Partial dependence plot - how the prediction responds in relation to one of the variables
pd <- partial(object = rf,
              pred.var = 'Spirochaetota',
              plot = TRUE)
# the yhat (metabolism) goes up steeply with any presence of Spiro
pd <- partial(object = rf,
              pred.var = 'Cyanobacteria',
              plot = TRUE)
pd <- partial(object = rf,
              pred.var = 'Bacteroidota',
              plot = TRUE)
pd <- partial(object = rf,
              pred.var = 'Planctomycetota',
              plot = TRUE)

# Run actual regression and see what the r-squared is
model <- lm(data = rf_test, formula = predict ~ EpiNPP_mgC_L)

summary(model)
# pretty good!



#--------------------------------------------------------------
# Plot trends in four desired phyla over time
#-------------------------------------------------------------

# Subset
fourphyla <- subset(data_summed, phylum %in% c("Cyanobacteria", "Spirochaetota", "Planctomycetota", "Bacteroidota"))

fourphyla$year <- as.character(fourphyla$year)
fourphyla$phylum <- factor(fourphyla$phylum,
                           levels = c("Spirochaetota",
                                      "Cyanobacteria",
                                      "Planctomycetota",
                                      "Bacteroidota"))

# Plot with a different line for each year, faceted by phylum
fourphyla %>%
  ggplot(aes(x = doy, y = abund, col = year)) + 
  geom_line(stat = "identity") +
  facet_wrap(~phylum, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(x = "Day of year", y = "Phylum relative abundance (%)", title = "")
