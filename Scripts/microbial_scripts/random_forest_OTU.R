# Random forest analysis of 16S microbial data predicting metabolism output
# Running the RF analysis on top OTUs
# April 2023
# AGS

pacman::p_load(tidyverse, devtools, lubridate, limony, randomForest, pdp, patchwork)

# Load the data file generated in the microbial_cleaning.R script
OTU_relabund <- read.csv(file = "Data/Microbial/OTU_relabund.csv")
OTU_relabund <- OTU_relabund[,-1]
OTU_relabund$sampledate <- as.Date(OTU_relabund$sampledate)

#----------------------------------------------------------------------
# Creating subsets of the data that we want to look at
#----------------------------------------------------------------------

# Subset to only keep top 50 OTUs in terms of max abundance
subset_max <- OTU_relabund %>%
  filter(seqID %in% max_subset$seqID)

# Subset to only keep top 50 OTUs in terms of sum abundance
subset_sum <- OTU_relabund %>%
  filter(seqID %in% sum_subset$seqID)

# Subset to only keep top 50 OTUs in terms of mean abundance
subset_mean <- OTU_relabund %>%
  filter(seqID %in% mean_subset$seqID)

# Subset to only keep top 50 OTUs in terms of presence
subset_pres <- OTU_relabund %>%
  filter(seqID %in% pres_subset$seqID)

#----------------------------------------------------------------------
# Random forest analysis
#----------------------------------------------------------------------

# Choose which version we're running the analysis on
subset <- subset_max
subset <- subset_sum
subset <- subset_mean
subset <- subset_pres

subset <- subset %>%
  select(-c(date, day, month, year, doy)) %>%
  pivot_wider(names_from = seqID, values_from = abund)

metabolism <- read.csv("Data/Metabolism_Model_Output/SimResultsMatrix_MetabData_run12jan23.csv") # load metabolism data

# Remove NPP outliers - cut off anything with NPP > 0.5
metabolism <- subset(metabolism, EpiNPP_mgC_L < 0.5)

metabolism <- metabolism %>% rename("sampledate" = "SimDate")
metabolism$sampledate <- as.Date(metabolism$sampledate)

#----------------------------------------------------------------------
# Predicting Epi NPP
#----------------------------------------------------------------------

rf_input <- merge(x = subset, y = metabolism, by.x = "sampledate")

rf_input <- rf_input %>%
  select(1:51, 56) %>%
  select(-contains("-")) # columns with dashes were causing us trouble

# Create training and test sets
rf_training <- sample_frac(tbl=rf_input,size=0.7,replace=FALSE)
rf_test <- rf_input[!(rf_input$sampledate %in% rf_training$sampledate),]

# Remove the date from the training set
rf_training <- rf_training %>% select(-c(sampledate))

# Run the random forest
rf <- randomForest(data=rf_training, EpiNPP_mgC_L~., localImp=TRUE, err.rate=TRUE)

rf
plot(rf)

# Plot predicted vs. observed
rf_test$predict <- predict(rf, rf_test)
rf_test %>%
  ggplot(aes(x= EpiNPP_mgC_L, y = predict)) +
  geom_point() +
  geom_abline()

# Variable importance plot
varImpPlot(rf)

# Run actual regression and see what the r-squared is
model <- lm(data = rf_test, formula = predict ~ EpiNPP_mgC_L)
summary(model)

#----------------------------------------------------------------------
# Predicting Epi DO
#----------------------------------------------------------------------

rf_input <- merge(x = subset, y = metabolism, by.x = "sampledate")

rf_input <- rf_input %>%
  select(1:51, 54) %>%
  select(-contains("-")) # columns with dashes were causing us trouble

# Create training and test sets
rf_training <- sample_frac(tbl=rf_input,size=0.7,replace=FALSE)
rf_test <- rf_input[!(rf_input$sampledate %in% rf_training$sampledate),]

# Remove the date from the training set
rf_training <- rf_training %>% select(-c(sampledate))

# Run the random forest
rf <- randomForest(data=rf_training, EpiDO_mgO2_L~., localImp=TRUE, err.rate=TRUE)

rf
plot(rf)

# Plot predicted vs. observed
rf_test$predict <- predict(rf, rf_test)
rf_test %>%
  ggplot(aes(x= EpiDO_mgO2_L, y = predict)) +
  geom_point() +
  geom_abline()

# Variable importance plot
varImpPlot(rf)

# Run actual regression and see what the r-squared is
model <- lm(data = rf_test, formula = predict ~ EpiDO_mgO2_L)
summary(model)

#--------------------------------------------------------------
# Further analysis of important OTUs
#-------------------------------------------------------------

# Subset to only keep top 30 important OTUs 
important_otus <- data.frame(seqID = c("otu_21", "otu_8", "otu_47", "otu_7", "otu_42", "otu_12", "otu_9", "otu_182", "otu_229", "otu_27", "otu_46",
                    "otu_53", "otu_39", "otu_181", "otu_123", "otu_56", "otu_118", "otu_63", "otu_139"))

important_subset <- OTU_relabund %>%
  filter(seqID %in% important_otus$seqID)

important_subset$year <- as.character(important_subset$year)

# Plot with a different line for each year, faceted by OTU
important_subset %>%
  ggplot(aes(x = doy, y = abund, col =  year)) + 
  geom_line(stat = "identity") +
  facet_wrap(~seqID, ncol = 4, scales = "free_y") +
  theme_classic() +
  labs(x = "Day of year", y = "Phylum relative abundance (%)", title = "")

# What are these OTUs??
# the data_long data frame has columns for all of the taxonomic levels, so let's subset that to only have our important OTUs
data_long_2 <- data_long %>%
  group_by(seqID) %>%
  mutate(test = sum(abund)) %>%
  select(-c(date, abund, test)) %>%
  filter(seqID %in% important_otus$seqID) %>%
  ungroup()

important_tax <- data_long_2[1:19,]

important_tax %>%
  ggplot(aes(x = phylum)) + 
  geom_histogram(stat = "count") 

ggplot(df, aes(x=weight)) + geom_histogram()