# Random forest analysis of 16S microbial data predicting metabolism output
# Running the RF analysis on top OTUs
# April 2023
# AGS

pacman::p_load(tidyverse, devtools, lubridate, limony, randomForest, pdp, patchwork)

# Load the data file generated in the microbial_cleaning.R script
OTU_relabund <- read.csv(file = "Data/OTU_relabund.csv")
OTU_relabund <- OTU_relabund[,-1]
OTU_relabund$sampledate <- as.Date(OTU_relabund$sampledate)

#----------------------------------------------------------------------
# Creating subsets of the data that we want to look at
#----------------------------------------------------------------------

df <- data_long

max <- df %>%
  group_by(seqID) %>%
  summarise(max = max(abund)) %>%
  ungroup() %>% arrange(desc(max))
max_subset <- max %>% slice_max(max, n = 50)

sum <- df %>%
  group_by(seqID) %>%
  summarise(sum = sum(abund)) %>%
  ungroup() %>% arrange(desc(sum))
sum_subset <- sum %>% slice_max(sum, n = 50)

mean <- df %>%
  group_by(seqID) %>%
  summarise(mean = mean(abund)) %>%
  ungroup() %>% arrange(desc(mean))
mean_subset <- mean %>% slice_max(mean, n = 50)

pres <- df %>%
  group_by(seqID) %>%
  summarise(pres = sum(abund > 0.0)) %>%
  ungroup() %>% arrange(desc(pres))
pres_subset <- pres %>% slice_max(pres, n = 50)


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

metabolism <- read.csv("Data/SimResultsMatrix_MetabData_run12jan23.csv") # load metabolism data

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
rf_training <- sample_frac(tbl=rf_input,size=0.8,replace=FALSE)
rf_test <- rf_input[!(rf_input$sampledate %in% rf_training$sampledate),]

# Remove the date from the training set
rf_training <- rf_training %>% select(-c(sampledate))

# search across many hyperparameter combinations and evaluate using ranger package (faster than randomForest package)
library(ranger)
# hyperparameter grid search 
hyper_grid <- expand.grid(
  mtry       = seq(5, 50, by = 2), # 
  node_size  = seq(3, 11, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)
# total number of combinations
nrow(hyper_grid)

# loop across hyperparameter combinations 
for(i in 1:nrow(hyper_grid)) {
  
  # train model
  model <- ranger(
    formula         = EpiNPP_mgC_L ~ ., 
    data            = rf_training, 
    num.trees       = 500,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sample_size[i],
    seed            = 123
  )
  
  # add OOB error to grid
  hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
}

hyper_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10) # top row signifies combination of mtry, node_size, and sample_size that minimizes RMSE

# Run the random forest
rf <- randomForest(data=rf_training, EpiNPP_mgC_L~., localImp=TRUE, err.rate=TRUE, mtry = 15, nodesize = 3) # choose optimal mtry and nodesize
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

# extract all OTUs and associated % increased MSE values
NPP_sum_OTUs <- randomForest::importance(rf, type = 1) # change name of object based on "mean", "sum", "presence", or "max"

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
rf_training <- sample_frac(tbl=rf_input, size=0.8, replace=FALSE)
rf_test <- rf_input[!(rf_input$sampledate %in% rf_training$sampledate),]

# Remove the date from the training set
rf_training <- rf_training %>% select(-c(sampledate))

# search across many hyperparameter combinations and evaluate using ranger package (faster than randomForest package)
library(ranger)
# hyperparameter grid search 
hyper_grid <- expand.grid(
  mtry       = seq(5, 50, by = 2), # 
  node_size  = seq(3, 11, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)
# total number of combinations
nrow(hyper_grid)

# loop across hyperparameter combinations 
for(i in 1:nrow(hyper_grid)) {
  
  # train model
  model <- ranger(
    formula         = EpiDO_mgO2_L ~ ., 
    data            = rf_training, 
    num.trees       = 500,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sample_size[i],
    seed            = 123
  )
  
  # add OOB error to grid
  hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
}

hyper_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10) # top row signifies combination of mtry, node_size, and sample_size that minimizes RMSE

# Run the random forest
rf <- randomForest(data=rf_training, EpiDO_mgO2_L~., localImp=TRUE, err.rate=TRUE, mtry = 9, nodesize = 3)
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

# extract all OTUs and associated % increased MSE values
DO_sum_OTUs <- randomForest::importance(rf, type = 1) # change name of object based on "mean", "sum", "presence", or "max"

# Run actual regression and see what the r-squared is
model <- lm(data = rf_test, formula = predict ~ EpiDO_mgO2_L)
summary(model)

#--------------------------------------------------------------
# Further analysis of important OTUs
#-------------------------------------------------------------

colnames(NPP_max_OTUs) <- c("OTU", "IncMSE")
colnames(NPP_mean_OTUs) <- c("OTU", "IncMSE")
colnames(NPP_presence_OTUs) <- c("OTU", "IncMSE")
colnames(NPP_sum_OTUs) <- c("OTU", "IncMSE")
colnames(DO_max_OTUs) <- c("OTU", "IncMSE")
colnames(DO_mean_OTUs) <- c("OTU", "IncMSE")
colnames(DO_presence_OTUs) <- c("OTU", "IncMSE")
colnames(DO_sum_OTUs) <- c("OTU", "IncMSE")

NPP_max_OTUs <- data.frame(NPP_max_OTUs)
NPP_mean_OTUs <- data.frame(NPP_mean_OTUs)
NPP_presence_OTUs <- data.frame(NPP_presence_OTUs)
NPP_sum_OTUs <- data.frame(NPP_sum_OTUs)
DO_max_OTUs <- data.frame(DO_max_OTUs)
DO_mean_OTUs <- data.frame(DO_mean_OTUs)
DO_presence_OTUs <- data.frame(DO_presence_OTUs)
DO_sum_OTUs <- data.frame(DO_sum_OTUs)

library(tibble)
NPP_max_OTUs <- tibble::rownames_to_column(NPP_max_OTUs, "OTU")
NPP_mean_OTUs <- tibble::rownames_to_column(NPP_mean_OTUs, "OTU")
NPP_presence_OTUs <- tibble::rownames_to_column(NPP_presence_OTUs, "OTU")
NPP_sum_OTUs <- tibble::rownames_to_column(NPP_sum_OTUs, "OTU")
DO_max_OTUs <- tibble::rownames_to_column(DO_max_OTUs, "OTU")
DO_mean_OTUs <- tibble::rownames_to_column(DO_mean_OTUs, "OTU")
DO_presence_OTUs <- tibble::rownames_to_column(DO_presence_OTUs, "OTU")
DO_sum_OTUs <- tibble::rownames_to_column(DO_sum_OTUs, "OTU")

write.csv(NPP_max_OTUs, file = "Data/NPP_max_OTUs.csv")
write.csv(NPP_mean_OTUs, file = "Data/NPP_mean_OTUs.csv")
write.csv(NPP_presence_OTUs, file = "Data/NPP_presence_OTUs.csv")
write.csv(NPP_sum_OTUs, file = "Data/NPP_sum_OTUs.csv")
write.csv(DO_max_OTUs, file = "Data/NDO_max_OTUs.csv")
write.csv(DO_mean_OTUs, file = "Data/DO_mean_OTUs.csv")
write.csv(DO_presence_OTUs, file = "Data/DO_presence_OTUs.csv")
write.csv(DO_sum_OTUs, file = "Data/DO_sum_OTUs.csv")

df_list <- list(NPP_max_OTUs, NPP_mean_OTUs, NPP_presence_OTUs, NPP_sum_OTUs, 
                  DO_max_OTUs, DO_mean_OTUs, DO_presence_OTUs, DO_sum_OTUs)

max_values <- do.call(rbind, lapply(df_list, function(df) {
  grouped_df <- df %>% group_by(OTU)
  max_values_df <- grouped_df %>% summarize(max_IncMSE = max(IncMSE))
  return(max_values_df)
}))

top_OTUs <- max_values %>%
  group_by(OTU) %>%
  filter(max_IncMSE == max(max_IncMSE))

top_OTUs <- top_OTUs %>% arrange(desc(max_IncMSE))

write.csv(top_OTUs, file = "Data/Top_OTUs.csv")

top_OTUs <- read.csv(file = "Data/Top_OTUS.csv")

# Subset to only keep top 30 important OTUs 
OTU_relabund <- read.csv(file = "Data/OTU_relabund.csv")

library(limony)
data("limony")
flat.list <- flatten.to.single.level(my.list = limony)
names <- data.frame(flat.list$names)

important_OTUs <- names %>% 
  filter(top_OTUs$OTU %in% names$seqID)

important_OTUs <- subset(names, names$seqID %in% top_OTUs$OTU)

write.csv(important_OTUs, "Data/important_OTUs.csv")

# important_otus <- data.frame(seqID = c("otu_21", "otu_8", "otu_47", "otu_7", "otu_42", "otu_12", "otu_9", "otu_182", "otu_229", "otu_27", "otu_46", "otu_53", "otu_39", "otu_181", "otu_123", "otu_56", "otu_118", "otu_63", "otu_139"))

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
