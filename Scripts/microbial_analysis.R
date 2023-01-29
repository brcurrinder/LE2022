# Analysis of microbial data in limony dataset
# November 2022, cleaned up January 2023
# AGS

pacman::p_load(tidyverse, devtools, lubridate, limony)

#devtools::install_github("rrohwer/limony", build_vignettes = TRUE)
#browseVignettes(package = "limony")

# Download data
data("limony")

# Remove Chloroplasts
lim.flat <- flatten.to.single.level(my.list = limony)
index <- get.taxon.indexes(my.list = lim.flat, taxa = "Chloroplast")
lim.flat <- subset.by.taxa(my.list = lim.flat, keep.index = -index, renormalize = TRUE) # note that the negative removes that index
limony <- group.by.tax(my.list = lim.flat)

#index <- get.sample.indexes(my.list = limony, remove.prefiltered = TRUE)
#no.pf <- subset.by.sample(my.list = limony, keep.index = index)


# Define date range of samples to keep
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "15-8-1", end.YY.MM.DD = "15-8-02")
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "15-9-1", end.YY.MM.DD = "15-9-02")
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "15-7-1", end.YY.MM.DD = "15-8-02")
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "15-1-1", end.YY.MM.DD = "15-12-31") # 2015 only
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "16-1-1", end.YY.MM.DD = "16-12-31") # 2016 only
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "15-1-1", end.YY.MM.DD = "18-12-31") # 2015 through 2018 only

index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "09-1-1", end.YY.MM.DD = "22-12-31") # 2009 through 2022

# 8/22/16 - exact match
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "16-8-22", end.YY.MM.DD = "16-8-22")

# 8/31/16 - day before (also have day after, but 8-30 was the same one for metabolism data)
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "16-8-30", end.YY.MM.DD = "16-8-30")

# all dates
index <- get.sample.indexes(my.list = limony, start.YY.MM.DD = "00-01-01", end.YY.MM.DD = "22-12-31")

subset <- subset.by.sample(my.list = limony, keep.index = index)
limony <- subset

# Get OTU table
otu.table <- limony$av$seqID
otu.taxonomy <- limony$names$seqID
sample.dates <- convert.sample.names.to.dates(sample.names = colnames(otu.table))
#otu.table[1:3,1:3] # numeric matrix
#otu.taxonomy[1:3,1:3] # character matrix
#sample.dates[1:3] # dates vector

# Questions:
#   Where is the 18S data? Is that included in here, or in a totally separate location? - found some Eukaryota
#   Can we get any sort of abundance data out of this? Or is it just relative abundance? Could look at proportion of community that is cyanobacteria. Or look at differential abundances.

# Playing around with making some figures
otu.table <- as.data.frame(otu.table)
otu.taxomony <- as.data.frame(otu.taxonomy)
otu.table <- tibble::rownames_to_column(otu.table, "seqID")
data <- merge(otu.taxonomy, otu.table, by="seqID")

# Plotting one date - jump below for looking across multiple dates
names(data)[names(data) == "22Aug2016"] <- "date"
names(data)[names(data) == "30Aug2016"] <- "date"

data_summed <- data %>% group_by(phylum) %>%
  summarise(abund =  sum(date, na.rm = TRUE)) %>%
  ungroup()

subset <- subset(data_summed, abund>1.4) 

subset %>% # bar plot
  ggplot(aes(x = phylum, y = abund, fill = phylum)) +
  geom_bar(stat = "identity", position="dodge") +
  labs(x = "Phylum", y = "Relative abundance (%)",
       title = "Top 7 most abundant phyla on 30Aug2016") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(fill = "Phylum") +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        title=element_text(size=14),
        legend.title=element_text(size=14,face="bold"))

#### Looking across multiple dates
data_long <- gather(data, date, abund, "06Apr2016":"16Nov2016")

data_summed <- data_long %>% group_by(date, phylum) %>%
  summarise(abund =  sum(abund, na.rm = TRUE)) %>%
  ungroup()

#subset <- subset(data_summed, abund>2) 

subset %>% # bar plot
  ggplot(aes(x = phylum, y = abund, fill = phylum)) +
  facet_wrap(vars(date), nrow=4) +
  geom_bar(stat = "identity", position="dodge") +
  labs(x = "Phylum", y = "Relative abundance (%)",
       title = "Phyla relative abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(fill = "Phylum") +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        title=element_text(size=14),
        legend.title=element_text(size=14,face="bold"))

# Looking at all 2015 dates at once
data_long <- gather(data, date, abund, "15Feb2015":"13Aug2015")
data_long <- gather(data, date, abund, "06Apr2016":"16Nov2016") # all 2016 dates
data_long <- gather(data, date, abund, "15Feb2015":"08Nov2018") # 2015-2018 dates

data_summed <- data_long %>% group_by(date, phylum) %>%
  summarise(abund =  sum(abund, na.rm = TRUE)) %>%
  ungroup()

subset <- data_summed

#subset <- subset(data_summed, abund>2) 

subset %>% # bar plot
  ggplot(aes(x = phylum, y = abund, fill = phylum)) +
  facet_wrap(vars(date), nrow=4) +
  geom_bar(stat = "identity", position="dodge") +
  labs(x = "Phylum", y = "Relative abundance (%)",
       title = "Phyla relative abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(fill = "Phylum") +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        title=element_text(size=14),
        legend.title=element_text(size=14,face="bold"))

# Getting a value for % Cyanobacteria and % heterotrophs for each sample

# collapse everything that's not Cyanobacteria into "Heterotrophs" and then summarize
summary <- subset %>% 
  mutate(group = case_when(phylum != "Cyanobacteria" ~ "Heterotroph",
                           phylum == "Cyanobacteria" ~ "Cyanobacteria")) %>% 
  group_by(date, group) %>%
  summarise(relabund =  sum(abund, na.rm = TRUE)) %>%
  ungroup()

summary %>% # bar plot
  ggplot(aes(x = date, y = relabund, fill = group)) +
  geom_bar(stat = "identity", position="stack") +
  labs(x = "Date", y = "Relative abundance (%)",
       title = "") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(fill = "Group") +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        title=element_text(size=14),
        legend.title=element_text(size=14,face="bold"))

# Next steps:
#   figure out how to sort the x-axis by time
#   also, shouldn't all of the bars add up to 100%??

summary <- summary %>% mutate(datenew = dmy(date),
                              day = day(datenew),
                              month = month(datenew),
                              year = year(datenew))


summary %>% # bar plot
  ggplot(aes(x = month, y = relabund, fill = group)) +
  geom_bar(stat = "identity", position="stack") +
  #scale_x_date(limits = as.Date(c('2015-01-01','2018-12-31'))) +
  facet_wrap(vars(year), ncol=1, scales="free") +
  labs(x = "Date", y = "Relative abundance (%)",
       title = "") +
  #theme(axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank()) +
  labs(fill = "Group") +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        title=element_text(size=14),
        legend.title=element_text(size=14,face="bold"))










# Plotting just Cyanobacterial relative abundance 
cyano <- subset(data_summed, phylum == "Cyanobacteria") 

month = c("Aug", "July", "June", "May", "Aug", "July", "June", "Aug", "June",
          "Aug", "July", "March", "June","Aug", "May", "July", "June", "April",
          "Aug", "June", "Aug",
          "Aug", "July", "June", "May", "Aug", "July", "June", "Aug", "June",
          "Aug", "July", "March", "June","Aug", "May", "July", "June", "April",
          "Aug", "June", "Aug")

cyano$month <- month

cyano %>% # bar plot
  ggplot(aes(x = date, y = abund, fill = month)) +
  geom_bar(stat = "identity", position="dodge")
