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
  select(-c(date, phylum, sampledate, day, month)) %>%
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