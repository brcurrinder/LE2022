# Analysis of 16S microbial data in limony dataset
# November 2022, cleaned up March 2023
# AGS

pacman::p_load(tidyverse, devtools, lubridate, limony)

# Load the data file generated in the microbial_cleaning.R script
phyla_relabund <- read.csv(file = "Data/Microbial/phyla_relabund.csv")
phyla_relabund <- phyla_relabund[,-1]
phyla_relabund$sampledate <- as.Date(phyla_relabund$sampledate)

#--------------------------------------------------------------
# Plot % cyanobacteria over time, make csv for time series fig
#-------------------------------------------------------------

# Subset to just look at Cyanobacteria
cyano <- subset(phyla_relabund, phylum == "Cyanobacteria") 

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
# write.csv(rel_abund_cyano, "Data/prelim_time_series_figure/Cyanobacteria_rel_abund_16S.csv", row.names = F, quote=T)

#-------------------------------------------------------------
# Looking at relative abundances of top phyla across all dates
#-------------------------------------------------------------

# Subset to only keep phyla that have a % relative abundance >2 at some point
phyla_relabund_subset <- subset(phyla_relabund, abund > 2)

# Subset by year (if desired)
phyla_relabund_subset <- subset(phyla_relabund_subset, year > 2015)

# Stacked bar plot
phyla_relabund_subset %>%
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