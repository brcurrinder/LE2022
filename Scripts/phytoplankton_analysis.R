# Analysis of Mendota LTER cyanobacterial and algal count data
# January 2023
# AGS

pacman::p_load(tidyverse, lubridate)

# Download the data
# Package ID: knb-lter-ntl.88.31 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Phytoplankton - Madison Lakes Area 1995 - current.
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/88/31/f2de15b2fff6ae962a04c150c0a1c510"
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")
phyto <- read_csv(infile1)

# Filter for Lake Mendota and add day, month, year, and doy columns
phyto <- phyto %>% 
  filter(lakeid == 'ME') %>%
  mutate(sampledate = ymd(sampledate),
         day = day(sampledate),
         month = month(sampledate),
         year = year(sampledate),
         doy = yday(sampledate))
phyto$year <- as.character(phyto$year)

# "phyto" is our dataset which will be used, and never re-written, in the next three sections

table(phyto$division) # there are 10 divisions (phyla)
table(phyto$genus) # genera that are present in the dataset

# Plot biovolume vs. biomass to see if they correlate - they do! (except for one point -- look into)
phyto %>% 
  ggplot(aes(x = biovolume_conc, y = biomass_conc)) +
  geom_point()

#---------------------------------------------------
# Plot trends in Cyanophyta and all groups over time
#---------------------------------------------------

# Plot Cyanophyta for each date
cyanophyta <- phyto %>%
  subset(division == "Cyanophyta") %>%
  group_by(sampledate) %>%
  summarise(summed = sum(biomass_conc),
            day = day[1], month = month[1], year = year[1]) %>%
  ungroup()

# Plot all groups for each date
all_groups <- phyto %>% 
  group_by(sampledate, division) %>%
  summarise(summed_biomass_conc = sum(biomass_conc),
            day = day[1], month = month[1], year = year[1]) %>%
  ungroup()

# Choose 1 to plot
since2014 <- subset(cyanophyta, year > 2013)
since2014 <- subset(all_groups, year > 2013)

# Plot (can clean this up and get rid of the fake date thing)
since2014 %>%
  mutate(FakeYeardate = ymd(paste("2024", month(sampledate), day(sampledate), sep = "-"))) %>% # make a Fake date we'll use to plot by
  ggplot(aes(x = FakeYeardate, y = summed)) + # add fill = division within the aes() for plotting all groups
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(~year, ncol = 1) +
  theme_classic() +
  labs(x = "Date", y = "Summed Cyanophyta biomass_conc", title = "") # Change y axis title for plotting all groups

#---------------------------------------------------------------------------
# Plot total abundance of all groups over time, make csv for time series fig
#---------------------------------------------------------------------------

# Sum biomass_conc of all groups for each date
all_groups <- phyto %>% 
  group_by(sampledate) %>%
  summarise(summed_biomass_conc = sum(biomass_conc),
            sampledate = sampledate[1], year = year[1], doy = doy[1]) %>%
  ungroup() %>%
  filter(summed_biomass_conc < 50) # remove outlier date - look into this and make sure I am okay to remove this

# Plot with a different line for each year
all_groups %>%
  ggplot(aes(x = doy, y = summed_biomass_conc, col = year)) + 
  geom_line(stat = "identity") +
  theme_classic() +
  labs(x = "Day of year", y = "Total algal biomass_conc (mg/L)", title = "")

# Plot across all years
all_groups %>%
  ggplot(aes(x = sampledate, y = summed_biomass_conc)) + 
  geom_line(stat = "identity") +
  theme_classic() +
  labs(x = "Date", y = "Total algal biomass_conc (mg/L)", title = "")

# Write .csv for time series figure
all_groups <- all_groups %>%
  #select(-c(sampledate)) %>%
  rename("total_phytoplankton_biomass_conc_mgL" = "summed_biomass_conc")
write.csv(all_groups, "Data/prelim_time_series_figure/total_phytoplankton_biomass_conc.csv", row.names = F, quote=T)

#------------------------------------------------------------------------------
# Plot relative abundance of Cyanophyta over time, make csv for time series fig
#------------------------------------------------------------------------------

# Add column for total_biomass_conc for each sampledate
phyto2 <- phyto %>% 
  group_by(sampledate) %>%
  mutate(total_biomass_conc =  sum(biomass_conc, na.rm = TRUE)) %>%
  ungroup()

# Calculate relative_biomass_conc for each division for each sampledate
phyto_rel <- phyto2 %>%
  group_by(sampledate, division) %>%
  summarise(biomass_conc = sum(biomass_conc, na.rm = TRUE),
            total_biomass_conc = total_biomass_conc[1],
            year = year[1], doy = doy[1]) %>%
  mutate(rel_biomass_conc = biomass_conc/total_biomass_conc) %>%
  ungroup()

# Just interested in Cyanophyta
phyto_rel_cyanophyta <- filter(phyto_rel, division == "Cyanophyta")

# Plot with a different line for each year
phyto_rel_cyanophyta %>%
  ggplot(aes(x = doy, y = rel_biomass_conc, col = year)) + 
  geom_line(stat = "identity") +
  theme_classic() +
  labs(x = "Day of year", y = "Cyanophyta relative biomass_conc (mg/L)", title = "")

# Plot across all years
phyto_rel_cyanophyta %>%
  ggplot(aes(x = sampledate, y = rel_biomass_conc)) + 
  geom_line(stat = "identity") +
  theme_classic() +
  labs(x = "Date", y = "Cyanophyta relative biomass_conc (mg/L)", title = "")

# Write .csv for time series figure
phyto_rel_cyanophyta <- phyto_rel_cyanophyta %>%
  #select(-c(sampledate, division, biomass_conc, total_biomass_conc)) %>%
  select(-c(division, biomass_conc, total_biomass_conc)) %>%
  rename("Cyanophyta_rel_biomass_conc_mgL" = "rel_biomass_conc")
write.csv(phyto_rel_cyanophyta, "Data/prelim_time_series_figure/Cyanophyta_relative_biomass_conc.csv", row.names = F, quote=T)

# This is all of the cyanobacteria (Cyanophyta phylum) grouped together - if we are interested, we could look at specific genera within Cynanophyta such as Microcystis. This dataset has that genus-level information.