#Packages
library(tidyverse)

# Data set title: North Temperate Lakes LTER: Physical Limnology of Primary Study Lakes 1981 - current.
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/29/34/03e232a1b362900e0f059859abe8eb97" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
physical = read_csv(infile1) 
# Filter for Mendota and surface measurement (= 0)
physical = physical %>% filter(lakeid == 'ME' & depth == 0) 

# Data set title: North Temperate Lakes LTER: Chemical Limnology of Primary Study Lakes: Major Ions 1981 - current.
inUrl2  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/2/36/3f740d0b77b3caf6930a8ce9cca4306a" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
ions = read_csv(infile2) 
# Filter for Mendota and surface measurement (= 0)
ions = ions %>% filter(lakeid == 'ME' & depth == 0) 

# Data set title: North Temperate Lakes LTER: Chemical Limnology of Primary Study Lakes: Nutrients, pH and Carbon 1981 - current.
inUrl3  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/1/57/802d63a4c35050b09ef6d1e7da3efd3f" 
infile3 <- tempfile()
try(download.file(inUrl3,infile3,method="curl"))
nutrients = read_csv(infile3)
# Filter for Mendota and surface measurement (= 0)
nutrients = nutrients %>% filter(lakeid == 'ME'& depth == 0) 

# Data set title: North Temperate Lakes LTER: Secchi Disk Depth; Other Auxiliary Base Crew Sample Data 1981 - current.
inUrl4  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/31/31/5a5a5606737d760b61c43bc59460ccc9" 
infile4 <- tempfile()
try(download.file(inUrl4,infile4,method="curl"))
secchi = read_csv(infile4)
# Filter for Mendota and surface measurement (= 0)
secchi = secchi %>% filter(lakeid == 'ME') 

# Data set title: North Temperate Lakes LTER: Chlorophyll - Madison Lakes Area 1995 - current.
inUrl5  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/38/28/66796c3bc77617e7cc95c4b09d4995c5" 
infile5 <- tempfile()
try(download.file(inUrl5,infile5,method="curl"))
chlorophyll = read_csv(infile5)
# Filter for Mendota 
chlorophyll = chlorophyll %>% filter(lakeid == 'ME') 

# Data set title: North Temperate Lakes LTER: Phytoplankton - Madison Lakes Area 1995 - current.
inUrl6  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/88/31/f2de15b2fff6ae962a04c150c0a1c510" 
infile6 <- tempfile()
try(download.file(inUrl6,infile6,method="curl"))
phytoplankton = read_csv(infile6)
# Filter for Mendota 
phytoplankton = phytoplankton %>% filter(lakeid == 'ME')
