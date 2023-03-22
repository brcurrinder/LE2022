The PRISMA output CSVs, clipped to the whole lake, are too big to push to github.

Follow these instructions to get them in your project:
1. Make two new folders in Data/PRISMA/. Name these folders "lake" and "buoy". (These names are included in the .gitignore file, so you won't accidentally try to push all these files next time you commit).
2. Download the PRISMA CSVs from google drive: https://drive.google.com/drive/folders/1qEnqdriNecd9iz3ua1zZaKPVKU1D44Vg?usp=share_link
(if this link doesn't work you can find the files at "Earth Observations (satellite and airborne)/PRISMA/L2R")
3. Move those files into Data/PRISMA/lake and Data/PRISMA/buoy
4. Now you should be able to use these files in any scripts and analyses within this project, without having to push or pull them.

**NOTE**: if these files are ever updated you will need to re-download the new version and replace the copy on your computer. This version is current as of March 21, 2023.
