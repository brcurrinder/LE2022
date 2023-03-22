The DESIS output CSVs, clipped to the whole lake, are too big to push to github.

Follow these instructions to get them in your project:
1. Make a new folder in Data/DESIS/. Name this folder L2W_destriped. (This name is included in the .gitignore file, so you won't accidentally try to push all these files next time you commit).
2. Download the zip archive containing all the DESIS CSVs from google drive: https://drive.google.com/file/d/1A2H1fp0FmqnsfHZtG4WBwq-4WyJiXSpL/view?usp=share_link
(if this link doesn't work you can find the files at "Earth Observations (satellite and airborne)/DESIS/L2W_destriped.zip")
3. Extract those files into Data/DESIS/L2W_destriped
4. Now you should be able to use these files in any scripts and analyses within this project, without having to push or pull them.

**NOTE**: if these files are ever updated you will need to re-download the new version and replace the copy on your computer. This version is current as of March 21, 2023.
