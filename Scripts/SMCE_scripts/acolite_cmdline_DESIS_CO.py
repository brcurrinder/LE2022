from multiprocessing import Process, freeze_support
freeze_support()
import sys, os
import glob
import datetime
import argparse
import numpy as np
import zipfile

############## Edit this block for mission and region #################
comSatDir = '../'
zipListDir = '../DESIS/L1C_zips'
backupDir = '../aco_backup/'

mission='DESIS/'

#settingsCode = 'HSI_ChesBay_20220622' # dsf range 400:1000; percentile=5; Alt glint corr 930:960 thresh 0.025 at 945; add L2W back
settingsCode = 'HSI_default'

# Checks in MISSION/ and MISSION/LEVEL/REGION for L2W, then L2R, then L1R,
#   and SKIPS writing that level if found, unless this is set (for reprocessing)
replaceL1R =0 # Very rare
replaceL2R =0
replaceL2W =0


########################################################################

parentDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
acoDir = parentDir+'/acolite'
os.sys.path.insert(0,acoDir)

## import acolite source
try:
    import acolite as ac
except:
    print('Could not import ACOLITE source')
    print("Error:", sys.exc_info())

## ignore numpy errors
olderr = np.seterr(all='ignore')

## set up command line interface
parser = argparse.ArgumentParser(description='ACOLITE')
parser.add_argument('--settings', help='settings file', default=None)
parser.add_argument('--inputfile', help='list of images', default=None)
parser.add_argument('--output', help='output directory', default=None)
args, unknown = parser.parse_known_args()
#### Temporary for testing. Strip off test dir
miss = mission.strip('/')
args.settings = f'{miss}_{settingsCode}.txt' # Changed below for SWIR (WV-3)


if 'DESIS' in mission:
    # inputL1list = f'{comSatDir}{mission}/L1/{region}/DESIS-HSI-L1C*{version}'
    # inputL1list = f'{comSatDir}{mission}/L1/{region}/DESIS-HSI-L1C*213'
    # inputL1list = f'{zipListDir}{mission}/L1/{region}/DESIS-HSI-L1C*.zip'
    inputL1list = f'{zipListDir}/DESIS-HSI-L1C*.zip'


    # Single File
    # inputL1list = f'{comSatDir}{mission}/L1/{region}/DESIS-HSI-L1C-DT0599335560_001-20210615T130432-V0213.zip'


else:
    # Setting file will depend on availability of SWIR.
    inputL1listSWIR = f'{zipListDir}{mission}/L1/WV3*S*.zip'
    Sfiles = sorted(glob.glob(inputL1listSWIR))

    # Create a list of SWIR files available
    dateTimeSWIR = []
    for Sfile in Sfiles:
        fName = Sfile.split('/')[-1]
        year = int(fName[3:7])
        mon = int(fName[7:9])
        day = int(fName[9:11])
        hr = int(fName[11:13])
        mn = int(fName[13:15])
        sec = int(fName[15:17])
 
        dateTimeSWIR.append(datetime.datetime(year,mon,day,hr,mn,sec,0,tzinfo=datetime.timezone.utc))
    inputL1list = f'{zipListDir}{mission}/L1/WV3*M*.zip'

args.output = f'{comSatDir}{mission}/'
outdirL1 = f'{comSatDir}{mission}L1/'
outdirL2R = f'{comSatDir}{mission}L2R/'
outdirL2RSWIR = f'{comSatDir}{mission}L2R/SWIR/'
outdirL2W = f'{comSatDir}{mission}L2W/'
outdirL2Wbackup = f'{backupDir}{mission}L2W/'
outdirRGB = f'{comSatDir}{mission}RGB/'


## Create folders where needed
if not os.path.exists(outdirL1):
    os.makedirs(outdirL1)
    print(f'Directory {outdirL1} created')
if not os.path.exists(outdirL2R):
    os.makedirs(outdirL2R)
    print(f'Directory {outdirL2R} created')
if not os.path.exists(outdirL2RSWIR) and 'WorldView3' in mission:
    os.makedirs(outdirL2RSWIR)
    print(f'Directory {outdirL2RSWIR} created')
if not os.path.exists(outdirL2W):
    os.makedirs(outdirL2W)
    print(f'Directory {outdirL2W} created')
if not os.path.exists(outdirRGB):
    os.makedirs(outdirRGB)
    print(f'Directory {outdirRGB} created')

# This is where processing begins
#### Batch imagery based on inputL1list directory files
time_start = datetime.datetime.now()
for inputZip in sorted(glob.glob(inputL1list)):
    start = datetime.datetime.now()

    print(inputZip)
    iFile = os.path.split(inputZip)[1]
    if 'DESIS' in mission:
        tile = iFile.split('-')[3]
        tile = tile.split('_')[1]
        fDate = iFile.split('-')[4]
        fDate = fDate[0:8]+fDate[9:]
        inputZip2 = inputZip
    else:
        fDate = iFile[3:17]
        # Target the raw GeoTIFF directory
        inputZip2 = glob.glob(f'{inputZip}/*')[0]

    year = fDate[0:4]
    mon = fDate[4:6]
    day = fDate[6:8]
    hr = fDate[8:10]
    mn = fDate[10:12]
    sec = fDate[12:14]

    # These are scripted for Acolite's output file name format
    if 'DESIS' in mission:
        # DESIS L1 names differ from imagery collection times by a couple of minutes. Acolite
        # names will be a couple minutes later. This could rollover the day, month, or year...
        # This only protects against day rollovers
        if 60 - float(mn) < 3:
            fbase1R = f'{args.output}DESIS_HSI_{tile}_{year}_{mon}_{day}_*L1R.nc'
            fbase1R2 = f'{outdirL1}DESIS_HSI_{tile}_{year}_{mon}_{day}_*L1R.nc'
            fbase2R = f'{args.output}DESIS_HSI_{tile}_{year}_{mon}_{day}_*L2R.nc'
            fbase2R2 = f'{outdirL2R}DESIS_HSI_{tile}_{year}_{mon}_{day}_*L2R.nc'
            fbase2W = f'{args.output}DESIS_HSI_{tile}_{year}_{mon}_{day}_*L2W.nc'
            fbase2W2 = f'{outdirL2W}DESIS_HSI_{tile}_{year}_{mon}_{day}_*L2W.nc'
            fbase2W3 = f'{outdirL2Wbackup}DESIS_HSI_{tile}_{year}_{mon}_{day}_*L2W.nc'
        else:
            fbase1R = f'{args.output}DESIS_HSI_{tile}_{year}_{mon}_{day}_{hr}_*_L1R.nc'
            fbase1R2 = f'{outdirL1}DESIS_HSI_{tile}_{year}_{mon}_{day}_{hr}_*_L1R.nc'
            fbase2R = f'{args.output}DESIS_HSI_{tile}_{year}_{mon}_{day}_{hr}_*_L2R.nc'
            fbase2R2 = f'{outdirL2R}DESIS_HSI_{tile}_{year}_{mon}_{day}_{hr}_*_L2R.nc'
            fbase2W = f'{args.output}DESIS_HSI_{tile}_{year}_{mon}_{day}_{hr}_*_L2W.nc'
            fbase2W2 = f'{outdirL2W}DESIS_HSI_{tile}_{year}_{mon}_{day}_{hr}_*_L2W.nc'
            fbase2W3 = f'{outdirL2Wbackup}DESIS_HSI_{tile}_{year}_{mon}_{day}_{hr}_*_L2W.nc'
        fbase2RSWIR = ''

    else:
        # WorldView folder names appear to match image times
        fbase1R = f'{args.output}{mission}_{year}_{mon}_{day}_{hr}_{mn}_{sec}_L1R.nc'
        fbase1R2 = f'{outdirL1}{mission}_{year}_{mon}_{day}_{hr}_{mn}_{sec}_L1R.nc'
        fbase2R = f'{args.output}{mission}_{year}_{mon}_{day}_{hr}_{mn}_{sec}_L2R.nc'
        fbase2R2 = f'{outdirL2R}{mission}_{year}_{mon}_{day}_{hr}_{mn}_{sec}_L2R.nc'
        fbase2RSWIR = f'{outdirL2RSWIR}{mission}_{year}_{mon}_{day}_{hr}_{mn}_{sec}_L2R.nc'
        fbase2W = f'{args.output}{mission}_{year}_{mon}_{day}_{hr}_{mn}_{sec}_L2W.nc'
        fbase2W2 = f'{outdirL2W}{mission}_{year}_{mon}_{day}_{hr}_{mn}_{sec}_L2W.nc'

        dateTimeM = datetime.datetime(int(year),int(mon),int(day),int(hr),
            int(mn),int(sec),0,tzinfo=datetime.timezone.utc)

    # Test for L2W. Skip if found; do not overwrite unless told otherwise.
    if ( glob.glob(fbase2W) or glob.glob(fbase2W2) or glob.glob(fbase2W3) ) and not replaceL2W:
        print(f'{fbase2W} or {fbase2W2} found on disk. Skip.')
        continue
    else:
        # Test for L2R
        # Start with SWIR (WV-3)
        if glob.glob(fbase2RSWIR) and not replaceL2R:
            print(f'{glob.glob(fbase2RSWIR)[0]} found. Using')
            args.inputfile = glob.glob(fbase2RSWIR)[0]
            level = 'L2W'
        elif glob.glob(fbase2R2) and not replaceL2R:
            print(f'{glob.glob(fbase2R2)[0]} found. Using')
            args.inputfile = glob.glob(fbase2R2)[0]
            level = 'L2W'

        # Test for L1R.
        elif glob.glob(fbase1R) and not replaceL1R:
            print(f'{glob.glob(fbase1R)[0]} found. Using')
            args.inputfile = glob.glob(fbase1R)[0]
            level = 'L2R'
        elif glob.glob(fbase1R2) and not replaceL1R:
            print(f'{glob.glob(fbase1R2)[0]} found. Using')
            args.inputfile = glob.glob(fbase1R2)[0]
            level = 'L2R'

        # New image. Run from GeoTIFF
        else:
            print(f'{fbase1R2} not found. Running L1R')
            # If the directory is not present already, the file needs to be unzipped
            fpList = inputZip2.split('/')
            basename = fpList[-1]
            testDir = basename.split('.zip')[0]
            testPath = '/'.join(fpList[:-1]) + '/' + testDir
            if not os.path.isdir(testPath):
                with zipfile.ZipFile(inputZip2, 'r') as zip_ref:
                    zip_ref.extractall(testPath)
            args.inputfile = testPath
            level = 'L1R'

        inputfile = args.inputfile.split(',')

        # Match SWIR for WorldView-3
        SWIR=0
        if 'WorldView3' in mission:
            timeDeltas = [abs(dateTimeSWIR[i] - dateTimeM) for i in range(0, len(dateTimeSWIR))]
            # If SWIR image found within 5 minutes of this Multi image...
            if min(timeDeltas) < datetime.timedelta(minutes = 5):
                indexSWIR = timeDeltas.index(min(timeDeltas))
                inputfile_swir = glob.glob(f'{Sfiles[indexSWIR]}/*')[0]
                print('SWIR file found to match VNIR. Using SWIR atmospheric correction.')
                print(inputfile_swir)
                SWIR=1

                # Edit settings file
                args.settings = f'{miss}_{settingsCode}_SWIR.txt'
                with open(args.settings, 'r') as f:
                    lines = f.readlines()
                    for i, line in enumerate(lines):
                        if line.find('inputfile_swir') != -1:
                            lines[i] =f'inputfile_swir={inputfile_swir}\n'
                with open(args.settings, 'w') as f:
                    for line in lines:
                        f.write(line)
        else:
            args.settings = f'{miss}_{settingsCode}.txt'


        #######
        print(f'Launching {level} for {inputfile[0]}')
        ac.acolite.acolite_run(args.settings, inputfile=inputfile, output=args.output)
        #######


        # Move L1R, L2R, and L2W to appropriate folders
        if glob.glob(fbase1R):
            os.replace(glob.glob(fbase1R)[0], f"{outdirL1}{glob.glob(fbase1R)[0].split('/')[-1]}")
        if glob.glob(fbase2R):
            if SWIR:
                os.replace(glob.glob(fbase2R)[0], f"{outdirL2RSWIR}/{glob.glob(fbase2R)[0].split('/')[-1]}")
            else:
                os.replace(glob.glob(fbase2R)[0], f"{outdirL2R}{glob.glob(fbase2R)[0].split('/')[-1]}")
        if glob.glob(fbase2W):
            os.replace(glob.glob(fbase2W)[0], f"{outdirL2W}{glob.glob(fbase2W)[0].split('/')[-1]}")

        stop = datetime.datetime.now()
        diff = stop - start
        print(f'Time elapsed for this file: {diff}')

    # Move RGBs and settings records
    for RGB in glob.glob(f'{args.output}*rgb*.png'):
        os.replace(RGB, f"{outdirRGB}{glob.glob(RGB)[0].split('/')[-1]}")
        # shutil.move(RGB, f"{outdirRGB}{glob.glob(RGB)[0].split('/')[-1]}")
    for l1r in glob.glob(f'{args.output}*l1r_settings.txt'):
        os.replace(l1r, f"{outdirL1}{glob.glob(l1r)[0].split('/')[-1]}")
    for l2r in glob.glob(f'{args.output}*l2r_settings.txt'):
        os.replace(l2r, f"{outdirL2R}{glob.glob(l2r)[0].split('/')[-1]}")

time_stop = datetime.datetime.now() ## time of processing stop
diff = time_stop - time_start
print(f'Time elapsed for this batch: {diff}')