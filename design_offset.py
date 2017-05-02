import yanny as y
import numpy as n
import os
from glob import glob
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
import sys
import fitsio
import matplotlib.pyplot as plt
import urllib2
"""
USAGE: python design_offset.py

PROCEDURE
(1) FROM INDIVIDUAL spec FILES: CALCULATE THE OBSERVED HOUR ANGLE & ALTITUDE FOR A SINGLE OBJECT
(2) FROM plateHoles FILES: FOR A GIVEN OBJECT, READ THE DESIGN HOUR ANGLE & ALTITUDE
"""
def RadToDeg(rad):
	""" Convert radians to degrees"""
	deg = rad*180./n.pi
	return deg
def DegToRad(deg):
	""" Convert degrees to radians"""
	rad = deg*n.pi/180.
	return rad

# READ THE PATHS TO THE SPECTRA
# This needs to be changed by the user to reflect the local path to the spectra
# In this case "spectralDir" has 2 subfolders, "QSO" and "FQSO", containing the quasar and failed quasar spectra
spectralDir   = os.path.normpath('/Volumes/T2TB/Isotropy/spectra/SDSS_DR12')
# QSO DIRECTORY
QSODir        = os.path.join(spectralDir,'QSO')
QSOList       = glob(os.path.join(QSODir,'*.fits'))
print QSODir
# FAILED QSO DIRECTORY
FQSODir       = os.path.join(spectralDir,'FQSO')
FQSOList      = glob(os.path.join(FQSODir,'*.fits'))
print FQSODir
# ALL SPECTRA
# Combine the two lists into a single list that is used by the program
spectralList  = QSOList + FQSOList
print "{0:>8d} spectra have been found.".format(len(spectralList))

# CREATE THE LIST OF PLATES
plateList = []
for file in spectralList:
    plate = file.split('-')[1]
    plateList.append(plate)
plateList = n.unique(plateList)
print "{0:>4d} plates are unique".format(len(plateList))

# LEAVE ONLY A SINGLE FILE FOR EACH PLATE
shortFileList = []
for plate in plateList:
    # print plate
    for file in spectralList:
        fileplate = file.split('-')[1]
        if fileplate == plate:
            shortFileList.append(file)
            break
    spectralList = [spec for spec in spectralList if plate not in spec]
    #print len(spectralList)
#print len(shortFileList)
spectralList = shortFileList

# DEFINE PATHS TO plateHoles FILES (to be used in calculations below)
# This needs to be changed by the user to reflect the local path to the "plateHole" files. 
# spCFrame 
#spCFrameDir  = os.path.normpath('/Volumes/T2TB/Isotropy/spectra/SDSS_DR12/spCFrame/v5_7_0/')
#spCFrameList = glob(os.path.join(spCFrameDir,'*.fits.gz'))
# plateHoles
plateHoleDir = os.path.normpath('/Volumes/T2TB/Isotropy/spectra/SDSS_DR12/plateHoles/')

# DEFINE OBSERVATORY COORDINATES (to be used in calculations below)
# APO Coordinates
APO = EarthLocation.of_site('Apache Point Observatory')
lat = APO.latitude.degree
lon = APO.longitude.degree
print 'APO Location: Lat = {0:.3}, Lon = {1:.3}'.format(lat, lon)

# ARRAYS WHICH WILL CONTAIN HOUR ANGLES (HA) AND ALTITUDES (Alt)
deltaHAArray  = []
deltaAltArray = []

for file in spectralList:
	# for each file in the list:
	# (1) calculate HA
	# (2) calculate Alt

    	# read the header of the first HDU in the fits file
    	h0           = fitsio.read_header(file,0)
    	# read the RA & Dec of telescope pointing at observation (raPointObs & decPointObs)
    	raPointObs   = h0['RA']
    	decPointObs  = h0['DEC']
	# read the PLATE number and number of exposures this co-added spectrum contains
	plate        = h0['PLATEID']
	nexp         = h0['NEXP']

	# read the plateHole file for the plate if it exists. if the file does not exist, download it from SDSS servers
	plateHoleFilename = os.path.join(plateHoleDir,"plateHoles-{0:06d}.par".format(plate))
	if os.path.isfile(plateHoleFilename)==False:
		print "plateHole file for plate {0} does not exist, downloading...".format(plate)
		url               = "https://svn.sdss.org/public/data/sdss/platelist/trunk/plates/00{id}XX/{plate:06d}/plateHoles-{plate:06d}.par".format(id=str(plate)[0:2],plate=plate)
		#print url
		plateHoleURL      = urllib2.urlopen(url)
		plateHoleFile     = open(plateHoleFilename,"w")
		data              = plateHoleURL.read(40)
		plateHoleFile.write(data)
		HAData            = data.split('\n')[1]
		HADes             = n.float(HAData[3:8])
	elif os.path.isfile(plateHoleFilename)==True:
		print "plateHole file for plate {0} exists, opening...".format(plate)
		plateHoleFile     = y.yanny(plateHoleFilename,np=True)
		# read the designed hourangle (in degrees)
		HADes             = n.float(plateHoleFile["ha"][0:5])
	# change HA to be in range [-180,180] if needed:
	if   HADes>180:
		HADes = HADes - 360.
	elif HADes<-180:
		HADes = HAdes + 360.
	# calculate design altitude
	AltDes       = RadToDeg(n.arcsin(n.sin(DegToRad(h0["PLUG_DEC"]))*n.sin(DegToRad(lat))+n.cos(DegToRad(h0["PLUG_DEC"]))*n.cos(DegToRad(lat))*n.cos(DegToRad(HADes))))
	print "Des hourangle: {0:-10.2f}   Des altitude: {1:-10.2f}".format(HADes, AltDes)
	# now calculate the median of the midpoints of individual exposures 
	# (this follows Margala's procedure in Sec 2, page 2, above figure 1)
	# this is possibly buggy, as the procedure is not explained in details
	taibegArr    = n.zeros(nexp)
	taiendArr    = n.zeros(nexp)
	for exposure in xrange(nexp):
		# read the header of each exposure and get the observation begin and end times
		h                   = fitsio.read_header(file,4+exposure)
		taibegArr[exposure] = h["TAI-BEG"]
		taiendArr[exposure] = h["TAI-END"]
	# get the median of the observing time (in MJD) from all exposures comprising the co-added spectrum
	taiobs       = n.median([beg+(end-beg)/2. for beg,end in zip(taibegArr,taiendArr)]) / 86400.
	# factor 86400 converts seconds into MJD
	# calculate the Local Sidereal Time at APO for the median of the midpoints
	# astropy can calculate "apparent" and "median" sidereal time. The difference is in decimals, but may be relevant. Needs further inquiry.
	LST          = Time(taiobs,scale='tai',format='mjd').sidereal_time("apparent",(lon,lat))
	# calculate the actual observing hour angle (in degrees)
	HAObs        = LST.deg - h0['PLUG_RA']
	# change HA to be in range [-180,180] if needed:
	if   HAObs>180:
		HAObs = HAObs - 360.
	elif HAObs<-180:
		HAObs = HAObs + 360.
	AltObs       = RadToDeg(n.arcsin(n.sin(DegToRad(h0["PLUG_DEC"]))*n.sin(DegToRad(lat))+n.cos(DegToRad(h0["PLUG_DEC"]))*n.cos(DegToRad(lat))*n.cos(DegToRad(HAObs))))
	print "Obs hourangle: {0:-10.2f}   Obs altitude: {1:-10.2f}".format(HAObs, AltObs)
	print "{:=>70}".format('')
	# calculate the difference in HA (observed-designed)
	deltaHA      = HAObs - HADes
    	deltaHAArray.append(deltaHA)
    	#deltaAltArray.append(deltaDec)
#sys.exit()
deltaHAArray = n.array(deltaHAArray)
#deltaAltArray = n.array(deltaAltArray)

# slice the data to make a histogram analogous to Margala et al. figure 2
slice = n.where(abs(deltaHAArray)<=45)[0]
data  = deltaHAArray[slice]

flag = 1
nbins = 45
if flag == 0:
    fig, ax = plt.subplots(1, 2)
    ax[0].hist(deltaRAarray,  bins = nbins, histtype = 'bar', normed = True)
    #ax[1].hist(deltaDecarray, bins = nbins, histtype = 'bar', normed = True)
elif flag == 1:
    	fig, ax = plt.subplots(1, 1)
    	ax.hist(data,   bins = nbins, histtype = 'stepfilled', normed = False)
	ax.set_xticks(n.arange(-45,46,15))
	ax.set_xlim(-45,45)
    	ax.set_xlabel('Degrees')
    	ax.set_ylabel('Number of plates')
plt.grid(True)
plt.savefig('/Volumes/T2TB/Isotropy/offset/plots/deltaOffset.pdf')
plt.show()
sys.exit()
