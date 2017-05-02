import numpy as np
import SDSSmodules.SDSSpaths as path
import os
from glob import glob
import fitsio
import shutil


#######################################################################################
# I want to identify  failed quasars in  DR12 that were observed by both  modified  and
# unmodified fibres in order to produce  the analogue of Fig 13 in Margala et al. 2015.
# The ouput of this script are two new folders containing the spectra of those objects.
#######################################################################################

# ============= READ THE spALL-DR12 CATALOG =============
columns = 'LAMBDA_EFF PLATE MJD FIBERID ANCILLARY_TARGET2 OBJTYPE CLASS ZWARNING THING_ID'.split()

# input catalog (spAll)
spAll  = fitsio.read(path.to_spAll_small(), 1, columns=columns )
nspAll = np.size(spAll)

#################################################################
###### DESCRIPTION OF COLUMNES USED TO DEFINE OBJECT TYPES ######
# LAMBDA_EFF      == effective wavelength
# ZWARNING bit 7  == the UNPLUGGED bit
# OBJTYPE         == the object is targeted as $OBJTYPE
# CLASS           == the object is classified as $CLASS
#################################################################

# ============= DEFINE TYPES OF OBJECTS =============
is_offset        = np.where((spAll['LAMBDA_EFF']==4000) & ((spAll['ZWARNING'] & 2**7) == 0))[0]
is_qso           = np.where((spAll['LAMBDA_EFF']==4000) & (spAll['OBJTYPE']=="QSO             ")
                          & (spAll['CLASS']=="QSO   ")  & ((spAll['ZWARNING'] & 2**7) == 0))[0]
is_failed_qso    = np.where((spAll['LAMBDA_EFF']==4000) & (spAll['OBJTYPE']=="QSO             ")
                          & (spAll['CLASS']=="STAR  ")  & ((spAll['ZWARNING'] & 2**7) == 0))[0]
is_spec_stand    = np.where((spAll['LAMBDA_EFF']==5400) & (spAll['OBJTYPE']=="SPECTROPHOTO_STD")
                          & (spAll['CLASS']=="STAR  ")  & ((spAll['ZWARNING'] & 2**7) == 0))[0]
is_offset_stand  = np.where((spAll['LAMBDA_EFF']==4000) & ((spAll['ANCILLARY_TARGET2'] & 2**20) != 0)
                          & (spAll['CLASS']=="STAR  ")  & ((spAll['ZWARNING'] & 2**7) == 0))[0]

print "OBJECT TYPES DEFINED                        UNIQUE"
print '---------------  OBJECTS IN DR12   ---------------'
# read a list of objects with no weather data
no_weather_obj = open(path.to_no_weather_dr12()).readlines()
#####################################################################
#####                OBJECT-RELATED CALCULATIONS               ######
#####################################################################
# =============  QUASARS  =============
# define QSO plates, mjd, fiberid
QSO_plates   = spAll['PLATE'][is_qso]
QSO_mjds     = spAll['MJD'][is_qso]
QSO_fiberids = spAll['FIBERID'][is_qso]
QSO_obj      = zip(QSO_plates, QSO_mjds, QSO_fiberids)

# put names of the failed QSO files in a list
QSO_list  = []

for obj in QSO_obj:
    plate, mjd, fiberid = [int(field) for field in obj]
    obj_name = 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)
    QSO_list.append(obj_name)
print '{0:<43} {1:6d}'.format('QUASARS:', np.size(QSO_list))

# =============  FAILED QUASARS  =============
# define failed QSO plates, mjd, fiberid
FailedQSO_plates   = spAll['PLATE'][is_failed_qso]
FailedQSO_mjds     = spAll['MJD'][is_failed_qso]
FailedQSO_fiberids = spAll['FIBERID'][is_failed_qso]
FailedQSO_obj      = zip(FailedQSO_plates, FailedQSO_mjds, FailedQSO_fiberids)

# put names of the failed QSO files in a list
FailedQSO_list  = []

for obj in FailedQSO_obj:
    plate, mjd, fiberid = [int(field) for field in obj]
    obj_name = 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)
    FailedQSO_list.append(obj_name)
print '{0:<43} {1:6d}'.format('FAILED QUASARS:', np.size(FailedQSO_list))

# =============  OFFSET STANDARDS  =============
# define offset photometric standards' plates, mjd, fiberid
OffsetStand_plates   = spAll['PLATE'][is_offset_stand]
OffsetStand_mjds     = spAll['MJD'][is_offset_stand]
OffsetStand_fiberids = spAll['FIBERID'][is_offset_stand]
OffsetStand_obj      = zip(OffsetStand_plates, OffsetStand_mjds, OffsetStand_fiberids)
# put names of the offset photometric standards files in a list
OffsetStand_list  = []
for obj in OffsetStand_obj:
    plate, mjd, fiberid = [int(field) for field in obj]
    obj_name = 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)
    OffsetStand_list.append(obj_name)
print '{0:<43} {1:6d}'.format('OFFSET SPEC STANDS:', np.size(OffsetStand_list))

# =============  NO-WEATHER OBJECTS  =============
# put names of the no-weather files in a list
no_weather_list = []
for obj in no_weather_obj:
    obj = os.path.basename(obj.rstrip())
    plate, mjd, fiberid = [int(field) for field in os.path.splitext(obj)[0].split('-')[1:]]
    obj_name = 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)
    no_weather_list.append(obj_name)
print '----------  WITH WEATHER DATA AVAILABLE ----------'
############################################################
#####                CROSS - CORRELATION               #####
############################################################
# =============  X-CORRELATE TO FIND OBJECTS W/ WEATHER =============
# cross-correlate to find all failed QSOs that have weather data 
# assuming every element in both lists is unique

# SET DEFINITIONS
# calculates set(OBJECTS_W_WEATHER) = set(ALL_OBJECTS) - set(OBJECTS_W/O_WEATHER)
QSO_set        = set(QSO_list)
FailedQSO_set  = set(FailedQSO_list)
no_weather_set = set(no_weather_list)

# QUASARS
QSO      = list(QSO_set.difference(no_weather_set))
nQSO     = np.size(QSO)
print '{0:<43} {1:6d}'.format('QUASARS:', nQSO)

# FAILED QUASARS
FailedQSO      = list(FailedQSO_set.difference(no_weather_set))
nFailedQSO     = np.size(FailedQSO)
print '{0:<43} {1:6d}'.format('FAILED QUASARS:', nFailedQSO)

# OFFSET STANDARDS
OffsetStand_set  = set(OffsetStand_list)
OffsetStand     = list(OffsetStand_set.difference(no_weather_set))
nOffsetStand     = np.size(OffsetStand)
print '{0:<43} {1:6d}'.format('OFFSET SPEC STANDS:', nOffsetStand)
############################################################
#####                VALIDATION SAMPLE                 #####
############################################################
# =============  X-CORRELATE TO FIND OBJECTS IN VALIDATION SET =============
# a new array flags those objects with plate number in the validation set
is_in_validation = np.full(nspAll, fill_value = False, dtype=np.dtype('bool'))
validation_set = set([6130, 6131, 6135, 6136, 6155, 6157, 6290, 6293, 6296, 6297, 6298, 6307, 6506, 6509, 6590, 6681, 6734, 6816, 6986])

# find validation failed QSO
i = 0
for plate in spAll['PLATE']:
    if plate in validation_set:
        is_in_validation[i] = True
    i += 1

is_validation_qso           = np.where((spAll['LAMBDA_EFF']==4000) & (spAll['OBJTYPE']=="QSO             ")
                                     & (spAll['CLASS']=="QSO   ")  & (is_in_validation == True)  & (spAll['ZWARNING']==0))[0]
is_validation_failed_qso    = np.where((spAll['LAMBDA_EFF']==4000) & (spAll['OBJTYPE']=="QSO             ")
                                     & (spAll['CLASS']=="STAR  ")  & (is_in_validation == True)  & (spAll['ZWARNING']==0))[0]
is_validation_offset_stand  = np.where((spAll['LAMBDA_EFF']==4000) & ((spAll['ANCILLARY_TARGET2'] & 2**20) != 0)
                                     & (spAll['CLASS']=="STAR  ")  & (is_in_validation == True)  & (spAll['ZWARNING']==0))[0]

nValidationQSO         = is_validation_qso.size
nValidationFailedQSO   = is_validation_failed_qso.size
nValidationOffsetStand = is_validation_offset_stand.size

print '---------------  IN VALIDATION SET ---------------'
print '{0:<43} {1:6d}'.format('QUASARS:', nValidationQSO)
print '{0:<43} {1:6d}'.format('FAILED QUASARS:', nValidationFailedQSO)
print '{0:<43} {1:6d}'.format('OFFSET SPEC STANDS:', nValidationOffsetStand)

#######################################################################################
# Find all observations of objects in the failed QSO validation data set
# (both with modified and unmodified SDSS fibres). THING_ID is unique to every object
# and can be used to identify those objects.
#######################################################################################
# defining THING_ID arrays for each object type

QSO_thingids     = spAll['THING_ID'][is_validation_qso]
QSO_thingids_set = set(np.unique(QSO_thingids))

FailedQSO_thingids     = spAll['THING_ID'][is_validation_failed_qso]
FailedQSO_thingids_set = set(np.unique(FailedQSO_thingids))

OffsetStand_thingids     = spAll['THING_ID'][is_validation_offset_stand]
OffsetStand_thingids_set = set(np.unique(OffsetStand_thingids))

# find all observations of objects with THING_IDs that are in the appropriate set
thingid_is_in_QSO         = np.full(nspAll, fill_value = False, dtype=np.dtype('bool'))
thingid_is_in_FailedQSO   = np.full(nspAll, fill_value = False, dtype=np.dtype('bool'))
thingid_is_in_OffsetStand = np.full(nspAll, fill_value = False, dtype=np.dtype('bool'))
for i, thingid in enumerate(spAll['THING_ID']):
    if thingid in QSO_thingids_set:
        thingid_is_in_QSO[i] = True
    if thingid in FailedQSO_thingids_set:
        thingid_is_in_FailedQSO[i] = True
    if thingid in OffsetStand_thingids_set:
        thingid_is_in_OffsetStand[i] = True

# select observations of validation objects made by (un)modified fibres
# some offset standards might have been observed by the unmodified fibres previously
print '--------------------------------------------------'
print 'USING "THING_ID", DR12 objects validation objects '
print '{0:<43}'.format('with ZWARNING = 0 were observed by fibres:   ')
print '----------------------  UNMODIFIED -- MODIFIED  --'
# FQ now stands for FailedQSO
# OS now stands for OffsetStand

# QSO 
QSO_unmodified  = np.where((spAll['LAMBDA_EFF']==5400) & (thingid_is_in_QSO == True) & (spAll['ZWARNING']==0))[0]
QSO_modified    = np.where((spAll['LAMBDA_EFF']==4000) & (thingid_is_in_QSO == True) & (spAll['ZWARNING']==0))[0]
nQSO_unmodified = np.size(QSO_unmodified)
nQSO_modified   = np.size(QSO_modified)
print '{0:<24} {1:9d} {2:11d}'.format('QUASARS:', nQSO_unmodified, nQSO_modified)
QSO_unmod_tuple    = zip(spAll['PLATE'][QSO_unmodified],spAll['MJD'][QSO_unmodified],spAll['FIBERID'][QSO_unmodified])
QSO_mod_tuple      = zip(spAll['PLATE'][QSO_modified],  spAll['MJD'][QSO_modified],  spAll['FIBERID'][QSO_modified])

# FQ
FQ_unmodified  = np.where((spAll['LAMBDA_EFF']==5400) & (thingid_is_in_FailedQSO == True) & (spAll['ZWARNING']==0))[0]
FQ_modified    = np.where((spAll['LAMBDA_EFF']==4000) & (thingid_is_in_FailedQSO == True) & (spAll['ZWARNING']==0))[0]
nFQ_unmodified = np.size(FQ_unmodified)
nFQ_modified   = np.size(FQ_modified)
print '{0:<24} {1:9d} {2:11d}'.format('FAILED QUASARS:', nFQ_unmodified, nFQ_modified)
FQ_unmod_tuple    = zip(spAll['PLATE'][FQ_unmodified],spAll['MJD'][FQ_unmodified],spAll['FIBERID'][FQ_unmodified])
FQ_mod_tuple      = zip(spAll['PLATE'][FQ_modified],  spAll['MJD'][FQ_modified],  spAll['FIBERID'][FQ_modified])

# OS
OS_unmodified  = np.where((spAll['LAMBDA_EFF']==5400) & (thingid_is_in_OffsetStand == True) & (spAll['ZWARNING']==0))[0]
OS_modified    = np.where((spAll['LAMBDA_EFF']==4000) & (thingid_is_in_OffsetStand == True) & (spAll['ZWARNING']==0))[0]
nOS_unmodified = OS_unmodified.size
nOS_modified   = OS_modified.size
print '{0:<24} {1:9d} {2:11d}'.format('OFFSET SPEC STANDS:', nOS_unmodified, nOS_modified)
OS_unmod_tuple    = zip(spAll['PLATE'][OS_unmodified],spAll['MJD'][OS_unmodified],spAll['FIBERID'][OS_unmodified])
OS_mod_tuple      = zip(spAll['PLATE'][OS_modified],  spAll['MJD'][OS_modified],  spAll['FIBERID'][OS_modified])

# are there Offset Spec standards that were observed by both modified and unmodified fibres?
OS_unmod_set          = set(spAll['THING_ID'][OS_unmodified])
OS_mod_set            = set(spAll['THING_ID'][OS_modified])
OS_is_obs_both = OS_unmod_set.intersection(OS_mod_set)
print 'OS unmodified',len(OS_unmod_set)
print 'OS modified',len(OS_mod_set)
print 'OS intersection',len(OS_is_obs_both)
print 'OS union',len(OS_unmod_set.union(OS_mod_set))

print 'OFFSET SPEC STAND observations by unmodified fibres' 
for tuple in OS_unmod_tuple:
    plate, mjd, fibre = [field for field in tuple]
    #print '{0:4d} | {1:5d} | {2:04d}'.format(plate,mjd,fibre)

# create a list of names for unmodified and modified targets to be copied
FQ_unmodified_names = []
FQ_modified_names   = []
for obj in xrange(FQ_unmodified.size):
    plate,mjd,fibre = [field for field in FQ_unmod_tuple[obj]]
    obj_name  = 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)
    spec_name = os.path.join(path.to_DR12(),'failedQSO',obj_name)
    #print spec_name
    FQ_unmodified_names.append(spec_name)
    
for obj in xrange(FQ_modified.size):
    plate,mjd,fibre = [field for field in FQ_mod_tuple[obj]]
    obj_name  = 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)
    spec_name = os.path.join(path.to_DR12(),'failedQSO',obj_name)
    #print spec_name
    FQ_modified_names.append(spec_name)

print np.size(FQ_unmodified_names)
print np.size(FQ_modified_names)



