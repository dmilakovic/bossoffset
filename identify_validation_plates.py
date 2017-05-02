#################################################################
# Use the DR12 spAll catalog to identify validation plates
#################################################################
# use: "python identify_validation_plates.py"
#################################################################
import numpy as np
import fitsio
import argparse
import os
import sys
from glob import glob
import SDSSmodules.SDSSpaths as path

# -------- READ THE spALL-DR12 CATALOG ---------

print "READING FROM spAll FITS FILE"
columns = 'LAMBDA_EFF PLATE MJD FIBERID ANCILLARY_TARGET2 OBJTYPE CLASS ZWARNING'.split()
spAll = fitsio.read(path.to_spAll_small(), 1, columns=columns )

#################################################################
###### DESCRIPTION OF COLUMNES USED TO DEFINE OBJECT TYPES ######
# LAMBDA_EFF      == effective wavelength
# ZWARNING bit 7  == the UNPLUGGED bit
# OBJTYPE         == the object is targeted as $OBJTYPE
# CLASS           == the object is classified as $CLASS
#################################################################

# ================= DEFINE TYPES OF OBJECTS =================
is_offset        = np.where((spAll['LAMBDA_EFF']==4000) & ((spAll['ZWARNING'] & 2**7) == 0))[0]
is_qso           = np.where((spAll['LAMBDA_EFF']==4000) & (spAll['OBJTYPE']=="QSO             ")
                          & (spAll['CLASS']=="QSO   ")  & ((spAll['ZWARNING'] & 2**7) == 0))[0]
is_failed_qso    = np.where((spAll['LAMBDA_EFF']==4000) & (spAll['OBJTYPE']=="QSO             ")
                          & (spAll['CLASS']=="STAR  ")  & ((spAll['ZWARNING'] & 2**7) == 0))[0]
is_spec_stand    = np.where((spAll['LAMBDA_EFF']==5400) & (spAll['OBJTYPE']=="SPECTROPHOTO_STD")
                          & (spAll['CLASS']=="STAR  ")  & ((spAll['ZWARNING'] & 2**7) == 0))[0]
is_offset_stand  = np.where((spAll['LAMBDA_EFF']==4000) & ((spAll['ANCILLARY_TARGET2'] & 2**20) != 0)
                          & (spAll['CLASS']=="STAR  ")  & ((spAll['ZWARNING'] & 2**7) == 0))[0]
print "OBJECTS DEFINED"

# ================= VALIDATION SET OF PLATES =================
# validation plates are defined as those plates with at least 10 offset standards per spectrograph
# COMMENT: obs is short for observation
offset_stands_fibres = spAll['FIBERID'][is_offset_stand]
offset_stands_plates = spAll['PLATE'][is_offset_stand]
offset_stands_mjds   = spAll['MJD'][is_offset_stand]

offset_stands_obs = zip(offset_stands_plates,offset_stands_mjds)

unique_plates = np.unique(offset_stands_plates)
unique_obs = set(zip(spAll['PLATE'][is_offset_stand], spAll['MJD'][is_offset_stand]))

# for all unique plates, check how many spectroscopic standards are observed by spectrographs #1 and #2,
# return only those that satisfy the 10 objects per spectrograph condition
unique_obs_at_least_10_per_spec = []
unique_plates_at_least_10_per_spec = []

j=0
for obs in unique_obs:
    spec1 = 0
    spec2 = 0
    plate, mjd = [field for field in obs]
    ii    = np.where((offset_stands_plates == plate) & (offset_stands_mjds == mjd))[0]
    fibres = offset_stands_fibres[ii]
    for fibre in fibres:
        if fibre <= 500:
            spec1 += 1
        if fibre > 500:
            spec2 += 1
    j += 1
    if spec1>=10 and spec2>=10:
        unique_obs_at_least_10_per_spec.append((plate,mjd))
        unique_plates_at_least_10_per_spec.append(plate)

validation_obs    = set(unique_obs_at_least_10_per_spec)
validation_plates = set(unique_plates_at_least_10_per_spec)

# ================= CROSS-MATCH WITH THE OBJECT TYPES =================

# locate those objects whose plate number is in the validation set
is_in_validation = np.full(np.size(spAll['PLATE']), fill_value = False, dtype=np.dtype('bool'))
i = 0
observations = zip(spAll['PLATE'],spAll['MJD'])
for obs in observations:
    plate, mjd = [field for field in obs]
    if plate in validation_plates:
        is_in_validation[i] = True
    #print i, obs, is_in_validation[i]
    i += 1

# cross-match to find objects with PLATE that is in the validation set
is_validation_offset        = np.where((spAll['LAMBDA_EFF']==4000) & (is_in_validation == True)  & (spAll['ZWARNING']==0))[0]
is_validation_qso           = np.where((spAll['LAMBDA_EFF']==4000) & (spAll['OBJTYPE']=="QSO             ")
                                     & (spAll['CLASS']=="QSO   ")  & (is_in_validation == True)  & (spAll['ZWARNING']==0))[0]
is_validation_failed_qso    = np.where((spAll['LAMBDA_EFF']==4000) & (spAll['OBJTYPE']=="QSO             ")
                                     & (spAll['CLASS']=="STAR  ")  & (is_in_validation == True)  & (spAll['ZWARNING']==0))[0]
is_validation_spec_stand    = np.where((spAll['LAMBDA_EFF']==5400) & (spAll['OBJTYPE']=="SPECTROPHOTO_STD")
                                     & (spAll['CLASS']=="STAR  ")  & (is_in_validation == True)  & (spAll['ZWARNING']==0))[0]
is_validation_offset_stand  = np.where((spAll['LAMBDA_EFF']==4000) & ((spAll['ANCILLARY_TARGET2'] & 2**20) != 0)
                                     & (spAll['CLASS']=="STAR  ")  & (is_in_validation == True)  & (spAll['ZWARNING']==0))[0]

# ================= PRINT OUT THE RESULTS =================
print '{0:>29s} {1:>10s}'.format('N_DR12', 'N_valid')
print '-------------------------------------------------'
print '{0:>20s} {1:8d} {2:8d}'.format('OFFSET TARGETS:',  np.size(is_offset),      np.size(is_validation_offset))
print '{0:>20s} {1:8d} {2:8d}'.format('QUASARS:',         np.size(is_qso),         np.size(is_validation_qso))
print '{0:>20s} {1:8d} {2:8d}'.format('FAILED QUASARS:',  np.size(is_failed_qso),  np.size(is_validation_failed_qso))
print '{0:>20s} {1:8d} {2:8d}'.format('SPEC STANDARDS:',  np.size(is_spec_stand),  np.size(is_validation_spec_stand))
print '{0:>20s} {1:8d} {2:8d}'.format('OFFSET STANDARDS:',np.size(is_offset_stand),np.size(is_validation_offset_stand))
print '-------------------------------------------------'
print 'VALIDATION PLATES:'
print len(validation_plates), '{0:10s}'.format(validation_plates)

