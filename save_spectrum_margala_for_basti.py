#!/usr/bin/env python
"""
Demonstrates how to apply throughput corrections to individual BOSS spectrum.

There a few different ways to specify a spectrum. In order of precedence:
    1) using the spec filename
    2) a target id string (must also provide top-level directory containing spec files)
    3) individual plate, mjd, and fiberid (must also provide top-level directory containing spec files)

Usage: 

python example_usage.py --tpcorr tpcorr-rc2.hdf5 --spec-file ~/data/boss/v5_7_0/spectra/lite/3615/spec-3615-55445-0008.fits 

python example_usage.py --tpcorr tpcorr-rc2.hdf5 --spec-dir ~/data/boss/v5_7_0/spectra/lite --target 3615-55445-8

python example_usage.py --tpcorr tpcorr-rc2.hdf5 --spec-dir ~/data/boss/v5_7_0/spectra/lite --plate 3615 --mjd 55445 --fiberid 8

"""
import argparse
import os
import sys

import h5py
import numpy as np

from astropy.io import fits
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from datetime import datetime
from corrutilis import *
from SDSSmodules.SDSSclasses import spectrum


def main():

    # PARSER 
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default='/Volumes/Transcend/Isotropy/spectra/BOSScorrections/tpcorr-0.1/tpcorr.hdf5',
        help='throughput correction filename, required')
    parser.add_argument('--all', action='store_true',
        help='save corrected spectra for all available QSOs')
    args = parser.parse_args()

    # MARGALA CORRECTION FILE
    tpcorr = h5py.File(args.tpcorr, 'r')
    tpcorr_wave = tpcorr['wave'].value

    # DIRECTORY TO SAVE THE CORRECTED SPECTRA TO, USER-DEFINED
    dirpath = os.path.normpath('/Volumes/Transcend/Isotropy/spectra/SDSS_DR12/QSO_Margala')
    # -------------------------------------------------------------------------------------
    # Parse arguments: extract plate, mjd, fiberid as integers and the spec filename
    # ids = list of tuples   : (plate,mjd,fiberid)
    # sid = list of strings  : string(plate+mjd+fiberid)  ---> used to check already corrected spectra
    ids = []
    sid_tot = []
    if args.all:
        print "CORRECTING ALL AVAILABLE SPECTRA"
        # CREATE A LIST OF ALL FILES TO BE CORRECTED : total_files
        PathToSpectra = os.path.normpath("/Volumes/Transcend/Isotropy/spectra/SDSS_DR12/QSO")
        total_files = glob(os.path.join(PathToSpectra,"*.fits"))
        print "Spectra available: {0:8d}".format(len(total_files))
        for file in total_files:
            # APPEND THE UNIQUE IDENTIFIER TO sid_tot
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            ids.append((plate, mjd, fiberid))
            sid_tot.append('%04d%5d%04d' %(plate, mjd, fiberid))

    # CREATE A PYTHON SET OBJECT FOR ALL SIDs (ALLOWS INTERSECTION, DIFERENCE)
    sid_tot_array = np.array(sid_tot)
    sid_tot = set(sid_tot)
    
       
    # CREATE A LIST OF ALL FILES ALREADY CORRECTED AND PRESENT IN THE OUTPUT DIRECTORY: corrected_files
    corrected_files = glob(os.path.join(dirpath,"*.fits"))
    print "Spectra corrected: {0:8d}".format(len(corrected_files))
    sid_corr = []
    for file in corrected_files:
        # APPEND THE UNIQUE IDENTIFIERS TO sid_corr
        spec_basename = os.path.basename(file).split('.fits')[0]
        plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
        sid_corr.append('%04d%5d%04d' %(plate, mjd, fiberid))
    # CREATE A PYTHON SET OBJECT FOR CORRECTED SIDs
    sid_corr = set(sid_corr)
    sid_unc = sid_tot.difference(sid_corr)

    # FIND THE DIFFERENCE OF THE TWO SETS (sid_tot - sid_corr)
    # APPEND SIDs TO uncorrected_files ---> used later
    uncorrected_files = []
    for i,file in enumerate(total_files):
        if sid_tot_array[i] in sid_unc:
            uncorrected_files.append(file)

    # ITERATE OVER FILES IN uncorrected_files AND APPLY THE CORRECTION
    print "To be corrected:  {0:8d}".format(len(uncorrected_files))
    file_array = np.array(uncorrected_files)
    for file in uncorrected_files:
        HDU = fitsio.FITS(file)
        flux = HDU[1]['flux'][:]
        ivar = HDU[1]['ivar'][:]
        wavelength = np.power(10, HDU[1]['loglam'][:])
        ra   = HDU[2]['ra'][0][0]
        dec  = HDU[2]['dec'][0][0]
        z    = HDU[2]['z'][0][0]
        HDU.close()

        plate,mjd,fiberid = ids[np.argwhere(file_array==file)[0][0]]
        # Read the target's throughput correction vector
        tpcorr_key = '%s/%s/%s' % (plate, mjd, fiberid)
        try:
            correction = tpcorr[tpcorr_key].value
        except:
            print "No data on object found"
            continue
        # Create an interpolated correction function
        correction_interp = interp1d(tpcorr_wave, correction, kind='linear')

        # Sample the interpolated correction using the observation's wavelength grid
        resampled_correction = correction_interp(wavelength)
    
        # Apply the correction to the observed flux and ivar
        corrected_flux = flux*resampled_correction
        corrected_ivar = ivar/resampled_correction**2
        # Save the data into a fits file
        spec = spectrum()

        spec.plateid = plate
        spec.mjd     = mjd
        spec.fiberid = fiberid
        spec.ra      = ra
        spec.dec     = dec
        spec.z       = z

        spec.npix       = np.size(wavelength)
        spec.wave       = wavelength
        spec.flux       = flux
        spec.ivar       = ivar
        spec.flux_corr  = corrected_flux
        spec.flux_error = np.sqrt(1./corrected_ivar)
        spec.corr       = resampled_correction

        spec.basename  = os.path.basename(file)
        spec.dirname   = os.path.dirname(file)

        spec.__save_fitsio__(dirpath=dirpath)
        del(spec)
        
    tpcorr.close()


if __name__ == '__main__':
    main()
