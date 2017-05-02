#!/usr/bin/env python
"""
Demonstrates how to apply throughput corrections to individual BOSS spectrum.
Saves the corrected spectrum into a FITS file.
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
from SDSSmodules.SDSSfiles import read_resid_corr
from SDSSmodules.SDSScorrections import apply_resid_corr

def main():
    np.set_printoptions(threshold=6000)
    # PARSER 
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default='/Volumes/T2TB/Isotropy/spectra/BOSScorrections/tpcorr-0.1/tpcorr.hdf5',
        help='Throughput correction filename, required')
    parser.add_argument('--spec-dir', type=str, default=None,
        help='Path to directory containing individual spec files')
    parser.add_argument('--all', action='store_true',
        help='Save corrected spectra for all available objects of the given data sample')
    parser.add_argument('--data', choices=['QSO','FQSO','STAR'], nargs=1, type=str, default='QSO',
        help='Object type of the chosen spectrum')
    parser.add_argument('--nobalmer', action='store_true', default=True,
        help='Balmer corrections flag. If the flag is true, no Balmer corrections are applied.')
    args = parser.parse_args()

    # data sample name
    sample = args.data[0]
    # MARGALA CORRECTION FILE
    tpcorr = h5py.File(args.tpcorr, 'r')
    tpcorr_wave = tpcorr['wave'].value
    # BALMER CORRECTION FILE
    if args.nobalmer == True:
        resid_corr = None
    elif args.nobalmer == False:
        resid_corr = read_resid_corr('/Volumes/T2TB/Isotropy/PySDSS/residcorr_v5_4_45.dat')
    print resid_corr

    # Parse arguments: extract plate, mjd, fiberid as integers and the spec filename
    ids = []
    sid_tot = []
    if args.spec_dir:
        spectral_list       = return_spectral_list(args.spec_dir)
        uncorrected_files   = spectral_list[0]
        total_files = uncorrected_files
        for file in uncorrected_files:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            ids.append((plate, mjd, fiberid))
            sid_tot.append('%04d%5d%04d' %(plate, mjd, fiberid))
    dtypedict = {'QSO':'QUASARS', 'FQSO':'FAILED QUASARS', 'STAR':'STARS'}
    if args.all:
        print "CORRECTING ALL AVAILABLE SPECTRA OF {data}".format(data=dtypedict[sample])
        PathToSpectra = os.path.join("/Volumes/T2TB/Isotropy/spectra/SDSS_DR12",sample)
        total_files = glob(os.path.join(PathToSpectra,"*.fits"))
        print "Spectra available: {0:8d}".format(len(total_files))
        for file in total_files:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            ids.append((plate, mjd, fiberid))
            sid_tot.append('%04d%5d%04d' %(plate, mjd, fiberid))
    #else:
    #    parser.error('Must specify either a spec file using the --spec-file option or path to directory containing spec files and a target.')
    #sid_tot = np.array(sid_tot)
    sid_tot_array = np.array(sid_tot)
    sid_tot = set(sid_tot)
    
    if args.all:
        if args.nobalmer == False:
            dirpath = os.path.join('/Volumes/T2TB/Isotropy/spectra/SDSS_DR12','{data}_{name}'.format(data=sample,name='Margala'))
        elif args.nobalmer == True:
            dirpath = os.path.join('/Volumes/T2TB/Isotropy/spectra/SDSS_DR12','{data}_{name}_nobalmer'.format(data=sample,name='Margala'))
    else:
        dirpath = os.path.join(os.path.normpath(args.spec_dir),'corrected','margala')
    # remove already corrected spectra from the list
    corrected_files = glob(os.path.join(dirpath,"*.fits"))
    print "Spectra previously corrected: {0:8d}".format(len(corrected_files))
    sid_corr = []
    for file in corrected_files:
        spec_basename = os.path.basename(file).split('.fits')[0]
        plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
        sid_corr.append('%04d%5d%04d' %(plate, mjd, fiberid))
    sid_corr = set(sid_corr)
    sid_unc = sid_tot.difference(sid_corr)
    uncorrected_files = []
    for i,file in enumerate(total_files):
        if sid_tot_array[i] in sid_unc:
            uncorrected_files.append(file)

    print "Spectra to be corrected:     {0:8d}".format(len(uncorrected_files))
    file_array = np.array(uncorrected_files)
    for file in uncorrected_files:
        print file
        HDU = fitsio.FITS(file)
        h    = HDU[0].read_header()
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
        
        # Save the data into a fits file
        spec = spectrum()

        spec.plateid = plate
        spec.mjd     = mjd
        spec.fiberid = fiberid
        spec.ra      = ra
        spec.dec     = dec
        spec.z       = z
        spec.beginwl = h['COEFF0']
        
        spec.npix       = np.size(wavelength)
        spec.wave       = wavelength
        spec.flux       = flux
        spec.ivar       = ivar
        spec.flux_error = np.sqrt(1./ivar)
        apply_resid_corr(spec, resid_corr)

        corrected_flux = spec.flux*resampled_correction
        corrected_ivar = spec.ivar/resampled_correction**2
        
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
