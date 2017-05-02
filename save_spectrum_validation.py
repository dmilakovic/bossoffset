#
# Example script for looking at BOSS spectra and redshift fits via Python.
#
# Copyright Adam S. Bolton, Oct. 2009.
# Licensed for free and unencumbered use in public domain.
# No warranty express or implied.
#

# Imports:
import numpy as n
import pyfits as pf
import matplotlib as mpl
from matplotlib import pyplot as p
import sys
import os
import argparse
from SDSSmodules.SDSSclasses import spectrum
from SDSSmodules.SDSScorrections import apply_resid_corr
""" 
Program to read and save the spectra of objects which were reduced using standard stars observed using modified fibres. 
"""
# Read the CSV file
def main(args):
    # files
    pathToFile = args.list
    # read Balmer corrections
    from SDSSmodules.SDSSfiles import read_resid_corr
    resid_corr = read_resid_corr('/Volumes/T2TB/Isotropy/PySDSS/residcorr_v5_4_45.dat')
    objects = []
    with open(pathToFile,'r') as file:
        for line in file:
            plate,mjd,fiber = [field for field in line.rstrip('\n').split(',')]
            objects.append((plate, mjd, fiber))
    k=0
    for object in objects:
        # define a spectrum instance:
        spec = spectrum()
        # read plate, mjd, fiber:
        plate, mjd, fiber = [int(value) for value in object]
        spec.plate   = plate
        spec.MJD     = mjd
        spec.fiberid = fiber
        spec.dirname = savedir
        spec.basename = 'spec-{plate:4d}-{mjd:5d}-{fiber:04d}.fits'.format(plate=plate,mjd=mjd,fiber=fiber)
        i = fiber - 1 #systematic offset
        # Pick your plate/mjd and read the data:
        spfile = '{top}{plate}/spPlate-{plate}-{mjd}.fits'.format(top=topdir,plate=plate,mjd=mjd)
        zbfile = '{top}{plate}/v5_7_0/spZbest-{plate}-{mjd}.fits'.format(top=topdir,plate=plate,mjd=mjd)
        zafile = '{top}{plate}/v5_7_0/spZall-{plate}-{mjd}.fits'.format(top=topdir,plate=plate,mjd=mjd)
        
        hdulist = pf.open(spfile)
        c0 = hdulist[0].header['coeff0']
        c1 = hdulist[0].header['coeff1']
        npix = hdulist[0].header['naxis1']
        wave = 10.**(c0 + c1 * n.arange(npix))
        spec.beginwl    = c0
        spec.npix       = npix
        spec.wave       = wave
        spec.lambda_eff = hdulist[5].data['LAMBDA_EFF'][i]
        
        # Following commented-out bit was needed for some of the early redux:
        #bzero = hdulist[0].header['bzero']
        bunit = hdulist[0].header['bunit']
        flux = hdulist[0].data
        ivar = hdulist[1].data
        # put Balmer _uncorrected_ flux and errors into the spectrum
        spec.flux       = flux[i,:] * (ivar[i,:] > 0); spec.flux[spec.flux==0] = n.nan
        spec.ivar       = ivar[i,:] * (ivar[i,:] > 0); spec.ivar[spec.ivar==0] = n.nan
        spec.flux_error = 1/n.sqrt(spec.ivar) * (ivar[i,:] > 0)
        # apply Balmer corrections to the spectrum
        if args.balmer == True:
            apply_resid_corr(spec, resid_corr)
        elif args.balmer == False:
            pass
        hdulist.close()
        hdulist = 0
        
        hdulist = pf.open(zbfile)
        metadata = hdulist[0].header
        synflux = hdulist[2].data
        zstruc = hdulist[1].data
        hdulist.close()
        hdulist = 0

        k+=1
        
        spec.alt      = metadata['ALT']
        spec.airtemp  = metadata['AIRTEMP']
        spec.seeing50 = metadata['SEEING50']
        spec.humidity = metadata['HUMIDITY']
        spec.pressure = metadata['PRESSURE']
        spec.az       = metadata['AZ']
        spec.z        = zstruc[i].field('z')
        spec.ra       = zstruc[i].field('plug_ra')
        spec.dec      = zstruc[i].field('plug_dec')

        #for key, val in vars(spec).items():
        #    print '{key} = {val}'.format(key=key, val=val)
        #p.figure()
        # Set starting fiber point (above), then copy and paste
        # the following repeatedly to loop over spectra:
        #i+=1

        # Following commented-out bit was needed for some of the early redux:
        #p.plot(wave, (flux[i,:]-bzero) * (ivar[i,:] > 0), 'k', hold=False)
        #auxflx = flux[i,:] * (ivar[i,:] > 0)
        #p.plot(wave, spec.flux, 'k', hold=False)
        #p.plot(wave, synflux[i,:], 'g', hold=True)
        #p.xlabel('Angstroms')
        #p.ylabel(bunit)
        #p.title(zstruc[i].field('class') + ', z = ' + str(zstruc[i].field('z')))
        #p.ylim(n.percentile(auxflx,1), 1.2*n.percentile(auxflx, 99))
        #p.show()
        spec.__save_fitsio__(dirpath=savedir, filetype='valid')
        #if k>10:
        #    sys.exit()
        #p.savefig('/Volumes/T2TB/Isotropy/offset/python_ancillary/margala/plot.pdf')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--list', type=str, default=None, 
        help = 'full path to the CSV file containing the Plate, MJD, Fiber of objects, required')
    # LISTS are in ../spectra/spAll/
    #parser.add_argument('--save', type=str, default = None,
    #    help = 'path to the save file')
    parser.add_argument('--data', choices=['QSO','FQSO','STAR'], nargs=1, type=str, default='QSO',
        help='data to use')
    parser.add_argument('--balmer', action='store_true', default=False,
        help='apply Balmer corrections')
    args = parser.parse_args()
    print '{0:=^80}'.format(' PROGRAM EXTRACT ')
    # Set topdir: (Folder which contains the spPlate, spZline, spZall files; provided by Daniel Margala & David Kirkby. Proprietary data, do not share.)
    topdir = '/Volumes/T2TB/Isotropy/spectra/SDSS_DR12/Margala/portal.nersc.gov/project/boss/temp/sjbailey/dmargala/reduxtest/v5_7_0/'
    # Set savedir:
    if args.balmer == False:
        savedir = os.path.join('/Volumes/T2TB/Isotropy/spectra/SDSS_DR12/Margala/validation/',args.data[0])
    if args.balmer == True:
        savedir = os.path.join('/Volumes/T2TB/Isotropy/spectra/SDSS_DR12/Margala/validation/',args.data[0]+'_balmer')
    main(args)
