#!/usr/bin/env python

# This program plots the spectrum of a FITS file, fits the continuum
# using two different fitting algorithms
# 1: using medians of fitting intervals
# 2: using all data points

from copy import deepcopy
import sys
#from SDSSmodules.SDSSfiles import check_filetype
#from SDSSmodules.SDSSclasses import spectrum
from SDSSmodules.SDSSfiles import *
from SDSSmodules.SDSSfitting import *
from SDSSmodules.SDSScorrections import *
from SDSSmodules.SDSSoutput import *
from SDSSmodules.SDSScspec import *
from SDSSmodules.SDSSutilities import *
from SDSSmodules.SDSSclasses import *
import SDSSmodules.SDSSpaths as path

import numpy as np
#from astropy.io import FITS
import fitsio
from pylab import *
from matplotlib.backend_bases import Event
from matplotlib import pyplot as plt

import argparse
import os
from glob import glob

# class which contains all necessary functions to work on the matplotlib graph
class WorkOnSpectrum:
    def __init__(self, filelist, filetype, dustmap, settings, resid_corr):
        self.filelist  = filelist
        self.filetype  = filetype
        self.dustmap   = dustmap
        self.nspec     = len(filelist)
        self.i         = 0
        self.settings  = settings
        self.resid_corr= resid_corr
    
    def start(self):
        #self.cidpress = self.fig.canvas.mpl_connect('key_press_event', self.press)
        while self.i < self.nspec:
            self.work_on_spectrum(self.filelist[self.i])
            self.i += 1

    def work_on_spectrum(self, filename):
        np.set_printoptions(threshold=6000)
        # create spectrum object and determine properties of spectrum
        global spec
        spec = spectrum()
        
        z_min = self.settings.z_min
        z_delta = self.settings.z_delta
        if self.filetype == 1:
            read_spSpec_fitsio(self.filelist[self.i].rstrip(), spec, None)
        if self.filetype == 2:
            spec.__read_BOSS__(self.filelist[self.i])
            #read_spec_fitsio(self.filelist[self.i].rstrip(), spec, None)
        if self.filetype == 3:
            read_speclya_fitsio(self.filelist[self.i].rstrip(), spec, None)

        print '================='
        spec.Ebv = obstools.get_SFD_dust(spec.coordinates.galactic.l.deg, spec.coordinates.galactic.b.deg, self.dustmap, interpolate=0)
        print filename
        spec.filename = filename.rstrip()
        #Gal_extinction_correction(spec)
        # perform the fit. 1st with medians 2nd with individual points

        spec_median = deepcopy(spec)

        emfree = np.array([1280.0,1312.0,1345.0,1440.0,1610])
        emfree_end = np.array([1292.0,1328.0,1365.0,1475.0,1790])
        spec.emfree = (1.0 + spec.z)*emfree
        spec.emfree_end = (1.0 + spec.z)*emfree_end

        zem_index = calc_zem_index(spec.z, z_min, z_delta)
        if zem_index == -1:
            zem_index = 0
        x_ind, y_ind, y_err_ind = fit_powerlaw_individual(spec, self.settings, 1, zem_index = zem_index, deviation_factor=5.0)
        # COORDINATES & METADATA
        #spec.ra  = 0.0
        #spec.dec = 0.0
        size = np.size(spec_median.powerlaw)
#        print np.shape(spec.wave), np.shape(spec_median.powerlaw[0:size])
        # ly alpha:
        ind = np.where(spec.wave/(1.0+spec.z) < 1200)[0]
        #flux_lya   = smooth_array(spec.flux[ind], 100, spec.flux_error[ind])
        #        flux_lya   = smooth_array(flux_lya[ind], 5, spec.flux_error[ind])

        # ADD ZENITH POSITION

        # spec.zenith   = input("Zenith = ")

        
        flux_corr = perform_flux_correction_adaptive(spec, self.settings, self.resid_corr)
        if flux_corr == True or flux_corr == False:
            print 'Spectrum unusable'
            return
        corr_4000 = flux_corr[0]
        corr_5400 = flux_corr[1]
        flux_corr = corr_5400 / corr_4000
        flux_corrected = spec.flux * flux_corr

        #######################################
        ### ADD THE CORRECTIONS TO THE SPECTRUM OBJECT
        #######################################
        spec.flux_corr  = flux_corrected
        spec.corr4000   = corr_4000
        spec.corr5400   = corr_5400
        spec.corr       = flux_corr
        spec.ivar_corr  = abs(spec.ivar/(flux_corr)**2)
        spec.flux_error = 1./np.sqrt(spec.ivar_corr)
        spec.dirname    = os.path.dirname(filename)
        spec.npix       = np.size(spec.flux)


        ### SAVE SPECTRUM TO FILE
        #save_spectrum(spec,settings = self.settings)
        dirpath = os.path.join(parser_args.spec_dir,'corrected')
        spec.__save_fitsio__(dirpath)


def main(args):

    if args.spec_list:
        if len(args) > 0:
            filelist = open(args[1], 'r').readlines()
    elif args.spec_file:
        filelist=[parser_args.spec_file]
        
    elif args.spec_dir:
        filelist = glob(os.path.join(parser_args.spec_dir,'uncorrected','*.fits'))
        numfiles = len(filelist)

    print filelist
    # determine if it's a DR7 or DR10 file
    filetype = check_filetype(filelist[0])

    dustmap = path.to_dustmaps()

    #print 'This program plots the spectra from a list of FITS files'
    #print 'Red:   spectrum'
    #print 'Blue:  error of spectrum'
    #print 'Pink:  fit_powerlaw_individual (based on individual pixels in intervals'
    #print 'Green: fit_powerlaw            (based on medians in intervals)'
    #print 'Cyan:  data points used for fit_powerlaw_individual'
    #print 'Black: data points used for fit_powerlaw'
    #print ''
    #print 'Pressing N (NEXT)   displays the next spectrum'
    #print 'Pressing B (BOTTOM) displays the last spectrum'
    #print 'Pressing P (PRINT)  save the current spectrum into a FITS file'

    # number of files in the list
    file_num = len(filelist)
    settings = program_settings()
    # create WorkOnSpectrum object
    if args.mock:
        resid_corr = None
    else:
        resid_corr = read_resid_corr(path.to_resid_corr())
    spectra = WorkOnSpectrum(filelist, filetype, dustmap, settings, resid_corr)
    spectra.start()

if __name__ == "__main__":
    import sys
    global parser_args
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spec-file', type=str, default=None,
        help='full path of a spec file')
    parser.add_argument('--spec-dir', type=str, default=None,
        help='path to directory containing individual spec files')
    parser.add_argument('--spec-list', type=str, default=None,
        help='path to a TXT containing paths to individual spec files')
    parser.add_argument('--mock', action='store_true',
        help='create mock spectra with uncorrected fluxes all equal to one')
    parser_args = parser.parse_args()
    main(parser_args)

