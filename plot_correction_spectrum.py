'''
PURPOSE:
--------
        Plots the BOSS, CORR, MARG spectrum of an object in a single panel,
        the correction function C(lambda) in the second panel, and the
        estimated continuum emission in the third panel.
Input:
------
       spec-file : path to the FITS file of the spectrum (fileclass = CORR)
Output:
-------
       Figure
        _____________________________
       |                             |
       |          SPECTRUM           |
       |_____________________________|
        _____________   _____________
       | Correction  | |    Power    |
       |  function   | |     law     |
       |_____________| |_____________| 
       
'''
import fitsio
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import argparse
import os
import sys
from glob import glob
from SDSSutilities import smooth_array
from SDSSmodules.SDSSclasses import spectrum
import corrutilis
from congrid import congrid

def get_limits(array):
    # however, this might be a spike
    arraysize = np.size(array)
    # calculate the average flux of the 4 neighbouring pixels
    avrflx = 0.
    # values of the min and max in the array
    minval = np.min(array) ; maxval = np.max(array)
    # position of the min and max in the array
    minpos = np.argmin(array) ; maxpos = np.argmax(array)
    # make sure not to go over the array boundaries
    lbound = maxpos-2  
    tbound = maxpos+2
    if lbound<0:
        while lbound<0:
            lbound = lbound + 1
    if tbound>=arraysize:
        while tbound>arraysize:
            tbound = tbound - 1
    # calculate the average flux
    for j in xrange(lbound,tbound,1):
        avrflx=float(avrflx)+array[j]
    avrflx = avrflx/5
    # if the current maximum value is more than 2x the average flux, we have a spike
    # define a new maximum value
    if maxval>=2*avrflx:
        maxval = 1.5* avrflx
    # repeat for minimum value (different condition: absolute value is more than 2x the average flux)
    avrflx = 0.
    lbound = minpos-2
    while lbound<0:
        lbound = lbound + 1
    if tbound>arraysize:
        while tbound>arraysize:
            tbound = tbound - 1
    for j in xrange(lbound,tbound,1):
        avrflx=avrflx+array[j]
    avrflx = float(avrflx)/5
    if abs(minval)>=2*abs(avrflx):
        minval = -1.3*abs(avrflx)

    extremes = (minval,maxval)
    return extremes
def main():
    # ======================= PARSER =======================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spec', type=str, default=None, 
        help='full path to the corrected FITS file, required')
    args = parser.parse_args()

    specpath = args.spec
    save_dirname = os.path.dirname(os.path.dirname(specpath))
    sdss_flag = os.path.isdir(os.path.join(save_dirname,'sdss'))
    if sdss_flag:
        sdss_specpath = glob(os.path.join(save_dirname,'sdss','*.fits'))
    #    spec_sdss     = spectrum()
    #    spec_sdss.__read_BOSS__(sdss_specpath[0])
    #    spec_sdss.__fit_powerlaw__('BOSS')
    # ======================= UNIVERSAL WAVELENGTH GRID =======================
    loglam_grid = np.arange(start=3.55, stop = 4.02, step = 0.0001)
    wave_grid   = np.power(10, loglam_grid)
    wave_npix   = wave_grid.size
    wave_shape  = wave_grid.shape
    npix        = wave_npix
    # ======================= ======================== ======================
    keys         = ['SDSS','BOSS','CORR']
    row_flux     = np.dtype([('BOSS','f8',1),('CORR','f8',1),('SDSS','f8',1)])
    row_data     = np.dtype([('WAVE','f8',1),('FLUX',row_flux,1),('FLUX_ERR',row_flux,1),('POWERLAW',row_flux,1), ('CORR',row_flux,1)])
    row_si       = np.dtype([('BOSS','f8',3),('CORR','f8',3),('SDSS','f8',3)])

    data         = np.zeros(shape = (1,wave_npix), dtype = row_data); data.fill(np.nan)
    data['WAVE'] = wave_grid
    si_data      = np.ones(shape = (1,), dtype=row_si)
    for key in keys:
        spec = spectrum()
        if key=='SDSS':
            spec.__read_BOSS__(sdss_specpath[0])
            filetype = 'BOSS'
        else:
            spec.__read_CORR__(specpath)
            filetype = key
        if key == 'CORR':
            spec.flux = spec.flux_corr
        spec.__fit_powerlaw__(filetype)
        # ======================= REBINNING  =======================
        minwav = spec.wave.min() ; maxwav = spec.wave.max()
        # find where in the rebinned wavelength grid to put the rebinned flux values
        index_start = np.argmin(abs(wave_grid-minwav))
        index_stop  = np.argmin(abs(wave_grid-maxwav))
        # the spectra is rebinned to 'newpix' pixels
        newpix = index_stop - index_start
        print key, data['FLUX'][key]
        # put new flux values & the powerlaw into the data structure
        data[0]['FLUX'][key][index_start:index_stop] = congrid.rebin_1d(spec.flux, newpix)
        si_data[0][key][0] = spec.alpha
        si_data[0][key][1] = spec.alpha_error
        si_data[0][key][2] = spec.chisq
     
    print save_dirname
    # READ SPECTRUM FROM FILE & FIT POWERLAWS TO THE UNCORRECTED AND CORRECTED SPECTRA
    
    # PRINT SI DATA
    print '{0:-^80s}'.format(' POWER-LAW FIT ')
    print '{0:>8}{1:^14s}{2:>10s}'.format('','SI', 'CHI2')
    for key in keys:
        print '{key} : {0:+6.4f}+-{1:6.4f}{2:>10.4f}'.format(key = key, *si_data[0][key])
    # CALCULATE THE CHI SQUARE OF BOSS AND CORR DATA (USING SDSS DATA AS EXPECTED VALUES)
    if sdss_flag:
        print '{0:-^80s}'.format(' AGREEMENT WITH SDSS DATA ')
        print '{0:>7s}{1:>10s}'.format('','CHI2')
        for key in keys:
            if key != 'SDSS':
                chisq = np.nansum((data[0]['FLUX'][key] - data[0]['FLUX']['SDSS'])**2 / data[0]['FLUX']['SDSS']) / npix
                print '{key} : {0:>10.4f}'.format(chisq, key = key)
    sys.exit()

    # ======================= BASTI'S SMOOTHING FUNCTION =======================
    basti = True
    if basti:
        box = 5
        flx      = smooth_array(array=spec.flux,smooth_width=box)
        flx_sdss = smooth_array(array=spec_sdss.flux,smooth_width=box)
        flx_corr = smooth_array(array=spec.flux_corr,smooth_width=box)
        

    # ======================= DEFINING THE CANVAS =======================
    fsize = (12,5)
    fig = plt.figure(figsize=fsize)
    plt.suptitle(r'PLATE = %04d, MJD = %5d, FIBRE = %04d,    RA = %6.3f, Dec = %6.3f , z = %5.3f' % (spec.plateid, spec.MJD, spec.fiberid, spec.ra, spec.dec, spec.z))
    fig.subplots_adjust(hspace=0.15)
    
    # ======================= LABELS =======================
    #fig.text(0.5,0.0, r'Observed wavelength $[\AA]$', ha='center')
    #fig.text(0.09, 0.5, r'Flux $[10^{-17} \rm{erg/s/cm^2/\AA}]$', va='center', rotation='vertical')
    # ======================= FIRST SUBPLOT =======================

    majorLocator   = MultipleLocator(500)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(250)

    #### PLOT
    ax2 = plt.subplot(223)
    ax3 = plt.subplot(224)
    ax1 = plt.subplot(211)
    ax1.plot(spec_sdss.wave,flx_sdss,color='orange',label='SDSS',linewidth=1.4)
    ax1.plot(spec.wave,flx,     color='b',label='BOSS',linewidth=1.4)
    ax1.plot(spec.wave,flx_corr,color='g',label='Corrected BOSS',linewidth=1.4)
    
    #### SET LIMITS

    xmin = np.min(spec.wave) ; xmax = np.max(spec.wave) ; xhalf = (xmin+xmax)/2
    ymin = min(0, np.percentile(spec.flux,  1))
    ymax = 1.1*np.percentile(spec.flux, 99.7)
    
    xplotlim = (xmin,xmax)
    yplotlim = (ymin,ymax)
    
    ax1.set_xlim(xplotlim[0],xplotlim[1])
    ax1.set_ylim(yplotlim[0],yplotlim[1])
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_minor_locator(minorLocator)
    fig.text(0.5,0.0,r'Observed wavelength $[\AA]$', ha='center', va='bottom', rotation='horizontal')
    ax1.set_ylabel(r'Flux $[10^{-17} \rm{erg/s/cm^2/\AA}]$')
    
    ax1.legend(loc='upper right',prop={'size':10})
 
    # ======================= SECOND SUBPLOT =======================
    
    majorLocator   = MultipleLocator(2000)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(1000)
    
    #### REDEFINE PLOTTING ARRAYS
    #flxarray             = np.zeros(np.size(flx))
    #flxcorrarray_margala = flx_corr_margala - flx
    #flxcorrarray         = flx_corr - flx

    #### PLOT

    ax2.plot(spec.wave,spec.corr,color='r',label='Correction function',linewidth=2)
    ax2.set_ylabel(r'$C(\lambda)$')
    ax2.set_xlim(xplotlim[0],xplotlim[1])
    #ax2.set_xlabel(r'Observed wavelength $[\AA]$')
    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_minor_locator(minorLocator)
    
    #axes[1].legend(loc='lower right',prop={'size':12})

    # ======================= THIRD SUBPLOT =======================
    
    majorLocator   = MultipleLocator(2000)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(1000)
    
    #### REDEFINE PLOTTING ARRAYS
    #flxarray             = np.zeros(np.size(flx))
    #flxcorrarray_margala = flx_corr_margala - flx
    #flxcorrarray         = flx_corr - flx

    #### PLOT
    ax3.plot(spec_sdss.wave,spec_sdss.powerlaw,color='orange',linestyle='-',linewidth=2.)
    ax3.plot(spec.wave,spec.powerlaw_BOSS,color='b',linestyle='-',linewidth=2.)
    ax3.plot(spec.wave,spec.powerlaw_CORR,color='g',linestyle='-',linewidth=2.)
    colors = ['orange','blue','green']
    alpha     = [spec_sdss.alpha_BOSS, spec.alpha_BOSS, spec.alpha_CORR]
    alpha_err = [spec_sdss.alpha_error_BOSS, spec.alpha_error_BOSS, spec.alpha_error_CORR]

    ymin, ymax = ax3.get_ylim()
    dy   = (ymax-ymin)/10
    for i in xrange(3):
        text = r'$\alpha_{{\nu}} = {0:=+{w}.{p}f}\pm{1:{w}.{p}f} $'.format(alpha[i], alpha_err[i], w=7, p=4)
        ax3.text(7700, 0.8*ymax-dy*i,  text, color=colors[i])
    ax3.set_xlim(xplotlim[0],xplotlim[1])
    #ax3.set_xlabel(r'Observed wavelength $[\AA]$')
    ax3.xaxis.set_major_locator(majorLocator)
    ax3.xaxis.set_minor_locator(minorLocator)
    ax3.set_ylabel(r'Flux $[10^{-17} \rm{erg/s/cm^2/\AA}]$')
    #ax3.legend(loc='upper right',prop={'size':9})
    
    plt.savefig(os.path.join(save_dirname,'%04d-%05d-%04d-corrected_spectrum_2.pdf' %(spec.plateid,spec.MJD,spec.fiberid)), format = 'pdf')
    plt.show()
    
if __name__=='__main__':
    main()
