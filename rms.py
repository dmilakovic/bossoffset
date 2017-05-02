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
from congrid import congrid
from corrutilis import *
from SDSSmodules.SDSSclasses import spectrum
from scipy.interpolate import interp1d
import offset
import progressbar

    
def main():
    # ======================= SETTINGS =======================
    settings = program_settings()
    settings.plot = True
    settings.smoothness = 5
    settings.resample = True
    keys = ['BOSS','CORR','MARG']
    # ======================= PARSER =======================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spec-dir', type=str, default=None, 
        help='full path to the directory containing the FITS files, required')
    parser.add_argument('--panels', type=int, default=3,
        help='number of panels, default = 3')
    args = parser.parse_args()

    data, medianSpec = read_spectra(args.spec_dir)
    
    # ======================= PRINT THE SPECTRAL INDEX =======================
    print '{0:>45s}'.format('SPECTRAL INDEX')
    print '{0:>33s} '.format('SI_medspec')

    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('UNCORRECTED:', np.nanmedian(medianSpec.alpha_BOSS), np.nanstd(medianSpec.alpha_BOSS))
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('CORRECTED (Ours):', np.nanmedian(medianSpec.alpha_CORR), np.nanstd(medianSpec.alpha_CORR))
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('CORRECTED (Margala):', np.nanmedian(medianSpec.alpha_MARG), np.nanstd(medianSpec.alpha_MARG))
    
    # ======================= CALCULATE THE VARIANCE =======================
    print 'Calculating the variance'
    # TRANSPOSE THE FLUXES TO CALCULATE RMS
    flux = np.transpose(data['FLUX'])
    
    row = data['FLUX'].dtype
    npix = data['WAVE'][0].size
    flux_rms  = np.zeros(shape = (npix,), dtype = row)
    flux_std  = np.zeros(shape = (npix,), dtype = row)
    for pix in xrange(npix):
        flux_rms['BOSS'][pix] = np.sqrt(np.mean(np.square(flux['BOSS'][pix])))
        flux_std['BOSS'][pix] = np.std(flux['BOSS'][pix])
        flux_rms['CORR'][pix] = np.sqrt(np.mean(np.square(flux['CORR'][pix])))
        flux_std['CORR'][pix] = np.std(flux['CORR'][pix])
        flux_rms['MARG'][pix] = np.sqrt(np.mean(np.square(flux['MARG'][pix])))
        flux_std['MARG'][pix] = np.std(flux['MARG'][pix])

    
    #medianSpec = median_spectrum(z, wave_grid, flux_rms, flux_std)
    print '{0:>50s}'.format('STANDARD DEVIATION OF FLUX')
    print '{0:>31s} {1:>15s}'.format('RMS','STD')
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('UNCORRECTED:',         np.nanmedian(flux_std['BOSS']), np.nanstd(flux_std['BOSS']))
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('CORRECTED (Ours):',    np.nanmedian(flux_std['CORR']), np.nanstd(flux_std['CORR']))
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('CORRECTED (Margala):', np.nanmedian(flux_std['MARG']), np.nanstd(flux_std['MARG']))



    settings.plot=False
    if settings.plot==True:
        npanels = args.panels
    # ======================= PLOT =======================
        majorLocator   = MultipleLocator(500)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(100)
        # PLOT NAME
        thingid = medianSpec.thingid
        if   npanels == 1:
            outfile = os.path.join(args.spec_dir,'%s_rms_single_panel.pdf' %thingid)
            fsize = (12,3)
        elif npanels == 2:
            outfile = os.path.join(args.spec_dir,'%s_rms_double_panel.pdf' %thingid)
            fsize = (12,5)
        elif npanels == 3:
            outfile = os.path.join(args.spec_dir,'%s_rms_triple_panel.pdf' %thingid)
            fsize = (12,7)
        # PLOT PARAMETERS
        cmap = plt.get_cmap('Paired')
        line_colors = cmap(np.linspace(0,1,medianSpec.nspectra))
        # ======================= DEFINING THE CANVAS =======================
        # PLOT SIZE
        fsize = (12,5)
        # MAKE SUBPLOTS
        fig, axs = plt.subplots(npanels,1,figsize=fsize,sharey=True)
        ax      = np.array(axs)
        fig.subplots_adjust(hspace=0.02)
        # SET TITLE
        plt.suptitle('THING ID = %9s' %(thingid))
       # plt.suptitle('RA = %10.6f, DEC = %10.6f , z = %5.3f' %(ra,dec,z))    
        # ======================= AXES LABELS =======================
        #fig.text(0.5,0.0, r'Observed wavelength $[\AA]$', ha='center')
        fig.text(0.08, 0.5, r'Flux $[10^{-17} \rm{erg/s/cm^2/\AA}]$', va='center', rotation='vertical')

        xmins = []
        xmaxs = []
        ymins = []
        ymaxs = []

        # PLOTTING DATA
        box=settings.smoothness
        auxwav = data['WAVE'][0]
        for i, key in zip(xrange(npanels),keys[0:npanels]):
            print 'PLOTTING {0:s}'.format(key)
            # i goes over the uncorrected, corrected, corrected_margala
            k = 0
            # ======================= SMOOTHING =======================
            auxmedian = smooth_array(flux_rms[key],box)
            auxdev  = smooth_array(flux_std[key],box)
            
            for j in xrange(medianSpec.nspectra):
                # j goes over the spectra
                auxflx = smooth_array(data['FLUX'][j][key], box)
            
                if np.size(auxflx)!=np.size(auxwav):
                    print 'error'
                    print auxwav.shape
                    print auxflx.shape
                    sys.exit()
                    #auxwav = auxwav[:np.size(auxflx)]
            
                # get the limits in y-axis
                ymin = np.percentile(auxflx,  0.5); ymins.append(ymin)
                ymax = 1.2*np.percentile(auxflx, 99.5);   ymaxs.append(ymax)             
                ax[i].plot(auxwav,auxflx,color=line_colors[k],linestyle=':', alpha=0.6)
                k = k+1
            # PLOT MEDIAN FLUX AND RMS
            sigma = flux_std[key]
            ax[i].plot(auxwav,auxmedian,color='Green',linewidth=1.3)
            ax[i].plot(auxwav,auxdev,color='Red',linewidth=1)
            # FILL AREA OF RMS
            ax[i].fill_between(auxwav, auxmedian - sigma, auxmedian + sigma, color="Green", alpha=0.3)
            # FILL AREA OF LY A FOREST
            #ax[i].axvspan(wave_grid.min(), 1215.67*(1+z), alpha=0.3, color="SkyBlue")
        xplotlim = (auxwav.min(), auxwav.max())
        yplotlim = (min(ymins),max(ymaxs))
        for i in xrange(npanels):
            ax[i].set_xlim(xplotlim[0],xplotlim[1])
            ax[i].set_ylim(yplotlim[0],yplotlim[1])
        ax[npanels-1].set_xlabel(r'Observed wavelength $[\AA]$')

        #0.08, 0.5, r'Flux $[10^{-17} \rm{erg/s/cm^2/\AA}]$', va='center', rotation='vertical')
        ax[0].text(0.96,0.5,'UNCORRECTED spectra', va='center', rotation='vertical')
        fig.text(6000,0.85*yplotlim[1],'CORRECTED spectra (this paper)')
        fig.text(6000,0.85*yplotlim[1],'CORRECTED spectra (Margala)')

        plt.savefig(outfile, format = 'pdf')
        print 'PLOT SAVED TO {0:s}'.format(outfile)
        plt.show()
        

if __name__=='__main__':
    main()
