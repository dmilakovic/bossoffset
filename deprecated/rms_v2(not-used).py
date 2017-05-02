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
    # ======================= PARSER =======================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spec-dir', type=str, default=None, 
        help='full path to the directory containing the FITS files, required')
    parser.add_argument('--panels', type=int, default=3,
        help='number of panels, default = 3')
    args = parser.parse_args()

    # ======================= FILES =======================
    # from the input argument, fetch file names
    if args.spec_dir:
        spectral_list       = return_spectral_list(args.spec_dir)
        numspec             = len(spectral_list[0])
        uncorrected_files   = spectral_list[0]
        corrected_files     = spectral_list[1]
        corrected_files_m   = spectral_list[2]
        # list that saves the number of files in each subfolder
        numfiles = [len(uncorrected_files),len(corrected_files),len(corrected_files_m)]
        labels = []
        dir_basename = os.path.dirname(args.spec_dir)
        thingid = dir_basename.split('/')[7]
        # read the data from the uncorrected FITS files
        ids = []
        for file in uncorrected_files:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            ids.append('%04d-%5d-%04d' % (plate, mjd, fiberid))
        labels.append(ids)
        # read the data from the corrected FITS files
        ids = []
        for file in corrected_files:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            ids.append('%04d-%5d-%04d' % (plate, mjd, fiberid))
        labels.append(ids)
        # read the data from the corrected FITS files (Margala)
        ids = []
        for file in corrected_files_m:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            ids.append('%04d-%5d-%04d' % (plate, mjd, fiberid))
        labels.append(ids)
        
    #data, header = fitsio.read(file,ext=0,header=True)
    # ======================= UNIVERSAL WAVELENGTH GRID =======================
    loglam_grid = np.arange(start=3.55, stop = 4.02, step = 0.0001)
    wave_grid   = np.power(10, loglam_grid)
    wave_npix   = wave_grid.size
    wave_shape  = wave_grid.shape
    # ======================= ======================== ======================
    flux_boss_list = []; flux_corr_list = []; flux_marg_list = []
    keys         = ['BOSS','CORR','MARG']
    row          = np.dtype([('BOSS','f8',1),('CORR','f8',1),('MARG','f8',1)])

    flux_data    = np.zeros(shape = (max(numfiles),wave_npix), dtype=row)
    si_data      = np.zeros(shape = (max(numfiles)), dtype=row)
    pwrlaw_data  = np.zeros(shape = (max(numfiles),wave_npix), dtype=row)
    # ======================= READ DATA FROM FILES =======================
    for i in xrange(len(spectral_list)):
        for j,file in enumerate(spectral_list[i]):
            # ======================= READ SPECTRA =======================
            spec = spectrum()
            spec_dirname  = os.path.dirname(file)
            spec_basename = os.path.basename(file)
            if 'corr' not in spec_basename:
                spec.__read_BOSS__(file)
                flux_boss = np.zeros(shape=wave_shape)
            elif 'corr' in spec_basename:
                spec.__read_CORR__(file)
                flux_corr = np.zeros(shape=wave_shape)
            # ======================= REBINNING  =======================
            minwav = spec.wave.min() ; maxwav = spec.wave.max()
            index_start = np.argmin(abs(wave_grid-minwav))
            index_stop  = np.argmin(abs(wave_grid-maxwav))
            newpix = index_stop - index_start
            
            if 'corr' not in spec_basename:
                spec.__fit_powerlaw__(type='BOSS')
                si_data[j]['BOSS'] = spec.alpha
                flux_data[j]['BOSS'][index_start:index_stop] = congrid.rebin_1d(spec.flux, newpix)
                pwrlaw_data[j]['BOSS'][index_start:index_stop] = congrid.rebin_1d(spec.powerlaw, newpix)
            if 'corr' in spec_basename:
                spec.__fit_powerlaw__(type='CORR')
                if 'margala' not in spec_dirname:
                    si_data[j]['CORR'] = spec.alpha
                    flux_data[j]['CORR'][index_start:index_stop] = congrid.rebin_1d(spec.flux_corr, newpix)
                    pwrlaw_data[j]['CORR'][index_start:index_stop] = congrid.rebin_1d(spec.powerlaw, newpix)
                elif 'margala' in spec_dirname:
                    si_data[j]['MARG'] = spec.alpha
                    flux_data[j]['MARG'][index_start:index_stop] = congrid.rebin_1d(spec.flux_corr, newpix)
                    pwrlaw_data[j]['MARG'][index_start:index_stop] = congrid.rebin_1d(spec.powerlaw, newpix)
    ra = spec.ra ; dec = spec.dec; z = spec.z
    del(spec)
    # TRANSPOSE THE FLUXES TO CALCULATE RMS 
    flux = np.transpose(flux_data)
    
    # ======================= CALCULATE THE VARIANCE =======================
    print 'Calculating the variance'

    
    flux_rms  = np.zeros(shape = wave_shape, dtype = row)
    flux_std  = np.zeros(shape = wave_shape, dtype = row)

    for pix in xrange(wave_npix):
        flux_rms[pix]['BOSS'] = np.sqrt(np.mean(np.square(flux['BOSS'][pix]))); flux_std[pix]['BOSS'] = np.std(flux['BOSS'][pix])
        flux_rms[pix]['CORR'] = np.sqrt(np.mean(np.square(flux['CORR'][pix]))); flux_std[pix]['CORR'] = np.std(flux['CORR'][pix])
        flux_rms[pix]['MARG'] = np.sqrt(np.mean(np.square(flux['MARG'][pix]))); flux_std[pix]['MARG'] = np.std(flux['MARG'][pix])

    medianSpec = median_spectrum(z, wave_grid, flux_rms, flux_std)
    print '{0:>50s}'.format('STANDARD DEVIATION OF FLUX')
    print '{0:>31s} {1:>15s}'.format('RMS','STD')
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('UNCORRECTED:',         np.nanmedian(flux_rms['BOSS']), np.nanstd(flux_rms['BOSS']))
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('CORRECTED (Ours):',    np.nanmedian(flux_rms['CORR']), np.nanstd(flux_rms['CORR']))
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('CORRECTED (Margala):', np.nanmedian(flux_rms['MARG']), np.nanstd(flux_rms['MARG']))



    settings.plot=True
    if settings.plot==True:
        npanels = args.panels
    # ======================= PLOT =======================
        majorLocator   = MultipleLocator(500)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(100)
        # PLOT NAME
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
        line_colors = cmap(np.linspace(0,1,len(uncorrected_files)))
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
        auxwav = wave_grid
        for i, key in zip(xrange(npanels),keys[0:npanels]):
            print 'PLOTTING {0:s}'.format(key)
            # i goes over the uncorrected, corrected, corrected_margala
            k = 0
            # ======================= SMOOTHING =======================
            auxmedian = smooth_array(flux_rms[key],box)
            auxdev  = smooth_array(flux_std[key],box)
            
            for j in xrange(numfiles[i]):
                # j goes over the spectra
                auxflx = smooth_array(flux_data[j][key], box)
            
                if np.size(auxflx)!=np.size(auxwav):
                    print 'error'
                    print auxwav.shape
                    print auxflx.shape
                    sys.exit()
                    #auxwav = auxwav[:np.size(auxflx)]
            
                # get the limits in y-axis
                ymin = np.percentile(auxflx,  0.5); ymins.append(ymin)
                ymax = 1.2*np.percentile(auxflx, 99.5);   ymaxs.append(ymax)             
                ax[i].plot(auxwav,auxflx,color=line_colors[k],label=labels[i][j],linestyle=':', alpha=0.6)
                k = k+1
            # PLOT MEDIAN FLUX AND RMS
            sigma = flux_std[key]
            ax[i].plot(auxwav,auxmedian,color='Green',linewidth=1.3)
            ax[i].plot(auxwav,auxdev,color='Red',linewidth=1)
            # FILL AREA OF RMS
            ax[i].fill_between(auxwav, auxmedian - sigma, auxmedian + sigma, color="Green", alpha=0.3)
            # FILL AREA OF LY A FOREST
            #ax[i].axvspan(wave_grid.min(), 1215.67*(1+z), alpha=0.3, color="SkyBlue")
        xplotlim = (wave_grid.min(), wave_grid.max())
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
