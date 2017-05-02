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
from corrutilis import *
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
    args = parser.parse_args()

    # ======================= FILES =======================
    # from the input argument, fetch file names
    if args.spec_dir:
        spectral_list       = return_spectral_list(args.spec_dir)
        uncorrected_files   = spectral_list[0]
        corrected_files     = spectral_list[1]
        corrected_files_m   = spectral_list[2]
        # list that saves the number of files in each subfolder
        numfiles = [len(uncorrected_files),len(corrected_files)]
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
        
    #data, header = fitsio.read(file,ext=0,header=True)
    # ======================= ======================== ======================
    wav_all  = []; wav_ucorr  = []; wav_corr  = []; wav_corr_m  = []
    flx_all  = []; flx_ucorr  = []; flx_corr  = []; flx_corr_m  = []
    npix_all = []; npix_ucorr = []; npix_corr = []; npix_corr_m = []
    loglam_start_all = []; loglam_start_ucorr = []; loglam_start_corr = []; loglam_start_corr_m = []
    loglam_end_all = []; loglam_end_ucorr = []; loglam_end_corr = []; loglam_end_corr_m = []
    # ======================= READ DATA FROM FILES =======================
    # uncorrected spectra
    npix = []
    for file in uncorrected_files:
        spec = read_spec(file)
        wav_ucorr.append(spec.wav)
        flx_ucorr.append(spec.flx)
        npix_ucorr.append(spec.npix)
        loglam_start_ucorr.append(spec.loglam_start)
        loglam_end_ucorr.append(spec.loglam_end)
    ra = spec.ra ; dec = spec.dec; z = spec.z
    wav_all.append(wav_ucorr)
    flx_all.append(flx_ucorr)
    npix_all.append(npix_ucorr)
    loglam_start_all.append(loglam_start_ucorr)
    loglam_end_all.append(loglam_end_ucorr)
    del spec
    # corrected spectra
    for file in corrected_files:
        spec = read_spec(file)
        wav_corr.append(spec.wav)
        #flx_corr.append(spec.corr_flux)
        flx_corr.append(spec.flx * spec.corr5400 / spec.corr4000)
        npix_corr.append(spec.npix)
        loglam_start_corr.append(spec.loglam_start)
        loglam_end_corr.append(spec.loglam_end)
        del spec
    wav_all.append(wav_corr)
    flx_all.append(flx_corr)
    npix_all.append(npix_corr)
    loglam_start_all.append(loglam_start_corr)
    loglam_end_all.append(loglam_end_corr)
    
    # ======================= RESAMPLING THE FLUX ARRAYS =======================
    # Each spectrum has a different number of pixels, 'starting' and 'end' wavelengths. We need to resample them a common pixel grid.
    if settings.resample==True:
        flx_all_resampled = []
        # reshaped = transposed
        flx_all_resampled_reshaped = []
        # create a resampled array in wavelengths: take minimum of 'starting' and maximum of 'end' wavelengths
        loglam_start = np.min(np.array(loglam_start_all))
        loglam_end   = np.max(np.array(loglam_end_all))
        loglam_step  = 0.0001
        loglam_resampled = np.arange(loglam_start,loglam_end,loglam_step)
        wav_resampled    = ma.power(10, loglam_resampled)
        ## # RESAMPLING THE FLUXES
        for i in xrange(np.size(numfiles)):
        # i goes over 2 indices: uncorrected, corrected
            flx_kind = []
            for j in xrange(numfiles[i]):
             # j goes over the spectra
                # locate the minimum and maximum of the wavelength array for this spectrum
                minwav = np.min(wav_all[i][j])
                maxwav = np.max(wav_all[i][j])
                # define auxiliarry arrays
                index = np.where((wav_resampled>=minwav)&(wav_resampled<=maxwav))[0]

                # find the interpolation function from data in the FITS file
                intfunc = interp1d(wav_all[i][j],flx_all[i][j])
                # apply the interpolation (but only within the proper boundaries)
                # defining a masked array flx_resampled with the same mask as the wavelength array, filled with NaNs
                flx_resampled = ma.array(np.full(wav_resampled.size,np.nan))#,mask = wav_mask)
                # changing the values of the flx_resampled array, but only for indices that are not masked
                flx_resampled[index] = intfunc(wav_resampled[index])

                # define a mask in wavelengths not covered by the j-th spectrum
                wav_mask = ma.masked_outside(wav_resampled,minwav,maxwav).mask
                # define a mask of outliers (further than nsigma)
                nsigma = 8
                outliers_boolean = is_outlier(flx_resampled,nsigma)
                outliers_mask    = ma.make_mask(outliers_boolean)

                # combine the masks into a new mask and attach it to the flx_resampled
                mask = ma.mask_or(wav_mask,outliers_mask)
                flx_resampled.mask = mask

                #for k in xrange(np.size(wav_resampled)):
                 #   print i,j,k, wav_resampled[k], flx_resampled[k]

                flx_kind.append(flx_resampled)
                
            # append the arrays
            flx_kind = ma.array(flx_kind)
            flx_all_resampled.append(flx_kind)
            flx_kind = ma.array(flx_kind).T
            flx_all_resampled_reshaped.append(flx_kind)

    # ======================= CALCULATE THE VARIANCE =======================
    flx_resampled_array = ma.array(flx_all_resampled_reshaped)
    print 'Calculating the variance'

    medFlux = []
    devFlux = []
    for i in xrange(np.size(numfiles)):
        # i goes over the uncorrected, corrected
        median_i = []
        dev_i    = []
        for j in xrange(np.size(flx_resampled_array[i])/numfiles[i]):
            # j goes over the pixels
            fluxes = []
            for k in xrange(numfiles[i]):
                # k goes over individual exposures (spectra)
                fluxes.append(flx_resampled_array[i][j][k])
            median = np.nanmedian(fluxes)
            dev    = np.nanstd(fluxes)
            median_i.append(median)
            dev_i.append(dev)
        medFlux.append(median_i)
        devFlux.append(dev_i)


    medFluxArray = np.array(medFlux)
    devFluxArray = np.array(devFlux)
    print '{0:>50s}'.format('STANDARD DEVIATION OF FLUX')
    print '{0:>31s} {1:>15s}'.format('MEDIAN','STD')
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('UNCORRECTED:',         np.nanmedian(devFluxArray[0]), np.nanstd(devFluxArray[0]))
    print '{0:>20s} {1:15.12f} {2:15.12f}'.format('CORRECTED (Ours):',    np.nanmedian(devFluxArray[1]), np.nanstd(devFluxArray[1]))

    # ======================= RMS IN WAVELENGTH BINS =======================
    # BINS
    LyA_range = [np.nanmin(wav_resampled),1215.67*(1+z[0])]
    range1    = [(1+z[0])*1650,(1+z[0])*1750]
    range2    = [(1+z[0])*2950,(1+z[0])*3050]
    range3    = [(1+z[0])*5100,(1+z[0])*5250]
    print range1,range2,range3
    # ARRAY SLICES
    indexLyA = np.where((wav_resampled >= LyA_range[0]) & (wav_resampled <= LyA_range[1]))[0]
    index1   = np.where((wav_resampled >= range1[0]) & (wav_resampled <=range1[1]))[0]
    index2   = np.where((wav_resampled >= range2[0]) & (wav_resampled <=range2[1]))[0]
    index3   = np.where((wav_resampled >= range3[0]) & (wav_resampled <=range3[1]))[0]

    # RMS IN LYMAN ALPHA FOREST
    medFluxArray_LyA = medFluxArray[:,indexLyA]
    devFluxArray_LyA = devFluxArray[:,indexLyA]
    print "LyA forest = ",LyA_range
    if LyA_range[0]<LyA_range[1]:
        print np.shape(medFluxArray_LyA), np.shape(devFluxArray_LyA)

        print '{0:>50s}'.format('STANDARD DEVIATION OF FLUX IN THE LYMAN ALPHA FOREST')
        print '{0:>31s} {1:>15s} {2:>15s} {3:>15s}'.format('MEDIAN','STD', 'MIN', 'MAX')
        print '{0:>20s} {1:15.12f} {2:15.12f} {3:15.12f} {4:15.12f}'.format('UNCORRECTED:',         np.nanmedian(devFluxArray_LyA[0]), np.nanstd(devFluxArray_LyA[0]), np.nanmin(devFluxArray_LyA[0]), np.nanmax(devFluxArray_LyA[0]))
        print '{0:>20s} {1:15.12f} {2:15.12f} {3:15.12f} {4:15.12f}'.format('CORRECTED (Ours):',    np.nanmedian(devFluxArray_LyA[1]), np.nanstd(devFluxArray_LyA[1]), np.nanmin(devFluxArray_LyA[1]), np.nanmax(devFluxArray_LyA[1]))

    # RMS IN RANGE 1  
    medFluxArray_r1 = medFluxArray[:,index1]
    devFluxArray_r1 = devFluxArray[:,index1]

    # RMS IN RANGE 2
    medFluxArray_r2 = medFluxArray[:,index2]
    devFluxArray_r2 = devFluxArray[:,index2]

    # RMS IN RANGE 3
    medFluxArray_r3 = medFluxArray[:,index3]
    devFluxArray_r3 = devFluxArray[:,index3]
    settings.plot=True
    if settings.plot==True:
    # ======================= PLOT =======================
        majorLocator   = MultipleLocator(500)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(100)
        # PLOT NAME
        outfile = os.path.join(args.spec_dir,'%s_rms_double_panel.pdf' %thingid)
        # PLOT PARAMETERS
        cmap = plt.cm.Spectral
        line_colors = cmap(np.linspace(0,1,len(uncorrected_files)))
        #colors = ['Black','Crimson','Gold','Green','Navy','DarkMagenta','CornflowerBlue','Sienna','Green','Purple','Teal','MediumSpringGreen','Turquoise','FireBrick','Brown','Bisque','Yellow','Coral','Linen']
        # ======================= DEFINING THE CANVAS =======================
        # PLOT SIZE
        fsize = (12,8)
        # MAKE SUBPLOTS
        fig, axes =plt.subplots(2,1,figsize=fsize,sharey=True)
        fig.subplots_adjust(hspace=0.02)
        ax1, ax2 = axes.flat
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
        for i in xrange(np.size(numfiles)):
            # i goes over the uncorrected, corrected
            k = 0
            auxmedian = smooth_array(np.array(medFlux[i]),box)
            auxdev  = smooth_array(np.array(devFlux[i]),box)
            
            for j in xrange(numfiles[i]):
                # j goes over the spectra
                print i,j,k, labels[i][j]
                auxwav = ma.array(wav_resampled,mask = flx_all_resampled[i][j].mask)
                auxflx = flx_all_resampled[i][j]
            
                # ======================= SMOOTHING =======================

                auxflx = smooth_array(auxflx,box)
                #if np.size(auxflx)!=np.size(auxwav):
                    
                    #auxwav = auxwav[:np.size(auxflx)]
            
                # get the limits in y-axis

                xmin = np.min(auxwav) ; xmax = np.max(auxwav)
                xmins.append(xmin)    ; xmaxs.append(xmax)    

                ymin = np.percentile(auxflx,  0.5); ymins.append(ymin)
                ymax = 1.2*np.percentile(auxflx, 99);   ymaxs.append(ymax)

                if i == 0:
                    ax1.plot(auxwav,auxflx,color=line_colors[k],label=labels[i][j],linestyle=':', alpha=0.3)
                elif i == 1:
                    ax2.plot(auxwav,auxflx,color=line_colors[k],label=labels[i][j],linestyle=':', alpha=0.3)
                k = k+1
            if i == 0:
                sigma = ma.array(devFluxArray[0],mask = flx_all_resampled[i][j].mask)
                # PLOT MEDIAN FLUX
                ax1.plot(auxwav,auxmedian,color='Green',linewidth=1.3)
                # FILL AREA OF RMS (along y-axis)
                ax1.fill_between(auxwav, auxmedian - sigma, auxmedian + sigma, color="Green", alpha=0.5)
                # PLOT RMS
                ax1.plot(auxwav,auxdev,color='Red',linewidth=1)
                # FILL WAVELENGTH BINS AREA (along x-axis)
                ax1.axvspan(LyA_range[0], LyA_range[1], alpha=0.2, color="SkyBlue")
                ax1.axvspan(range1[0], range1[1], alpha=0.2, color="Crimson")
                ax1.axvspan(range2[0], range2[1], alpha=0.2, color="Gold")
                ax1.axvspan(range3[0], range3[1], alpha=0.2, color="Lavender")
            elif i == 1:
                sigma = ma.array(devFluxArray[1],mask = flx_all_resampled[i][j].mask)
                ax2.plot(auxwav,auxmedian,color='Green',linewidth=1.3)
                ax2.fill_between(auxwav, auxmedian - sigma, auxmedian + sigma, color="Green", alpha=0.5)
                ax2.plot(auxwav,auxdev,color='Red',linewidth=1)
                ax2.axvspan(LyA_range[0], LyA_range[1], alpha=0.2, color="SkyBlue")
                ax2.axvspan(range1[0], range1[1], alpha=0.2, color="Crimson")
                ax2.axvspan(range2[0], range2[1], alpha=0.2, color="Gold")
                ax2.axvspan(range3[0], range3[1], alpha=0.2, color="Lavender")

        xplotlim = (min(xmins),max(xmaxs),0.5*(min(xmins)+max(xmaxs)))
        yplotlim = (0.0,max(ymaxs))

        # TEXT IN PANELS

        xpos = 0.8*xplotlim[1]
        ax1.text(xpos,0.85*yplotlim[1],'UNCORRECTED spectra')
        ax1.text(xpos,0.75*yplotlim[1],'{0:>6s}{1:5.2f}{2:2s}{3:5.2f}'.format(r'$\sigma_{total}$ = ',np.nanmedian(devFluxArray[0]),'$\pm$', np.nanstd(devFluxArray[0])))
        if LyA_range[0]<LyA_range[1]:
            ax1.text(xpos,0.70*yplotlim[1],'{0:>5s}{1:5.2f} {2:2s} {3:5.2f}'.format('$\sigma_{Ly}$ = ',np.nanmedian(devFluxArray_LyA[0]),'$\pm$', np.nanstd(devFluxArray_LyA[0])))
        ax1.text(xpos,0.65*yplotlim[1],'{0:>6s}{1:5.2f}{2:2s}{3:5.2f}'.format(r'$\sigma_{1650-1750}$ = ',np.nanmedian(devFluxArray_r1[0]),'$\pm$',np.nanstd(devFluxArray_r1[0])))
        #ax1.text(xpos,0.60*yplotlim[1],'{0:>6s}{1:5.2f}{2:2s}{3:5.2f}'.format(r'$\sigma_{2950-3050}$ = ',np.nanmedian(devFluxArray_r2[0]),'$\pm$',np.nanstd(devFluxArray_r2[0])))
        #ax1.text(xpos,0.55*yplotlim[1],'{0:>5s}{1:5.2f}{2:2s}{3:5.2f}'.format(r'$\sigma_{5100-5250}$ = ',np.nanmedian(devFluxArray_r3[0]),'$\pm$',np.nanstd(devFluxArray_r3[0])))
        
        ax2.text(xpos,0.85*yplotlim[1],'CORRECTED spectra (this paper)')
        ax2.text(xpos,0.75*yplotlim[1],'{0:>5s}{1:5.2f}{2:2s}{3:5.2f}'.format(r'$\sigma_{total}$ = ',np.nanmedian(devFluxArray[1]),'$\pm$', np.nanstd(devFluxArray[1])))
        if LyA_range[0]<LyA_range[1]:
            ax2.text(xpos,0.70*yplotlim[1],'{0:>5s}{1:5.2f} {2:2s} {3:5.2f}'.format('$\sigma_{Ly}$ = ',np.nanmedian(devFluxArray_LyA[1]),'$\pm$',np.nanstd(devFluxArray_LyA[1])))
        ax2.text(xpos,0.65*yplotlim[1],'{0:>6s}{1:5.2f}{2:2s}{3:5.2f}'.format(r'$\sigma_{1650-1750}$ = ',np.nanmedian(devFluxArray_r1[1]),'$\pm$',np.nanstd(devFluxArray_r1[1])))
        #ax2.text(xpos,0.60*yplotlim[1],'{0:>6s}{1:5.2f}{2:2s}{3:5.2f}'.format(r'$\sigma_{2950-3050}$ = ',np.nanmedian(devFluxArray_r2[1]),'$\pm$',np.nanstd(devFluxArray_r2[1])))
        #ax2.text(xpos,0.55*yplotlim[1],'{0:>5s}{1:5.2f}{2:2s}{3:5.2f}'.format(r'$\sigma_{5100-5250}$ = ',np.nanmedian(devFluxArray_r3[1]),'$\pm$',np.nanstd(devFluxArray_r3[1])))
        

        # AXES LIMITS
        ax1.set_xlim(xplotlim[0],xplotlim[1])
        ax2.set_xlim(xplotlim[0],xplotlim[1])

        ax1.set_ylim(yplotlim[0],yplotlim[1])
        ax2.set_xlabel(r'Observed wavelength $[\AA]$')

        # LEGEND
        #ax1.legend(loc='upper right',prop={'size':9})

        plt.savefig(outfile, format = 'pdf')
        plt.show()

if __name__=='__main__':
    main()
