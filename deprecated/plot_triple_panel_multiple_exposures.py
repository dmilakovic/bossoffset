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


def main():
    # ======================= PARSER =======================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spec-dir', type=str, default=None, 
        help='full path to the directory containing the FITS files, required')
    args = parser.parse_args()

    # ======================= FILES =======================
    # from the input argument, fetch file names
    if args.spec_dir:
        uncorrected_files = glob(os.path.join(args.spec_dir,'uncorrected','*.fits'))
        corrected_files   = glob(os.path.join(args.spec_dir,'corrected','*.fits'))
        corrected_files_m = glob(os.path.join(args.spec_dir,'corrected','margala','*.fits'))
        # list that saves the number of files in each subfolder
        numfiles = [len(uncorrected_files),len(corrected_files),len(corrected_files_m)]
        labels = []
        dir_basename = os.path.dirname(args.spec_dir)
        thingid = dir_basename.split('/')[6]
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
    # ======================= =======================
    wav_all  = []; wav_ucorr = []; wav_corr = []; wav_corr_m = []
    flx_all  = []; flx_ucorr = []; flx_corr = []; flx_corr_m = []
    # ======================= READ DATA FROM FILES =======================
    # uncorrected spectra
    for file in uncorrected_files:
        spec = read_spec(file)
        wav_ucorr.append(spec.wav)
        flx_ucorr.append(spec.flx)
    wav_all.append(wav_ucorr)
    ra = spec.ra ; dec = spec.dec; z = spec.z
    flx_all.append(flx_ucorr)
    # corrected spectra
    for file in corrected_files:
        spec = read_spec(file)
        wav_corr.append(spec.wav)
        flx_corr.append(spec.flx_corr)
    wav_all.append(wav_corr)
    flx_all.append(flx_corr)
    # corrected margala spectra
    for file in corrected_files_m:
        spec = read_spec(file)
        wav_corr_m.append(spec.wav)
        flx_corr_m.append(spec.flx_corr)
    wav_all.append(wav_corr_m)
    flx_all.append(flx_corr_m)
    
    majorLocator   = MultipleLocator(500)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(100)

    # ======================= PLOT =======================
    # PLOT NAME
    outfile = outfile = os.path.join(args.spec_dir,'%s_triple_panel_spectrum.pdf' %thingid)
    # PLOT PARAMETERS
    colors = ['Black','Crimson','Gold','Green','Navy','DarkMagenta','CornflowerBlue','MediumSpringGreen','Sienna','Green','Purple','Teal','Turquoise','FireBrick','Brown','Bisque','Yellow','Coral','Linen']
    # ======================= DEFINING THE CANVAS =======================
    # PLOT SIZE
    fsize = (12,8)
    # MAKE SUBPLOTS
    fig, axes =plt.subplots(3,1,figsize=fsize,sharey=True)
    fig.subplots_adjust(hspace=0.02)
    ax1, ax2, ax3 = axes.flat
    # SET TITLE
    plt.suptitle('RA = %10.6f, DEC = %10.6f , z = %5.3f' %(ra,dec,z))    
    # ======================= AXES LABELS =======================
    #fig.text(0.5,0.0, r'Observed wavelength $[\AA]$', ha='center')
    fig.text(0.09, 0.5, r'Flux $[10^{-17} \rm{erg/s/cm^2/\AA}]$', va='center', rotation='vertical')
    
    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []

    # PLOTTING DATA
    for i in xrange(np.size(numfiles)):
        k = 1
        for j in xrange(numfiles[i]):
            auxwav = np.array(wav_all[i][j])
            auxflx = np.array(flx_all[i][j])
            # ======================= SMOOTHING =======================
            box=5
            auxflx = smooth_array(auxflx,box)
            if np.size(auxflx)!=np.size(auxwav):
                auxwav = auxwav[:np.size(auxflx)]
            #print np.size(auxwav), np.size(auxflx)
            # get the limits in y-axis

            xmin = np.min(auxwav) ; xmax = np.max(auxwav)
            xmins.append(xmin)    ; xmaxs.append(xmax)    
        
            ymin = np.percentile(auxflx,  0.5); ymins.append(ymin)
            ymax = 1.2*np.percentile(auxflx, 99);   ymaxs.append(ymax)

            if i == 0:
                ax1.plot(auxwav,auxflx,color=colors[k],label=labels[i][j])
            elif i == 1:
                ax2.plot(auxwav,auxflx,color=colors[k],label=labels[i][j])
            elif i == 2:
                ax3.plot(auxwav,auxflx,color=colors[k],label=labels[i][j])
            k = k+1
        
    xplotlim = (min(xmins),max(xmaxs),0.5*(min(xmins)+max(xmaxs)))
    yplotlim = (min(ymins),max(ymaxs))
    # AXES FINE TUNING

    ax1.text(6000,0.85*yplotlim[1],'Uncorrected spectra')
    ax2.text(6000,0.85*yplotlim[1],'Corrected spectra (this paper)')
    ax3.text(6000,0.85*yplotlim[1],'Corrected spectra (Margala)')
    
    ax1.set_xlim(xplotlim[0],xplotlim[1])
    ax2.set_xlim(xplotlim[0],xplotlim[1])
    ax3.set_xlim(xplotlim[0],xplotlim[1])

    ax1.set_ylim(yplotlim[0],yplotlim[1])
    ax3.set_xlabel(r'Observed wavelength $[\AA]$')

    # LEGEND
    ax1.legend(loc='upper right',prop={'size':10})

    plt.savefig(outfile, format = 'pdf')
    plt.show()

if __name__=='__main__':
    main()
