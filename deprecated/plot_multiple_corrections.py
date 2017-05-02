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

def main():
    # ======================= PARSER =======================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spec-dir', type=str, default=None, 
        help='full path to the directory containing the FITS files, required')
    args = parser.parse_args()

    # ======================= FILES =======================
    if args.spec_dir:
        spec_allfiles = glob(os.path.join(args.spec_dir,'*.fits'))
        numfiles = len(spec_allfiles)
        labels = []
        dir_basename = os.path.dirname(args.spec_dir)
        thingid = dir_basename.split('/')[6]
        for file in spec_allfiles:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            #spec_filename.append(os.path.join(args.dir, 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)))
            labels.append('%04d-%5d-%04d' % (plate, mjd, fiberid))
    #data, header = fitsio.read(file,ext=0,header=True)
    # ======================= =======================
    wav_all  = []
    flx_all  = []
    ivar_all = []
    dim_all  = []
    # ======================= READ DATA FROM FILES =======================
    for i in xrange(numfiles):
        wav  = []
        flx  = []
        err = []
        spec_dirname  = os.path.dirname(spec_allfiles[i])
        spec_basename = os.path.basename(spec_allfiles[i])
        print spec_basename
        print spec_dirname
        if 'corr' not in spec_basename:
            #print 'uncorrected'
            info, head = fitsio.read(spec_allfiles[i],ext=2,columns=['Z','RA','DEC'],header=True)
            data = fitsio.read(spec_allfiles[i],ext=1,columns=['FLUX','LOGLAM','IVAR'])
            outfile = '../paper/plots/'+thingid+'_spectrum.pdf'
        else:
            #print 'corrected'
            info, head = fitsio.read(spec_allfiles[i],ext=0,columns=['Z','RA','DEC'],header=True)
            data = fitsio.read(spec_allfiles[i],ext=1,columns=['WAV','FLUX_CORR','FLUX_ERR'])
            if 'margala' not in spec_dirname:
                outfile = '../paper/plots/'+thingid+'_corr_spectrum.pdf'
            else:
                outfile = '../paper/plots/'+thingid+'_corr_spectrum_margala.pdf'
        s = data.size
        dim_all.append(s)
        for j in xrange(s):
            if 'corr' not in spec_basename:
                z    = info[0][0]
                ra   = info[0][1]
                dec  = info[0][2]
                flx.append(data[j][0])
                wav.append(10**data[j][1])
                err.append(data[j][2])
            else:
                z    = head['Z']
                ra   = head['RA']
                dec  = head['DEC']
                flx.append(data[j][2])
                wav.append(data[j][0])
                err.append(data[j][1])
            # print wav[j], flx[j], err[j]
        flx_all.append(flx)
        wav_all.append(wav)
        err.append(err)
    #l = min(dim_all)
    
    

    majorLocator   = MultipleLocator(500)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(100)

    # ======================= PLOT =======================

    # PLOT
    xmin = np.min(wav_all[0]) ; xmax = np.max(wav_all[0])
    xplotlim = (xmin, xmax)
    colors = ['Black','Crimson','Blue','Green','Gold','Purple','Brown','Sienna','DarkMagenta','Teal','Turquoise','FireBrick','MediumSpringGreen','CornflowerBlue','Bisque','Yellow','Coral','Linen']
     # ======================= DEFINING THE CANVAS =======================
    fsize = (20,5)
    fig, axes =plt.subplots(2,1,figsize=fsize,sharey=True)
    #plt.suptitle('PLATEID=%04d, MJD=%5d, FIBREID=%04d,   RA = %10.6f, DEC = %10.6f , z = %5.3f' % (plate, mjd, fiberid, ra, dec, z))
    fig.subplots_adjust(hspace=0.15)
    plt.suptitle('RA = %10.6f, DEC = %10.6f , z = %5.3f' %(ra,dec,z))
    # ======================= LABELS =======================
    #fig.text(0.5,0.0, r'Observed wavelength $[\AA]$', ha='center')
    fig.text(0.09, 0.5, r'Flux $[10^{-17} \rm{erg/s/cm^2/\AA}]$', va='center', rotation='vertical')
    ax1, ax2 = axes.flat
    
    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []
    j = 0
    for i in xrange(numfiles):
        auxwav = np.array(wav_all[i])
        auxflx = np.array(flx_all[i])
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
        ymax = 2.0*np.percentile(auxflx, 99);   ymaxs.append(ymax)

        ax1.plot(auxwav,auxflx,color=colors[j],label=labels[i])
        ax2.plot(auxwav,auxflx,color=colors[j],label=labels[i])
        j = j+1

    xplotlim = (min(xmins),max(xmaxs),0.5*(min(xmins)+max(xmaxs)))
    yplotlim = (min(ymins),max(ymaxs))
    #print 'NEW GLOBAL MIN & MAX = ', yplotlim[0], yplotlim[1]
    ax1.set_xlim(xplotlim[0],xplotlim[2])
    ax2.set_xlim(xplotlim[2],xplotlim[1])

    ax1.set_ylim(yplotlim[0],yplotlim[1])
    ax2.set_xlabel(r'Observed wavelength $[\AA]$')
    ax1.legend(loc='upper right',prop={'size':10})

    plt.savefig(outfile, format = 'pdf')
    plt.show()

if __name__=='__main__':
    main()
