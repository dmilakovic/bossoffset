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
    parser.add_argument('--dir', type=str, default=None, 
        help='full path to the directory containing the FITS files, required')
    args = parser.parse_args()

    # ======================= FILES =======================
    if args.dir:
        spec_allfiles = glob(os.path.join(args.dir,'*.fits'))
        numfiles = len(spec_allfiles)
        labels = []
        dir_basename = os.path.dirname(args.dir)
        thingid = dir_basename.split('/')[6]
        for file in spec_allfiles:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid, zenith = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            #spec_filename.append(os.path.join(args.dir, 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)))
            labels.append('z = %2d' % zenith)
    #data, header = fitsio.read(file,ext=0,header=True)
    # ======================= =======================
    wav_all  = []
    flx_all  = []
    ivar_all = []
    dim_all  = []
    corr_all = []
    # ======================= READ DATA FROM FILES =======================
    for i in xrange(numfiles):
        wav  = []
        flx  = []
        err = []
        corr4000 = []
        corr5400 = []
        corr = []
        spec_dirname  = os.path.dirname(spec_allfiles[i])
        spec_basename = os.path.basename(spec_allfiles[i])
        #print spec_basename
        #print spec_dirname
        #print 'corrected'
        info, head = fitsio.read(spec_allfiles[i],ext=0,columns=['Z','RA','DEC'],header=True)
        data = fitsio.read(spec_allfiles[i],ext=1,columns=['WAV','FLUX_CORR','FLUX_ERR','CORR4000','CORR5400'])
        #if 'margala' not in spec_dirname:
        #    outfile = '../paper/plots/'+thingid+'_corr_spectrum.pdf'
        #else:
        outfile = '../paper/plots/'+thingid+'_corr_zenith.pdf'
        s = data.size
        dim_all.append(s)
        for j in xrange(s):
            z    = head['Z']
            ra   = head['RA']
            dec  = head['DEC']
            flx.append(data[j][2])
            wav.append(data[j][0])
            err.append(data[j][1])
            corr4000.append(data[j][3])
            corr5400.append(data[j][4])
            #print i, j, wav[j], corr4000[j], corr5400[j]
            if corr4000[j]!=0.0:
                corr.append(corr5400[j]/corr4000[j])
            else:
                corr.append(0.0)
        flx_all.append(flx)
        wav_all.append(wav)
        corr_all.append(corr)
        err.append(err)
    

    majorLocator   = MultipleLocator(1e-1)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(1e-2)

    # ======================= PLOT =======================

    # PLOT
    xmin = np.min(wav_all[0]) ; xmax = np.max(wav_all[0])
    xplotlim = (xmin, xmax)
    colors = ['Blue','Crimson','Green','DarkOrange','Skyblue','DarkMagenta','Teal','Burlywood','MediumSpringGreen','CornflowerBlue','Bisque','Gold','Coral','Linen']
     # ======================= DEFINING THE CANVAS =======================
    fsize = (9,6)
    fig , ax1 = plt.subplots(1,1,figsize=fsize,sharey=True)
    #plt.suptitle('PLATEID=%04d, MJD=%5d, FIBREID=%04d,   RA = %10.6f, DEC = %10.6f , z = %5.3f' % (plate, mjd, fiberid, ra, dec, z))
    #fig.subplots_adjust(hspace=0.15)
    #plt.suptitle('RA = %10.6f, DEC = %10.6f , z = %5.3f' %(ra,dec,z))
    # ======================= LABELS =======================
    fig.text(0.5,0.0, r'Observed wavelength $[\AA]$', ha='center')
    fig.text(0.03, 0.5, r'Correction function', va='center', rotation='vertical')
    
    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []
    for i in xrange(numfiles):
        auxwav = np.array(wav_all[i])
        auxflx = np.array(corr_all[i])
        #for i in xrange(np.size(auxflx)):
        #    print i, auxwav[i], auxflx[i]
        # ======================= SMOOTHING =======================
        box=5
        auxflx = smooth_array(auxflx,box)
        if np.size(auxflx)!=np.size(auxwav):
            auxwav = auxwav[:np.size(auxflx)]
            #print np.size(auxwav), np.size(auxflx)
        # get the limits in y-axis

        xmin = np.min(auxwav) ; xmax = np.max(auxwav)
        xmins.append(xmin)    ; xmaxs.append(xmax)
        
        ymin = min(0, np.percentile(auxflx,  1)); ymins.append(ymin)
        ymax = 1.1*np.percentile(auxflx, 99.9);   ymaxs.append(ymax)

        ax1.plot(auxwav,auxflx,color=colors[i],label=labels[i])
        #ax2.plot(auxwav,auxflx,color=colors[i],label=labels[i])

    
    xplotlim = (min(xmins),max(xmaxs))
    xplotlim = (min(xmins),5500)
    yplotlim = (3e-1,2e0)
    #print 'NEW GLOBAL MIN & MAX = ', yplotlim[0], yplotlim[1]
    #ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xplotlim[0],xplotlim[1])
    #ax2.set_xlim(xplotlim[2],xplotlim[1])

    ax1.set_ylim(yplotlim[0],yplotlim[1])
    #ax1.set_xlabel(r'Observed wavelength $[\AA]$')
    ax1.legend(loc='lower right',prop={'size':10})

    plt.savefig(outfile, format = 'pdf')
    plt.show()

if __name__=='__main__':
    main()
