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
from congrid import congrid

def main():
    # ======================= PARSER =======================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spec-dir', type=str, default=None, 
        help='full path to the directory containing the FITS files, required')
    args = parser.parse_args()

    # ======================= FILES =======================
    if args.spec_dir:
        specList = glob(os.path.join(args.spec_dir,'*.fits'))
        numfiles = len(specList)
        labels = []
        dir_basename = os.path.dirname(args.spec_dir)
        thingid = dir_basename.split('/')[6]
        for file in specList:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            #spec_filename.append(os.path.join(args.dir, 'spec-%04d-%5d-%04d.fits' % (plate, mjd, fiberid)))
            labels.append('%04d-%5d-%04d' % (plate, mjd, fiberid))
    #data, header = fitsio.read(file,ext=0,header=True)
    # ======================= UNIVERSAL WAVELENGTH GRID =======================
    loglam_grid = np.arange(start=3.55, stop = 4.02, step = 0.0001)
    wave_grid   = np.power(10, loglam_grid)
    wave_npix   = wave_grid.size
    wave_shape  = wave_grid.shape
    # ======================= ======================== ======================
    flux_boss_list = []; flux_corr_list = []; flux_marg_list = []
    keys         = ['BOSS','CORR','MARG']
    row          = np.dtype([('WAVE','f8',1),('FLUX','f8',1),('FLUX_ERR','f8',1)])

    data         = np.zeros(shape = (numfiles,wave_npix), dtype=row)
    # ======================= READ DATA FROM FILES =======================
    for i,file in enumerate(specList):
        spec = spectrum()
        spec_dirname  = os.path.dirname(file)
        spec_basename = os.path.basename(file)
        if 'corr' not in spec_basename:
            #print 'uncorrected'
            spec.__read_BOSS__(file)
            key = 'BOSS'
            outfile = '../paper/plots/'+thingid+'_spectrum.pdf'
        else:
            #print 'corrected'
            spec.__read_CORR__(file)
            if 'margala' not in spec_dirname:
                key = 'CORR'
                outfile = '../paper/plots/'+thingid+'_corr_spectrum.pdf'
            else:
                key = 'MARG'
                outfile = '../paper/plots/'+thingid+'_corr_spectrum_margala.pdf'

        # ======================= REBINNING  =======================
        minwav = spec.wave.min() ; maxwav = spec.wave.max()
        index_start = np.argmin(abs(wave_grid-minwav))
        index_stop  = np.argmin(abs(wave_grid-maxwav))
        newpix = index_stop - index_start
        data[i]['WAVE'][index_start:index_stop]     = congrid.rebin_1d(spec.wave, newpix)
        if key=='BOSS':
            data[i]['FLUX'][index_start:index_stop] = congrid.rebin_1d(spec.flux, newpix)
        elif (key == 'CORR' or key == 'MARG'):
            data[i]['FLUX'][index_start:index_stop] = congrid.rebin_1d(spec.flux, newpix)
    
    ra = spec.ra ; dec = spec.dec ; z = spec.z

    majorLocator   = MultipleLocator(500)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(100)

    # ======================= PLOT =======================

    # PLOT
    xmin = wave_grid.min ; xmax = wave_grid.max
    xplotlim = (xmin, xmax)
    cmap = plt.get_cmap('Paired')
    line_colors = cmap(np.linspace(0,1,numfiles))
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
    for i in xrange(0,numfiles):
        auxwav = wave_grid
        auxflx = data['FLUX'][i]
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

        ax1.plot(auxwav,auxflx,color=line_colors[j],label=labels[i])
        ax2.plot(auxwav,auxflx,color=line_colors[j],label=labels[i])
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
