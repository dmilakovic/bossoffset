import fitsio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import argparse
import os
import sys
from glob import glob
from SDSSutilities import smooth_array
from SDSSmodules.SDSSclasses import spectrum
from congrid import congrid
import corrutilis
import itertools
import pandas as pd

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])


def main():
    # ======================= SETTINGS =======================
    settings = corrutilis.program_settings()
    settings.plot = True
    settings.smoothness = 5
    settings.resample = True
    # ======================= PARSER =======================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spec-dir', type=str, default=None, 
        help='full path to the directory containing the FITS files, required')
    parser.add_argument('--data', choices=['BOSS','CORR','MARG'], nargs='+', type=str, default=['BOSS','CORR'],
        help='data to use')
    parser.add_argument('--noplot', action='store_true',
        help='number of panels')
    parser.add_argument('--rms', action='store_true',
        help='plot the RMS')
    parser.add_argument('--pl', action='store_true',
        help='plot the median power-law fit')
    parser.add_argument('--corr', action='store_true',
        help='plot the mean of the correction functions')
    parser.add_argument('--indpl', action='store_true',
        help='plot the median power-law fit')
    parser.add_argument('--show', action='store_true',
        help='show the figure instead of saving it')
    parser.add_argument('--step', type=int, default=1,
        help='plots every x-th spectrum, where x = step')
    parser.add_argument('--smooth', type=int, default=None,
        help='smooth the spectra with a boxcar b, where b = smooth')
    parser.add_argument('--xmax', type=float, default=None,
        help='define the maximum of xrange')
    parser.add_argument('--xmin', type=float, default=None,
        help='define the minimum of xrange')
    args = parser.parse_args()

    keys          = args.data
    # READ DATA FROM INDIVIDUAL SPECTRA
    data, rmsSpec = corrutilis.read_spectra(args.spec_dir)

    # CREATE FILENAME & OUTFILE NAMES
    filename      = '{}'.format(rmsSpec.thingid)
    if args.rms:  filename = filename+'_RMS'
    if args.pl :  filename = filename+'_PL'
    if args.corr: filename = filename+'_CORRECTION'
    for key in args.data:
        filename = filename+'_{}'.format(key)
    outfile       = os.path.join(args.spec_dir,filename)
    # OPEN A FILE TO SAVE THE RESULTS
    outdb         = outfile + '.hdf5'
    if os.path.isfile(outdb):
        os.remove(outdb)
    outhdf        = pd.HDFStore(outdb,'w')
    wave = pd.Series(data['WAVE'][0])
    outhdf.put('WAVE',wave)
    datakeys = [item[0] for item in data['FLUX'].__array_interface__['descr']]
    for key in datakeys:
        df = pd.DataFrame(data['FLUX'][key].T)
        outhdf.put(key,df)

    outhdf.close()
    
    corrutilis.print_info(rmsSpec)
    corrutilis.print_si_data(rmsSpec)
    corrutilis.print_rms_data(rmsSpec)
    corrutilis.print_pl_data(rmsSpec)

    if args.noplot==False:
        ndata = len(args.data)
        # number of panels corresponds to the number of data used
        if not args.corr:
            npanels = ndata
        if args.corr:
            npanels = 1
    # ======================= PLOT =======================
        xmajorLocator   = MultipleLocator(1000)
        xmajorFormatter = FormatStrFormatter('%d')
        xminorLocator   = MultipleLocator(250)
        # PLOT NAME
        outfig = outfile+'.pdf'
        if not args.corr:
            fsize = (12, 2+2*ndata)
        if args.corr:
            fsize = (6,4)
        # PLOT PARAMETERS
        cmap = plt.get_cmap('Paired')
        line_colors = cmap(np.linspace(1,0,rmsSpec.nspectra))
        c     = {'BOSS':'b','CORR':'g','MARG':'purple'}
        lab   = {'BOSS':'BOSS','CORR':'Corrected BOSS','MARG':'Margala et al. 2015'}
        ls    = {'BOSS':':','CORR':'-','MARG':'--'}
        # ======================= DEFINING THE CANVAS =======================
        fig = plt.figure(figsize=fsize)
        delta  = 0.01
        left   = 0.1 ; right = 0.9 ; width  = right-left 
        bottom = 0.15 ; top   = 0.95
        # if plotting RMS, add aditional panel of height g, adjust other panels
        if args.rms:
            g = 0.15
            x = delta/g
        if args.rms and args.pl:
            g = 0.15
            x = delta/g
        if not args.rms:
            g = 0.0
            x = 0.0
        # height = total height for plotting data panles;        h = height of a panel (not including SNR)
        height = (top - bottom) - (npanels-1)*delta - g*(1+x) ;  h = height/npanels
        # add axes for data panels
        for i in xrange(npanels):
            b = top - (i+1)*h - i*delta
            fig.add_axes([left,b,width,h])
        # add panel for Q, if plotting RMS
        if args.rms:
            b = b - g - 2*delta
            fig.add_axes([left,b,width,g])
        # attach axes to figure
        ax = fig.axes
        # SET TITLE
        #plt.suptitle('THING ID = %9s' %(rmsSpec.thingid))
        # ======================= AXES LABELS =======================
        if not args.corr:
            fig.text(0.06, (g+delta+bottom+top)/2, r'Flux $[10^{-17} \rm{erg/s/cm^2/\AA}]$',va='center',  rotation='vertical')
        if args.corr:
            fig.text(0.0, (g+delta+bottom+top)/2, r'Correction',va='center',  rotation='vertical')
        xmins = [] ; ymins = [] ; qmin = []
        xmaxs = [] ; ymaxs = [] ; qmax = []

        # PLOTTING DATA
        box=settings.smoothness
        auxwav = rmsSpec.wave
        print '{:-^80}'.format('PLOTTING')
        
        
        alpha_BOSS, alpha_error_BOSS, delta_BOSS = (np.mean(rmsSpec.alpha_BOSS.T[0]), np.std(rmsSpec.alpha_BOSS.T[0]), np.mean(rmsSpec.alpha_BOSS.T[2]))
        alpha_CORR, alpha_error_CORR, delta_CORR = (np.mean(rmsSpec.alpha_CORR.T[0]), np.std(rmsSpec.alpha_CORR.T[0]), np.mean(rmsSpec.alpha_CORR.T[2]))
        alpha_MARG, alpha_error_MARG, delta_MARG = (np.mean(rmsSpec.alpha_MARG.T[0]), np.std(rmsSpec.alpha_MARG.T[0]), np.mean(rmsSpec.alpha_MARG.T[2]))
        alpha       = {'BOSS':alpha_BOSS,'CORR':alpha_CORR,'MARG':alpha_MARG}
        alpha_error = {'BOSS':alpha_error_BOSS,'CORR':alpha_error_CORR,'MARG':alpha_error_MARG}
        delta       = {'BOSS':delta_BOSS,'CORR':delta_CORR,'MARG':delta_MARG}
        
        for i, key in zip(xrange(ndata),keys):
            print '{0:>20s}'.format(key) , 
            # i goes over the uncorrected, corrected, corrected_margala
            k = 0
            # ----------------------------------------------------------------------------------------------------------
            # PLOT FLUXES
            if not args.rms and not args.pl and not args.corr:
                for j in xrange(0,rmsSpec.nspectra,args.step):
                    # j goes over the spectra
                    if args.smooth:
                        auxflx = smooth_array(data[j]['FLUX'][key], args.smooth)
                    else:
                        auxflx = data[j]['FLUX'][key]
                        auxpl  = data[j]['POWERLAW'][key]
                            

                    if np.size(auxflx)!=np.size(auxwav):
                        print np.shape(auxwav)
                        print np.shape(auxflx)
                        sys.exit('ERROR: Flux and wavelength grids not equal in size!')
                        #auxwav = auxwav[:np.size(auxflx)]

                    # get the limits in y-axis
                    ymin = np.nanpercentile(auxflx,  0.5);       ymins.append(ymin)
                    ymax = 1.2*np.nanpercentile(auxflx, 99.5);     ymaxs.append(ymax)
                    if args.indpl:
                        alpha=0.3
                    elif not args.indpl:
                        alpha=1
                    ax[i].plot(auxwav,auxflx,color=line_colors[j],label=rmsSpec.labels[j], linestyle='-', alpha=alpha)
                    if args.indpl:
                        ax[i].plot(auxwav,auxpl,color=line_colors[j],label=rmsSpec.labels[j], linestyle='--', linewidth=2, alpha=1)
                    k = k+1
            # ----------------------------------------------------------------------------------------------------------
            # PLOT RMS and Q
            if args.rms and not args.pl:
                lw = 1.5
                if args.smooth:
                    auxrms  = smooth_array(rmsSpec.flux[key], args.smooth)
                    auxdev  = smooth_array(rmsSpec.flux_error[key], args.smooth)
                    bossdev = smooth_array(rmsSpec.flux_error['BOSS'], args.smooth)
                    auxcorr = smooth_array(rmsSpec.corr[key], args.smooth) -1
                    R2      = smooth_array(rmsSpec.R2[key], args.smooth)
                else:
                    auxrms  = rmsSpec.flux[key]
                    auxdev  = rmsSpec.flux_error[key]
                    bossdev = rmsSpec.flux_error['BOSS']
                    auxcorr = rmsSpec.corr[key] - 1
                    R2      = rmsSpec.R2[key]
                qval = auxdev/bossdev - 1
                #absqval = np.absolute(qval[~np.isnan(qval)])
                #n = np.argmin(absqval)
                #print wave[n]
                
                ymin = np.nanpercentile(auxrms,  0.5);         ymins.append(ymin)
                ymax = 1.2*np.nanpercentile(auxrms, 99.5);     ymaxs.append(ymax)
                # plot RMS and RMSD
                ax[i].plot(auxwav, auxrms, color='g', linewidth=lw, label = 'Mean flux')
                ax[i].fill_between(auxwav, auxrms - auxdev, auxrms + auxdev, color='g', alpha=0.3)
                ax[i].plot(auxwav, R2, color='r', linewidth=lw, label = 'Residuals')
                # plot the Q RMS value
                ax[-1].plot(auxwav, qval, color=c[key], linewidth=lw, label = lab[key], linestyle = '-')
                #ax[-1].plot(auxwav, auxcorr, color=c[key], linewidth=2, label = lab[key], linestyle = '--')
                qmin.append(0.5*np.nanpercentile(qval, 1))
                qmax.append(1.5*np.nanpercentile(qval, 99))
            # ----------------------------------------------------------------------------------------------------------
            # PLOT POWER-LAWS
            #alpha = {'BOSS':rmsSpec.alpha_BOSS,'CORR':rmsSpec.alpha_CORR,'MARG':rmsSpec.alpha_MARG}
            if args.pl and not args.rms:
                lw = 1.5
                auxpl  = rmsSpec.powerlaw[key]
                auxple = rmsSpec.powerlaw_error[key]
                for j in xrange(0,rmsSpec.nspectra,args.step):
                    # j goes over the spectra
                    if args.smooth:
                        auxflx = smooth_array(data[j]['FLUX'][key], args.smooth)
                    else:
                        auxflx = data[j]['FLUX'][key]

                    if np.size(auxflx)!=np.size(auxwav):
                        print np.shape(auxwav)
                        print np.shape(auxflx)
                        sys.exit('ERROR: Flux and wavelength grids not equal in size!')
                        #auxwav = auxwav[:np.size(auxflx)]

                    # get the limits in y-axis
                    ymin = np.nanpercentile(auxflx,  0.5);       ymins.append(ymin)
                    ymax = 1.2*np.nanpercentile(auxflx, 99.5);     ymaxs.append(ymax)
                    ax[i].plot(auxwav,auxflx,color=line_colors[j],label=rmsSpec.labels[j], linestyle='-', alpha=1)
                    k = k+1
                #ax[i].plot(auxwav,auxpl, color = 'r', linewidth = lw, label = 'Powerlaw', linestyle = '--')
                ax[i].plot(auxwav,auxpl, color = 'b', linewidth = lw, label = 'Power-law', linestyle = '--')
                ax[i].fill_between(auxwav,auxpl - auxple, auxpl + auxple, color = 'b', alpha=0.6 )
            # ----------------------------------------------------------------------------------------------------------
            # PLOT RMS AND POWER-LAW
            if args.pl and args.rms:
                rmsAlpha    = {'BOSS': rmsSpec.alpha[0], 'CORR': rmsSpec.alpha[1],'MARG':rmsSpec.alpha[2]}
                rmsAlphaErr = {'BOSS': rmsSpec.alpha_error[0], 'CORR': rmsSpec.alpha_error[1],'MARG':rmsSpec.alpha_error[2]}
                rmsDelta    = {'BOSS': rmsSpec.delta[0], 'CORR': rmsSpec.delta[1],'MARG':rmsSpec.delta[2]}
                # emission-free regions:
                emfree_regions    = [[1280, 1292],[1312, 1328],[1345, 1365],[1440, 1475], 
                                [1685, 1715],[1730, 1742],[1805, 1837],[2020, 2055], [2190, 2210]]
                wave_regions = []
                for region in emfree_regions:
                    wave_regions.append([region[0]*(1.+rmsSpec.z),region[1]*(1.+rmsSpec.z)])
                # spectral indices
                lw = 1.5
                # spectral shapes
                auxpl  = rmsSpec.powerlaw[key]
                auxple = rmsSpec.powerlaw_error[key]
                if args.smooth:
                    auxrms = smooth_array(rmsSpec.flux[key], args.smooth)
                    auxdev = smooth_array(rmsSpec.flux_error[key], args.smooth)
                    bossdev = smooth_array(rmsSpec.flux_error['BOSS'], args.smooth)
                    auxcorr = smooth_array(rmsSpec.corr[key], args.smooth) -1 
                else:
                    auxrms = rmsSpec.flux[key]
                    auxdev = rmsSpec.flux_error[key]
                    bossdev = rmsSpec.flux_error['BOSS']
                    auxcorr = rmsSpec.corr[key] - 1
                qval = rmsSpec.powerlaw_error[key]/rmsSpec.powerlaw_error['BOSS'] - 1
                ymin = np.nanpercentile(auxrms,  0.5);         ymins.append(ymin)
                ymax = 1.2*np.nanpercentile(auxrms, 99.5);     ymaxs.append(ymax)
                qmax.append(1.2*np.nanpercentile(auxrms, 90))
                # plot RMS and RMSD
                rmspl = 10**(rmsDelta[key]+(-2-rmsAlpha[key])*np.log10(auxwav)); print 
                ax[i].plot(auxwav, auxrms, color='g', linewidth=lw, label = 'Mean flux')
                #ax[i].plot(auxwav, rmspl, color='r', linewidth=lw, label = 'RMS PL')
                ax[i].plot(auxwav,auxpl, color = 'orange', linewidth = lw, label = 'Power-law', linestyle = '--')
                ax[i].fill_between(auxwav,auxpl - auxple, auxpl + auxple, color = 'orange', alpha=0.6 )
                print key, 'MEAN of individual:', alpha[key], alpha_error[key], delta[key]
                print key, 'RMS', rmsAlpha[key], rmsAlphaErr[key], rmsDelta[key]
                
                # plot the Q value
                if key!='BOSS':
                    ax[-1].plot(auxwav, qval, color=c[key], linewidth=lw, label = lab[key], linestyle = '-')
                qmin.append(0.5*np.nanpercentile(qval, 1))
                qmax.append(1.5*np.nanpercentile(qval, 99))
                # plot vertical lines for regions
                for region in wave_regions:
                    ax[i].axvspan(region[0], region[1], alpha=0.2, color='grey')
            # ----------------------------------------------------------------------------------------------------------
            # PLOT CORRECTION FUNCTION
            if args.corr:
                lw = 1.5
                corr = rmsSpec.corr[key].T
                if args.smooth:
                    auxcorr    = smooth_array(corr[0], args.smooth)
                    auxdevlow  = smooth_array(corr[1], args.smooth)
                    auxdevhigh = smooth_array(corr[2], args.smooth)
                else:
                    auxcorr    = corr[0]
                    auxdevlow  = corr[1]
                    auxdevhigh = corr[2]
                ymin = np.nanpercentile(auxcorr,  0.5);         ymins.append(ymin)
                ymax = 1.2*np.nanpercentile(auxcorr, 99.5);     ymaxs.append(ymax)
                # plot CORR and limits
                print auxwav.shape, auxcorr.shape
                ax[-1].plot(auxwav, auxcorr, color=c[key], linewidth=lw, linestyle = ls[key], label = lab[key])
                ax[-1].fill_between(auxwav, auxdevlow, auxdevhigh, color=c[key], alpha=0.3)
                ax[-1].xaxis.set_major_locator(xmajorLocator)
                ax[-1].xaxis.set_minor_locator(xminorLocator)
            if not args.corr: 
                ax[i].xaxis.set_major_locator(xmajorLocator)
                ax[i].xaxis.set_minor_locator(xminorLocator)
        # XRANGE 
        if not args.xmin:
            xmin = rmsSpec.wave.min()
        if args.xmin:
            xmin = args.xmin
        if not args.xmax:
            xmax = rmsSpec.wave.max()
        if args.xmax:
            xmax = args.xmax
        xplotlim = (xmin, xmax)
        # YRANGE
        yplotlim = (min(ymins),max(ymaxs))
        for i in xrange(len(ax)):
            ax[i].set_xlim(xplotlim[0],xplotlim[1])
            ax[i].set_ylim(yplotlim[0],yplotlim[1])
        # X-AXIS LABELS
        for i in xrange(len(ax)-1):
            ax[i].xaxis.set_ticklabels([])
        ax[-1].set_xlabel(r'Observed wavelength $[\AA]$')
        # LEGEND
        if args.corr and ndata>1:
            ax[0].legend(loc='upper left', prop={'size':9})
        if (rmsSpec.nspectra <= 10) and not args.rms and not args.corr:
            ncol = int(rmsSpec.nspectra//2.)
            handles, labels = ax[0].get_legend_handles_labels()
            ax[0].legend(flip(handles,2), flip(labels,2), loc='upper right',prop={'size':9},ncol=ncol)
        if args.rms and not args.pl and not args.corr:
            ymajorLocator   = MultipleLocator(1)
            ymajorFormatter = FormatStrFormatter('%d')
            yminorLocator   = MultipleLocator(.25)
            ax[0].legend(loc='upper right',prop={'size':12})
            fig.text(0.06, (2*bottom + g)/2, r'$Q_{Flux}$', va='center', rotation='vertical')
            handles, labels = ax[-1].get_legend_handles_labels()
            handles = handles[1::2] ; labels = labels[1::2]
            ax[-1].legend(handles, labels, loc='upper left',prop={'size':12}, ncol=len(args.data))
            ax[-1].set_ylim(-0.85, max(qmax))
            ax[-1].xaxis.set_major_locator(xmajorLocator)
            ax[-1].xaxis.set_minor_locator(xminorLocator)
            ax[-1].yaxis.set_major_locator(ymajorLocator)
            ax[-1].yaxis.set_minor_locator(yminorLocator)
            print
        if args.rms and args.pl and not args.corr:
            ymajorLocator   = MultipleLocator(1)
            ymajorFormatter = FormatStrFormatter('%d')
            yminorLocator   = MultipleLocator(.25)
            ax[-1].plot(auxwav, np.zeros(shape=auxwav.shape), color='k', linewidth=lw, label = '', linestyle = '--')
            ax[0].legend(loc='upper right',prop={'size':12})
            fig.text(0.06, (2*bottom + g)/2, r'$Q_{PL}$', va='center', rotation='vertical')
            handles, labels = ax[-1].get_legend_handles_labels()
            #handles = handles[1::2] ; labels = labels[1::2]
            ax[-1].legend(handles, labels, loc='upper left',prop={'size':12}, ncol=len(args.data))
            ax[-1].set_ylim(-0.85, 1.75)
            ax[-1].xaxis.set_major_locator(xmajorLocator)
            ax[-1].xaxis.set_minor_locator(xminorLocator)
            ax[-1].yaxis.set_major_locator(ymajorLocator)
            ax[-1].yaxis.set_minor_locator(yminorLocator)
            print
        if args.show:
            plt.show()
        if not args.show:
            plt.savefig(outfig, format = 'pdf')
            print
            print 'PLOT SAVED TO {0:s}'.format(outfig)
        
if __name__=='__main__':
    main()
