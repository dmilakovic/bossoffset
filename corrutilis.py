import fitsio
import numpy as np
import argparse
import os
import sys
from glob import glob
from SDSSmodules.SDSSclasses import spectrum
from congrid import congrid
import copy
import warnings

class exposure:
    def __init__(self):
        self.filename     =  ""
        self.cameras      =  ""
        self.number       = int(0)
        self.mjd          = int(0)
        self.plateid      = int(0)
        self.RA           = 0.0
        self.DEC          = 0.0
        self.az           = 0.0
        self.alt          = 0.0
        self.airmass      = 0.0
        self.pressure     = 0.0
        self.airtemp      = 0.0
        self.humidity     = 0.0
        self.seeing50     = 0.0
        self.extname      = ""

        self.npix         = int(0)
        self.loglam       = []
        self.flux         = []
        self.ivar         = []
        

    
# Class for the program settings

class program_settings:
    def __init__(self):
        self.debug            = 0
        self.program_name     = 'RMS'
        self.plot             = False
        self.smoothness       = 5

def read_spec(filepath):
    """
    Reads a FITS file containing the spectrum of a Quasar from the input path
    and returns a 'spectrum' class object into the main program.

    Parameters:
    -----------
          filepath : Path to the FITS file

    Returns:
    --------
          spectrum : A 'spectrum' class object
    """
    
    spec_dirname  = os.path.dirname(filepath)
    spec_basename = os.path.basename(filepath)
    #initialise spectrum
    spec = spectrum()
    spec.filename = spec_basename
    if 'corr' not in spec_basename:
        #print 'uncorrected'
        fits = fitsio.FITS(filepath)

        
        spec.ra     = fits[2]['RA'][0]
        spec.dec    = fits[2]['DEC'][0]
        spec.z      = fits[2]['Z'][0]
        
        spec.loglam = fits[1]['LOGLAM'][:]
        spec.wav    = np.power(10,spec.loglam)
        spec.flx    = fits[1]['FLUX'][:]
        spec.err    = fits[1]['IVAR'][:]
        spec.ivar   = fits[1]['IVAR'][:]
        spec.npix   = np.size(spec.loglam)

        spec.numexp = fits[-1].get_extnum() - 3 #(-4 because of the 4 standard HDUs + 1 because we start at index 0, i.e. fits[0])

        spec.loglam_start = np.min(spec.loglam)
        spec.loglam_end   = np.max(spec.loglam)
    else:
        #print 'corrected'
        fits = fitsio.FITS(filepath)
        header = fits[0].read_header()
    
        spec.ra     = header['RA']
        spec.dec    = header['DEC']
        spec.z      = header['Z']
        spec.npix   = header['NPIX']

        spec.wav      = fits[1]['WAV'][:]
        spec.flx      = fits[1]['FLUX'][:]
        spec.flx_corr = fits[1]['FLUX_CORR'][:]
        spec.err      = fits[1]['FLUX_ERR'][:]
        if 'margala' not in spec_dirname:
            spec.corr5400 = fits[1]['CORR5400'][:]
            spec.corr4000 = fits[1]['CORR4000'][:]
            spec.corr     = spec.corr5400/spec.corr4000

        spec.loglam_start = np.log10(spec.wav[0])
        spec.loglam_end   = np.log10(spec.wav[spec.npix-1])

    return spec

def return_spectral_list(directory):
    """
    Returns a tuple with paths to the FITS files in the specified directory.

    Parameters:
    -----------
        directory : A directory with FITS files. Needs to have subdirectories:
                    'uncorrected', 'corrected', and 'corrected/margala'

    Returns:
    --------
             list : A tuple containing the paths to the FITS files in the subfolders.
    """
    from glob import glob
    import os
    uncorrected_files = glob(os.path.join(directory,'uncorrected','*.fits'))
    corrected_files   = glob(os.path.join(directory,'corrected','*.fits'))
    corrected_files_m = glob(os.path.join(directory,'corrected','margala','*.fits'))

    if not (uncorrected_files):
        print 'Could not find uncorrected spectra'
    if not (corrected_files):
        print 'Could not find corrected spectra'
    if not (corrected_files_m):
        print 'Could not find Margala corrected spectra'
    list = [uncorrected_files,corrected_files,corrected_files_m]
    return list

def return_exposure(HDU):
    """
    Returns an 'exposure' class object with the data on the individual SDSS exposure and the spectrum.

    Parameters:
    -----------
           HDU : A Header-and-Data-Unit containing the data on the individual exposure. Object returned by FITSIO.
    Returns:
    --------
           exp : Exposure class object.
    """
    import fitsio
    import numpy as np

    # initialise exposure
    exp = exposure()

    # read data from the header
    header = HDU.read_header()
    exp.filename     = header["FILENAME"]
    exp.cameras      = header["CAMERAS"]
    exp.number       = header["EXPOSURE"]
    exp.mjd          = header["MJD"]
    exp.plateid      = header["PLATEID"]
    exp.RA           = header["RA"]
    exp.DEC          = header["DEC"]
    exp.az           = header["AZ"]
    exp.alt          = header["ALT"]
    exp.airmass      = header["AIRMASS"]
    try:
        exp.pressure = header["PRESSURE"]
    except ValueError:
        exp.pressure = None
    try:
        exp.airtemp  = header["AIRTEMP"]
    except ValueError:
        try:
            exp.airtemp = header["TEMP"]
        except:
            exp.airtemp  = None
            raise ValueError('No temperature found!')
    exp.humidity     = header["HUMIDITY"]
    exp.seeing50     = header["SEEING50"]
    exp.extname      = header["EXTNAME"]

    exp.npix         = np.size(HDU["LOGLAM"][:])
    exp.loglam       = HDU["LOGLAM"][:]
    exp.flux         = HDU["FLUX"][:]
    exp.ivar         = HDU["IVAR"][:]

    return exp

def pause():
    """
    Pauses the execution of the program until a keyboard key is pressed.
    """
    x = raw_input('Press any key to continue')
    return

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.nanmedian(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.nanmedian(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
def median_spectrum(redshift, wave_grid, flux_rms, flux_std):

    medianSpec = spectrum()
    medianSpec.z         = redshift
    medianSpec.npix      = wave_grid.size
    medianSpec.wave      = wave_grid
    
    medianSpec.flux      = flux_rms['BOSS']
    medianSpec.flux_error= flux_std['BOSS']
    medianSpec.__fit_powerlaw__(type='BOSS')
    medianSpec.alpha_BOSS = medianSpec.alpha
    

    medianSpec.flux_corr = flux_rms['CORR']
    medianSpec.flux_error= flux_std['CORR']
    medianSpec.__fit_powerlaw__(type='CORR')
    medianSpec.alpha_CORR = medianSpec.alpha

    medianSpec.flux_corr = flux_rms['MARG']
    medianSpec.flux_error= flux_std['MARG']
    medianSpec.__fit_powerlaw__(type='CORR')
    medianSpec.alpha_MARG = medianSpec.alpha

    print '{0:>45s}'.format('SPECTRAL INDEX')
    print '{0:>33s} '.format('SI_medspec')

    print '{0:>20s} {1:15.12f}'.format('UNCORRECTED:', medianSpec.alpha_BOSS)
    print '{0:>20s} {1:15.12f}'.format('CORRECTED (Ours):', medianSpec.alpha_CORR)
    print '{0:>20s} {1:15.12f}'.format('CORRECTED (Margala):',medianSpec.alpha_MARG)

    return medianSpec
def rms(data,variance=None):
    '''
    Return the variance weighted mean of the array
    '''
    data     = data[~np.isnan(data)]
    if variance!=None:
        variance = variance[~np.isnan(data)]
    
    if data.size>0:
        rms = np.average(data, weights=variance)
    else:
        rms = np.nan
    #rms = np.sqrt(np.nanmean(np.square(data)))
    return rms
def read_spectra(path):
    np.set_printoptions(threshold=6000)
    outspec = spectrum()
    
    if os.path.isdir(path):
        filetype = 'dir'
    elif os.path.isfile(path):
        filetype = 'file'

    labels = []
    if os.path.isdir(path):
        spectral_list       = return_spectral_list(path)
        numspec             = len(spectral_list[0])
        uncorrected_files   = spectral_list[0]
        corrected_files     = spectral_list[1]
        corrected_files_m   = spectral_list[2]
        # list that saves the number of files in each subfolder
        numfiles = [len(uncorrected_files),len(corrected_files),len(corrected_files_m)]
        dir_basename = os.path.dirname(path)
        thingid = dir_basename.split('/')[7]
        # read the data from the uncorrected FITS files
        for file in uncorrected_files:
            spec_basename = os.path.basename(file).split('.fits')[0]
            plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
            labels.append('%04d-%5d-%04d' % (plate, mjd, fiberid))
    if os.path.isfile(path):
        spectral_list = []
        spectral_list.append([path])
        numspec              = 1
        dir_basename = os.path.dirname(path)
        thingid = dir_basename.split('/')[7]
        spec_basename = os.path.basename(path).split('.fits')[0]
        plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
        labels.append('%04d-%5d-%04d' % (plate, mjd, fiberid))

    # ======================= UNIVERSAL WAVELENGTH GRID =======================
    loglam_grid = np.arange(start=3.55, stop = 4.02, step = 0.0001)
    wave_grid   = np.power(10, loglam_grid)
    wave_npix   = wave_grid.size
    wave_shape  = wave_grid.shape
    npix        = wave_npix
    # ======================= ======================== ======================
    keys         = ['BOSS','CORR','MARG']
    row_flux     = np.dtype([('BOSS',np.float32,1),('CORR',np.float32,1),('MARG',np.float32,1)])
    row_data     = np.dtype([('WAVE',np.float32,1),('FLUX',row_flux,1),('FLUX_ERR',row_flux,1),('POWERLAW',row_flux,1),('POWERLAW_ERR',row_flux,1), ('CORR',row_flux,1)])
    row_si       = np.dtype([('BOSS',np.float32,3),('CORR',np.float32,3),('MARG',np.float32,3)])

    data         = np.zeros(shape = (numspec,wave_npix), dtype = row_data); data.fill(np.nan)
    data['WAVE'] = wave_grid
    si_data      = np.ones(shape = (numspec,), dtype=row_si)
    # ======================= ======================== ======================
    outspec.keys = keys
    outspec.npix = wave_npix
    outspec.wave = wave_grid
    # ======================= READ DATA FROM FILES =======================
    for i in xrange(len(spectral_list)):
        for j,file in enumerate(spectral_list[i]):
            # ======================= READ SPECTRA =======================
            spec = spectrum()
            spec_dirname        = os.path.dirname(file)
            spec_basename       = os.path.basename(file)
            filetype, fileclass = check_filetype(file)
            if filetype == 'BOSS':
                spec.__read_BOSS__(file)
                flux_boss = np.zeros(shape=wave_shape)
            elif filetype == 'CORR':
                spec.__read_CORR__(file)
                flux_corr = np.zeros(shape=wave_shape)
                # __IMPORTANT__: use corrected flux as 'flux'
                spec.flux = spec.flux_corr
            
            # fit a powerlaw to the spectrum
            spec.__fit_powerlaw__(type=filetype)
            # ======================= REBINNING  =======================
            minwav = spec.wave.min() ; maxwav = spec.wave.max()
            # find where in the rebinned wavelength grid to put the rebinned flux values
            index_start = np.argmin(abs(wave_grid-minwav))
            index_stop  = np.argmin(abs(wave_grid-maxwav))
            # the spectra is rebinned to 'newpix' pixels
            newpix = index_stop - index_start
            # put new flux values & the powerlaw into the data structure
            data[j]['FLUX']        [fileclass][index_start:index_stop] = congrid.rebin_1d(spec.flux, newpix)
            data[j]['FLUX_ERR']    [fileclass][index_start:index_stop] = congrid.rebin_1d(spec.flux_error, newpix)
            data[j]['POWERLAW']    [fileclass][index_start:index_stop] = congrid.rebin_1d(spec.powerlaw, newpix)
            data[j]['POWERLAW_ERR'][fileclass][index_start:index_stop] = congrid.rebin_1d(spec.powerlaw_error, newpix)
            if len(spec.corr)>0:
                data[j]['CORR'][fileclass][index_start:index_stop] = congrid.rebin_1d(spec.corr, newpix)
            #for k in xrange(spec.npix):
            #    print j, k, fileclass, spec.flux_error[k], data[j]['FLUX_ERR'][fileclass][k]
            # save the data on the spectral index & error into the si_data structure
            si_data[j][fileclass][0] = spec.alpha
            si_data[j][fileclass][1] = spec.alpha_error
            si_data[j][fileclass][2] = spec.delta

            outspec.ra = spec.ra ; outspec.dec = spec.dec; outspec.z = spec.z
            
            del(spec)
    # ======================= ======================== ======================      
    
    # ======================= ======================== ======================      
    # calculate the RMS of flux for each spectrum class
    dataT = data.transpose()
    rmsFlux  = np.zeros(shape=wave_shape, dtype=row_flux); rmsFlux.fill(np.nan)
    rmsR2    = np.zeros(shape=wave_shape, dtype=row_flux); rmsR2.fill(np.nan)
    rmsErr   = np.zeros(shape=wave_shape, dtype=row_flux); rmsErr.fill(np.nan)
    rmsCorr  = np.zeros(shape=wave_shape, dtype=row_si); rmsCorr.fill(np.nan)
    #rmsCorrE = np.zeros(shape=wave_shape, dtype=row_flux); rmsCorrE.fill(np.nan)
    rmsPL    = np.zeros(shape=wave_shape, dtype=row_flux); rmsPL.fill(np.nan)
    rmsPLe   = np.zeros(shape=wave_shape, dtype=row_flux); rmsPLe.fill(np.nan)
    for key in keys:
        for pix in xrange(npix):
            # calculate the variance weigthed mean and deviation from numspec values of flux
            rmsFlux [pix][key]    = rms(dataT[pix]['FLUX'][key],variance=dataT[pix]['FLUX_ERR'][key]**2)
            rmsR2   [pix][key]    = np.sum((dataT[pix]['FLUX'][key] - rmsFlux[pix][key])**2/(np.nanstd(dataT[pix]['FLUX'][key]))**2)
            rmsErr  [pix][key]    = np.nanstd(dataT[pix]['FLUX'][key])
            rmsCorr [pix][key][0] = np.nanpercentile(dataT[pix]['CORR'][key], 50)  # MEDIAN
            rmsCorr [pix][key][1] = np.nanpercentile(dataT[pix]['CORR'][key], 36)  # LOWER RANGE
            rmsCorr [pix][key][2] = np.nanpercentile(dataT[pix]['CORR'][key], 84)  # UPPER RANGE
            rmsPL   [pix][key]    = rms(dataT[pix]['POWERLAW'][key],variance=dataT[pix]['POWERLAW_ERR'][key]**2)
            rmsPLe  [pix][key]    = np.nanstd(dataT[pix]['POWERLAW'][key])
    outspec.flux           = rmsFlux  # has three keys: ['BOSS', 'CORR', 'MARG']
    outspec.flux_error     = rmsErr   # has three keys: ['BOSS', 'CORR', 'MARG']
    outspec.corr           = rmsCorr  # has three keys: ['BOSS', 'CORR', 'MARG']
    outspec.R2             = rmsR2    # has three keys: ['BOSS', 'CORR', 'MARG']
    outspec.powerlaw       = rmsPL    # has three keys: ['BOSS', 'CORR', 'MARG']
    outspec.powerlaw_error = rmsPLe   # has three keys: ['BOSS', 'CORR', 'MARG']
    # ======================= ======================== ======================
    auxspec = copy.deepcopy(outspec)
    rmsSIalpha   = np.zeros(shape = (1,),      dtype=row_si)
    rmsSIdelta   = np.zeros(shape = (1,),      dtype=row_si)
    for key in keys:
        auxspec.flux       = rmsFlux[key]
        auxspec.flux_error = rmsErr[key]
        auxspec.__fit_powerlaw__(type='BOSS')
        rmsSIalpha[key][0][0] = auxspec.alpha
        rmsSIalpha[key][0][1] = auxspec.alpha_error
        rmsSIdelta[key][0][0] = auxspec.delta
    del auxspec
    
    outspec.alpha       = np.array([rmsSIalpha[key][0][0] for key in keys])
    outspec.alpha_error = np.array([rmsSIalpha[key][0][1] for key in keys])
    outspec.delta       = np.array([rmsSIdelta[key][0][0] for key in keys])
    
    # spectrum metadata to return to the main program
    outspec.nspectra   = numspec
    # spectral indices of _INDIVIDUAL_(!) spectra
    outspec.alpha_BOSS = si_data['BOSS'] # has three values: [alpha, alpha_error, delta]
    outspec.alpha_CORR = si_data['CORR'] # has three values: [alpha, alpha_error, delta]
    outspec.alpha_MARG = si_data['MARG'] # has three values: [alpha, alpha_error, delta]

    outspec.labels  = labels
    outspec.thingid = thingid
    #del(spec)
    return data, outspec

def check_filetype(filepath):
    '''
    Function that checks for file type and class 

    Input:
    ------
          filepath : path to the spectrum file
    Output:
    -------
          filetype : 'BOSS' or 'CORR'
         fileclass : 'BOSS' or 'CORR' or 'MARG'
    '''
    file_basename = os.path.basename(filepath)
    file_dirname  = os.path.dirname(filepath)
    if 'corr' not in file_basename:
        filetype = 'BOSS'
        fileclass = 'BOSS'
    if 'corr' in file_basename:
        filetype = 'CORR'
        if 'margala' not in file_dirname:
            fileclass = 'CORR'
        elif 'margala' in file_dirname:
            fileclass = 'MARG'
    return filetype, fileclass

def print_si_data(outspec):
    '''
    Print data on spectral index of data read by the function 'read_spectra'
    '''
    print '{0:-^80}'.format('')
    print '{0:^80}'.format('SPECTRAL INDEX')
    print '{0:^80}'.format('----------------------------')
    print '{0:^26}{1:^26}{2:^26}'.format('BOSS', 'CORR', 'MARG')
    # print the spectral index for each spectrum and each spectrum type
    for s in xrange(outspec.nspectra):
        si, err, delta = (np.transpose(val) for val in zip(outspec.alpha_BOSS[s], outspec.alpha_CORR[s], outspec.alpha_MARG[s]))
        h = np.hstack((si, err))
        print '{0:>+12.3f}+-{3:<12.4f}{1:>+12.3f}+-{4:<12.4f}{2:>+12.3f}+-{5:<12.4f}'.format(*h)
    # calculate and print the RMS of the spectral index for each spectrum type
    itemlist = [outspec.alpha_BOSS.transpose(), outspec.alpha_CORR.transpose(), outspec.alpha_MARG.transpose()]
    si_list  = []
    err_list = []
    for item in itemlist:
        #si, err = ((np.nanmedian(val), np.nanstd(val)) for val in item)
        si, err, delta = ((np.nanmedian(val), abs(np.nanmedian(val)-np.nanpercentile(val,16)), abs(np.nanmedian(val)-np.nanpercentile(val,84))) for val in item)
        si_list.append(si)
        err_list.append(err)
    k = np.hstack(si_list)
    print '{0:-^80}'.format(' QUANTILES: 50%, -(50%-16%), +(84%-50%)')
    print '{0:>+12.3f}-{1:<4.2f}+{2:<4.2f}{3:>+16.3f}-{4:<4.2f}+{5:<4.2f}{6:>+16.3f}-{7:<4.2f}+{8:<4.2f}'.format(*k)

    return

def print_rms_data(outspec):
    '''
    Print data on the root-mean-square of signal from data read by the function 'read_spectra'
    '''
    row_flux = np.dtype([('BOSS',np.float32,1),('CORR',np.float32,1),('MARG',np.float32,1)])
    #row_data = np.dtype([(])
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'Degrees of freedom <= 0 for (slice|axis).')
        warnings.filterwarnings('ignore', r'Mean of empty (slice|axis).')
        print '{0:-^80}'.format('')
        print '{0:^80}'.format('ROOT-MEAN-SQUARE OF FLUX')
        print '{0:^80}'.format('----------------------------'); print
        print '{0:-^10}|{1:-^42}|{2:-^26}'.format(' REGION ',' RMS +- RMSD ',' Q ')
        print '{0:^10}|{1:^14}{2:^14}{3:^14}|{1:^12}{3:^12}'.format('','BOSS', 'CORR', 'MARG')
        # print the SNR in specific regions (restframe)
        #regions   = [[912,1215],[1215,1550],[1650,1750],[1780,1910],[2790,2950],[2950,3050],[5100,5250]]
        regions   = [[912,1215],[1215,1550],[1550,1750],[1750,2000],[2000,2500],[2790,2950],[2950,3050],[5100,5250],
                     [1280, 1292],[1312, 1328],[1345, 1365],[1440, 1475],[1685, 1715],[1730, 1742],[1805, 1837],[2020, 2055], [2190, 2210]]
        waves      = []
        for region in regions:
            waves.append([region[0]*(1.+outspec.z),region[1]*(1.+outspec.z)])
        wave_rest = outspec.wave/(1.+outspec.z)
        for region in regions:
            start_pix  = np.argmin(abs(region[0] - wave_rest))
            stop_pix   = np.argmin(abs(region[1] - wave_rest))
            n=[] ; s = []
            boss_flux_error = outspec.flux_error['BOSS'][start_pix:stop_pix]
            for key in outspec.keys:                
                flux       = outspec.flux[key][start_pix:stop_pix]
                flux_error = outspec.flux_error[key][start_pix:stop_pix]
                corr       = outspec.corr[key][start_pix:stop_pix]
                # remove outliers over 5-sigma threshold
                thresh = 5
                good_flux  = ~is_outlier(flux, thresh)       ; flux = flux[good_flux]
                good_error = ~is_outlier(flux_error, thresh) ; flux_error = flux_error[good_error]
                signal     = np.nanmedian(flux)
                noise      = np.nanmedian(flux_error)
                bossnoise  = np.nanmedian(boss_flux_error[good_error])
                # Q = 1 - RMSD_[key]/RMSD_BOSS
                qval       = - 1 + np.nanmedian(flux_error/boss_flux_error[good_error])
                n.append((signal,noise, qval))
            h = np.hstack(n)
            f = np.vstack(n) #; print f
            print '{s[0]:>04d}-{s[1]:>4d}{0:>7.2f}+-{1:<6.2f}{3:>7.2f}+-{4:<6.2f}{6:>7.2f}+-{7:<6.2f} {5:>12.3%}{8:>12.3%}'.format(s=region, *h)
            
        print '{0:-^80}'.format(' ENTIRE SPECTRUM ')
        l = []
        boss_flux_error = outspec.flux_error['BOSS']
        for key in outspec.keys:
            flux       = outspec.flux[key]
            flux_error = outspec.flux_error[key]
            # remove outliers over 5-sigma threshold
            thresh = 5
            good_flux  = ~is_outlier(flux, thresh)         ; flux = flux[good_flux]
            good_error = ~is_outlier(flux_error, thresh)   ; flux_error = flux_error[good_error]
            signal     = np.nanmedian(flux)
            noise      = np.nanmedian(flux_error)
            bossnoise  = np.nanmedian(boss_flux_error[good_error])
            qval       = - 1 + np.nanmedian(flux_error/boss_flux_error[good_error])
            pval       = noise/bossnoise
            l.append((signal,noise, qval))
        k = np.hstack(l)
        print '{s:^9s}{0:>7.2f}+-{1:<6.2f}{3:>7.2f}+-{4:<6.2f}{6:>7.2f}+-{7:<6.2f} {5:12.3%}{8:12.3%}'.format(s='entire', *k)
    return

def print_pl_data(outspec):
    '''
    Print data on the root-mean-square of signal from data read by the function 'read_spectra'
    '''
    row_flux = np.dtype([('BOSS',np.float32,1),('CORR',np.float32,1),('MARG',np.float32,1)])
    #row_data = np.dtype([(])
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'Degrees of freedom <= 0 for (slice|axis).')
        warnings.filterwarnings('ignore', r'Mean of empty (slice|axis).')
        print '{0:-^80}'.format('')
        print '{0:^80}'.format('POWERLAW OF FLUX')
        print '{0:^80}'.format('----------------------------'); print
        print '{0:-^10}|{1:-^42}|{2:-^26}'.format(' REGION ',' PL FLUX +- SIGMA ',' Q ')
        print '{0:^10}|{1:^14}{2:^14}{3:^14}|{1:^12}{3:^12}'.format('','BOSS', 'CORR', 'MARG')
        # print the SNR in specific regions (restframe)
        #regions   = [[912,1215],[1215,1550],[1650,1750],[1780,1910],[2790,2950],[2950,3050],[5100,5250]]
        regions   = [[912,1215],[1215,1550],[1550,1750],[1750,2000],[2000,2500],[2790,2950],[2950,3050],[5100,5250],
                     [1280, 1292],[1312, 1328],[1345, 1365],[1440, 1475],[1685, 1715],[1730, 1742],[1805, 1837],[2020, 2055], [2190, 2210]]
        waves      = []
        for region in regions:
            waves.append([region[0]*(1.+outspec.z),region[1]*(1.+outspec.z)])
        wave_rest = outspec.wave/(1.+outspec.z)
        for region in regions:
            start_pix  = np.argmin(abs(region[0] - wave_rest))
            stop_pix   = np.argmin(abs(region[1] - wave_rest))
            n=[] ; s = []
            boss_powerlaw_error = outspec.powerlaw_error['BOSS'][start_pix:stop_pix]
            for key in outspec.keys:                
                powerlaw       = outspec.powerlaw[key][start_pix:stop_pix]
                powerlaw_error = outspec.powerlaw_error[key][start_pix:stop_pix]
                corr       = outspec.corr[key][start_pix:stop_pix]
                # remove outliers over 5-sigma threshold
                thresh = 5
                good_flux  = ~is_outlier(powerlaw, thresh)       ; powerlaw = powerlaw[good_flux]
                good_error = ~is_outlier(powerlaw_error, thresh) ; powerlaw_error = powerlaw_error[good_error]
                signal     = np.nanmedian(powerlaw)
                noise      = np.nanmedian(powerlaw_error)
                bossnoise  = np.nanmedian(boss_powerlaw_error[good_error])
                # Q = 1 - RMSD_[key]/RMSD_BOSS
                qval       = - 1 + np.nanmedian(powerlaw_error/boss_powerlaw_error[good_error])
                n.append((signal,noise, qval))
            h = np.hstack(n)
            f = np.vstack(n) #; print f
            print '{s[0]:>04d}-{s[1]:>4d}{0:>7.2f}+-{1:<6.2f}{3:>7.2f}+-{4:<6.2f}{6:>7.2f}+-{7:<6.2f} {5:>12.3%}{8:>12.3%}'.format(s=region, *h)
            
        print '{0:-^80}'.format(' ENTIRE SPECTRUM ')
        l = []
        boss_powerlaw_error = outspec.powerlaw_error['BOSS']
        for key in outspec.keys:
            powerlaw       = outspec.powerlaw[key]
            powerlaw_error = outspec.powerlaw_error[key]
            # remove outliers over 5-sigma threshold
            thresh = 5
            good_flux  = ~is_outlier(powerlaw, thresh)         ; powerlaw = powerlaw[good_flux]
            good_error = ~is_outlier(powerlaw_error, thresh)   ; powerlaw_error = powerlaw_error[good_error]
            signal     = np.nanmedian(powerlaw)
            noise      = np.nanmedian(powerlaw_error)
            bossnoise  = np.nanmedian(boss_powerlaw_error[good_error])
            qval       = - 1 + np.nanmedian(powerlaw_error/boss_powerlaw_error[good_error])
            pval       = noise/bossnoise
            l.append((signal,noise, qval))
        k = np.hstack(l)
        print '{s:^9s}{0:>7.2f}+-{1:<6.2f}{3:>7.2f}+-{4:<6.2f}{6:>7.2f}+-{7:<6.2f} {5:12.3%}{8:12.3%}'.format(s='entire', *k)
    return

def print_info(outspec):
    '''
    Print information about the spectrum
    '''
    print '{:-^80}'.format('INFO')
    position = {'RA':outspec.ra, 'Dec':outspec.dec, 'z':outspec.z}
    info     = {'Thing ID':outspec.thingid, 'Observations': outspec.nspectra}
    for name, value in position.iteritems():
        print '{name:>12s}: {value:<+10.4f}'.format(name=name, value=value),
    print ''
    for name, value in info.iteritems():
        print '{name:>25s}: {value:<10}'.format(name=name, value=value),
    print ''
    print '{:-^80}'.format('')
