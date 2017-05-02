import fitsio
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import argparse
import os
import sys
from glob import glob
from SDSSmodules.SDSSutilities import smooth_array
from SDSSmodules.SDSSfitting import fit_powerlaw_individual
from SDSSmodules.SDSSclasses import spectrum
from corrutilis import program_settings, return_spectral_list, is_outlier
from congrid import congrid


def main():
    # ======================= SETTINGS =======================
    settings = program_settings()
    settings.plot = True
    settings.smoothness = 5
    settings.resample = True
    settings.z_min = 0
    settings.z_max = 4
    settings.z_delta = 0.1
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
        for fileList in spectral_list:
            ids = []
            for file in fileList:
                spec_basename = os.path.basename(file).split('.fits')[0]
                plate, mjd, fiberid = [int(field) for field in os.path.splitext(spec_basename)[0].split('-')[1:]]
                ids.append('%04d-%5d-%04d' % (plate, mjd, fiberid))
            labels.append(ids)
    #data, header = fitsio.read(file,ext=0,header=True)
    # ======================= ======================== ======================
    flux_boss_list = []; flux_corr_list = []; flux_marg_list = []
    # ======================= UNIVERSAL WAVELENGTH GRID =======================
    loglam_grid = np.arange(start=3.55, stop = 4.02, step = 0.0001)
    wave_grid   = np.power(10, loglam_grid)
    wave_npix   = wave_grid.size
    wave_shape  = wave_grid.shape

    row  = np.dtype([('BOSS','d8',1),('CORR','d8',1),('MARG','d8',1)])
    si_data = np.zeros(shape = (max(numfiles)), dtype=row)
    pl_flux = np.zeros(shape = (max(numfiles),wave_npix), dtype=row)
    fl_data = np.zeros(shape = (max(numfiles),wave_npix), dtype=row)
    print pl_flux[0].size
    # ======================= CALCULATIONS =======================
    for i in xrange(len(spectral_list)):
        for j,file in enumerate(spectral_list[i]):
            # ======================= READ SPECTRA =======================
            spec = spectrum()
            spec_dirname  = os.path.dirname(file)
            spec_basename = os.path.basename(file)
            if 'corr' not in spec_basename:
                spec.__read_BOSS__(file)
            elif 'corr' in spec_basename:
                spec.__read_CORR__(file)

            minwav = spec.wave.min() ; maxwav = spec.wave.max()
            index_start = np.argmin(abs(wave_grid-minwav))
            index_stop  = np.argmin(abs(wave_grid-maxwav))
            newpix = index_stop - index_start
            
            # ======================= GET SPECTRAL INDEX =======================
            if 'corr' not in spec_basename:
                spec.__fit_powerlaw__(type='BOSS')
                flux_old = spec.flux
                si_data[j]['BOSS'] = spec.alpha
                fl_data[j]['BOSS'][index_start:index_stop] = congrid.rebin_1d(spec.flux, newpix)
                pl_flux[j]['BOSS'][index_start:index_stop] = congrid.rebin_1d(spec.powerlaw, newpix)#10.0**(spec.delta + spec.beta*loglam_grid)
            if 'corr' in spec_basename:
                spec.__fit_powerlaw__(type='CORR')
                flux_old = spec.flux_corr
                if 'margala' not in spec_dirname:
                    si_data[j]['CORR'] = spec.alpha
                    fl_data[j]['CORR'][index_start:index_stop] = congrid.rebin_1d(spec.flux_corr, newpix)
                    pl_flux[j]['CORR'][index_start:index_stop] = congrid.rebin_1d(spec.powerlaw, newpix)
                elif 'margala' in spec_dirname:
                    si_data[j]['MARG'] = spec.alpha
                    fl_data[j]['MARG'][index_start:index_stop] = congrid.rebin_1d(spec.flux_corr, newpix)
                    pl_flux[j]['MARG'][index_start:index_stop] = congrid.rebin_1d(spec.powerlaw, newpix)
    ra = spec.ra ; dec = spec.dec; z = spec.z
    del(spec)

    print np.median(si_data['BOSS']), np.std(si_data['BOSS'])
    print np.median(si_data['CORR']), np.std(si_data['CORR'])
    print np.median(si_data['MARG']), np.std(si_data['MARG'])

if __name__ == '__main__':
    main()
