import fitsio
import numpy as np
from corrutilis import *
from SDSSmodules.SDSSfiles import get_array_from_ind_exposures
import sys
import os
from astropy.table import Table,Column
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons

class catalogue:
    def __init__(self):
        self.plateid  = []
        self.mjd      = []
        self.fiberid  = []
        self.alt      = []
        self.az       = []
        self.seeing50 = []
        self.pressure = []
        self.airtemp  = []
        self.humidity = []

def main():
    spec_dir            = "/Volumes/Transcend/Isotropy/offset/spectra/qso/97176486/"
    spectral_list       = return_spectral_list(spec_dir)
    uncorrected_files   = spectral_list[0]
    corrected_files     = spectral_list[1]
    corrected_files_m   = spectral_list[2]

    nfiles = len(uncorrected_files)
    columns = ("PLATE","MJD","FIBER","ALT", "AZ", "SEE50", "PRESS", "TEMP", "HUM")
    types   = ("i4", "i5", "i4", 'f8', 'f8', 'f4', 'f4', 'f4', 'f4')
    

    cat = Table(names=columns,dtype=types)
    for file in uncorrected_files:
        
        fits      = fitsio.FITS(file)
        nHDU      = fits[-1].get_extnum()+1
        main_head = fits[0].read_header()

        spec       = spectrum()

        # spAll HDU
        HDU2 = fits[2]
        H2   = HDU2.read_header()

        spec.plateid    = HDU2["PLATE"][:][0]
        spec.mjd        = HDU2["MJD"][:][0]
        spec.fiberid    = HDU2["FIBERID"][:][0]
        
        spec.ra         = HDU2["RA"][:][0]
        spec.dec        = HDU2["DEC"][:][0]
        spec.z          = HDU2["Z"][:][0]
        spec.airmass    = HDU2["AIRMASS"][:][0]

        # INDIVIDUAL EXPOSURE HDUs
        altarr     = np.zeros(nHDU-4)
        azarr      = np.zeros(nHDU-4)
        see50arr   = np.zeros(nHDU-4)
        pressarr   = np.zeros(nHDU-4)
        airtemparr = np.zeros(nHDU-4)
        humarr     = np.zeros(nHDU-4)
        for i in xrange(4,nHDU):
            j=i-4
            HDUx          = fits[i]
            exp           = return_exposure(HDUx)
            altarr[j]     = exp.alt
            azarr[j]      = exp.az
            see50arr[j]   = exp.seeing50
            pressarr[j]   = exp.pressure
            airtemparr[j] = exp.airtemp
            humarr[j]     = exp.humidity
            exp = None

        spec.filename = os.path.basename(file)
        spec.alt      = np.mean(altarr)
        spec.az       = np.mean(azarr)
        spec.seeing50 = np.mean(see50arr)
        spec.pressure = np.mean(pressarr)
        spec.airtemp  = np.mean(airtemparr)
        spec.humidity = np.mean(humarr)

        cat.add_row((spec.plateid, spec.mjd, spec.fiberid, spec.alt, spec.az, spec.seeing50, spec.pressure, spec.airtemp, spec.humidity))
        #print "{0:6.4f} {1:6.4f} {2:6.4f} {3:6.4f} {4:6.4f} {5:6.4f}".format(spec.alt, spec.az, spec.seeing50, spec.pressure, spec.airtemp, spec.humidity)

    cat.write(os.path.join(spec_dir,"cat.fits"),format='fits')

    mjdArray   = cat.columns["MJD"]
    plateArray = cat.columns["PLATE"]
    fiberArray = cat.columns["FIBER"]
    altArray   = cat.columns["ALT"]
    azArray    = cat.columns["AZ"]
    see50Array = cat.columns["SEE50"]
    pressArray = cat.columns["PRESS"]
    tempArray  = cat.columns["TEMP"]
    humidArray = cat.columns["HUM"]

    fig, ax = plt.subplots()
    l0, = ax.plot(mjdArray, altArray, '.', visible=True, lw=2)
    l1, = ax.plot(mjdArray, azArray, '.', lw=2)
    l2, = ax.plot(mjdArray, tempArray, '.', lw=2)
    plt.subplots_adjust(left=0.2)

    rax = plt.axes([0.05, 0.4, 0.1, 0.15])
    check = CheckButtons(rax, ('Alt', 'Az', 'Temp'), (True, True, True))


    def func(label):
        if label == 'Alt':
            l0.set_visible(not l0.get_visible())
        elif label == 'Az':
            l1.set_visible(not l1.get_visible())
        elif label == 'Temp':
            l2.set_visible(not l2.get_visible())
        plt.draw()
    check.on_clicked(func)

    plt.show()

    
if __name__=="__main__":
    main()
