import fitsio
import numpy as np
import matplotlib.pyplot as plt
import os

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator

specPath = '/Volumes/Transcend/Isotropy/spectra/SDSS_DR12/'
fitsFilePath = os.path.join(specPath,'Margala/portal.nersc.gov/project/boss/temp/sjbailey/dmargala/reduxtest/v5_7_0/6130/spFluxcorr-b1-00148896.fits.gz')
hduList = fitsio.FITS(fitsFilePath,'r')

print hduList

imgData1 = hduList[0].read()
imgData2 = hduList[1].read()

fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

im1 = ax1.imshow(imgData1, cmap = 'gray', aspect='auto')
divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right",size="20%", pad=0.05)
ax1.set_ylim(0,500)
cbar1 = plt.colorbar(im1, cax = cax1, ticks=MultipleLocator(0.5), format="%.2f")

im2 = ax2.imshow(imgData2, cmap = 'gray', aspect='auto')
divider2 = make_axes_locatable(ax2)
cax2 = divider2.append_axes("right",size="20%", pad=0.05)
ax2.set_ylim(0,500)
cbar2 = plt.colorbar(im2, cax = cax2, ticks=MultipleLocator(.01), format="%.2f")

plt.show()


