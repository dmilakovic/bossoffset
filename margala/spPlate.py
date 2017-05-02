#
# Example script for looking at BOSS spectra and redshift fits via Python.
#
# Copyright Adam S. Bolton, Oct. 2009.
# Licensed for free and unencumbered use in public domain.
# No warranty express or implied.
#

# Imports:
import numpy as n
import pyfits as pf
import matplotlib as mpl
#mpl.use('TkAgg')
#mpl.interactive(True)
from matplotlib import pyplot as p

# Set topdir:
topdir = '/Volumes/Transcend/Isotropy/spectra/SDSS_DR12/Margala/portal.nersc.gov/project/boss/temp/sjbailey/dmargala/reduxtest/v5_7_0/'

# Pick your plate/mjd and read the data:
plate = '6130'
mjd = '56214'
spfile = topdir + plate + '/spPlate-' + plate + '-' + mjd + '.fits'
zbfile = topdir + plate + '/v5_7_0' + '/spZbest-' + plate + '-' + mjd + '.fits'
hdulist = pf.open(spfile)
c0 = hdulist[0].header['coeff0']
c1 = hdulist[0].header['coeff1']
npix = hdulist[0].header['naxis1']
wave = 10.**(c0 + c1 * n.arange(npix))
# Following commented-out bit was needed for some of the early redux:
#bzero = hdulist[0].header['bzero']
bunit = hdulist[0].header['bunit']
flux = hdulist[0].data
ivar = hdulist[1].data
hdulist.close()
hdulist = 0
hdulist = pf.open(zbfile)
print hdulist
synflux = hdulist[2].data
zstruc = hdulist[1].data
hdulist.close()
hdulist = 0

i = 711
p.figure()
# Set starting fiber point (above), then copy and paste
# the following repeatedly to loop over spectra:
#i+=1

# Following commented-out bit was needed for some of the early redux:
#p.plot(wave, (flux[i,:]-bzero) * (ivar[i,:] > 0), 'k', hold=False)
auxflx = flux[i,:] * (ivar[i,:] > 0)
p.plot(wave, auxflx, 'k', hold=False)
p.plot(wave, synflux[i,:], 'g', hold=True)
p.xlabel('Angstroms')
p.ylabel(bunit)
p.title(zstruc[i].field('class') + ', z = ' + str(zstruc[i].field('z')))
p.ylim(n.percentile(auxflx,1), 1.2*n.percentile(auxflx, 99))
p.show()
#p.savefig('/Volumes/Transcend/Isotropy/offset/python_ancillary/margala/plot.pdf')

