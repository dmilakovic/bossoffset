import os
import numpy as np

topdir  = os.path.normpath('/Volumes/Transcend/Isotropy/spectra/SDSS_DR12/FQSO')
infile  = '/Volumes/Transcend/Isotropy/spectra/spAll/failedQSO_validation.csv'
outfile = '/Volumes/Transcend/Isotropy/spectra/spAll/failedQSO_validation_filepaths.txt'

g = open(infile,'r')
f = open(outfile,'w')
for line in g:
    line = line.rstrip('\n')
    plate,mjd,fiber = [int(item) for item in line.split(',')]
    filename = os.path.join(topdir,'spec-{plate:4d}-{mjd:5d}-{fiber:04d}.fits\n'.format(plate=plate,mjd=mjd,fiber=fiber))
    f.write(filename)
f.close
