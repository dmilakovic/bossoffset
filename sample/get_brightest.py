import fitsio

def getKey(item):
    return item[6][2]

filename = '/Volumes/Transcend/QSOSIM10/DR10Q_r.fits'
data,head = fitsio.read(filename,ext=1,header=True, columns=['RA','DEC','SDSS_NAME','PSFMAG','PLATE','MJD','FIBER'])
print head
l = sorted(data, key=getKey, reverse=False)
for i in xrange(5):
    print i,'th brightest QSO in r band'
    print 'PLATE = ', l[i][3], 'MJD = ', l[i][4], 'FIBER = ', l[i][5] 
    print 'MAGNITUDES = ', l[i][6]
    print 'RA = ', l[i][1]
    print 'DEC = ', l[i][2]
    print 'SDSS_NAME = ', l[i][0]
    print '--------------------------------'

