import fitsio
import argparse
import os
import numpy as np
from astropy.io import fits
import sys
from astropy.coordinates import SkyCoord
import astropy.units as u


def getKey(item):
    return item[6][2]
def list_duplicates(seq,item):
    start_at = -1
    locs = []
    index = start_at
    while True:
        try:
            loc = np.where(seq,item)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs
def print_info_verbose(arg):
    print 'PLATE    = ',plate[arg]
    print 'MJD      = ',mjd[arg]
    print 'FIBER    = ',fiber[arg]
    print 'RA       = ',ra[arg]
    print 'DEC      = ',dec[arg]
    print 'THING_ID = ',thing_id[arg]
    print '----------------------------'
def print_info_short(arg):
    print '%4i|%5i|%4i' %(plate[arg],mjd[arg],fiber[arg])
    
def main():

    global plate
    global mjd
    global fiber
    global ra
    global dec
    global thing_id 
    # ================= PARSER =================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--cat', type=str, default=None,
        help='full path to the QUASAR catalogue FITS file, required')
    parser.add_argument('--nobs', type=int, default=None,
        help='whith this option, the program will show only objects that were observed a specified number of times')
    parser.add_argument('--id', type=int, default=None,
        help='the unique identifier of the object')
    parser.add_argument('--ra', type=float, default=None,
        help='target right ascension')
    parser.add_argument('--dec', type=float, default=None,
        help='target declination')
    parser.add_argument('--rad', type=float, default=None,
        help='search radius, must specify --ra and --dec')
    parser.add_argument('-v', action='store_true',
        help='verbose')
    parser.add_argument('--index', type=int, default=None, nargs='+',
        help='position of the object in the catalog (row number)')
    parser.add_argument('--multi', action='store_true',
        help='identify duplicate QSOs by comparing unique object identifiers')
    args = parser.parse_args()

    if args.cat:
        catalog_path = args.cat
        columns = ['PLATE', 'MJD', 'FIBERID', 'RA', 'DEC','THING_ID']
        catdata, cathead = fitsio.read(catalog_path,ext=1, header=True,columns = columns)
    else:
        parser.error('Must specify the path to the QSO catalog using --cat option.')

    # CODE TO REDUCE THE spAll-DR12.fits file to contain only neccessary data
    #reduced_cat = fitsio.FITS('/Volumes/Transcend/Isotropy/spectra/spAll-DR12-reduced.fits','rw')
    #hdict = {'PLATE', 'MJD', 'FIBERID', 'RA', 'DEC','THING_ID'}
    #nrows = catdata.size
    #reduced_cat.write(catdata)
    #reduced_cat.close()
    #sys.exit()

    
    # ================= READ DATA FROM CATALOG =================
    s = catdata.size
    plate     = np.zeros(s,dtype='int')
    mjd       = np.zeros(s,dtype='int')
    fiber     = np.zeros(s,dtype='int')
    ra        = np.zeros(s,dtype='float')
    dec       = np.zeros(s,dtype='float')
    thing_id  = np.zeros(s,dtype='int')
    for i in xrange(s):
        plate[i]    = catdata[i][0]
        mjd[i]      = catdata[i][1]
        fiber[i]    = catdata[i][2]
        ra[i]       = catdata[i][3]
        dec[i]      = catdata[i][4]
        thing_id[i] = catdata[i][5]
    # ================= INDEX tag =================
    if args.index:
        index = args.index
        print 'PRINTING ALL OBSERVATIONS ON THE OBJECT IN ROW', index
        for i in index:
            if args.v:
                print_info_verbose(i)
            else:
                print_info_short(i)
    # ================= MULTI tag =================
    if args.multi:
        # the search criterion is the unique 'thing_id' identifier
        # return the positions of objects with the same thing_id value
        uniq_src = np.unique(thing_id)
        s = uniq_src.size
        if s == thing_id.size:
            print 'All QSOs in the catalog are unique, quitting program'
            sys.exit()
#        index = np.zeros(np.size(thing_id),dtype=('i8,i6'))
        k=0
        X = []
        for i in uniq_src:
            ind = np.nonzero(thing_id==i)[0]
            #index[count] = (i,np.nonzero(thing_id==i)[0])
            entry = (i,ind.size,ind)
            if args.nobs:
                if ind.size==args.nobs:
                    X.append(entry)
                    print k, X[k]
                    k = k+1
            else:
                X.append(entry)
                print k, X[k]
                k = k+1
    # ================= ID tag =================
    if args.id:
        pos = np.nonzero(thing_id==args.id)[0]
        print 'PRINTING DATA ON ALL OBSERVATIONS OF OBJECT ', args.id
        print '======================================================'
        for i in pos:
            if args.v:
                print_info_verbose(i)
            else:
                print_info_short(i)
    # ================= COORDINATES tag =================
    if (args.ra and args.dec and not args.rad):
        print 'Radius not given'
    elif (not args.ra and args.dec and args.rad):
        print 'RA not given'
    elif (args.ra and not args.dec and args.rad):
        print 'DEC not given'
    elif (args.ra and args.dec and args.rad):
        ra0  = args.ra
        dec0 = args.dec
        rad  = args.rad

        coord_0 = SkyCoord(ra0,dec0,unit='deg')
        dist_list = []
        coord_array = SkyCoord(ra,dec,unit='deg')
        dist_array = coord_0.separation(coord_array)
        src_list = []
        for src in xrange(catdata.size):
            inside_radius = dist_array[src].is_within_bounds(0*u.deg,rad*u.arcsec)
            if inside_radius:
                src_list.append(src)
                print src
        print src_list
    #elif (not args.ra and not args.dec and not args.rad):
    #    print 'Incorrect coordinate input, please enter RA (--ra ##), DEC (--dec ##) and the search radius (--rad ##) in arcsec'

        
if __name__=='__main__':
    main()
