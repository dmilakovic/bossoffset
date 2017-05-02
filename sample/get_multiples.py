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
def main():

    # ================= PARSER =================
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--catalog', type=str, default=None,
        help='full path to the QUASAR catalogue FITS file, required')
    parser.add_argument('--sort', type=str, default=None,
        help='sorting key, e.g. brightest, position')
    parser.add_argument('--ra', type=float, default=None,
        help='target right ascension')
    parser.add_argument('--dec', type=float, default=None,
        help='target declination')
    parser.add_argument('--rad', type=float, default=None,
        help='search radius, must specify --ra and --dec')
    parser.add_argument('--id', action='store_true',
        help='look for duplicate target IDs')
    args = parser.parse_args()

    if args.catalog:
        catalog_path = args.catalog
        catalog_hdulist = fits.open(catalog_path)
        catdata = catalog_hdulist[1].data
        print catdata.size
    else:
        parser.error('Must specify the path to the QSO catalog using --cat option.')
    # ================= READ DATA FROM CATALOG =================
    SDSS_name = catdata.field('sdss_name')
    ra        = catdata.field('ra')
    dec       = catdata.field('dec')
    thing_id  = catdata.field('thing_id')
    plate     = catdata.field('plate')
    mjd       = catdata.field('mjd')
    fiber     = catdata.field('fiber')
    psfmag    = catdata.field('psfmag')
    # ================= THING ID tag =================
    if args.id:
        # the search criterion is the unique 'thing_id' identifier
        # return the positions of objects with the same thing_id value
        uniq_src = np.unique(thing_id)
        print uniq_src.size, thing_id.size
        if uniq_src.size == thing_id.size:
            print 'All QSOs in the catalog are unique, quitting program'
            sys.exit()
        index = np.zeros(np.size(thing_id),dtype=('i8,i6,'))
        count=-1
        for i in uniq_src:
            count = count+1
            index[count] = (i,np.nonzero(thing_id==i)[0])
            print count, index[count]
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
    else:
        print 'Incorrect coordinate input, please enter RA (--ra ##), DEC (--dec ##) and the search radius (--rad ##) in arcsec'

        
if __name__=='__main__':
    main()
