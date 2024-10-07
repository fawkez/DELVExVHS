""" script to download the entire DELVE catalog (or the useful sections)

    The idea is to organice in a HEALPix by HEALPix basis the download. 
    For each N-side 2048 HEALPix a TAP+ query will be performed to
    download the point sources that lie in the correct g - i range. 

    Then a crossmatch will be made with VHS to obtain the Ks band, after
    which the catalog for that HEALPix will be saved.

    After all the catalogs are constructed a master catalog will be made
    by using Dask to join all the fits files into a large parquet.

    Then color selections will be made using the parquet file.
"""
# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
#from astroquery import TAP
import healpy as hp
from astroquery.utils.tap.core import TapPlus
from astroquery.xmatch import XMatch
from tqdm import tqdm
import astropy.units as u
import os

# Set up astroquery
delve = TapPlus(url='https://datalab.noirlab.edu/tap')


def query_delve_256(healpix):
    query = f"""SELECT o.*
            FROM delve_dr2.objects AS o
            WHERE o.extended_class_g = 0
            AND o.ring256 = {healpix} 
            AND o.extended_class_i = 0
            AND o.extended_class_z = 0
            AND (o.mag_psf_g - o.mag_psf_i) > 1.2
            AND (o.mag_psf_g - o.mag_psf_i) < 2.6
            AND o.flags_g = 0
            AND o.flags_i = 0
            AND o.flags_r = 0
            AND o.flags_z = 0"""
    
    # Run query using the TAP service 
    job = delve.launch_job(query= query)

    return job.get_results()


if __name__ == '__main__':

    output_delve =  '/fast_scratch3/mncavieres/DELVE'
    output_matched =  '/fast_scratch3/mncavieres/delveXvhs'

    print('Running query for NSIDE 256')
    # Constants for Nside ring 256
    Nside = 256
    npix = hp.nside2npix(Nside)

    # Generate the array of HEALPix pixel indices
    healpix_indices = np.arange(npix)

    # Get the angular coordinates of the center of each pixel
    l, b = hp.pix2ang(Nside, healpix_indices, lonlat= True)

    # Filter indices based on galactic latitude
    southern_hemisphere_indices = healpix_indices[b < -20]

    # Now, `southern_hemisphere_indices` contains the HEALPix indices in the southern galactic hemisphere

    # Iterate over healpix
    for healpix in tqdm(southern_hemisphere_indices[6517:]):

        # query data
        delve_data = query_delve_256(healpix)


        # if healpix is empty pass: 
        if len(delve_data) < 1: 
            print('HEALPix', healpix, 'empty')
            continue
    
        # save data
        delve_data.write(os.path.join(output_delve, f'delve_{healpix}.dat'), format = 'ascii', overwrite= True)

        
        # xmatch with VHS data
        xmatch_result = XMatch.query(cat1=delve_data, cat2='vizier:II/367/vhs_dr5', max_distance=1 * u.arcsec, colRA1='ra', colDec1='dec', colRA2='RA2000', colDec2= 'DEC2000')
        
        # save the xmatch result
        xmatch_result.write(os.path.join(output_matched,  f'delve_x_vhs_{healpix}.dat'), format = 'ascii', overwrite= True)
