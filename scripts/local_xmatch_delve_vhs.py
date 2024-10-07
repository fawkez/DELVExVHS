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
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, hstack

# Set up astroquery
delve = TapPlus(url='https://datalab.noirlab.edu/tap')


# def query_delve_ra_dec(ra_min, ra_max, dec_min, dec_max):
#     query = f"""SELECT o.*
#             FROM delve_dr2.objects AS o
#             WHERE o.extended_class_g = 0
#             AND o.extended_class_i = 0
#             AND o.extended_class_z = 0
#             AND (o.mag_psf_g - o.mag_psf_i) > 1.2
#             AND (o.mag_psf_g - o.mag_psf_i) < 2.6
#             AND o.flags_g = 0
#             AND o.flags_i = 0
#             AND o.flags_r = 0
#             AND o.flags_z = 0
#             AND o.ra > {ra_min}
#             AND o.ra < {ra_max}
#             AND o.dec > {dec_min}
#             AND o.dec < {dec_max}
#             """
    
#     # Run query using the TAP service 
#     job = delve.launch_job(query= query)

#     return job.get_results()

def query_delve_ra_dec(ra_min, ra_max, dec_min, dec_max):
    query = f"""SELECT o.*
            FROM delve_dr2.objects AS o
            WHERE o.extended_class_g = 0
            AND o.extended_class_i = 0
            AND o.extended_class_z = 0
            AND o.flags_g = 0
            AND o.flags_i = 0
            AND o.flags_r = 0
            AND o.flags_z = 0
            AND o.ra > {ra_min}
            AND o.ra < {ra_max}
            AND o.dec > {dec_min}
            AND o.dec < {dec_max}
            """
    
    # Run query using the TAP service 
    job = delve.launch_job(query= query)

    return job.get_results()

def sky_crossmatch(table1, table2, tolerance=1*u.arcsec):
    """
    Perform a sky crossmatch between two astropy tables with a 1 arcsecond tolerance
    and include the on-sky distances for matched objects.

    Parameters:
    - table1, table2: astropy.table.Table
        The input tables with columns for right ascension and declination.
    - tolerance: astropy.units.Quantity
        The tolerance for the crossmatch, default is 1 arcsecond.

    Returns:
    - crossmatched_table : astropy.table.Table
        An astropy table containing matched rows from both tables with an additional
        column for the on-sky distances.
    """

    # Ensure right ascension and declination columns are present
    if 'ra' not in table1.colnames or 'dec' not in table1.colnames:
        raise ValueError("Table 1 must have 'ra' and 'dec' columns.")
    if 'RA2000' not in table2.colnames or 'DEC2000' not in table2.colnames:
        raise ValueError("Table 2 must have 'RA2000' and 'DEC2000' columns.")

    # Creating SkyCoord objects for both tables
    coords1 = SkyCoord(ra=table1['ra'], dec=table1['dec'])
    coords2 = SkyCoord(ra=table2['RA2000'], dec=table2['DEC2000'])

    # Performing the crossmatch
    idx, d2d, _ = coords1.match_to_catalog_sky(coords2)

    # Selecting matches within the tolerance
    match_mask = d2d < tolerance
    matched_rows_table1 = table1[match_mask]
    matched_rows_table2 = table2[idx[match_mask]]
    matched_distances = d2d[match_mask]

    # Adding the distance column to one of the matched tables
    matched_rows_table1['on_sky_distance_arcsec'] = matched_distances

    # Combine matched rows into a new table
    crossmatched_table = hstack([matched_rows_table1, matched_rows_table2])

    return crossmatched_table


if __name__ == '__main__':

    vhs_dir = '/fast_scratch3/mncavieres/VHS'
    output_delve =  '/fast_scratch3/mncavieres/DELVE_radec'
    output_matched =  '/fast_scratch3/mncavieres/xmatch_delve_large'

    # Iterate over each file
    for file in tqdm(os.listdir(vhs_dir)):
        if not file.endswith('.fits'):
            continue

        # open file
        vhs_data = Table.read(os.path.join(vhs_dir, file))

        # get min max coordinates
        ra_min = vhs_data['RA2000'].min() 
        ra_max = vhs_data['RA2000'].max()
        dec_min = vhs_data['DEC2000'].min()
        dec_max = vhs_data['DEC2000'].max()

        # query data
        delve_data = query_delve_ra_dec(ra_max=ra_max, ra_min=ra_min, dec_max=dec_max, dec_min=dec_min)


        # if region is empty pass: 
        if len(delve_data) < 1: 
            print('Region empty')
            continue
    
        # save data
        delve_data.write(os.path.join(output_delve, f'delve_{ra_min}_{ra_max}_{dec_min}_{dec_max}.dat'), format = 'ascii', overwrite= True)

        
        # xmatch with VHS data
        xmatch_result = sky_crossmatch(delve_data, vhs_data)

        # save the xmatch result
        xmatch_result.write(os.path.join(output_matched,  f'delve_x_vhs_{ra_min}_{ra_max}_{dec_min}_{dec_max}.dat'), format = 'ascii', overwrite= True)
