
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
import json
import time



# Set up astroquery
delve = TapPlus(url='https://datalab.noirlab.edu/tap')


def query_delve_ra_dec(ra_min, ra_max, dec_min, dec_max):
    query = f"""SELECT o.*
            FROM delve_dr2.objects AS o
            WHERE o.extended_class_g = 0
            AND o.extended_class_i = 0
            AND o.flags_g = 0
            AND o.flags_i = 0
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


def load_last_processed(filename='/home/mcavieres/halo_map/last_processed.json'):
    """ Load the last processed filename from a JSON file. """
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
            return data.get('last_processed', None)
        
    return None

def save_last_processed(filename, state_file='/home/mcavieres/halo_map/last_processed.json'):
    """ Save the last processed filename to a JSON file. """
    with open(state_file, 'w') as f:
        json.dump({'last_processed': filename}, f)

if __name__ == '__main__':
    vhs_dir = '/fast_scratch3/mncavieres/VHS'
    output_delve = '/fast_scratch3/mncavieres/DELVE_radec_large'
    output_matched = '/fast_scratch3/mncavieres/xmatch_delve_large'
    last_processed = load_last_processed()

    # Sort files to ensure processing order
    files = sorted(os.listdir(vhs_dir))

    # Iterate over each file
    for file in tqdm(files):
        if not file.endswith('.fits') or (last_processed and file <= last_processed):
            continue  # Skip non-FITS files and already processed files

        # Open file

        vhs_data = Table.read(os.path.join(vhs_dir, file))
        # Ensure right ascension and declination columns are present
        if 'RA2000' not in vhs_data.colnames or 'DEC2000' not in vhs_data.colnames:
            continue
        if len(vhs_data['RA2000']) < 1:
            print('Skipping empty file')
            continue
        # Get min max coordinates
        ra_min, ra_max = vhs_data['RA2000'].min(), vhs_data['RA2000'].max()
        dec_min, dec_max = vhs_data['DEC2000'].min(), vhs_data['DEC2000'].max()

        # Check if it has been processed through filename
        filename_xmatch = f'delve_x_vhs_{ra_min}_{ra_max}_{dec_min}_{dec_max}.dat'
        if os.path.isfile(os.path.join(output_matched, filename_xmatch)):
            print('File already processed')
            continue # Skip already processed files before the json implementation
        
        # Query data with two attempts
        try:
            delve_data = query_delve_ra_dec(ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max)
        except:
            print('Query Failed waiting 5 min')
            time.sleep(300)
            try:
                delve_data = query_delve_ra_dec(ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max)
            except:
                print('Query failed again file will be skipped')
                continue

        # If region is empty, continue to next file
        if len(delve_data) < 1:
            print('Region empty')
            continue

        # Save DELVE data
        filename_delve = f'delve_{ra_min}_{ra_max}_{dec_min}_{dec_max}.dat'
        delve_data.write(os.path.join(output_delve, filename_delve), format='ascii', overwrite=True)

        # Crossmatch with VHS data
        xmatch_result = sky_crossmatch(delve_data, vhs_data)

        # Save the crossmatch result
        filename_xmatch = f'delve_x_vhs_{ra_min}_{ra_max}_{dec_min}_{dec_max}.dat'
        xmatch_result.write(os.path.join(output_matched, filename_xmatch), format='ascii', overwrite=True)

        # Save the name of the last processed file
        save_last_processed(file)
