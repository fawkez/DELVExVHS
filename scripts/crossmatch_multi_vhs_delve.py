import numpy as np
import os
import healpy as hp
from astroquery.utils.tap.core import TapPlus
from astroquery.xmatch import XMatch
from tqdm import tqdm
import astropy.units as u
from concurrent.futures import ThreadPoolExecutor, as_completed

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
    job = delve.launch_job(query=query)
    return job.get_results()

def process_healpix(healpix):
    output_delve =  '/Users/mncavieres/Documents/2024-1/Delve/Data/test/delve' #'/fast_scratch3/mncavieres/DELVE'
    output_matched =  '/Users/mncavieres/Documents/2024-1/Delve/Data/test/delve/vhs_delve' #'/fast_scratch3/mncavieres/delveXvhs'

    # Query data
    delve_data = query_delve_256(healpix)

    # Save data
    filename_delve = os.path.join(output_delve, f'delve_{healpix}.dat')
    delve_data.write(filename_delve, format='ascii')

    # XMatch with VHS data
    xmatch_result = XMatch.query(cat1=delve_data, cat2='vizier:II/367/vhs_dr5', max_distance=1 * u.arcsec, colRA1='ra', colDec1='dec', colRA2='RA2000', colDec2='DEC2000')

    # Save the xmatch result
    filename_xmatch = os.path.join(output_matched, f'delve_x_vhs_{healpix}.dat')
    xmatch_result.write(filename_xmatch, format='ascii')

if __name__ == '__main__':
    print('Running query for NSIDE 256')
    # Constants for Nside ring 256
    Nside = 256
    npix = hp.nside2npix(Nside)

    # Generate the array of HEALPix pixel indices and filter by galactic latitude
    healpix_indices = np.arange(npix)
    l, b = hp.pix2ang(Nside, healpix_indices, lonlat=True)
    southern_hemisphere_indices = healpix_indices[b < -20]

    # Set up ThreadPoolExecutor to manage parallel execution
    max_workers = 6  # Adjust based on your system and service limits

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks to the executor
        futures = {executor.submit(process_healpix, healpix): healpix for healpix in southern_hemisphere_indices}

        # Progress bar setup
        progress = tqdm(as_completed(futures), total=len(southern_hemisphere_indices), desc='Processing HEALPix')

        # Wait for and track the completion of all submitted tasks
        for future in progress:
            healpix_processed = futures[future]
            try:
                future.result()  # You can handle results or exceptions here
            except Exception as exc:
                print(f'HEALPix {healpix_processed} generated an exception: {exc}')
