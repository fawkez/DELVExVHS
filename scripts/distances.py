import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import BarycentricInterpolator
from scipy.interpolate import KroghInterpolator
from scipy.interpolate import interp1d

import itertools

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

def calculate_distance(apparent_mag, absolute_mag):
    distance_pc =  10 ** (0.2 * (apparent_mag - absolute_mag + 5))
    return distance_pc




# compute distance for every star in each evolutionary stage

def photometric_distance(isochrone, catalog):
    # assuming a parsec isochrone with decam bands
    # Compute g_i color
    isochrone['g_i'] = isochrone['gmag']- isochrone['imag']
    if isinstance(isochrone, Table):
        isochrone = isochrone.to_pandas()

    catalog = catalog.to_pandas()

    

    for sequence in range(1,5):


        # Select evolutionary stange from 1 main sequence to 3
        isochrone_segment = isochrone.loc[isochrone.label == sequence ]

        #  build interpolator
        interp_parsec = interp1d(isochrone_segment['g_i'], isochrone_segment['gmag']) 

        # set up catalog
        
        catalog[f'dist_{sequence}'] = np.nan
        catalog[f'g_ABS_{sequence}'] = np.nan

        # select objects within the interpolation range

        within_interp_parsec = catalog.loc[(catalog['g_i']<max(isochrone_segment['g_i'])) & (catalog['g_i']>min(isochrone_segment['g_i']))]

        # compute distances
        g_abs = interp_parsec(within_interp_parsec['g_i'])
        dists = calculate_distance(within_interp_parsec['gmag'], g_abs)

        # save to catalog
        catalog.loc[(catalog['g_i']<max(isochrone_segment['g_i'])) & (catalog['g_i']>min(isochrone_segment['g_i'])), f'g_ABS_{sequence}'] = g_abs
        catalog.loc[(catalog['g_i']<max(isochrone_segment['g_i'])) & (catalog['g_i']>min(isochrone_segment['g_i'])), f'dist_{sequence}'] = dists

    return catalog
    
if __name__ == '__main__':
    # import isochrone for distance computation

    isochrone = pd.read_csv('/Users/mncavieres/Documents/2023-1/Investigacion2/Data/Isochrones/PARSEC/1.2_1.5_10Gyrs_SDSS_2Mass.csv')#pd.read_csv('/Users/mncavieres/Documents/2023-1/Investigacion2/Data/Isochrones/PARSEC/10gyr1.5dex.csv')
    isochrone['g_i'] = isochrone['gmag']- isochrone['imag']


    # import data for stars
    # It must contain g_i colors

    stars_data = Table.read('/Users/mncavieres/Documents/2024-1/Delve/Data/giants_polygon.fits')#Table.read('/Users/mncavieres/Documents/2023-1/Investigacion2/Data/Catalogs/Stars_sigma2.fits') #Table.read('/Users/mncavieres/Documents/2023-1/Investigacion2/Data/Catalogs/Stars_sigma.fits')
    #stars_data['g_i'] = stars_data['WAVG_MAG_PSF_G'] - stars_data['WAVG_MAG_PSF_I'] # uncomment for our data

    stars_data['gmag'] = stars_data['mag_psf_g']   
    

    distance_cat = photometric_distance(isochrone, stars_data)

    Table.from_pandas(distance_cat).write('/Users/mncavieres/Documents/2024-1/Delve/Data/giants_polygon_distances.fits', format = 'fits', overwrite = True)