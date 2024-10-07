import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import os
from astroquery import xmatch
# define path for VHS data

VHS_directory = 'dir_geryon'

for vhs_file in os.listdir(VHS_directory):
    if vhs_file.endswith('.fits') == False:
        # only work with fits files
        pass
    # load vhs data
    vhs_table = Table.read(os.path.join(VHS_directory, vhs_file))

    # select only sources with no flags
    vhs_table = vhs_table[vhs_table['Ks_flags'] == 0]

    # select sources with magntiudes equivalent to a K-giant beyond 40 kpc
    # and at less than 200 kpc, based on a PARSEC 2.0, 10 Gyr  -1.5 dex
    # tip of the RGB K-giant.

    mag_low = -6.79 + 5*np.log10(4e4) - 5
    mag_max = -6.79 + 5*np.log10(2e5) - 5
    vhs_table = vhs_table[(vhs_table['Ks_mag'] > mag_low) & (vhs_table['Ks_mag'] < mag_max) ] 
    
    # crossmatch with DELVE

    result = xmatch(vhs_table['ra'], vhs_table['dec'], 'DELVE - CODE') # no se cual es el acceso a esta tabla desde DELVE.


    # QUERY TAP

    query = '''SELECT TOP 100000 o.*
                FROM delve_dr2.objects AS o
                WHERE o.extended_class_g = 0
                AND o.extended_class_i = 0
                AND o.extended_class_z = 0
                AND (o.mag_psf_g - o.mag_psf_i) > 1.2
                AND (o.mag_psf_g - o.mag_psf_i) < 2.6
                AND o.flags_g = 0
                AND o.flags_i = 0
                AND o.flags_r = 0
                AND o.flags_z = 0
                AND o.glat < -30
                AND o.glon > 0
                AND o.glon < 120'''