from mpi4py import MPI
import os
import json
import numpy as np
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
from tqdm import tqdm

# MPI Initialization
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Define the directories (assuming they are accessible from all MPI nodes)
vhs_dir = '/fast_scratch3/mncavieres/VHS'
output_delve = '/fast_scratch3/mncavieres/DELVE_radec'
output_matched = '/fast_scratch3/mncavieres/xmatch_delve_large'
state_file = '/home/mcavieres/halo_map/last_processed_2.json'



def sky_crossmatch(table1, table2, tolerance=1*u.arcsec):
    """ Crossmatch function as defined earlier. """
    coords1 = SkyCoord(ra=table1['ra']*u.degree, dec=table1['dec']*u.degree)
    coords2 = SkyCoord(ra=table2['RA2000'], dec=table2['DEC2000'])
    idx, d2d, _ = coords1.match_to_catalog_sky(coords2)
    match_mask = d2d < tolerance
    matched_rows_table1 = table1[match_mask]
    matched_rows_table2 = table2[idx[match_mask]]
    matched_rows_table1['on_sky_distance_arcsec'] = d2d[match_mask]
    crossmatched_table = hstack([matched_rows_table1, matched_rows_table2])
    return crossmatched_table

def process_file(file):
    """Process each file."""
    if not file.endswith('.fits'):
        return  # Skip non-FITS files
    
    print(f'Working on {file}')

    vhs_data = Table.read(os.path.join(vhs_dir, file))
    if 'RA2000' not in vhs_data.colnames or 'DEC2000' not in vhs_data.colnames:
        return  # Ensure necessary columns are present
    if len(vhs_data)<1:
        return # Skip empty regions

    ra_min, ra_max = vhs_data['RA2000'].min(), vhs_data['RA2000'].max()
    dec_min, dec_max = vhs_data['DEC2000'].min(), vhs_data['DEC2000'].max()
    filename_xmatch = f'delve_x_vhs_{ra_min}_{ra_max}_{dec_min}_{dec_max}.dat'

    if os.path.isfile(os.path.join(output_matched, filename_xmatch)):
        return  # Skip already processed files

    filename_delve = f'delve_{ra_min}_{ra_max}_{dec_min}_{dec_max}.dat'
    if os.path.isfile(os.path.join(output_delve, filename_delve)):
        delve_data = Table.read(os.path.join(output_delve, filename_delve), format='ascii')
    else:
        print('Skipping Missing Delve data')
        return

    if len(delve_data) < 1:
        return  # Skip empty regions

    xmatch_result = sky_crossmatch(delve_data, vhs_data)
    xmatch_result.write(os.path.join(output_matched, filename_xmatch), format='ascii', overwrite = True)


# Main execution block
if __name__ == '__main__':
    all_files = sorted(os.listdir(vhs_dir))
    files_to_process = np.array_split(all_files, size)[rank]

    for file in files_to_process:
        if rank == 0:
            print(f"Process {rank} processing {file}")
        process_file(file)

    comm.Barrier()
    if rank == 0:
        print("All files processed.")

    MPI.Finalize()
