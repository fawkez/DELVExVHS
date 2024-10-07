from mpi4py import MPI
import numpy as np
import os
from astropy.table import Table
from matplotlib.path import Path

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Define the function to check if points are inside a given polygon
def isInside(teffs, loggs, *vertices):
    points = np.array(vertices).reshape(-1, 2)
    polygon = Path(points)
    test_points = np.column_stack((teffs, loggs))
    return polygon.contains_points(test_points)

# Define vertices
#vertices = (1.31, 2.15, 1.88, 2.43, 3.12, 3.86, 1.32, 2.47)
vertices =  (1.63,2.43, 1.08,2.17, 1.54,2.61, 2.00,2.95, 2.64,3.28, 2.04,2.65)
#(1.790,2.386, 1.364,2.250, 2.121,2.864, 2.098,2.714)

# Define the paths
path_delve_x_vhs = '/fast_scratch3/mncavieres/xmatch_delve_large' #/fast_scratch3/mncavieres/xmatch_delve'
output_giants = '/fast_scratch3/mncavieres/giants_polygon_v2'

# Function to process each file
def process_file(file):
    data = Table.read(os.path.join(path_delve_x_vhs, file), format = 'ascii')
    data['g_i'] = data['mag_psf_g'] - data['mag_psf_i']
    data['i_Ks'] = data['mag_psf_i'] - data['KSAPERMAG3']
    hc_predictions = isInside(data['g_i'], data['i_Ks'], vertices)
    data_giants = data[hc_predictions]
    data_giants.write(os.path.join(output_giants, file), overwrite=True, format = 'ascii')

if __name__ == '__main__':
    # List all FITS files
    all_files = [f for f in os.listdir(path_delve_x_vhs) if f.endswith('.dat')]
    # Split files among processes
    files_to_process = np.array_split(all_files, size)[rank]

    # Process each file assigned to this MPI process
    for file in files_to_process:
        process_file(file)

    # Sync all processes
    comm.Barrier()
    if rank == 0:
        print("All files processed.")


# small change to test git