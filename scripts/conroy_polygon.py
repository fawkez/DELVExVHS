import pandas as pd
import numpy as np
from astropy.table import Table
from tqdm import tqdm

from matplotlib.path import Path
import matplotlib.pyplot as plt
import os



# Define polygon selection function

def isInside(teffs, loggs, *vertices):
    # Assume teffs and loggs are iterable arrays of equal length
    points = np.array(vertices).reshape(-1, 2)  # Reshape flat list into a 2D array of (x, y) pairs
    polygon = Path(points)  # Create a Path object representing the polygon

    # Stack teffs and loggs into a 2D array of points
    test_points = np.column_stack((teffs, loggs))

    # Return a boolean array where each element tells if the corresponding point is inside the polygon
    return polygon.contains_points(test_points)

# Define vertices
vertices= ( 1.31,2.15, 1.88,2.43, 3.12,3.86, 1.32,2.47)



if __name__ == '__main__':

    path_delve_x_vhs = '/fast_scratch3/mncavieres/xmatch_delve'
    output_giants = '/fast_scratch3/mncavieres/giants_polygon'
    

    for file in tqdm(os.listdir(path_delve_x_vhs)):
        if not file.endswith('.fits'):
            continue
        data = Table.read(os.path.join(path_delve_x_vhs, file))
        # Predicting test set results
        data['g_i'] = data['mag_psf_g'] - data['mag_psf_i']
        data['i_Ks'] = data['mag_psf_i'] - data['KSAPERMAG3']
        hc_predictions = isInside(data['g_i'], data['i_Ks'], vertices)
        data_giants = data[hc_predictions]
        data_giants.writeto(os.path.join(output_giants, file))
