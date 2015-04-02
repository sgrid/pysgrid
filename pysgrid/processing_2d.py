'''
Created on Apr 2, 2015

@author: ayan
'''
import numpy as np


def determine_avg_axis(array_shape, dim_0_max, dim_1_max):
    """
    Only works where x and y are have padding values
    of 'both' right now.
    
    """
    try:
        avg_axis = array_shape.index(dim_0_max)
    except ValueError:
        avg_axis = array_shape.index(dim_1_max)
    return avg_axis


def vector_sum(x_arr, y_arr):
    vector_sum = np.sqrt(x_arr**2 + y_arr**2)
    return vector_sum


def rotate_vectors(x_arr, y_arr, angle_arr):
    x_rot = x_arr*np.cos(angle_arr) - y_arr*np.sin(angle_arr)
    y_rot = x_arr*np.sin(angle_arr) + y_arr*np.cos(angle_arr)
    return x_rot, y_rot


def avg_to_cell_center(data_array, avg_dim):
    if avg_dim == 0:
        da = np.transpose(data_array)
    else:
        da = data_array
    da_trim_low = da[:, 1:]
    da_trim_high = da[:, :-1]
    da_avg_raw = 0.5 * (da_trim_low + da_trim_high)
    if avg_dim == 0:
        da_avg = np.transpose(da_avg_raw)
    else:
        da_avg = da_avg_raw
    return da_avg