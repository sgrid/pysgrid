'''
Created on Mar 23, 2015

@author: ayan
'''
from collections import namedtuple

import numpy as np


GridPadding = namedtuple('GridPadding', ['mesh_topology_var',  # the variable containing the padding information
                                         'dim',  # the topology attribute
                                         'sub_dim',  # node dimension within the topology attribute
                                         'padding'  # padding type for the node dimension
                                         ]
                         )


def pair_arrays(x_array, y_array):
    """
    Given two arrays to equal dimensions,
    pair their values element-wise.
    
    For example given arrays [[1, 2], [3, 4]]
    and [[-1, -2], [-3, -4]], this function will
    return [[[1, -1], [2, -2]], [[3, -3], [4, -4]]].
    
    :param np.array x_array: a numpy array containing "x" coordinates
    :param np.array y_array: a numpy array containing "y" coordinates
    :return: array containing (x, y) arrays
    :rtype: np.array
    
    """
    x_shape = x_array.shape
    paired_array_shape = x_shape + (2,)
    paired_array = np.empty(paired_array_shape, dtype=np.float64)
    paired_array[..., 0] = x_array[:]
    paired_array[..., 1] = y_array[:]
    return paired_array


def check_element_equal(lst):
    """
    Check that all elements in an
    iterable are the same.
    
    :params lst: iterable object to be checked
    :type lst: np.array, list, tuple
    :return: result of element equality check
    :rtype: bool
    
    """
    return lst[1:] == lst[:-1]


def determine_variable_slicing(sgrid_obj, nc_dataset, variable, method='center'):
    """
    Figure out how to slice a variable. This function
    only knows who to figure out slices that would be
    used to trim data before averaging to grid cell
    centers; grid cell nodes will be supported later.
    
    :param sgrid_obj: an SGrid object derived from a netCDF file or netCDF4.Dataset object
    :type sgrid_obj: sgrid.SGrid
    :param nc_dataset: a netCDF4.Dataset object from which the sgrid_obj was derived
    :type nc_dataset: netCDF4.Dataset
    :param str variable: the name of a variable to be sliced
    :param str method: slice method for analysis at grid cell centers or grid cell nodes; accepts either 'center' or 'node'
    :return: the slice for the varible for the given method
    :rtype: tuple
    
    """
    var_obj = nc_dataset.variables[variable]
    grid_variables = sgrid_obj.grid_variables
    if grid_variables is None:
        grid_variables = []
    var_dims = var_obj.dimensions
    padding_summary = sgrid_obj._define_face_padding_summary()
    slice_indices = tuple()
    if method == 'center':
        for var_dim in var_dims:
            try:
                padding_info = next((info for info in padding_summary if info[0] == var_dim))
            except StopIteration:
                slice_index = np.s_[:]
                slice_indices += (slice_index,)
            else:
                padding_val = padding_info[-1]
                slice_datum = sgrid_obj.padding_slices[padding_val]
                lower_slice, upper_slice = slice_datum
                slice_index = np.s_[lower_slice:upper_slice]
                slice_indices += (slice_index,)
    else:
        pass
    return slice_indices                  