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


def does_intersection_exist(a, b):
    set_a = set(a)
    try:
        set_b = set(b)
    except TypeError:
        intersect_exists = False
    else:
        intersect = set_a.intersection(set_b)
        if len(intersect) > 0:
            intersect_exists = True
        else:
            intersect_exists = False
    return intersect_exists


def determine_variable_slicing(sgrid_obj, nc_variable, method='center'):
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
    grid_variables = sgrid_obj.grid_variables
    if grid_variables is None:
        grid_variables = []
    var_dims = nc_variable.dimensions
    node_dims = tuple(sgrid_obj.node_dimensions.split(' '))
    separate_edge_dim_exists = does_intersection_exist(var_dims, node_dims)
    slice_indices = tuple()
    if separate_edge_dim_exists:
        try:
            padding = sgrid_obj.face_padding  # try 2D sgrid
        except AttributeError:
            padding = sgrid_obj.volume_padding  # if not 2D, try 3D sgrid
    else:
        padding = sgrid_obj.all_padding()
    if method == 'center':
        for var_dim in var_dims:
            try:
                padding_info = next((info for info in padding if info.dim == var_dim))
            except StopIteration:
                slice_index = np.s_[:]
                slice_indices += (slice_index,)
            else:
                padding_val = padding_info[-1]
                slice_datum = sgrid_obj.padding_slices[padding_val]
                lower_slice, upper_slice = slice_datum
                slice_index = np.s_[lower_slice:upper_slice]
                slice_indices += (slice_index, )
    else:
        pass
    return slice_indices


def infer_avg_axes(sgrid_obj, nc_var_obj):
    """
    Infer which numpy axis to average over given
    the a variable defined on the grid. Works
    well for 2D. Not so sure about 3D.
    
    """
    fe_padding = []
    if hasattr(sgrid_obj, 'face_padding') and sgrid_obj.face_padding is not None:
        fe_padding += sgrid_obj.face_padding
    if hasattr(sgrid_obj, 'face1_padding') and sgrid_obj.face1_padding is not None:
        fe_padding += sgrid_obj.face1_padding
    if hasattr(sgrid_obj, 'face2_padding') and sgrid_obj.face2_padding is not None:
        fe_padding += sgrid_obj.face2_padding
    if hasattr(sgrid_obj, 'face3_padding') and sgrid_obj.face3_padding is not None:
        fe_padding += sgrid_obj.face3_padding
    if sgrid_obj.edge1_padding is not None:
        fe_padding += sgrid_obj.edge1_padding
    if sgrid_obj.edge2_padding is not None:
        fe_padding += sgrid_obj.edge2_padding
    if hasattr(sgrid_obj, 'edge3_padding') and sgrid_obj.edge3_padding is not None:
        fe_padding += sgrid_obj.edge3_padding
    var_dims = nc_var_obj.dimensions
    # define center averaging axis for a variable
    for var_dim in var_dims:
        try:
            padding_info = next((info for info in fe_padding if info.dim == var_dim))
        except StopIteration:
            padding_info = None
            avg_dim = None
            continue
        else:
            avg_dim = var_dim  # name of the dimension we're averaging over
            break  # exit the loop once it's found
    if padding_info is not None and avg_dim is not None:
        var_position = var_dims.index(avg_dim)
        center_avg_axis = len(var_dims) - var_position - 1
    else:
        center_avg_axis = None
    # define the node averaging axis for a variable
    if center_avg_axis == 1:
        node_avg_axis = 0
    elif center_avg_axis == 0:
        node_avg_axis = 1
    else:
        node_avg_axis = None
    return center_avg_axis, node_avg_axis