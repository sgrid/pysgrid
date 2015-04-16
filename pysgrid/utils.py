'''
Created on Mar 23, 2015

@author: ayan
'''
import re
from collections import namedtuple
import numpy as np
from .custom_exceptions import CannotFindPaddingError, DimensionMismatchError


GridPadding = namedtuple('GridPadding', ['mesh_topology_var',  # the variable containing the padding information
                                         'dim',  # the topology attribute
                                         'sub_dim',  # node dimension within the topology attribute
                                         'padding'  # padding type for the node dimension
                                         ]
                         )


def check_array_dims(*args):
    """
    Given an unspecified number of
    numpy arrays, make sure they all
    have the same dimensions. A 
    DimensionMismatch exception is raised
    if there is a mismatch.
    
    """
    array_shapes = [arr.shape for arr in args]
    arrays_match = check_element_equal(array_shapes)
    if not arrays_match:
        raise DimensionMismatchError(*array_shapes)
    else:
        pass


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
    check_array_dims(x_array, y_array)
    x_shape = x_array.shape
    paired_array_shape = x_shape + (2,)
    slices = (np.s_[:],) * len(x_shape)
    x_slices = slices + (0,)
    y_slices = slices + (1,)
    paired_array = np.empty(paired_array_shape, dtype=np.float64)
    paired_array[x_slices] = x_array[:]
    paired_array[y_slices] = y_array[:]
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
                padding_val = padding_info[-1]
                slice_datum = sgrid_obj.padding_slices[padding_val]
                lower_slice, upper_slice = slice_datum
                slice_index = np.s_[lower_slice:upper_slice]
                slice_indices += (slice_index,)
            except StopIteration:
                slice_index = np.s_[:]
                slice_indices += (slice_index,)
    else:
        pass
    return slice_indices


class ParsePadding(object):
    """
    Parse out the padding types from
    variables with a cf_role of 'grid_topology'.
    
    """
    padding_slices = {'both': (1, -1),
                      'none': (None, None),
                      'low': (1, None),
                      'high': (None, 1)
                      }
    
    def __init__(self, mesh_topology_var=None):
        self.mesh_topology_var = mesh_topology_var

    def parse_padding(self, padding_str):
        """
        Use regex expressions to break apart an
        attribute string containining padding types
        for each variable with a cf_role of 
        'grid_topology'.
        
        Padding information is returned within a named tuple
        for each node dimension of an edge, face, or vertical
        dimension. The named tuples have the following attributes:
        mesh_topology_var, dim_name, dim_var, and padding.
        Padding information is returned as a list
        of these named tuples.
        
        :param str padding_str: string containing padding types from a netCDF attribute
        :return: named tuples with padding information
        :rtype: list
        
        """
        p = re.compile('([a-zA-Z0-9_]+:) ([a-zA-Z0-9_]+) (\(padding: [a-zA-Z]+\))')
        padding_matches = p.findall(padding_str)
        padding_type_list = []
        for padding_match in padding_matches:
            raw_dim, raw_sub_dim, raw_padding_var = padding_match
            dim = raw_dim.split(':')[0]
            sub_dim = raw_sub_dim
            cleaned_padding_var = re.sub('[\(\)]', '', raw_padding_var)  # remove parentheses
            padding_type = cleaned_padding_var.split(':')[1].strip()  # get the padding value and remove spaces
            grid_padding = GridPadding(mesh_topology_var=self.mesh_topology_var,
                                       dim=dim,
                                       sub_dim=sub_dim,
                                       padding=padding_type
                                       )
            padding_type_list.append(grid_padding)
        if len(padding_type_list) > 0:
            final_padding_types = padding_type_list
        else:
            final_padding_types = None
            raise CannotFindPaddingError
        return final_padding_types        