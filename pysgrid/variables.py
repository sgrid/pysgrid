'''
Created on Apr 15, 2015

@author: ayan
'''

from __future__ import (absolute_import, division, print_function)

from .read_netcdf import parse_axes, parse_vector_axis
from .utils import determine_variable_slicing, infer_avg_axes, infer_variable_location


class SGridVariable(object):
    """
    Object of variables found and inferred
    from an SGRID compliant dataset.

    """
    def __init__(self,
                 center_axis=None,
                 center_slicing=None,
                 coordinates=None,
                 dimensions=None,
                 dtype=None,
                 grid=None,
                 location=None,
                 node_axis=None,
                 node_slicing=None,
                 standard_name=None,
                 variable=None,
                 vector_axis=None,
                 x_axis=None,
                 y_axis=None,
                 z_axis=None,):
        self.center_axis = center_axis
        self.center_slicing = center_slicing
        self.coordinates = coordinates
        self.dimensions = dimensions
        self.dtype = dtype
        self.grid = grid
        self.location = location
        self.node_axis = node_axis
        self.node_slicing = node_slicing
        self.standard_name = standard_name
        self.variable = variable
        self.vector_axis = vector_axis
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.z_axis = z_axis

    @classmethod
    def create_variable(cls, nc_var_obj, sgrid_obj):
        variable = nc_var_obj.name
        try:
            grid = nc_var_obj.grid
        except AttributeError:
            grid = None
            center_axis = None
            node_axis = None
        else:
            center_axis, node_axis = infer_avg_axes(sgrid_obj, nc_var_obj)
        center_slicing = determine_variable_slicing(sgrid_obj,
                                                    nc_var_obj,
                                                    method='center')
        dimensions = nc_var_obj.dimensions
        dtype = nc_var_obj.dtype
        try:
            location = nc_var_obj.location
        except AttributeError:
            location = infer_variable_location(sgrid_obj, nc_var_obj)
        if location == 'edge':
            if center_axis == 0:
                location = 'edge2'
            elif center_axis == 1:
                location = 'edge1'
            else:
                location = None
        try:
            axes = nc_var_obj.axes
        except AttributeError:
            x_axis = None
            y_axis = None
            z_axis = None
        else:
            x_axis, y_axis, z_axis = parse_axes(axes)
        try:
            standard_name = nc_var_obj.standard_name
        except AttributeError:
            standard_name = None
            vector_axis = None
        else:
            vector_axis = parse_vector_axis(standard_name)
        try:
            raw_coordinates = nc_var_obj.coordinates.strip()
        except AttributeError:
            coordinates = None
        else:
            coordinates = tuple(raw_coordinates.split())
        sgrid_var = cls(variable=variable,
                        grid=grid,
                        x_axis=x_axis,
                        y_axis=y_axis,
                        z_axis=z_axis,
                        center_slicing=center_slicing,
                        center_axis=center_axis,
                        node_slicing=None,
                        node_axis=node_axis,
                        dimensions=dimensions,
                        dtype=dtype,
                        location=location,
                        standard_name=standard_name,
                        vector_axis=vector_axis,
                        coordinates=coordinates)
        return sgrid_var
