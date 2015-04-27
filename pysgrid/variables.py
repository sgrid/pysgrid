'''
Created on Apr 15, 2015

@author: ayan
'''
from .read_netcdf import parse_axes
from .utils import determine_variable_slicing, infer_avg_axes


class SGridVariable(object):
    """
    Object of variables found and inferred 
    from an SGRID compliant dataset.
    
    """
    def __init__(self, 
                 variable=None, 
                 grid=None, 
                 x_axis=None,
                 y_axis=None,
                 z_axis=None, 
                 center_slicing=None,
                 center_axis=None,
                 node_slicing=None,
                 node_axis=None, 
                 dimensions=None, 
                 dtype=None, 
                 location=None):
        self.variable = variable
        self.grid = grid
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.z_axis = z_axis
        self.center_slicing = center_slicing
        self.center_axis = center_axis
        self.node_slicing = node_slicing
        self.node_axis = node_axis
        self.dimensions = dimensions
        self.dtype = dtype
        self.location = location
        
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
                                                    method='center'
                                                    )
        dimensions = nc_var_obj.dimensions
        dtype = nc_var_obj.dtype
        try:
            location = nc_var_obj.location
        except AttributeError:
            location = None
        try:
            axes = nc_var_obj.axes
        except AttributeError:
            x_axis = None
            y_axis = None
            z_axis = None
        else:
            x_axis, y_axis, z_axis = parse_axes(axes)
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
                        location=location
                        )
        return sgrid_var