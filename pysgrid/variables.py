'''
Created on Apr 15, 2015

@author: ayan
'''
from .read_netcdf import parse_axes
from .utils import determine_variable_slicing


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
                 node_slicing=None, 
                 dimensions=None, 
                 dtype=None, 
                 location=None):
        self._variable = variable
        self._grid = grid
        self._x_axis = x_axis
        self._y_axis = y_axis
        self._z_axis = z_axis
        self._center_slicing = center_slicing
        self._node_slicing = node_slicing
        self._dimensions = dimensions
        self._dtype = dtype
        self._location = location
        
    @classmethod
    def create_variable(cls, nc_var_obj, sgrid_obj):
        variable = nc_var_obj.name
        try:
            grid = nc_var_obj.grid
        except:
            grid = None
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
                        node_slicing=None,
                        dimensions=dimensions,
                        dtype=dtype,
                        location=location
                        )
        return sgrid_var
        
    @property
    def variable(self):
        return self._variable
        
    @property
    def grid(self):
        return self._grid
    
    @property
    def x_axis(self):
        return self._x_axis
    
    @property
    def y_axis(self):
        return self._y_axis
    
    @property
    def z_axis(self):
        return self._z_axis
        
    @property
    def center_slicing(self):
        """
        Get the slicing necessary
        when averaging to grid cell
        centers.
        
        """
        return self._center_slicing
        
    @property
    def node_slicing(self):
        """
        Get the slicing necessary
        when averaging to grid nodes
        
        """
        return self._node_slicing
        
    @property
    def dimensions(self):
        return self._dimensions
        
    @property
    def dtype(self):
        return self._dtype
        
    @property
    def location(self):
        return self._location