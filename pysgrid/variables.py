'''
Created on Apr 15, 2015

@author: ayan
'''


def load_variable(nc_var_obj, sgrid_var_obj):
    """
    Create a SGridVariable object from
    a netCDF4.Dataset.variable object.
    
    :param netCDF4.variable nc_var_obj: a netCDF4 variable
    :param SGridVariable sgrid_var_obj: the target SGridVariable object
    :return: SGridVariable object with attributes from nc_var_obj
    :rtype: SGridVariable
    
    """
    nc_var_name = nc_var_obj.name
    nc_var_dims = nc_var_obj.dimensions
    nc_var_dtype = nc_var_obj.dtype
    try:
        nc_var_grid = nc_var_obj.grid
    except AttributeError:
        nc_var_grid = None
    try:
        nc_var_location = nc_var_obj.location
    except AttributeError:
        nc_var_location = None
    sgrid_var_obj.variable = nc_var_name
    sgrid_var_obj.grid = nc_var_grid
    sgrid_var_obj.dimensions = nc_var_dims
    sgrid_var_obj.dtype = nc_var_dtype
    sgrid_var_obj.location = nc_var_location
    return sgrid_var_obj


class SGridVariable(object):
    """
    Object of variables found and inferred 
    from an SGRID compliant dataset.
    
    """
    def __init__(self, variable=None, grid=None, 
                 axes=None, slicing=None,
                 dimensions=None, dtype=None,
                 location=None):
        self._variable = variable
        self._grid = grid
        self._axes = axes
        self._slicing = slicing
        self._dimensions = dimensions
        self._dtype = dtype
        self._location = location
        
    @classmethod
    def create_variable(cls, nc_var_obj):
        sgrid_var = cls()
        load_variable(nc_var_obj, sgrid_var)
        return sgrid_var
        
    @property
    def variable(self):
        return self._variable
    
    @variable.setter
    def variable(self, variable_name):
        self._variable = variable_name
        
    @property
    def grid(self):
        return self._grid
    
    @grid.setter
    def grid(self, grid_name):
        self._grid = grid_name
        
    @property
    def axes(self):
        return self._axes
    
    @axes.setter
    def axes(self, variable_axes):
        self._axes = variable_axes
        
    @property
    def slicing(self):
        return self._slicing
    
    @slicing.setter
    def slicing(self, variable_slicing):
        self._slicing = variable_slicing
        
    @property
    def dimensions(self):
        return self._dimensions
    
    @dimensions.setter
    def dimensions(self, variable_dimensions):
        self._dimensions = variable_dimensions
        
    @property
    def dtype(self):
        return self._dtype
    
    @dtype.setter
    def dtype(self, variable_dtype):
        self._dtype = variable_dtype
        
    @property
    def location(self):
        return self._location
    
    @location.setter
    def location(self, variable_location):
        self._location = variable_location