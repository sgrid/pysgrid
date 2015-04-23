'''
Created on Apr 15, 2015

@author: ayan
'''
from .utils import determine_variable_slicing


class SGridVariable(object):
    """
    Object of variables found and inferred 
    from an SGRID compliant dataset.
    
    """
    def __init__(self, 
                 variable=None, 
                 grid=None, 
                 axes=None, 
                 center_slicing=None,
                 node_slicing=None, 
                 dimensions=None, 
                 dtype=None, 
                 location=None):
        self._variable = variable
        self._grid = grid
        self._axes = axes
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
        sgrid_var = cls(variable=variable,
                        grid=grid,
                        axes=None,
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
    def axes(self):
        return self._axes
        
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