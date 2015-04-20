'''
Created on Mar 19, 2015

@author: ayan
'''
import netCDF4 as nc4
from .read_netcdf import load_grid_from_nc_dataset, load_grid_from_nc_file


class SGrid(object):
    
    padding_slices = {'both': (1, -1),
                      'none': (None, None),
                      'low': (1, None),
                      'high': (None, 1)
                      }
    
    def __init__(self):
        # general attributes
        self._nodes = None
        self._centers = None
        self._edges = None
        self._node_padding = None
        self._edge1_padding = None
        self._edge2_padding = None
        self._grid_topology_vars = None
        self._grid_times = None
        self._variables = None
        self._grid_variables = None
        self._topology_dimension = None
        self._dimensions = None
        self._node_coordinates = None
        self._edge1_coordinates = None
        self._edge2_coordinates = None
        self._angles = None
        # attributes specific to topology_dimension 2
        self._faces = None
        self._face_padding = None
        self._face_coordinates = None
        self._vertical_padding = None
        # attributes specific to topology_dimension 3
        self._volume_padding = None
        self._volume_coordinates = None
        self._edge3_padding = None
        self._edge3_coordinates = None
        self._face1_padding = None
        self._face1_coordinates = None
        self._face2_padding = None
        self._face2_coordinates = None
        self._face3_padding = None
        self._face3_coordinates = None
        # attributes for the verbatim padding text
        # used for saving SGrid as a netCDF file
        self._node_dimensions = None
        self._face_dimensions = None
        self._volume_dimensions = None
        self._vertical_dimensions = None
        self._edge1_dimensions = None
        self._edge2_dimensions = None
        self._edge3_dimensions = None
        self._face1_dimensions = None
        self._face2_dimensions = None
        self._face3_dimensions = None
        
    @classmethod
    def from_nc_file(cls, nc_url, grid_topology_vars=None, load_data=False):
        grid = cls()
        load_grid_from_nc_file(nc_url, grid, 
                               grid_topology_vars, load_data
                               )
        return grid
    
    @classmethod
    def from_nc_dataset(cls, nc_dataset, grid_topology_vars=None, load_data=False):
        grid = cls()
        load_grid_from_nc_dataset(nc_dataset, grid, 
                                  grid_topology_vars, load_data
                                  )
        return grid
    
    @property
    def grid_topology_vars(self):
        return self._grid_topology_vars
    
    @grid_topology_vars.setter
    def grid_topology_vars(self, grid_topo_vars):
        self._grid_topology_vars = grid_topo_vars
        
    @property
    def variables(self):
        """
        Return a list of variables
        
        """
        return self._variables
    
    @variables.setter
    def variables(self, dataset_variables):
        self._variables = dataset_variables
        
    @property
    def grid_variables(self):
        """
        Return a list of variables
        with a grid attribute.
        
        """
        return self._grid_variables
    
    @grid_variables.setter
    def grid_variables(self, dataset_grid_variables):
        self._grid_variables = dataset_grid_variables
        
    @property
    def topology_dimension(self):
        return self._topology_dimension
    
    @topology_dimension.setter
    def topology_dimension(self, dataset_topology_dimension):
        self._topology_dimension = dataset_topology_dimension
        
    @property
    def dimensions(self):
        """
        Return list of tuples containing
        dimension name and size.
        
        """
        return self._dimensions
    
    @dimensions.setter
    def dimensions(self, dataset_dims):
        self._dimensions = dataset_dims
    
    # start topology_dimension = 2   
    @property
    def face_padding(self):
        return self._face_padding
    
    @face_padding.setter
    def face_padding(self, f_padding):
        self._face_padding = f_padding
        
    @face_padding.deleter
    def face_padding(self):
        del self._face_padding
        
    @property
    def face_coordinates(self):
        return self._face_coordinates
    
    @face_coordinates.setter
    def face_coordinates(self, dataset_face_coordinates):
        self._face_coordinates = dataset_face_coordinates
        
    @face_coordinates.deleter
    def face_coordinates(self):
        del self._face_coordinates
        
    @property
    def face_dimensions(self):
        return self._face_dimensions
    
    @face_dimensions.setter
    def face_dimensions(self, face_dim):
        self._face_dimensions = face_dim
        
    @face_dimensions.deleter
    def face_dimensions(self):
        del self._face_dimensions
        
    @property
    def vertical_padding(self):
        return self._vertical_padding
    
    @vertical_padding.setter
    def vertical_padding(self, vert_padding):
        self._vertical_padding = vert_padding
        
    @vertical_padding.deleter
    def vertical_padding(self):
        del self._vertical_padding
        
    @property
    def vertical_dimensions(self):
        return self._vertical_dimensions
    
    @vertical_dimensions.setter
    def vertical_dimensions(self, vertical_dim):
        self._vertical_dimensions = vertical_dim
        
    @vertical_dimensions.deleter
    def vertical_dimensions(self):
        del self._vertical_dimensions
    # end topology_dimension = 2
    
    # start topology_dimension = 3
    @property
    def volume_padding(self):
        return self._volume_padding
    
    @volume_padding.setter
    def volume_padding(self, vol_padding):
        self._volume_padding = vol_padding
        
    @volume_padding.deleter
    def volume_padding(self):
        del self._volume_padding
        
    @property
    def volume_dimensions(self):
        return self._volume_dimensions
    
    @volume_dimensions.setter
    def volume_dimensions(self, volume_dims):
        self._volume_dimensions = volume_dims
        
    @volume_dimensions.deleter
    def volume_dimensions(self):
        del self._volume_dimensions
        
    @property
    def volume_coordinates(self):
        return self._volume_coordinates
    
    @volume_coordinates.setter
    def volume_coordinates(self, vol_coordinates):
        self._volume_coordinates = vol_coordinates
        
    @volume_coordinates.deleter
    def volume_coordinates(self):
        del self._volume_coordinates
        
    @property
    def face1_padding(self):
        return self._face1_padding
    
    @face1_padding.setter
    def face1_padding(self, f1_padding):
        self._face1_padding = f1_padding
        
    @face1_padding.deleter
    def face1_padding(self):
        del self._face1_padding
        
    @property
    def face1_coordinates(self):
        return self._face1_coordinates
    
    @face1_coordinates.setter
    def face1_coordinates(self, f1_coordinates):
        self._face1_coordinates = f1_coordinates
        
    @face1_coordinates.deleter
    def face1_coordinates(self):
        del self._face1_coordinates
        
    @property
    def face1_dimensions(self):
        return self._face1_dimensions
    
    @face1_dimensions.setter
    def face1_dimensions(self, f1_dims):
        self._face1_dimensions = f1_dims
        
    @face1_dimensions.deleter
    def face1_dimensions(self):
        del self._face1_dimensions
        
    @property
    def face2_padding(self):
        return self._face2_padding
    
    @face2_padding.setter
    def face2_padding(self, f2_padding):
        self._face2_padding = f2_padding
        
    @face2_padding.deleter
    def face2_padding(self):
        del self._face2_padding
        
    @property
    def face2_coordinates(self):
        return self._face2_coordinates
    
    @face2_coordinates.setter
    def face2_coordinates(self, f2_coordinates):
        self._face2_coordinates = f2_coordinates
        
    @face2_coordinates.deleter
    def face2_coordinates(self):
        del self._face2_coordinates
    
    @property    
    def face2_dimensions(self):
        return self._face2_dimensions
    
    @face2_dimensions.setter
    def face2_dimensions(self, f2_dims):
        self._face2_dimensions = f2_dims
        
    @face2_dimensions.deleter
    def face2_dimensions(self):
        del self._face2_dimensions
        
    @property
    def face3_padding(self):
        return self._face3_padding
    
    @face3_padding.setter
    def face3_padding(self, f3_padding):
        self._face3_padding = f3_padding
        
    @face3_padding.deleter
    def face3_padding(self):
        del self._face3_padding
        
    @property
    def face3_coordinates(self):
        return self._face3_coordinates
    
    @face3_coordinates.setter
    def face3_coordinates(self, f3_coordinates):
        self._face3_coordinates = f3_coordinates
        
    @face3_coordinates.deleter
    def face3_coordinates(self):
        del self._face3_coordinates
        
    @property
    def face3_dimensions(self):
        return self._face3_dimensions
    
    @face3_dimensions.setter
    def face3_dimensions(self, f3_dims):
        self._face3_dimensions = f3_dims
        
    @face3_dimensions.deleter
    def face3_dimensions(self):
        del self._face3_dimensions
        
    @property
    def edge3_padding(self):
        return self._edge3_padding
    
    @edge3_padding.setter
    def edge3_padding(self, e3_padding):
        self._edge3_padding = e3_padding
        
    @edge3_padding.deleter
    def edge3_padding(self):
        del self._edge3_padding
    
    @property
    def edge3_coordinates(self):
        return self._edge3_coordinates
    
    @edge3_coordinates.setter
    def edge3_coordinates(self, e3_coordinates):
        self._edge3_coordinates = e3_coordinates
        
    @edge3_coordinates.deleter
    def edge3_coordinates(self):
        del self._edge3_coordinates
        
    @property
    def edge3_dimensions(self):
        return self._edge3_dimensions
    
    @edge3_dimensions.setter
    def edge3_dimensions(self, e3_dims):
        self._edge3_dimensions = e3_dims
        
    @edge3_dimensions.deleter
    def edge3_dimensions(self):
        del self._edge3_dimensions
    # end topology_dimension = 3
        
    @property
    def nodes(self):
        """
        return the vertices of the grid as arrays 
        of lon, lat pairs.
        
        """
        return self._nodes
        
    @nodes.setter
    def nodes(self, nodes):
        self._nodes = nodes
        
    @property
    def node_dimensions(self):
        return self._node_dimensions
    
    @node_dimensions.setter
    def node_dimensions(self, node_dim):
        self._node_dimensions = node_dim
        
    @property
    def node_coordinates(self):
        return self._node_coordinates
    
    @node_coordinates.setter
    def node_coordinates(self, dataset_node_coordinates):
        self._node_coordinates = dataset_node_coordinates
        
    @property
    def edge1_coordinates(self):
        return self._edge1_coordinates
    
    @edge1_coordinates.setter
    def edge1_coordinates(self, dataset_edge1_coordinates):
        self._edge1_coordinates = dataset_edge1_coordinates
        
    @property
    def edge1_padding(self):
        return self._edge1_padding
    
    @edge1_padding.setter
    def edge1_padding(self, e1_padding):
        self._edge1_padding = e1_padding
        
    @property
    def edge1_dimension(self):
        return self._edge1_dimensions
    
    @edge1_dimension.setter
    def edge1_dimension(self, e1_dim):
        self._edge1_dimensions = e1_dim
        
    @property
    def edge2_coordinates(self):
        return self._edge2_coordinates
    
    @edge2_coordinates.setter
    def edge2_coordinates(self, dataset_edge2_coordinates):
        self._edge2_coordinates = dataset_edge2_coordinates
        
    @property
    def edge2_padding(self):
        return self._edge2_padding
    
    @edge2_padding.setter
    def edge2_padding(self, e2_padding):
        self._edge2_padding = e2_padding
        
    @property
    def edge2_dimensions(self):
        return self._edge2_dimensions
    
    @edge2_dimensions.setter
    def edge2_dimensions(self, e2_dim):
        self._edge2_dimensions = e2_dim
        
    @property
    def angles(self):
        return self._angles
    
    @angles.setter
    def angles(self, dataset_angles):
        self._angles = dataset_angles
        
    @property
    def centers(self):
        """
        return the coordinates of the grid centers
        as arrays of lon, lat pairs.
        
        """
        return self._centers
    
    @centers.setter
    def centers(self, centers):
        self._centers = centers
        
    @property
    def grid_times(self):
        return self._grid_times
    
    @grid_times.setter
    def grid_times(self, grid_times):
        self._grid_times = grid_times
        
    def _define_face_padding_summary(self):
        all_padding = []
        if self._volume_padding is not None:
            all_padding += self._volume_padding
        if self._face_padding is not None:
            all_padding += self._face_padding
        if self._face1_padding is not None:
            all_padding += self._face1_padding
        if self._face2_padding is not None:
            all_padding += self._face2_padding
        if self._face3_padding is not None:
            all_padding += self._face3_padding
        if self._vertical_padding is not None:
            all_padding += self._vertical_padding
        if self._edge1_padding is not None:
            all_padding += self._edge1_padding
        if self._edge2_padding is not None:
            all_padding += self._edge2_padding
        if self._edge3_padding is not None:
            all_padding += self._edge3_padding
        padding_summary = []
        for padding_datum in all_padding:
            dim = padding_datum.dim
            sub_dim = padding_datum.sub_dim
            padding_val = padding_datum.padding
            pad_short = (dim, sub_dim, padding_val)
            padding_summary.append(pad_short)
        return padding_summary
    
    def _set_property(self, name, value):
        prop_name = '_{0}'.format(name)
        setattr(self, prop_name, value)
        
    def _get_property(self, name):
        prop_name = '_{0}'.format(name)
        return getattr(self, prop_name)
    
    def add_property(self, name, value):
        """
        Method to dynamically add attributes
        to a class.
        
        """
        fget = lambda self: self._get_property(name)
        fset = lambda self, value: self._set_property(name, value)
        # add property/attribute to the object with getting and setting functions
        setattr(self.__class__, name, property(fget, fset))
        prop_name = '_{0}'.format(name)
        # set the value of the property that was just created
        setattr(self, prop_name, value)
        
    def save_as_netcdf(self, filepath):
        with nc4.Dataset(filepath, 'w') as nclocal:
            grid_var = self._grid_topology_vars
            # create dimensions
            for grid_dim in self._dimensions:
                dim_name, dim_size = grid_dim
                nclocal.createDimension(dim_name, dim_size)
            # create variables
            center_lon, center_lat = self._face_coordinates
            center_lon_obj = getattr(self, center_lon)
            center_lat_obj = getattr(self, center_lat)
            node_lon, node_lat = self._node_coordinates
            node_lon_obj = getattr(self, node_lon)
            node_lat_obj = getattr(self, node_lat)
            grid_center_lon = nclocal.createVariable(center_lon, 
                                                     center_lon_obj.dtype, 
                                                     center_lon_obj.dimensions
                                                     )
            grid_center_lat = nclocal.createVariable(center_lat, 
                                                     center_lat_obj.dtype, 
                                                     center_lat_obj.dimensions
                                                     )
            grid_node_lon = nclocal.createVariable(node_lon, 
                                                   node_lon_obj.dtype, 
                                                   node_lon_obj.dimensions
                                                   )
            grid_node_lat = nclocal.createVariable(node_lat, 
                                                   node_lat_obj.dtype, 
                                                   node_lat_obj.dimensions
                                                   )
            grid_var_obj = getattr(self, grid_var)
            grid_vars = nclocal.createVariable(grid_var, grid_var_obj.dtype)
            # not the most robust here... time and angle are hard-coded
            # need to address this
            time_obj = getattr(self, 'time')
            grid_time = nclocal.createVariable('time', 
                                               time_obj.dtype, 
                                               time_obj.dimensions
                                               )
            
            if hasattr(self, 'angle'):
                angle_obj = getattr(self, 'angle', None)
                grid_angle = nclocal.createVariable('angle', 
                                                    angle_obj.dtype, 
                                                    angle_obj.dimensions
                                                    )
                if self._angles is not None:
                    grid_angle[:] = self._angles[:]
            # save the grid variables with attributes
            for dataset_variable in self._variables:
                dataset_var_obj = getattr(self, dataset_variable)
                if dataset_var_obj.grid is not None:
                    dataset_grid_var = nclocal.createVariable(dataset_variable, 
                                                              dataset_var_obj.dtype, 
                                                              dataset_var_obj.dimensions
                                                              )
                    dataset_grid_var.grid = grid_var
            # add attributes to the variables
            grid_vars.cf_role = 'grid_topology'
            grid_vars.topology_dimension = self._topology_dimension
            grid_vars.node_dimensions = self._node_dimensions
            if self._face_dimensions is not None:
                grid_vars.face_dimensions = self._face_dimensions
            if hasattr(self, 'volume_dimensions') and self._volume_dimensions is not None:
                grid_vars.volume_dimensions = self._volume_dimensions
            if self._edge1_dimensions is not None:
                grid_vars.edge1_dimensions = self._edge1_dimensions
            if self._edge2_dimensions is not None:
                grid_vars.edge2_dimensions = self._edge2_dimensions
            if self._vertical_dimensions is not None:
                grid_vars.vertical_dimensions = self._vertical_dimensions
            if self._face_coordinates is not None:
                grid_vars.face_coordinates = ' '.join(self._face_coordinates)
            if hasattr(self, 'volume_coordinates') and self._volume_coordinates is not None:
                grid_vars.volume_coordinates = ' '.join(self._volume_coordinates)
            if self._node_coordinates is not None:
                grid_vars.node_coordinates = ' '.join(self._node_coordinates)
            if self._edge1_coordinates is not None:
                grid_vars.edge1_coordinates = ' '.join(self._edge1_coordinates)
            if self._edge2_coordinates is not None:
                grid_vars.edge2_coordinates = ' '.join(self._edge2_coordinates)
            # populate variables with data
            grid_time[:] = self._grid_times[:]
            grid_center_lon[:, :] = self._centers[:, :, 0]
            grid_center_lat[:, :] = self._centers[:, :, 1] 
            grid_node_lon[:, :] = self._nodes[:, :, 0]
            grid_node_lat[:, :] = self._nodes[:, :, 1]