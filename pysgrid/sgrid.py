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
    
    def __init__(self, nodes=None, centers=None, faces=None, 
                 edges=None, node_padding=None, face_padding=None,
                 edge_1_padding=None, edge_2_padding=None,
                 vertical_padding=None, grid_topology_vars=None,
                 grid_cell_center_vars=None, grid_times=None, 
                 variables=None, dimensions=None, face_coordinates=None,
                 node_coordinates=None, edge_1_coordinates=None,
                 edge_2_coordinates=None, angles=None,
                 node_dim=None, face_dim=None,
                 vertical_dim=None, edge_1_dim=None,
                 edge_2_dim=None, grid_variables=None):
        self._nodes = nodes
        self._centers = centers
        self._faces = faces
        self._edges = edges
        self._node_padding = node_padding
        self._face_padding = face_padding
        self._edge_1_padding = edge_1_padding
        self._edge_2_padding = edge_2_padding
        self._vertical_padding = vertical_padding
        self._grid_topology_vars = grid_topology_vars
        self._grid_times = grid_times
        self._variables = variables
        self._grid_variables = grid_variables
        self._dimensions = dimensions
        self._face_coordinates = face_coordinates
        self._node_coordinates = node_coordinates
        self._edge_1_coordinates = edge_1_coordinates
        self._edge_2_coordinates = edge_2_coordinates
        self._angles = angles
        # attributes for the verbatim padding text
        self._node_dimensions = node_dim
        self._face_dimensions = face_dim
        self._vertical_dimensions = vertical_dim
        self._edge_1_dimensions = edge_1_dim
        self._edge_2_dimensions = edge_2_dim
        
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
    def dimensions(self):
        """
        Return list of tuples containing
        dimension name and size.
        """
        return self._dimensions
    
    @dimensions.setter
    def dimensions(self, dataset_dims):
        self._dimensions = dataset_dims
        
    @property
    def face_padding(self):
        return self._face_padding
    
    @face_padding.setter
    def face_padding(self, f_padding):
        self._face_padding = f_padding
        
    @property
    def face_coordinates(self):
        return self._face_coordinates
    
    @face_coordinates.setter
    def face_coordinates(self, dataset_face_coordinates):
        self._face_coordinates = dataset_face_coordinates
        
    @property
    def face_dimensions(self):
        return self._face_dimensions
    
    @face_dimensions.setter
    def face_dimensions(self, face_dim):
        self._face_dimensions = face_dim
        
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
    def edge_1_coordinates(self):
        return self._edge_1_coordinates
    
    @edge_1_coordinates.setter
    def edge_1_coordinates(self, dataset_edge_1_coordinates):
        self._edge_1_coordinates = dataset_edge_1_coordinates
        
    @property
    def edge_1_padding(self):
        return self._edge_1_padding
    
    @edge_1_padding.setter
    def edge_1_padding(self, e1_padding):
        self._edge_1_padding = e1_padding
        
    @property
    def edge_1_dimension(self):
        return self._edge_1_dimensions
    
    @edge_1_dimension.setter
    def edge_1_dimension(self, e1_dim):
        self._edge_1_dimensions = e1_dim
        
    @property
    def edge_2_coordinates(self):
        return self._edge_2_coordinates
    
    @edge_2_coordinates.setter
    def edge_2_coordinates(self, dataset_edge_2_coordinates):
        self._edge_2_coordinates = dataset_edge_2_coordinates
        
    @property
    def edge_2_padding(self):
        return self._edge_2_padding
    
    @edge_2_padding.setter
    def edge_2_padding(self, e2_padding):
        self._edge_2_padding = e2_padding
        
    @property
    def edge_2_dimensions(self):
        return self._edge_2_dimensions
    
    @edge_2_dimensions.setter
    def edge_2_dimensions(self, e2_dim):
        self._edge_2_dimensions = e2_dim
        
    @property
    def angles(self):
        return self._angles
    
    @angles.setter
    def angles(self, dataset_angles):
        self._angles = dataset_angles
    
    @property
    def vertical_padding(self):
        return self._vertical_padding
    
    @vertical_padding.setter
    def vertical_padding(self, vert_padding):
        self._vertical_padding = vert_padding
        
    @property
    def vertical_dimensions(self):
        return self._vertical_dimensions
    
    @vertical_dimensions.setter
    def vertical_dimensions(self, vertical_dim):
        self._vertical_dimensions = vertical_dim
        
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
        all_padding = self._face_padding + self._vertical_padding
        if self._edge_1_padding is not None:
            all_padding += self._edge_1_padding
        if self._edge_2_padding is not None:
            all_padding += self._edge_2_padding
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
            grid_var = self._grid_topology_vars[0]
            # create dimensions
            for grid_dim in self._dimensions:
                dim_name, dim_size = grid_dim
                nclocal.createDimension(dim_name, dim_size)
            # create variables
            center_lon, center_lat = self._face_coordinates
            center_lon_data = next((nc_var for nc_var in self._variables if nc_var[0] == center_lon))
            center_lat_data = next((nc_var for nc_var in self._variables if nc_var[0] == center_lat))
            node_lon, node_lat = self._node_coordinates
            node_lon_data = next((nc_var for nc_var in self._variables if nc_var[0] == node_lon))
            node_lat_data = next((nc_var for nc_var in self._variables if nc_var[0] == node_lat))
            grid_center_lon = nclocal.createVariable(center_lon, 'f4', center_lon_data[-1])
            grid_center_lat = nclocal.createVariable(center_lat, 'f4', center_lat_data[-1])
            grid_node_lon = nclocal.createVariable(node_lon, 'f4', node_lon_data[-1])
            grid_node_lat = nclocal.createVariable(node_lat, 'f4', node_lat_data[-1])
            grid_vars = nclocal.createVariable(grid_var, 'i2')
            grid_time = nclocal.createVariable('time', 'f8', ('time',))
            grid_angle = nclocal.createVariable('angle', 'f8', center_lat_data[-1])
            for dataset_variable in self._variables:
                dataset_variable_name = dataset_variable[0]
                if dataset_variable_name in self._grid_variables:
                    dataset_variable_dims = dataset_variable[2]
                    dataset_grid_var = nclocal.createVariable(dataset_variable_name, 'f4', dataset_variable_dims)
                    dataset_grid_var.grid = grid_var
            # add attributes to the variables
            grid_vars.cf_role = 'grid_topology'
            grid_vars.topology_dimension = 2
            grid_vars.node_dimensions = self._node_dimensions
            grid_vars.face_dimensions = self._face_dimensions
            if self._edge_1_dimensions is not None:
                grid_vars.edge1_dimensions = self._edge_1_dimensions
            if self._edge_2_dimensions is not None:
                grid_vars.edge2_dimensions = self._edge_2_dimensions
            if self._vertical_dimensions is not None:
                grid_vars.vertical_dimensions = self._vertical_dimensions
            if self._face_coordinates is not None:
                grid_vars.face_coordinates = ' '.join(self._face_coordinates)
            if self._node_coordinates is not None:
                grid_vars.node_coordinates = ' '.join(self._node_coordinates)
            if self._edge_1_coordinates is not None:
                grid_vars.edge1_coordinates = ' '.join(self._edge_1_coordinates)
            if self._edge_2_coordinates is not None:
                grid_vars.edge2_coordinates = ' '.join(self._edge_2_coordinates)
            # populate variables with data
            grid_time[:] = self._grid_times[:]
            grid_center_lon[:, :] = self._centers[:, :, 0]
            grid_center_lat[:, :] = self._centers[:, :, 1] 
            grid_node_lon[:, :] = self._nodes[:, :, 0]
            grid_node_lat[:, :] = self._nodes[:, :, 1]
            if self._angles is not None:
                grid_angle[:] = self._angles[:]