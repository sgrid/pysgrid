'''
Created on Mar 19, 2015

@author: ayan
'''
import netCDF4 as nc4
from read_netcdf import load_grid_from_nc_dataset, load_grid_from_nc_file


class SGrid(object):
    
    padding_slices = {'both': (1, -1),
                      'none': (None, None),
                      'low': (1, None),
                      'high': (None, 1)
                      }
    
    def __init__(self, nodes=None, faces=None, edges=None,
                 node_padding=None, face_padding=None,
                 edge_1_padding=None, edge_2_padding=None,
                 vertical_padding=None, grid_topology_vars=None,
                 grid_cell_center_vars=None, grid_cell_center_lon=None,
                 grid_cell_center_lat=None, grid_times=None, variables=None):
        self._nodes = nodes
        self._faces = faces
        self._edges = edges
        self._node_padding = node_padding
        self._face_padding = face_padding
        self._edge_1_padding = edge_1_padding
        self._edge_2_padding = edge_2_padding
        self._vertical_padding = vertical_padding
        self._grid_topology_vars = grid_topology_vars
        self._grid_cell_center_vars = grid_cell_center_vars  # (lat, lon)
        self._grid_cell_center_lon = grid_cell_center_lon
        self._grid_cell_center_lat = grid_cell_center_lat
        self._grid_times = grid_times
        self._variables = variables
        
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
        return self._variables
    @variables.setter
    def variables(self, dataset_variables):
        self._variables = dataset_variables
        
    @property
    def grid_cell_center_vars(self):
        return self._grid_cell_center_vars
    @grid_cell_center_vars.setter
    def grid_cell_center_vars(self, grid_cell_center_variables):
        self._grid_cell_center_vars = grid_cell_center_variables
        
    @property
    def face_padding(self):
        return self._face_padding
    @face_padding.setter
    def face_padding(self, f_padding):
        self._face_padding = f_padding
        
    @property
    def edge_1_padding(self):
        return self._edge_1_padding
    @edge_1_padding.setter
    def edge_1_padding(self, e1_padding):
        self._edge_1_padding = e1_padding
        
    @property
    def edge_2_padding(self):
        return self._edge_2_padding
    @edge_2_padding.setter
    def edge_2_padding(self, e2_padding):
        self._edge_2_padding = e2_padding
    
    @property
    def vertical_padding(self):
        return self._vertical_padding
    @vertical_padding.setter
    def vertical_padding(self, vert_padding):
        self._vertical_padding = vert_padding
        
    @property
    def grid_cell_center_lon(self):
        return self._grid_cell_center_lon
    @grid_cell_center_lon.setter
    def grid_cell_center_lon(self, grid_cell_center_lon):
        self._grid_cell_center_lon = grid_cell_center_lon
        
    @property
    def grid_cell_center_lat(self):
        return self._grid_cell_center_lat
    @grid_cell_center_lat.setter
    def grid_cell_center_lat(self, grid_cell_center_lat):
        self._grid_cell_center_lat = grid_cell_center_lat
        
    @property
    def nodes(self):
        """
        return the vertices of the grid
        """
        return self._nodes
    @nodes.setter
    def nodes(self, nodes):
        self._nodes = nodes
        
    @property
    def grid_times(self):
        return self._grid_times
    @grid_times.setter
    def grid_times(self, grid_times):
        self._grid_times = grid_times
        
    def _define_face_padding_summary(self):
        all_padding = self._face_padding + self._edge_1_padding + self._edge_2_padding + self._vertical_padding
        padding_summary = []
        for padding_datum in all_padding:
            dim = padding_datum.face_dim
            padding_val = padding_datum.padding
            pad_short = (dim, padding_val)
            padding_summary.append(pad_short)
        return padding_summary
    
    def _set_property(self, name, value):
        prop_name = '_{0}'.format(name)
        setattr(self, prop_name, value)
        
    def _get_property(self, name):
        prop_name = '_{0}'.format(name)
        return getattr(self, prop_name)
    
    def add_property(self, name, value):
        fget = lambda self: self._get_property(name)
        fset = lambda self, value: self._set_property(name, value)
        # add property/attribute to the object with getting and setting functions
        setattr(self.__class__, name, property(fget, fset))
        prop_name = '_{0}'.format(name)
        # set the value of the property that was just created
        setattr(self, prop_name, value)
        
    def save_as_netcdf(self, filepath):
        with nc4.Dataset(filepath, 'w') as nclocal:
            time_dim_size = self._grid_time.shape[0]
            nclocal.createDimension('time', time_dim_size)