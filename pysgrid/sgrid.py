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
                 grid_cell_center_lat=None, grid_times=None):
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
    def grid_times(self):
        return self._grid_times
    @grid_times.setter
    def grid_times(self, grid_times):
        self._grid_times = grid_times
        
    def save_as_netcdf(self, filepath):
        with nc4.Dataset(filepath, 'w') as nclocal:
            time_dim_size = self._grid_time.shape[0]
            nclocal.createDimension('time', time_dim_size)