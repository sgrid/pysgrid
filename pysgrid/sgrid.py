'''
Created on Mar 19, 2015

@author: ayan
'''
from read_netcdf import load_grid_from_nc_dataset, load_grid_from_nc_file


class SGrid(object):
    
    def __init__(self, nodes=None, faces=None, edges=None,
                 node_padding=None, face_padding=None,
                 edge_1_padding=None, edge_2_padding=None,
                 vertical_padding=None, grid_topology_vars=None):
        self._nodes = nodes
        self._faces = faces
        self._edges = edges
        self._node_padding = node_padding
        self._face_padding = face_padding
        self._edge_1_padding = edge_1_padding
        self._edge_2_padding = edge_2_padding
        self._vertical_padding = vertical_padding
        self._grid_topology_vars = grid_topology_vars
        
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