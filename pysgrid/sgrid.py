'''
Created on Mar 19, 2015

@author: ayan
'''
from read_netcdf import load_grid_from_nc_dataset, load_grid_from_nc_file


class SGrid(object):
    
    def __init__(self, nodes=None, faces=None, edges=None,
                 node_padding=None, face_padding=None,
                 edge_padding=None, grid_topology_vars=None):
        self._nodes = nodes
        self._faces = faces
        self._edges = edges
        self._node_padding = node_padding
        self._face_padding = face_padding
        self._edge_padding = edge_padding
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