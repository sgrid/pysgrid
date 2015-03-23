'''
Created on Mar 19, 2015

@author: ayan
'''
from read_netcdf import load_grid_from_nc_dataset, load_grid_from_nc_file


class SGrid(object):
    
    def __init__(self, nodes=None, faces=None, edges=None,
                 node_padding=None, face_padding=None,
                 edge_padding=None, grid_topology_vars=None):
        self.nodes = nodes
        self.faces = faces
        self.edges = edges
        self.node_padding = node_padding
        self.face_padding = face_padding
        self.edge_padding = edge_padding
        self.grid_topology_vars = grid_topology_vars
        
    @classmethod
    def from_nc_file(cls, nc_url, grid_topology_vars=None, load_data=False):
        grid = cls()
        sgrid_obj = load_grid_from_nc_file(nc_url, grid, grid_topology_vars, load_data)
        return sgrid_obj
    
    @classmethod
    def from_nc_dataset(cls, nc_dataset, grid_topology_vars=None, load_data=False):
        grid = cls()
        sgrid_obj = load_grid_from_nc_dataset(nc_dataset, grid, grid_topology_vars, load_data)
        return sgrid_obj
        