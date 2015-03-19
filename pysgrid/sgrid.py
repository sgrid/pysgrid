'''
Created on Mar 19, 2015

@author: ayan
'''


class SGrid(object):
    
    def __init__(self, nodes=None, faces=None, edges=None,
                 node_padding=None, face_padding=None,
                 edge_padding=None, grid_topology=None):
        self.nodes = nodes
        self.faces = faces
        self.edges = edges
        self.node_padding = node_padding
        self.face_padding = face_padding
        self.edge_padding = edge_padding
        self.grid_topology = grid_topology
        
    