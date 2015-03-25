'''
Created on Mar 23, 2015

@author: ayan
'''
import re
from collections import namedtuple
from custom_exceptions import CannotFindPadding

'''
    int grid;
      :long_name = "";
      :cf_role = "grid_topology";
      :topology_dimension = 2; // int
      :node_dimensions = "xi_psi eta_psi";
      :face_dimensions = "xi_rho: xi_psi (padding: both) eta_rho: eta_psi (padding: both)";
      :edge1_dimensions = "xi_u: xi_psi eta_u: eta_psi (padding: both)";
      :edge2_dimensions = "xi_v: xi_psi (padding: both) eta_v: eta_psi";
      :node_coordinates = "lon_psi lat_psi";
      :face_coordinates = "lon_rho lat_rho";
      :edge1_coordinates = "lon_u lat_u";
      :edge2_coordinates = "lon_v lat_v";
      :vertical_dimensions = "s_rho: s_w (padding: none)";
'''


GridPadding = namedtuple('GridPadding', ['mesh_topology_var',  # the variable containing the padding information
                                         'face_dim',  # the topology attribute
                                         'node_dim',  # node dimension within the topology attribute
                                         'padding'  # padding type for the node dimension
                                         ]
                         )


class ParsePadding(object):
    """
    Parse out the padding types from
    variables with a cf_role of 'grid_topology'.
    
    """
    padding_slices = {'both': (1, -1),
                      'none': (None, None),
                      'low': (1, None),
                      'high': (None, 1)
                      }
    
    def __init__(self, mesh_topology_var=None):
        self.mesh_topology_var = mesh_topology_var

    def parse_padding(self, padding_str):
        """
        Use regex expressions to break apart an
        attribute string containining padding types
        for each variable with a cf_role of 
        'grid_topology'.
        
        Padding information is returned within a named tuple
        for each node dimension of an edge, face, or vertical
        dimension. The named tuples have the following attributes:
        mesh_topology_var, dim_name, dim_var, and padding.
        Padding information is returned as a list
        of these named tuples.
        
        :param str padding_str: string containing padding types from a netCDF attribute
        :return: named tuples with padding information
        :rtype: list
        
        """
        p = re.compile('([a-zA-Z0-9_]+:) ([a-zA-Z0-9_]+) (\(padding: [a-zA-Z]+\))')
        padding_matches = p.findall(padding_str)
        padding_type_list = []
        for padding_match in padding_matches:
            raw_face_dim, raw_node_dim, raw_padding_var = padding_match
            face_dim = raw_face_dim.split(':')[0]
            node_dim = raw_node_dim
            cleaned_padding_var = re.sub('[\(\)]', '', raw_padding_var)  # remove parentheses
            padding_type = cleaned_padding_var.split(':')[1].strip()  # get the padding value and remove spaces
            grid_padding = GridPadding(mesh_topology_var=self.mesh_topology_var,
                                       face_dim=face_dim,
                                       node_dim=node_dim,
                                       padding=padding_type
                                       )
            padding_type_list.append(grid_padding)
        if len(padding_type_list) > 0:
            final_padding_types = padding_type_list
        else:
            final_padding_types = None
            raise CannotFindPadding
        return final_padding_types
    
    def define_recommended_slices(self, padding_str):
        if padding_str is not None:
            padding_data = self.parse_padding(padding_str)
            rec_slices = [] 
            for padding_datum in padding_data:
                face_dim = '{0}_padding_slice'.format(padding_datum.face_dim)
                padding_type = padding_datum.padding
                padding_slices = self.padding_slices[padding_type]
                rec_padding_slice = (face_dim, padding_slices)
                rec_slices.append(rec_padding_slice)
            return rec_slices
        else:
            raise CannotFindPadding
        