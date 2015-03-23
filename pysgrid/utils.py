'''
Created on Mar 23, 2015

@author: ayan
'''
import re
from collections import namedtuple


GridPadding = namedtuple('GridPadding', ['mesh_topology_var',  # the variable containing the padding information
                                         'dim_name',  # the topology attribute
                                         'dim_var',  # node dimension within the topology attribute
                                         'padding'  # padding type for the node dimension
                                         ]
                         )


class ParsePadding(object):
    """
    Parse out the padding types from
    variables with a cf_role of 'grid_topology'.
    
    """
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
            raw_dim_name, raw_dim_var, raw_padding_var = padding_match
            dim_name = raw_dim_name.split(':')[0]
            dim_var = raw_dim_var
            cleaned_padding_var = re.sub('[\(\)]', '', raw_padding_var)  # remove parentheses
            padding_type = cleaned_padding_var.split(':')[1].strip()  # get the padding value and remove spaces
            grid_padding = GridPadding(mesh_topology_var=self.mesh_topology_var,
                                       dim_name=dim_name,
                                       dim_var=dim_var,
                                       padding=padding_type
                                       )
            padding_type_list.append(grid_padding)
        if len(padding_type_list) > 0:
            final_padding_types = padding_type_list
        else:
            final_padding_types = None
        return final_padding_types
        
        