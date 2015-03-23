'''
Created on Mar 23, 2015

@author: ayan
'''
import re
from collections import namedtuple


GridPadding = namedtuple('GridPadding', ['mesh_topology_var', 
                                         'dim_name', 
                                         'dim_var', 
                                         'padding'
                                         ]
                         )


class ParsePadding(object):
    
    def __init__(self, mesh_topology_var=None):
        self.mesh_topology_var = mesh_topology_var

    def parse_padding(self, padding_str):
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
        
        