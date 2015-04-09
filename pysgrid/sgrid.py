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
                 variables=None, face_coordinates=None,
                 node_coordinates=None, edge_1_coordinates=None,
                 edge_2_coordinates=None):
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
        self._grid_cell_center_vars = grid_cell_center_vars  # (lat, lon)
        self._grid_times = grid_times
        self._variables = variables
        self._face_coordinates = face_coordinates
        self._node_coordinates = node_coordinates
        self._edge_1_coordinates = edge_1_coordinates
        self._edge_2_coordinates = edge_2_coordinates
        
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
    def face_coordinates(self):
        return self._face_coordinates
    
    @face_coordinates.setter
    def face_coordinates(self, dataset_face_coordinates):
        self._face_coordinates = dataset_face_coordinates
        
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
    def edge_2_coordinates(self):
        return self._edge_2_coordinates
    
    @edge_2_coordinates.setter
    def edge_2_coordinates(self, dataset_edge_2_coordinates):
        self._edge_2_coordinates = dataset_edge_2_coordinates
        
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
            nclocal.createDimension('time', len(self._grid_times))
            grid_x_center_dim = '{0}_x_center'.format(grid_var)
            nclocal.createDimension(grid_x_center_dim, self._centers.shape[1])
            grid_y_center_dim = '{0}_y_center'.format(grid_var)
            nclocal.createDimension(grid_y_center_dim, self._centers.shape[0])
            grid_x_node_dim = '{0}_x_node'.format(grid_var)
            nclocal.createDimension(grid_x_node_dim, self._nodes.shape[1])
            grid_y_node_dim = '{0}_y_node'.format(grid_var)
            nclocal.createDimension(grid_y_node_dim, self._nodes.shape[0])
            # create variables
            gc_lon_name = '{0}_center_lon'.format(grid_var)
            grid_center_lon = nclocal.createVariable(gc_lon_name, 'f4', (grid_y_center_dim, grid_x_center_dim))
            gc_lat_name = '{0}_center_lat'.format(grid_var)
            grid_center_lat = nclocal.createVariable(gc_lat_name, 'f4', (grid_y_center_dim, grid_x_center_dim))
            gn_lon_name = '{0}_node_lon'.format(grid_var)
            grid_node_lon = nclocal.createVariable(gn_lon_name, 'f4', (grid_y_node_dim, grid_x_node_dim))
            gn_lat_name = '{0}_node_lat'.format(grid_var)
            grid_node_lat = nclocal.createVariable(gn_lat_name, 'f4', (grid_y_node_dim, grid_x_node_dim))
            grid_var = nclocal.createVariable(grid_var, 'i2')
            grid_time = nclocal.createVariable('time', 'f8', ('time',))
            # add attributes to the variables
            grid_var.cf_role = 'grid_topology'
            grid_var.topology_dimension = 2
            grid_var.node_dimensions = '{0} {1}'.format(grid_x_node_dim, grid_y_node_dim)
            grid_var.face_dimensions = ('{x_center}: {x_node} (padding: {x_padding}) '
                                        '{y_center}: {y_node} (padding: {y_padding})').format(x_center=grid_x_center_dim,
                                                                                              x_node=grid_x_node_dim,
                                                                                              x_padding=self._face_padding[0].padding,
                                                                                              y_center=grid_y_center_dim,
                                                                                              y_node=grid_y_node_dim,
                                                                                              y_padding=self._face_padding[1].padding
                                                                                              )
            if self._edge_1_padding is not None and self._edge_2_padding is not None:
                dim_1 = self._edge_1_padding[0].dim
                dim_2 = self._edge_2_padding[0].dim
                dim_1_split = dim_1.split('_')
                dim_2_split = dim_2.split('_')
                dim_1_prefix = dim_1_split[0]
                dim_2_prefix = dim_2_split[0]
                dim_1_suffix = dim_1_split[1]
                dim_2_suffix = dim_2_split[1]
                sub_dim_1 = self._edge_1_padding[0].sub_dim
                sub_dim_2 = self._edge_2_padding[0].sub_dim
                padding_1 = self._edge_1_padding[0].padding
                padding_2 = self._edge_2_padding[0].padding
                sub_dim_suffix = sub_dim_1.split('_')[1]
                # edge 1 padding value
                edge1_str = '{dim_2}: {sub_dim_2} {dim_1}: {sub_dim_1} (padding: {padding})'
                e1_padding_str = edge1_str.format(dim_2='{0}_{1}'.format(dim_2_prefix, dim_1_suffix),
                                                  sub_dim_2='{0}_{1}'.format(dim_2_prefix, sub_dim_suffix),
                                                  dim_1=dim_1,
                                                  sub_dim_1=sub_dim_1,
                                                  padding=padding_1
                                                  )
                grid_var.edge1_padding = e1_padding_str
                # edge 2 padding value
                edge2_str = '{dim_2}: {sub_dim_2} (padding: {padding}) {dim_1}: {sub_dim_1}'
                e2_padding_str = edge2_str.format(dim_2=dim_2,
                                                  sub_dim_2=sub_dim_2,
                                                  padding=padding_2,
                                                  dim_1='{0}_{1}'.format(dim_1_prefix, dim_2_suffix),
                                                  sub_dim_1='{0}_{1}'.format(dim_1_prefix, sub_dim_suffix)
                                                  )
                grid_var.edge2_padding = e2_padding_str
            if self._face_coordinates is not None:
                grid_var.face_coordinates = ' '.join(self._face_coordinates)
            if self._node_coordinates is not None:
                grid_var.node_dimensions = ' '.join(self._node_coordinates)
            if self._edge_1_coordinates is not None:
                grid_var.edge1_coordinates = ' '.join(self._edge_1_coordinates)
            if self._edge_2_coordinates is not None:
                grid_var.edge2_coordinates = ' '.join(self._edge_2_coordinates)
            # populate variables with data
            grid_time[:] = self._grid_times[:]
            grid_center_lon[:, :] = self._centers[:, :, 0]
            grid_center_lat[:, :] = self._centers[:, :, 1] 
            grid_node_lon[:, :] = self._nodes[:, :, 0]
            grid_node_lat[:, :] = self._nodes[:, :, 1]