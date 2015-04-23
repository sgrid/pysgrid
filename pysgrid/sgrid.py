'''
Created on Apr 20, 2015

@author: ayan
'''
import abc

import netCDF4 as nc4

from .custom_exceptions import SGridNonCompliantError
from .read_netcdf import NetCDFDataset, parse_padding
from .utils import pair_arrays
from .variables import SGridVariable


class SGridND(object):
    
    __metaclass__ = abc.ABCMeta

    padding_slices = {'both': (1, -1),
                      'none': (None, None),
                      'low': (1, None),
                      'high': (None, 1)
                      }
    topology_dimensions = None
    
    def __init__(self, 
                 nodes=None,
                 centers=None,
                 edges=None,
                 node_padding=None,
                 edge1_padding=None,
                 edge2_padding=None,
                 grid_topology_vars=None,
                 grid_times=None,
                 variables=None,
                 grid_variables=None,
                 dimensions=None,
                 node_dimensions=None,
                 node_coordinates=None,
                 edge1_coordinates=None,
                 edge2_coordinates=None,
                 angles=None,
                 edge1_dimensions=None,
                 edge2_dimensions=None):
        # general attributes
        self._nodes = nodes
        self._centers = centers
        self._edges = edges
        self._node_padding = node_padding
        self._edge1_padding = edge1_padding
        self._edge2_padding = edge2_padding
        self._grid_topology_vars = grid_topology_vars
        self._grid_times = grid_times
        self._variables = variables
        self._grid_variables = grid_variables
        self._dimensions = dimensions
        self._node_dimensions = node_dimensions
        self._node_coordinates = node_coordinates
        self._edge1_coordinates = edge1_coordinates
        self._edge2_coordinates = edge2_coordinates
        self._angles = angles
        self._edge1_dimensions = edge1_dimensions
        self._edge2_dimensions = edge2_dimensions
        
    @property
    def grid_topology_vars(self):
        return self._grid_topology_vars
        
    @property
    def variables(self):
        """
        Return a list of variables
        
        """
        return self._variables
        
    @property
    def grid_variables(self):
        """
        Return a list of variables
        with a grid attribute.
        
        """
        return self._grid_variables
    
    @property
    def non_grid_variables(self):
        non_grid_variables = [variable for variable in self._variables if variable not in self._grid_variables]
        return non_grid_variables
        
    @property
    def dimensions(self):
        """
        Return list of tuples containing
        dimension name and size.
        
        """
        return self._dimensions

    @property
    def nodes(self):
        """
        return the vertices of the grid as arrays 
        of lon, lat pairs.
        
        """
        return self._nodes
    
    @property
    def node_dimensions(self):
        return self._node_dimensions
        
    @property
    def node_coordinates(self):
        return self._node_coordinates

    @property
    def edge1_coordinates(self):
        return self._edge1_coordinates

    @property
    def edge1_padding(self):
        return self._edge1_padding

    @property
    def edge1_dimension(self):
        return self._edge1_dimensions
        
    @property
    def edge2_coordinates(self):
        return self._edge2_coordinates
        
    @property
    def edge2_padding(self):
        return self._edge2_padding
        
    @property
    def edge2_dimensions(self):
        return self._edge2_dimensions
        
    @property
    def angles(self):
        return self._angles
        
    @property
    def centers(self):
        """
        return the coordinates of the grid centers
        as arrays of lon, lat pairs.
        
        """
        return self._centers

    @property
    def grid_times(self):
        return self._grid_times
    
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
        
    @abc.abstractmethod
    def _define_face_padding_summary(self):
        return
    
    @abc.abstractmethod
    def save_as_netcdf(self):
        return


class SGrid2D(SGridND):
    
    topology_dimension = 2
    
    def __init__(self,
                 faces=None,
                 face_padding=None,
                 face_coordinates=None,
                 face_dimensions=None,
                 vertical_padding=None,
                 vertical_dimensions=None,
                 *args,
                 **kwargs):
        self._faces = faces
        self._face_padding = face_padding
        self._face_coordinates = face_coordinates
        self._face_dimensions = face_dimensions
        self._vertical_padding = vertical_padding
        self._vertical_dimensions = vertical_dimensions
        super(SGrid2D, self).__init__(*args, **kwargs)
        
    @classmethod
    def sgrid_from_dataset(cls, nc_dataset, topology_variable=None):
        sa = SGridAttributes(nc_dataset, 2, topology_variable)
        dimensions = sa.get_dimensions()
        node_dimensions, node_coordinates = sa.get_node_coordinates()
        grid_topology_vars = sa.get_topology_vars()
        edge1_dimensions, edge1_padding = sa.get_attr_dimension('edge1_dimensions')
        edge2_dimensions, edge2_padding = sa.get_attr_dimension('edge2_dimensions')
        edge1_coordinates = sa.get_attr_coordinates('edge1_coordinates')
        edge2_coordinates = sa.get_attr_coordinates('edge2_coordinates')
        angles = sa.get_angles()
        grid_times = sa.get_time()
        vertical_dimensions, vertical_padding = sa.get_attr_dimension('vertical_dimensions')
        centers = sa.get_cell_center_lat_lon()
        face_dimensions, face_padding = sa.get_attr_dimension('face_dimensions')
        face_coordinates = sa.get_attr_coordinates('face_coordinates')
        nodes = sa.get_cell_node_lat_lon()
        sgrid = cls(angles=angles,
                    centers=centers,
                    dimensions=dimensions,
                    edge1_coordinates=edge1_coordinates,
                    edge1_dimensions=edge1_dimensions,
                    edge1_padding=edge1_padding,
                    edge2_coordinates=edge2_coordinates,
                    edge2_dimensions=edge2_dimensions,
                    edge2_padding=edge2_padding,
                    edges=None,
                    face_coordinates=face_coordinates,
                    face_dimensions=face_dimensions,
                    face_padding=face_padding,
                    faces=None,
                    grid_times=grid_times,
                    grid_topology_vars=grid_topology_vars,
                    grid_variables=None,
                    node_coordinates=node_coordinates,
                    node_dimensions=node_dimensions,
                    node_padding=None,
                    nodes=nodes,
                    variables=None,
                    vertical_dimensions=vertical_dimensions,
                    vertical_padding=vertical_padding
                    )
        sa.get_variable_attributes(sgrid)
        return sgrid
    
    @classmethod
    def sgrid_from_file(cls, nc_file_path, topology_variable=None):
        with nc4.Dataset(nc_file_path) as nc_dataset:
            sgrid = cls.sgrid_from_dataset(nc_dataset, topology_variable)
        return sgrid
        
    @property
    def face_padding(self):
        return self._face_padding

    @property
    def face_coordinates(self):
        return self._face_coordinates
    
    @property
    def face_dimensions(self):
        return self._face_dimensions

    @property
    def vertical_padding(self):
        return self._vertical_padding

    @property
    def vertical_dimensions(self):
        return self._vertical_dimensions

    def _define_face_padding_summary(self):
        all_padding = []
        if self._face_padding is not None:
            all_padding += self._face_padding
        if self._vertical_padding is not None:
            all_padding += self._vertical_padding
        if self._edge1_padding is not None:
            all_padding += self._edge1_padding
        if self._edge2_padding is not None:
            all_padding += self._edge2_padding
        padding_summary = []
        for padding_datum in all_padding:
            dim = padding_datum.dim
            sub_dim = padding_datum.sub_dim
            padding_val = padding_datum.padding
            pad_short = (dim, sub_dim, padding_val)
            padding_summary.append(pad_short)
        return padding_summary
        
    def save_as_netcdf(self, filepath):
        with nc4.Dataset(filepath, 'w') as nclocal:
            grid_var = self._grid_topology_vars
            # create dimensions
            for grid_dim in self._dimensions:
                dim_name, dim_size = grid_dim
                nclocal.createDimension(dim_name, dim_size)
            # create variables
            center_lon, center_lat = self._face_coordinates
            center_lon_obj = getattr(self, center_lon)
            center_lat_obj = getattr(self, center_lat)
            node_lon, node_lat = self._node_coordinates
            node_lon_obj = getattr(self, node_lon)
            node_lat_obj = getattr(self, node_lat)
            grid_center_lon = nclocal.createVariable(center_lon, 
                                                     center_lon_obj.dtype, 
                                                     center_lon_obj.dimensions
                                                     )
            grid_center_lat = nclocal.createVariable(center_lat, 
                                                     center_lat_obj.dtype, 
                                                     center_lat_obj.dimensions
                                                     )
            grid_node_lon = nclocal.createVariable(node_lon, 
                                                   node_lon_obj.dtype, 
                                                   node_lon_obj.dimensions
                                                   )
            grid_node_lat = nclocal.createVariable(node_lat, 
                                                   node_lat_obj.dtype, 
                                                   node_lat_obj.dimensions
                                                   )
            grid_var_obj = getattr(self, grid_var)
            grid_vars = nclocal.createVariable(grid_var, grid_var_obj.dtype)
            # not the most robust here... time and angle are hard-coded
            # need to address this
            time_obj = getattr(self, 'time')
            grid_time = nclocal.createVariable('time', 
                                               time_obj.dtype, 
                                               time_obj.dimensions
                                               )
            
            if hasattr(self, 'angle'):
                angle_obj = getattr(self, 'angle', None)
                grid_angle = nclocal.createVariable('angle', 
                                                    angle_obj.dtype, 
                                                    angle_obj.dimensions
                                                    )
                if self._angles is not None:
                    grid_angle[:] = self._angles[:]
            # save the grid variables with attributes
            for dataset_variable in self._variables:
                dataset_var_obj = getattr(self, dataset_variable)
                if dataset_var_obj.grid is not None:
                    dataset_grid_var = nclocal.createVariable(dataset_variable, 
                                                              dataset_var_obj.dtype, 
                                                              dataset_var_obj.dimensions
                                                              )
                    dataset_grid_var.grid = grid_var
            # add attributes to the variables
            grid_vars.cf_role = 'grid_topology'
            grid_vars.topology_dimension = self.topology_dimension
            grid_vars.node_dimensions = self._node_dimensions
            if self._face_dimensions is not None:
                grid_vars.face_dimensions = self._face_dimensions
            if self._edge1_dimensions is not None:
                grid_vars.edge1_dimensions = self._edge1_dimensions
            if self._edge2_dimensions is not None:
                grid_vars.edge2_dimensions = self._edge2_dimensions
            if self._vertical_dimensions is not None:
                grid_vars.vertical_dimensions = self._vertical_dimensions
            if self._face_coordinates is not None:
                grid_vars.face_coordinates = ' '.join(self._face_coordinates)
            if self._node_coordinates is not None:
                grid_vars.node_coordinates = ' '.join(self._node_coordinates)
            if self._edge1_coordinates is not None:
                grid_vars.edge1_coordinates = ' '.join(self._edge1_coordinates)
            if self._edge2_coordinates is not None:
                grid_vars.edge2_coordinates = ' '.join(self._edge2_coordinates)
            # populate variables with data
            grid_time[:] = self._grid_times[:]
            grid_center_lon[:, :] = self._centers[:, :, 0]
            grid_center_lat[:, :] = self._centers[:, :, 1] 
            grid_node_lon[:, :] = self._nodes[:, :, 0]
            grid_node_lat[:, :] = self._nodes[:, :, 1]
        
        
class SGrid3D(SGridND):
    
    topology_dimension = 3
    
    def __init__(self,
                 volume_padding=None,
                 volume_coordinates=None,
                 edge3_padding=None,
                 edge3_coordinates=None,
                 face1_padding=None,
                 face1_coordinates=None,
                 face2_padding=None,
                 face2_coordinates=None,
                 face3_padding=None,
                 face3_coordinates=None,
                 volume_dimensions=None,
                 edge1_dimensions=None,
                 edge2_dimensions=None,
                 edge3_dimensions=None,
                 face1_dimensions=None,
                 face2_dimensions=None,
                 face3_dimensions=None,
                 *args,
                 **kwargs):
        # attributes specific to topology_dimension 3
        self._volume_padding = volume_padding
        self._volume_coordinates = volume_coordinates
        self._edge3_padding = edge3_padding
        self._edge3_coordinates = edge3_coordinates
        self._face1_padding = face1_padding
        self._face1_coordinates = face1_coordinates
        self._face2_padding = face2_padding
        self._face2_coordinates = face2_coordinates
        self._face3_padding = face3_padding
        self._face3_coordinates = face3_coordinates
        self._volume_dimensions = volume_dimensions
        self._edge1_dimensions = edge1_dimensions
        self._edge2_dimensions = edge2_dimensions
        self._edge3_dimensions = edge3_dimensions
        self._face1_dimensions = face1_dimensions
        self._face2_dimensions = face2_dimensions
        self._face3_dimensions = face3_dimensions
        super(SGrid3D, self).__init__(*args, **kwargs)
        
    @classmethod
    def sgrid_from_dataset(cls, nc_dataset, topology_variable=None):
        sa = SGridAttributes(nc_dataset, 3, topology_variable)
        dimensions = sa.get_dimensions()
        node_dimensions, node_coordinates = sa.get_node_coordinates()
        grid_topology_vars = sa.get_topology_vars()
        edge1_dimensions, edge1_padding = sa.get_attr_dimension('edge1_dimensions')
        edge2_dimensions, edge2_padding = sa.get_attr_dimension('edge2_dimensions')
        edge1_coordinates = sa.get_attr_coordinates('edge1_coordinates')
        edge2_coordinates = sa.get_attr_coordinates('edge2_coordinates')
        grid_times = sa.get_time()
        edge3_dimensions, edge3_padding = sa.get_attr_dimension('edge3_dimensions')
        edge3_coordinates = sa.get_attr_coordinates('edge3_coordinates')
        face1_dimensions, face1_padding = sa.get_attr_dimension('face1_dimensions')
        face1_coordinates = sa.get_attr_coordinates('face1_coordinates')
        face2_dimensions, face2_padding = sa.get_attr_dimension('face2_dimensions')
        face2_coordinates = sa.get_attr_coordinates('face2_coordinates')
        face3_dimensions, face3_padding = sa.get_attr_dimension('face3_dimensions')
        face3_coordinates = sa.get_attr_coordinates('face3_coordinates')
        volume_dimensions, volume_padding = sa.get_attr_dimension('volume_dimensions')
        volume_coordinates = sa.get_attr_coordinates('volume_coordinates')
        centers = sa.get_cell_center_lat_lon_3d()
        nodes = sa.get_cell_node_lat_lon_3d()
        sgrid = cls(angles=None,
                    centers=centers,
                    dimensions=dimensions,
                    edge1_coordinates=edge1_coordinates,
                    edge1_dimensions=edge1_dimensions,
                    edge1_padding=edge1_padding,
                    edge2_coordinates=edge2_coordinates,
                    edge2_dimensions=edge2_dimensions,
                    edge2_padding=edge2_padding,
                    edge3_coordinates=edge3_coordinates,
                    edge3_dimensions=edge3_dimensions,
                    edge3_padding=edge3_padding,
                    edges=None,
                    face1_coordinates=face1_coordinates,
                    face1_dimensions=face1_dimensions,
                    face1_padding=face1_padding,
                    face2_coordinates=face2_coordinates,
                    face2_dimensions=face2_dimensions,
                    face2_padding=face2_padding,
                    face3_coordinates=face3_coordinates,
                    face3_dimensions=face3_dimensions,
                    face3_padding=face3_padding,
                    grid_times=grid_times,
                    grid_topology_vars=grid_topology_vars,
                    grid_variables=None,
                    node_coordinates=node_coordinates,
                    node_dimensions=node_dimensions,
                    node_padding=None,
                    nodes=nodes,
                    variables=None,
                    volume_coordinates=volume_coordinates,
                    volume_dimensions=volume_dimensions,
                    volume_padding=volume_padding
                    )
        sa.get_variable_attributes(sgrid)
        return sgrid
    
    @classmethod
    def sgrid_from_file(cls, nc_file_path, topology_variable=None):
        with nc4.Dataset(nc_file_path) as nc_dataset:
            sgrid = cls.sgrid_from_dataset(nc_dataset, topology_variable)
        return sgrid
        
    @property
    def volume_padding(self):
        return self._volume_padding

    @property
    def volume_dimensions(self):
        return self._volume_dimensions

    @property
    def volume_coordinates(self):
        return self._volume_coordinates

    @property
    def face1_padding(self):
        return self._face1_padding

    @property
    def face1_coordinates(self):
        return self._face1_coordinates

    @property
    def face1_dimensions(self):
        return self._face1_dimensions
        
    @property
    def face2_padding(self):
        return self._face2_padding
        
    @property
    def face2_coordinates(self):
        return self._face2_coordinates

    @property    
    def face2_dimensions(self):
        return self._face2_dimensions

    @property
    def face3_padding(self):
        return self._face3_padding

    @property
    def face3_coordinates(self):
        return self._face3_coordinates

    @property
    def face3_dimensions(self):
        return self._face3_dimensions

    @property
    def edge3_padding(self):
        return self._edge3_padding

    @property
    def edge3_coordinates(self):
        return self._edge3_coordinates

    @property
    def edge3_dimensions(self):
        return self._edge3_dimensions

    def _define_face_padding_summary(self):
        all_padding = []
        if self._volume_padding is not None:
            all_padding += self._volume_padding
        if self._face1_padding is not None:
            all_padding += self._face1_padding
        if self._face2_padding is not None:
            all_padding += self._face2_padding
        if self._face3_padding is not None:
            all_padding += self._face3_padding
        if self._edge1_padding is not None:
            all_padding += self._edge1_padding
        if self._edge2_padding is not None:
            all_padding += self._edge2_padding
        if self._edge3_padding is not None:
            all_padding += self._edge3_padding
        padding_summary = []
        for padding_datum in all_padding:
            dim = padding_datum.dim
            sub_dim = padding_datum.sub_dim
            padding_val = padding_datum.padding
            pad_short = (dim, sub_dim, padding_val)
            padding_summary.append(pad_short)
        return padding_summary
        
    def save_as_netcdf(self, filepath):
        with nc4.Dataset(filepath, 'w') as nclocal:
            grid_var = self._grid_topology_vars
            # create dimensions
            for grid_dim in self._dimensions:
                dim_name, dim_size = grid_dim
                nclocal.createDimension(dim_name, dim_size)
            # create variables
            center_lon, center_lat = self._face_coordinates
            center_lon_obj = getattr(self, center_lon)
            center_lat_obj = getattr(self, center_lat)
            node_lon, node_lat = self._node_coordinates
            node_lon_obj = getattr(self, node_lon)
            node_lat_obj = getattr(self, node_lat)
            grid_center_lon = nclocal.createVariable(center_lon, 
                                                     center_lon_obj.dtype, 
                                                     center_lon_obj.dimensions
                                                     )
            grid_center_lat = nclocal.createVariable(center_lat, 
                                                     center_lat_obj.dtype, 
                                                     center_lat_obj.dimensions
                                                     )
            grid_node_lon = nclocal.createVariable(node_lon, 
                                                   node_lon_obj.dtype, 
                                                   node_lon_obj.dimensions
                                                   )
            grid_node_lat = nclocal.createVariable(node_lat, 
                                                   node_lat_obj.dtype, 
                                                   node_lat_obj.dimensions
                                                   )
            grid_var_obj = getattr(self, grid_var)
            grid_vars = nclocal.createVariable(grid_var, grid_var_obj.dtype)
            # not the most robust here... time and angle are hard-coded
            # need to address this
            time_obj = getattr(self, 'time')
            grid_time = nclocal.createVariable('time', 
                                               time_obj.dtype, 
                                               time_obj.dimensions
                                               )
            
            if hasattr(self, 'angle'):
                angle_obj = getattr(self, 'angle', None)
                grid_angle = nclocal.createVariable('angle', 
                                                    angle_obj.dtype, 
                                                    angle_obj.dimensions
                                                    )
                if self._angles is not None:
                    grid_angle[:] = self._angles[:]
            # save the grid variables with attributes
            for dataset_variable in self._variables:
                dataset_var_obj = getattr(self, dataset_variable)
                if dataset_var_obj.grid is not None:
                    dataset_grid_var = nclocal.createVariable(dataset_variable, 
                                                              dataset_var_obj.dtype, 
                                                              dataset_var_obj.dimensions
                                                              )
                    dataset_grid_var.grid = grid_var
            # add attributes to the variables
            grid_vars.cf_role = 'grid_topology'
            grid_vars.topology_dimension = self.topology_dimension
            grid_vars.node_dimensions = self._node_dimensions
            grid_vars.volume_dimensions = self._volume_dimensions
            if self._volume_coordinates is not None:
                grid_vars.volume_coordinates = ' '.join(self._volume_coordinates)
            if self._node_coordinates is not None:
                grid_vars.node_coordinates = ' '.join(self._node_coordinates)
            if self._face_dimensions is not None:
                grid_vars.face_dimensions = self._face_dimensions
            if self._face1_coordinates is not None:
                grid_vars.face1_coordinates = ' '.join(self._face1_coordinates)
            if self._face2_coordinates is not None:
                grid_vars.face2_coordinates = self._face2_coordinates
            if self._face3_coordinates is not None:
                grid_vars.face3_coordinates = self._face3_coordinates
            if self._edge1_dimensions is not None:
                grid_vars.edge1_dimensions = self._edge1_dimensions
            if self._edge2_dimensions is not None:
                grid_vars.edge2_dimensions = self._edge2_dimensions
            if self._edge3_dimensions is not None:
                grid_vars.edge3_dimensions = self._edge3_dimensions
            if self._edge1_coordinates is not None:
                grid_vars.edge1_coordinates = ' '.join(self._edge1_coordinates)
            if self._edge2_coordinates is not None:
                grid_vars.edge2_coordinates = ' '.join(self._edge2_coordinates)
            if self._edge3_coordinates is not None:
                grid_vars.edge3_coordinates = ' '.join(self._edge3_coordinates)
            # populate variables with data
            grid_time[:] = self._grid_times[:]
            grid_center_lon[:] = self._centers[..., 0]
            grid_center_lat[:] = self._centers[..., 1] 
            grid_node_lon[:] = self._nodes[..., 0]
            grid_node_lat[:] = self._nodes[..., 1]
    

class SGridAttributes(object):
    """
    Class containg methods to help with getting the
    attributes for either a 2D or 3D SGrid.
    
    """
    def __init__(self, nc_dataset, topology_dim, topology_variable=None):
        self.nc_dataset = nc_dataset
        self.ncd = NetCDFDataset(self.nc_dataset)
        self.topology_dim = topology_dim
        if topology_variable is None:
            # the netCDF variable with a cf_role of 'grid_topology'
            self.topology_variable = self.ncd.find_grid_topology_var()
        else:
            self.topology_variable = topology_variable
        self.topology_var = self.nc_dataset.variables[self.topology_variable]
    
    def get_dimensions(self):
        ds_dims = self.nc_dataset.dimensions
        grid_dims = [(ds_dim, len(ds_dims[ds_dim])) for ds_dim in ds_dims]
        return grid_dims
        
    def get_topology_vars(self):
        grid_topology_vars = self.ncd.find_grid_topology_var()
        return grid_topology_vars
    
    def get_attr_dimension(self, attr_name):
        try:
            attr_dim = getattr(self.topology_var, attr_name)
        except AttributeError:
            attr_dim = None
            attr_padding = None
        else:
            attr_dim_padding = parse_padding(attr_dim, self.topology_variable)
            attr_padding = attr_dim_padding
        return attr_dim, attr_padding
    
    def get_attr_coordinates(self, attr_name):
        try:
            attr_coordinates_raw = getattr(self.topology_var, attr_name)
        except AttributeError:
            location_name = attr_name.split('_')[0]
            attr_coordinates = self.ncd.find_coordinates_by_location(location_name, self.topology_dim)
        else:
            attr_coordinates_val = attr_coordinates_raw.split(' ')
            attr_coordinates = tuple(attr_coordinates_val)
        return attr_coordinates

    def get_node_coordinates(self):
        node_dims = self.topology_var.node_dimensions
        node_dimensions = node_dims
        try:
            node_coordinates = self.topology_var.node_coordinates
        except AttributeError:
            grid_cell_node_vars = self.ncd.find_grid_cell_node_vars()
            node_coordinates = grid_cell_node_vars
        else:
            node_coordinate_val = node_coordinates.split(' ')
            node_coordinates = tuple(node_coordinate_val)
        return node_dimensions, node_coordinates
            
    def get_variable_attributes(self, sgrid):
        dataset_variables = []
        grid_variables = []
        nc_variables = self.nc_dataset.variables
        for nc_variable in nc_variables:
            nc_var = nc_variables[nc_variable]
            sgrid_var = SGridVariable.create_variable(nc_var, sgrid)
            sgrid.add_property(sgrid_var.variable, sgrid_var)
            dataset_variables.append(nc_var.name)
            if hasattr(nc_var, 'grid'):
                grid_variables.append(nc_var.name)
        sgrid._variables = dataset_variables
        sgrid._grid_variables = grid_variables
        
    def get_angles(self):
        try:
            # remove hard coding of variable name moving forward
            grid_angles = self.nc_dataset.variables['angle'][:]
            angles = grid_angles
        except KeyError:
            angles = None
        return angles
        
    def get_time(self):
        try:
            # hard coding the time variable is not the best way to go...
            # change this in the future
            grid_time = self.nc_dataset.variables['time'][:]
        except KeyError:
            grid_time = self.nc_dataset.variables['Times'][:]
        return grid_time
        
    def get_cell_center_lat_lon(self):
        grid_cell_center_lon_var, grid_cell_center_lat_var = self.get_attr_coordinates('face_coordinates')
        grid_cell_center_lat = self.nc_dataset.variables[grid_cell_center_lat_var][:]
        grid_cell_center_lon = self.nc_dataset.variables[grid_cell_center_lon_var][:]
        return pair_arrays(grid_cell_center_lon, grid_cell_center_lat)
        
    def get_cell_node_lat_lon(self):
        grid_cell_nodes_lon_var, grid_cell_nodes_lat_var = self.get_node_coordinates()[1]
        grid_cell_nodes_lat = self.nc_dataset.variables[grid_cell_nodes_lat_var][:]
        grid_cell_nodes_lon = self.nc_dataset.variables[grid_cell_nodes_lon_var][:]
        return pair_arrays(grid_cell_nodes_lon, grid_cell_nodes_lat)
        
    def get_cell_center_lat_lon_3d(self):
        volume_coordinates = self.get_attr_coordinates('volume_coordinates')
        grid_cell_center_lon_var = volume_coordinates[0]
        grid_cell_center_lat_var = volume_coordinates[1]
        grid_cell_center_lon = self.nc_dataset.variables[grid_cell_center_lon_var][:]
        grid_cell_center_lat = self.nc_dataset.variables[grid_cell_center_lat_var][:]
        return pair_arrays(grid_cell_center_lon, grid_cell_center_lat)
        
    def get_cell_node_lat_lon_3d(self):
        pass
        
        
def _load_grid_from_nc_dataset(nc_dataset,
                               topology_dim,
                               grid_topology_var=None
                               ):
    """
    Create an SGridND object from an SGRID
    compliant netCDF4.Dataset object. An
    exception is raised if the dataset is
    non-compliant. This function will introspect
    the datast to determine whether a 2D or 3D
    SGRID object is returned.
    
    :param nc_dataset: a netCDF resource read into a netCDF4.Dataset object
    :type nc_dataset: netCDF4.Dataset
    :param grid_topology_var: the name of the grid topology variable; defaults to None
    :type grid_topology_var: str
    :return: an SGrid object
    :rtype: sgrid.SGrid2D or sgrid.SGrid3D
    
    """
    if topology_dim == 2:
        grid = SGrid2D.sgrid_from_dataset(nc_dataset, grid_topology_var)
    elif topology_dim == 3:
        grid = SGrid3D.sgrid_from_dataset(nc_dataset, grid_topology_var)
    else:
        raise ValueError('Only topology dimensions of 2 or 3 are supported')
    return grid
    
    
def _return_grid_topology_dim(nc_dataset, grid_topology_var=None):
    """
    Given a netCDF dataset, determine the topology
    dimension.
    
    :param nc_dataset: a netCDF dataset
    :type nc_dataset: netCDF4.Dataset
    :param str grid_topology_vars: the name of the grid topology variable; defaults to None
    :return: topology dimension
    :rtype: int
    
    """
    ncd = NetCDFDataset(nc_dataset)
    if ncd.sgrid_compliant_file():
        if grid_topology_var is not None:
            topology_var = grid_topology_var
        else:
            topology_var = ncd.find_grid_topology_var()
        nc_grid_topology_var = nc_dataset.variables[topology_var]
        topology_dim = nc_grid_topology_var.topology_dimension
        if topology_dim == 2 or topology_dim == 3:
            return topology_dim, grid_topology_var
        else:
            raise ValueError('Only topology dimensions of 2 or 3 are supported')
    else:
        raise SGridNonCompliantError(nc_dataset)
    
    
def from_nc_file(nc_url, grid_topology_var=None):
    """
    Get a SGrid object from a file. There is no need
    to know the topology dimensions a priori.
    
    :param str nc_url: URL or filepath to the netCDF file
    :param str grid_topology_vars: the name of the grid topology variable; defaults to None
    :return: SGrid object
    :rtype: sgrid.SGrid2D or sgrid.SGrid3D
    
    """
    with nc4.Dataset(nc_url, 'r') as nc_dataset:
        topology_dim, introspected_grid_topology_var = _return_grid_topology_dim(nc_dataset, grid_topology_var)
        if grid_topology_var is not None:
            topology_var = grid_topology_var
        else:
            topology_var = introspected_grid_topology_var
        grid = _load_grid_from_nc_dataset(nc_dataset, 
                                          topology_dim, 
                                          topology_var
                                          )
    return grid


def from_nc_dataset(nc_dataset, grid_topology_var=None):
    """
    Get a SGrid object from a netCDF4.Dataset. There is no need
    to know the topology dimensions a priori.
    
    :param netCDF4.Dataset nc_dataset: a netCDF4 Dataset
    :param str grid_topology_vars: the name of the grid topology variable; defaults to None
    :return: SGrid object
    :rtype: sgrid.SGrid2D or sgrid.SGrid3D
    
    """
    topology_dim, introspected_grid_topology_var = _return_grid_topology_dim(nc_dataset, grid_topology_var)
    if grid_topology_var is not None:
        topology_var = grid_topology_var
    else:
        topology_var = introspected_grid_topology_var
    grid = _load_grid_from_nc_dataset(nc_dataset, 
                                      topology_dim, 
                                      topology_var
                                      )
    return grid