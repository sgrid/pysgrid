'''
Created on Apr 20, 2015

@author: ayan
'''
import abc

import netCDF4 as nc4

from .custom_exceptions import SGridNonCompliantError
from .read_netcdf import NetCDFDataset, parse_padding
from .utils import calculate_angle_from_true_east, pair_arrays
from .variables import SGridVariable


class SGridND(object):
    
    __metaclass__ = abc.ABCMeta

    padding_slices = {'both': (1, -1),
                      'none': (None, None),
                      'low': (1, None),
                      'high': (None, 1)
                      }
    topology_dimension = None
    
    def __init__(self, 
                 nodes=None,
                 centers=None,
                 edges=None,
                 node_padding=None,
                 edge1_padding=None,
                 edge2_padding=None,
                 grid_topology_var=None,
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
        self.nodes = nodes
        self.centers = centers
        self.edges = edges
        self.node_padding = node_padding
        self.edge1_padding = edge1_padding
        self.edge2_padding = edge2_padding
        self.grid_topology_var = grid_topology_var
        self.variables = variables
        self.grid_variables = grid_variables
        self.dimensions = dimensions
        self.node_dimensions = node_dimensions
        self.node_coordinates = node_coordinates
        self.edge1_coordinates = edge1_coordinates
        self.edge2_coordinates = edge2_coordinates
        self.angles = angles
        self.edge1_dimensions = edge1_dimensions
        self.edge2_dimensions = edge2_dimensions
        
    @classmethod
    def from_ncfile(cls, nc_file_path, topology_variable=None):
        with nc4.Dataset(nc_file_path) as nc_dataset:
            sgrid = cls.sgrid_from_dataset(nc_dataset, topology_variable)
        return sgrid
    
    @property
    def non_grid_variables(self):
        non_grid_variables = [variable for variable in self.variables if variable not in self.grid_variables]
        return non_grid_variables
    
    def _save_common_components(self, nc_file):
        grid_var = self.grid_topology_var
        # create dimensions
        for grid_dim in self.dimensions:
            dim_name, dim_size = grid_dim
            nc_file.createDimension(dim_name, dim_size)
        # create variables
        center_lon, center_lat = self.face_coordinates
        center_lon_obj = getattr(self, center_lon)
        center_lat_obj = getattr(self, center_lat)
        grid_center_lon = nc_file.createVariable(center_lon_obj.variable,
                                                 center_lon_obj.dtype,
                                                 center_lon_obj.dimensions
                                                 )
        grid_center_lat = nc_file.createVariable(center_lat_obj.variable,
                                                 center_lat_obj.dtype,
                                                 center_lat_obj.dimensions
                                                 )
        grid_center_lon[:] = self.centers[..., 0]
        grid_center_lat[:] = self.centers[..., 1]
        try:
            node_lon, node_lat = self.node_coordinates
        except TypeError:
            pass
        else:
            node_lon_obj = getattr(self, node_lon)
            grid_node_lon = nc_file.createVariable(node_lon_obj.variable,
                                                   node_lon_obj.dtype,
                                                   node_lon_obj.dimensions
                                                   )
            node_lat_obj = getattr(self, node_lat)
            grid_node_lat = nc_file.createVariable(node_lat_obj.variable,
                                                   node_lat_obj.dtype,
                                                   node_lat_obj.dimensions
                                                   )
            grid_node_lon[:] = self.nodes[..., 0]
            grid_node_lat[:] = self.nodes[..., 1]
        grid_var_obj = getattr(self, grid_var)
        grid_vars = nc_file.createVariable(grid_var_obj.variable, grid_var_obj.dtype)
        grid_vars.cf_role = 'grid_topology'
        grid_vars.topology_dimension = self.topology_dimension
        grid_vars.node_dimensions = self.node_dimensions
        if self.edge1_dimensions is not None:
            grid_vars.edge1_dimensions = self.edge1_dimensions
        if self.edge2_dimensions is not None:
            grid_vars.edge2_dimensions = self.edge2_dimensions
        if self.node_coordinates is not None:
            grid_vars.node_coordinates = ' '.join(self.node_coordinates)
        if self.edge1_coordinates is not None:
            grid_vars.edge1_coordinates = ' '.join(self.edge1_coordinates)
        if self.edge2_coordinates is not None:
            grid_vars.edge2_coordinates = ' '.join(self.edge2_coordinates)
        if hasattr(self, 'angle'):
            angle_obj = getattr(self, 'angle', None)
            grid_angle = nc_file.createVariable(angle_obj.variable,
                                                angle_obj.dtype,
                                                angle_obj.dimensions
                                                )
            if self.angles is not None:
                grid_angle[:] = self.angles[:]
        for dataset_variable in self.grid_variables:
            dataset_var_obj = getattr(self, dataset_variable)
            axes = []
            dataset_grid_var = nc_file.createVariable(dataset_var_obj.variable,
                                                      dataset_var_obj.dtype,
                                                      dataset_var_obj.dimensions
                                                      )
            if dataset_var_obj.grid is not None:
                dataset_grid_var.grid = grid_var
            if dataset_var_obj.standard_name is not None:
                dataset_grid_var.standard_name = dataset_var_obj.standard_name
            if dataset_var_obj.x_axis is not None:
                x_axis = 'X: {0}'.format(dataset_var_obj.x_axis)
                axes.append(x_axis)
            if dataset_var_obj.y_axis is not None:
                y_axis = 'Y: {0}'.format(dataset_var_obj.y_axis)
                axes.append(y_axis)
            if dataset_var_obj.z_axis is not None:
                z_axis = 'Z: {0}'.format(dataset_var_obj.z_axis)
                axes.append(z_axis)
            if len(axes) > 0:
                dataset_grid_var.axes = ' '.join(axes)
        for dataset_variable in self.non_grid_variables:
            dataset_var_obj = getattr(self, dataset_variable)
            try:
                nc_file.createVariable(dataset_var_obj.variable,
                                       dataset_var_obj.dtype,
                                       dataset_var_obj.dimensions
                                       )
            except RuntimeError:
                # lat/lon and grid variables will already exist
                continue
        return grid_vars
    
    @abc.abstractmethod
    def from_nc_dataset(self):
        return
    
    @abc.abstractmethod
    def get_all_face_padding(self):
        return
    
    @abc.abstractmethod
    def get_all_edge_padding(self):
        return
        
    @abc.abstractmethod
    def all_padding(self):
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
        self.faces = faces
        self.face_padding = face_padding
        self.face_coordinates = face_coordinates
        self.face_dimensions = face_dimensions
        self.vertical_padding = vertical_padding
        self.vertical_dimensions = vertical_dimensions
        super(SGrid2D, self).__init__(*args, **kwargs)
        
    @classmethod
    def from_nc_dataset(cls, nc_dataset, topology_variable=None):
        sa = SGridAttributes(nc_dataset, 2, topology_variable)
        dimensions = sa.get_dimensions()
        node_dimensions, node_coordinates = sa.get_node_coordinates()
        grid_topology_var = sa.get_topology_var()
        edge1_dimensions, edge1_padding = sa.get_attr_dimension('edge1_dimensions')
        edge2_dimensions, edge2_padding = sa.get_attr_dimension('edge2_dimensions')
        edge1_coordinates = sa.get_attr_coordinates('edge1_coordinates')
        edge2_coordinates = sa.get_attr_coordinates('edge2_coordinates')
        angles = sa.get_angles()
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
                    grid_topology_var=grid_topology_var,
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
    
    def get_all_face_padding(self):
        if self.face_padding is not None:
            all_face_padding = self.face_padding
        else:
            all_face_padding = []
        return all_face_padding
    
    def get_all_edge_padding(self):
        all_edge_padding = []
        if self.edge1_padding is not None:
            all_edge_padding += self.edge1_padding
        if self.edge2_padding is not None:
            all_edge_padding += self.edge2_padding
        return all_edge_padding
                    
    def all_padding(self):
        all_padding = self.get_all_face_padding() + self.get_all_edge_padding()
        if self.vertical_padding is not None:
            all_padding += self.vertical_padding
        return all_padding
        
    def save_as_netcdf(self, filepath):
        with nc4.Dataset(filepath, 'w') as nclocal:
            grid_vars = self._save_common_components(nclocal)
            # add attributes to the grid_topology variable
            grid_vars.face_dimensions = self.face_dimensions
            if self.vertical_dimensions is not None:
                grid_vars.vertical_dimensions = self.vertical_dimensions
            if self.face_coordinates is not None:
                grid_vars.face_coordinates = ' '.join(self.face_coordinates)

              
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
        self.volume_padding = volume_padding
        self.volume_coordinates = volume_coordinates
        self.edge3_padding = edge3_padding
        self.edge3_coordinates = edge3_coordinates
        self.face1_padding = face1_padding
        self.face1_coordinates = face1_coordinates
        self.face2_padding = face2_padding
        self.face2_coordinates = face2_coordinates
        self.face3_padding = face3_padding
        self.face3_coordinates = face3_coordinates
        self.volume_dimensions = volume_dimensions
        self.edge1_dimensions = edge1_dimensions
        self.edge2_dimensions = edge2_dimensions
        self.edge3_dimensions = edge3_dimensions
        self.face1_dimensions = face1_dimensions
        self.face2_dimensions = face2_dimensions
        self.face3_dimensions = face3_dimensions
        super(SGrid3D, self).__init__(*args, **kwargs)
        
    @classmethod
    def from_nc_dataset(cls, nc_dataset, topology_variable=None):
        sa = SGridAttributes(nc_dataset, 3, topology_variable)
        dimensions = sa.get_dimensions()
        node_dimensions, node_coordinates = sa.get_node_coordinates()
        grid_topology_var = sa.get_topology_var()
        edge1_dimensions, edge1_padding = sa.get_attr_dimension('edge1_dimensions')
        edge2_dimensions, edge2_padding = sa.get_attr_dimension('edge2_dimensions')
        edge1_coordinates = sa.get_attr_coordinates('edge1_coordinates')
        edge2_coordinates = sa.get_attr_coordinates('edge2_coordinates')
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
                    grid_topology_var=grid_topology_var,
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
    
    def get_all_face_padding(self):
        all_face_padding = []
        if self.face1_padding is not None:
            all_face_padding += self.face1_padding
        if self.face2_padding is not None:
            all_face_padding += self.face2_padding
        if self.face3_padding is not None:
            all_face_padding += self.face3_padding
        return all_face_padding
    
    def get_all_edge_padding(self):
        all_edge_padding = []
        if self.edge1_padding is not None:
            all_edge_padding += self.edge1_padding
        if self.edge2_padding is not None:
            all_edge_padding += self.edge2_padding
        if self.edge3_padding is not None:
            all_edge_padding += self.edge3_padding
        return all_edge_padding

    def all_padding(self):
        all_padding = self.volume_padding + self.get_all_face_padding() + self.get_all_edge_padding()
        return all_padding
    
    def save_as_netcdf(self, filepath):
        with nc4.Dataset(filepath, 'w') as nclocal:
            grid_vars = self._save_common_components(nclocal)
            # add attributes to the variables
            grid_vars.volume_dimensions = self.volume_dimensions
            if self.volume_coordinates is not None:
                grid_vars.volume_coordinates = ' '.join(self.volume_coordinates)
            if self.face3_coordinates is not None:
                grid_vars.face3_coordinates = self.face3_coordinates
                grid_vars.edge3_dimensions = self.edge3_dimensions
            if self.edge3_coordinates is not None:
                grid_vars.edge3_coordinates = ' '.join(self.edge3_coordinates)
    

class SGridAttributes(object):
    """
    Class containing methods to help with getting the
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
        
    def get_topology_var(self):
        grid_topology_var = self.ncd.find_grid_topology_var()
        return grid_topology_var
    
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
            grid_cell_node_vars = self.ncd.find_node_coordinates(node_dimensions)
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
            setattr(sgrid, sgrid_var.variable, sgrid_var)
            dataset_variables.append(nc_var.name)
            if hasattr(nc_var, 'grid'):
                grid_variables.append(nc_var.name)
        sgrid.variables = dataset_variables
        sgrid.grid_variables = grid_variables
        
    def get_angles(self):
        try:
            # get angles if they exist, otherwise calculate them
            grid_angles = self.nc_dataset.variables['angle'][:]
            angles = grid_angles
        except KeyError:
            cell_centers = self.get_cell_center_lat_lon()
            centers_start = cell_centers[..., :-1, :]
            centers_end = cell_centers[..., 1:, :]
            angles = calculate_angle_from_true_east(centers_start, centers_end)
        return angles
        
    def get_cell_center_lat_lon(self):
        grid_cell_center_lon_var, grid_cell_center_lat_var = self.get_attr_coordinates('face_coordinates')
        grid_cell_center_lat = self.nc_dataset.variables[grid_cell_center_lat_var][:]
        grid_cell_center_lon = self.nc_dataset.variables[grid_cell_center_lon_var][:]
        return pair_arrays(grid_cell_center_lon, grid_cell_center_lat)
        
    def get_cell_node_lat_lon(self):
        try:
            grid_cell_nodes_lon_var, grid_cell_nodes_lat_var = self.get_node_coordinates()[1]
        except TypeError:
            cell_nodes = None
        else:
            grid_cell_nodes_lat = self.nc_dataset.variables[grid_cell_nodes_lat_var][:]
            grid_cell_nodes_lon = self.nc_dataset.variables[grid_cell_nodes_lon_var][:]
            cell_nodes = pair_arrays(grid_cell_nodes_lon, grid_cell_nodes_lat)
        return cell_nodes
        
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
        grid = SGrid2D.from_nc_dataset(nc_dataset, grid_topology_var)
    elif topology_dim == 3:
        grid = SGrid3D.from_nc_dataset(nc_dataset, grid_topology_var)
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
    
    
def from_ncfile(nc_url, grid_topology_var=None):
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