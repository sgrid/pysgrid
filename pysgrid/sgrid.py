'''
Created on Apr 20, 2015

@author: ayan
'''

from __future__ import (absolute_import, division, print_function)

from netCDF4 import Dataset

from .read_netcdf import NetCDFDataset, parse_padding, find_grid_topology_var
from .utils import calculate_angle_from_true_east, pair_arrays
from .variables import SGridVariable


class SGrid(object):

    padding_slices = {'both': (1, -1),
                      'none': (None, None),
                      'low': (1, None),
                      'high': (None, 1)
                      }

    topology_dimension = 2

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
                 edge2_dimensions=None,
                 faces=None,
                 face_padding=None,
                 face_coordinates=None,
                 face_dimensions=None,
                 vertical_padding=None,
                 vertical_dimensions=None,
                 *args,
                 **kwargs):

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
        self.faces = faces
        self.face_padding = face_padding
        self.face_coordinates = face_coordinates
        self.face_dimensions = face_dimensions
        self.vertical_padding = vertical_padding
        self.vertical_dimensions = vertical_dimensions

    @classmethod
    def load_sgrid(cls, nc):
        if isinstance(nc, Dataset):
            pass
        else:
            nc = Dataset(nc, 'r')
        topology_var = find_grid_topology_var(nc)
        sa = SGridAttributes(nc, cls.topology_dimension, topology_var)
        dimensions = sa.get_dimensions()
        node_dimensions, node_coordinates = sa.get_node_coordinates()
        grid_topology_var = sa.get_topology_var()
        edge1_dimensions, edge1_padding = sa.get_attr_dimension('edge1_dimensions')  # noqa
        edge2_dimensions, edge2_padding = sa.get_attr_dimension('edge2_dimensions')  # noqa
        edge1_coordinates = sa.get_attr_coordinates('edge1_coordinates')
        edge2_coordinates = sa.get_attr_coordinates('edge2_coordinates')
        angles = sa.get_angles()
        vertical_dimensions, vertical_padding = sa.get_attr_dimension('vertical_dimensions')  # noqa
        centers = sa.get_cell_center_lat_lon()
        face_dimensions, face_padding = sa.get_attr_dimension('face_dimensions')  # noqa
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
        with Dataset(filepath, 'w') as nclocal:
            grid_vars = self._save_common_components(nclocal)
            # Add attributes to the grid_topology variable.
            grid_vars.face_dimensions = self.face_dimensions
            if self.vertical_dimensions is not None:
                grid_vars.vertical_dimensions = self.vertical_dimensions
            if self.face_coordinates is not None:
                grid_vars.face_coordinates = ' '.join(self.face_coordinates)

    @property
    def non_grid_variables(self):
        non_grid_variables = [variable for variable in self.variables if
                              variable not in self.grid_variables]
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
        grid_vars = nc_file.createVariable(grid_var_obj.variable,
                                           grid_var_obj.dtype)
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
        for dataset_variable in self.variables:
            dataset_var_obj = getattr(self, dataset_variable)
            try:
                dataset_grid_var = nc_file.createVariable(
                    dataset_var_obj.variable,
                    dataset_var_obj.dtype,
                    dataset_var_obj.dimensions
                )
            except RuntimeError:
                continue
            else:
                axes = []
                if dataset_var_obj.grid is not None:
                    dataset_grid_var.grid = grid_var
                if dataset_var_obj.standard_name is not None:
                    dataset_grid_var.standard_name = dataset_var_obj.standard_name  # noqa
                if dataset_var_obj.coordinates is not None:
                    dataset_grid_var.coordinates = ' '.join(dataset_var_obj.coordinates)  # noqa
                if dataset_var_obj.x_axis is not None:
                    x_axis = 'X: {0}'.format(dataset_var_obj.x_axis)
                    axes.append(x_axis)
                if dataset_var_obj.y_axis is not None:
                    y_axis = 'Y: {0}'.format(dataset_var_obj.y_axis)
                    axes.append(y_axis)
                if dataset_var_obj.z_axis is not None:
                    z_axis = 'Z: {0}'.format(dataset_var_obj.z_axis)
                    axes.append(z_axis)
                if axes:
                    dataset_grid_var.axes = ' '.join(axes)
        return grid_vars


class SGridAttributes(object):
    """
    Class containing methods to help with getting the
    attributes for either SGrid.

    """
    def __init__(self, nc, topology_dim, topology_variable):
        self.nc = nc
        self.ncd = NetCDFDataset(self.nc)
        self.topology_dim = topology_dim
        self.topology_variable = topology_variable
        self.topology_var = self.nc.variables[self.topology_variable]

    def get_dimensions(self):
        ds_dims = self.nc.dimensions
        grid_dims = [(ds_dim, len(ds_dims[ds_dim])) for ds_dim in ds_dims]
        return grid_dims

    def get_topology_var(self):
        grid_topology_var = find_grid_topology_var(self.nc)
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
            attr_coordinates = self.ncd.find_coordinates_by_location(location_name, self.topology_dim)  # noqa
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
            grid_cell_node_vars = self.ncd.find_node_coordinates(node_dimensions)  # noqa
            node_coordinates = grid_cell_node_vars
        else:
            node_coordinate_val = node_coordinates.split(' ')
            node_coordinates = tuple(node_coordinate_val)
        return node_dimensions, node_coordinates

    def get_variable_attributes(self, sgrid):
        dataset_variables = []
        grid_variables = []
        nc_variables = self.nc.variables
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
        angles = self.nc.variables.get('angle')
        if not angles:
            cell_centers = self.get_cell_center_lat_lon()
            centers_start = cell_centers[..., :-1, :]
            centers_end = cell_centers[..., 1:, :]
            angles = calculate_angle_from_true_east(centers_start, centers_end)
        return angles

    def get_cell_center_lat_lon(self):
        grid_cell_center_lon_var, grid_cell_center_lat_var = self.get_attr_coordinates('face_coordinates')  # noqa
        grid_cell_center_lat = self.nc[grid_cell_center_lat_var]
        grid_cell_center_lon = self.nc[grid_cell_center_lon_var]
        return pair_arrays(grid_cell_center_lon, grid_cell_center_lat)

    def get_cell_node_lat_lon(self):
        try:
            grid_cell_nodes_lon_var, grid_cell_nodes_lat_var = self.get_node_coordinates()[1]  # noqa
        except TypeError:
            cell_nodes = None
        else:
            grid_cell_nodes_lat = self.nc[grid_cell_nodes_lat_var]
            grid_cell_nodes_lon = self.nc[grid_cell_nodes_lon_var]
            cell_nodes = pair_arrays(grid_cell_nodes_lon, grid_cell_nodes_lat)
        return cell_nodes


def load_sgrid(nc):
    """
    Get a SGrid object from a netCDF4.Dataset or file/URL.

    :param str or netCDF4.Dataset nc: a netCDF4 Dataset or URL/filepath
                                       to the netCDF file
    :return: SGrid object
    :rtype: sgrid.SGrid

    """
    if isinstance(nc, Dataset):
        pass
    else:
        nc = Dataset(nc, 'r')

    return SGrid.load_sgrid(nc)
