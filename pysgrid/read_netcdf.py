'''
Created on Mar 19, 2015

@author: ayan
'''
import netCDF4 as nc4
from .custom_exceptions import SGridNonCompliant, deprecated
from .utils import ParsePadding, pair_arrays, determine_variable_slicing
from .lookup import (LAT_GRID_CELL_CENTER_LONG_NAME, LON_GRID_CELL_CENTER_LONG_NAME,
                     LAT_GRID_CELL_NODE_LONG_NAME, LON_GRID_CELL_NODE_LONG_NAME)


def read_netcdf_file(dataset_url):
    """
    Read a netCDF file into a dataset
    object.
    
    :param str dataset_url: path or URL to a netCDF file
    :return: netCDF dataset object
    :rtype: netCDF4.Dataset
    
    """
    nc_dataset = nc4.Dataset(dataset_url)
    return nc_dataset


class NetCDFDataset(object):
    
    def __init__(self, nc_dataset_obj):
        self.ncd = nc_dataset_obj
    
    @deprecated
    def find_grid_cell_center_vars(self):
        """
        Find the variables for the grid
        cell centers.
        
        """
        nc_vars = self.ncd.variables
        grid_cell_center_lon = None
        grid_cell_center_lat = None
        for nc_var in nc_vars.keys():
            try:
                nc_var_obj = nc_vars[nc_var]
                # need to revisit this... long_name is not a required attribute
                nc_var_long_name = nc_var_obj.long_name
                if nc_var_long_name in LON_GRID_CELL_CENTER_LONG_NAME:
                    grid_cell_center_lon = nc_var
                if nc_var_long_name in LAT_GRID_CELL_CENTER_LONG_NAME:
                    grid_cell_center_lat = nc_var
            except AttributeError:
                continue
        return grid_cell_center_lon, grid_cell_center_lat
    
    def find_grid_cell_node_vars(self):
        """
        Find the variables for the grid
        cell vertices.
        
        """
        nc_vars = self.ncd.variables
        grid_cell_node_lon = None
        grid_cell_node_lat = None
        for nc_var in nc_vars.keys():
            try:
                nc_var_obj = nc_vars[nc_var]
                nc_var_long_name = nc_var_obj.long_name
                if nc_var_long_name in LON_GRID_CELL_NODE_LONG_NAME:
                    grid_cell_node_lon = nc_var
                if nc_var_long_name in LAT_GRID_CELL_NODE_LONG_NAME:
                    grid_cell_node_lat = nc_var
            except AttributeError:
                continue
        return grid_cell_node_lon, grid_cell_node_lat
        
    def find_grid_topology_vars(self):
        """
        Get the variables from a netCDF dataset
        that have a cf_role attribute of 'grid_topology'.
        
        :params nc: netCDF dataset
        :type nc: netCDF4.Dataset
        :return: list of variables that contain grid topologies
        :rtype: list
        
        """
        nc_vars = self.ncd.variables
        grid_topology_vars = []
        for nc_var in nc_vars.keys():
            nc_var_obj = nc_vars[nc_var]
            try:
                cf_role = nc_var_obj.cf_role.strip()
                topology_dim = nc_var_obj.topology_dimension
            except AttributeError:
                cf_role = None
                topology_dim = None
            if cf_role == 'grid_topology' and topology_dim == 2:
                grid_topology_vars.append(nc_var)
        return grid_topology_vars
    
    def find_coordinations_by_location(self, location_str):
        """
        Find a variable with a location attribute equal
        to location_str.
        
        Location is a required attribute per SGRID conventions.
        
        :param str location_str: the location value to search for
        
        """
        nc_vars = self.ncd.variables
        for nc_var in nc_vars.keys():
            print(nc_var)
            nc_var_obj = nc_vars[nc_var]
            try:
                nc_var_location = nc_var_obj.location
                if nc_var_location == location_str:
                    print(nc_var)
                    nc_var_coordinates = nc_var_obj.coordinates
                    nc_var_coord_split = nc_var_coordinates.strip().split(' ')
                    print(nc_var_coord_split)
                    x_coordinate = None
                    y_coordinate = None
                    for nc_var_coord in nc_var_coord_split:
                        print(nc_var_coord)
                        var_coord = nc_vars[nc_var_coord]
                        try:
                            var_coord_standard_name = var_coord.standard_name
                            if var_coord_standard_name == 'longitude':
                                x_coordinate = nc_var_coord
                            elif var_coord_standard_name == 'latitude':
                                y_coordinate = nc_var_coord
                        except AttributeError:
                            continue
                    result = (x_coordinate, y_coordinate)
                    break
                else:
                    result = None
            except AttributeError:
                result = None
                continue
        return result

    def sgrid_compliant_file(self):
        """
        Determine whether a dataset is
        SGRID compliant.
        
        :param nc: netCDF dataset
        :type nc: netCDF4.Dataset
        :return: True if dataset is compliant, False if it is not
        :rtype: bool
        
        """
        grid_vars = self.find_grid_topology_vars()
        if len(grid_vars) > 0:
            sgrid_compliant = True
        else:
            sgrid_compliant = False
        return sgrid_compliant


def load_grid_from_nc_file(nc_path, grid, grid_topology_vars=None, load_data=True):
    """
    Create a SGRID object from a path to an
    SGRID compliant netCDF resource. An 
    exception is raised if the resource is
    found to be non-compliant.
    
    :param str nc_path: path to the resource; this can be a filepath or a URL
    :param grid: an SGRID object
    :type grid: sgrid.SGrid
    :return: an SGrid object
    :rtype: sgrid.SGrid
    
    """
    with nc4.Dataset(nc_path, 'r') as nc_dataset:
        grid = load_grid_from_nc_dataset(nc_dataset, grid, 
                                         grid_topology_vars=grid_topology_vars, 
                                         load_data=load_data
                                         )
    return grid


def load_grid_from_nc_dataset(nc_dataset, grid, 
                              grid_topology_vars=None, 
                              load_data=True):
    """
    Create an SGRID object from an SGRID
    compliant netCDF4.Dataset object. An
    exception is raised if the dataset is
    non-compliant.
    
    :param nc_dataset: a netCDF resource read into a netCDF4.Dataset object
    :type nc_dataset: netCDF4.Dataset
    :param grid: an SGRID object
    :type grid: sgrid.SGrid
    :return: an SGrid object
    :rtype: sgrid.SGrid
    
    """
    ncd = NetCDFDataset(nc_dataset)
    is_sgrid_compliant = ncd.sgrid_compliant_file()
    if is_sgrid_compliant:
        ds_dims = nc_dataset.dimensions
        grid_dims = [(ds_dim, len(ds_dims[ds_dim])) for ds_dim in ds_dims]
        grid.dimensions = grid_dims
        if grid_topology_vars is None:
            grid_topology_vars_attr = ncd.find_grid_topology_vars()
        else:
            grid_topology_vars_attr = grid_topology_vars
        grid.grid_topology_vars = grid_topology_vars_attr  # set grid variables 
        for topology_var in grid_topology_vars_attr:
            nc_grid_topology_var = nc_dataset.variables[topology_var]
            pp = ParsePadding(topology_var)
            try:
                face_dim = nc_grid_topology_var.face_dimensions
                face_dim_padding = pp.parse_padding(face_dim)
                grid.face_dimensions = face_dim
                grid.face_padding = face_dim_padding  # set face padding
            except AttributeError:
                pass
            try:
                edge1_dim = nc_grid_topology_var.edge1_dimensions
                edge1_dim_padding = pp.parse_padding(edge1_dim)
                grid.edge_1_dimension = edge1_dim
                grid.edge_1_padding = edge1_dim_padding  # set edge 1 padding
            except AttributeError:
                pass
            try:
                edge2_dim = nc_grid_topology_var.edge2_dimensions
                edge2_dim_padding = pp.parse_padding(edge2_dim)
                grid.edge_2_dimensions = edge2_dim
                grid.edge_2_padding = edge2_dim_padding  # set edge 2 padding
            except AttributeError:
                pass
            try:
                vertical_dim = nc_grid_topology_var.vertical_dimensions
                vertical_dim_padding = pp.parse_padding(vertical_dim)
                grid.vertical_dimensions = vertical_dim
                grid.vertical_padding = vertical_dim_padding  # set vertical padding
            except AttributeError:
                pass
            try:
                face_coordinates = nc_grid_topology_var.face_coordinates
                face_coordinate_val = face_coordinates.split(' ')
                grid.face_coordinates = tuple(face_coordinate_val)
            except AttributeError:
                grid_cell_center_vars = ncd.find_coordinations_by_location('face')
                grid.face_coordinates = grid_cell_center_vars
            try:
                node_coordinates = nc_grid_topology_var.node_coordinates
                node_coordinate_val = node_coordinates.split(' ')
                grid.node_coordinates = tuple(node_coordinate_val)
            except AttributeError:
                # grid_cell_node_vars = ncd.find_grid_cell_node_vars()
                # grid.node_coordinates = grid_cell_node_vars
                pass
            try:
                edge_1_coordinates = nc_grid_topology_var.edge1_coordinates
                edge_1_coordinates_val = edge_1_coordinates.split(' ')
                grid.edge_1_coordinates = tuple(edge_1_coordinates_val)
            except AttributeError:
                edge_1_coordinates_val = ncd.find_coordinations_by_location('edge1')
                grid.edge_1_coordinates = edge_1_coordinates_val
            try:
                edge_2_coordinates = nc_grid_topology_var.edge2_coordinates
                edge_2_coordinates_val = edge_2_coordinates.split(' ')
                grid.edge_2_coordinates = tuple(edge_2_coordinates_val)
            except AttributeError:
                edge_2_coordinates_val = ncd.find_coordinations_by_location('edge2')
                grid.edge_2_coordinates = edge_2_coordinates_val
        grid_cell_center_lon_var, grid_cell_center_lat_var = grid.face_coordinates
        grid_cell_center_lat = nc_dataset.variables[grid_cell_center_lat_var][:]
        grid_cell_center_lon = nc_dataset.variables[grid_cell_center_lon_var][:]
        grid.centers = pair_arrays(grid_cell_center_lon, grid_cell_center_lat)
        # get the variables names for the cell vertices
        grid_cell_nodes_lat_var, grid_cell_nodes_lon_var = ncd.find_grid_cell_node_vars()
        grid_cell_nodes_lat = nc_dataset.variables[grid_cell_nodes_lat_var][:]
        grid_cell_nodes_lon = nc_dataset.variables[grid_cell_nodes_lon_var][:]
        grid.nodes = pair_arrays(grid_cell_nodes_lon, grid_cell_nodes_lat)
        grid.node_dimensions = nc_grid_topology_var.node_dimensions
        # get time data
        grid_time = nc_dataset.variables['time'][:]
        nc_variables = nc_dataset.variables
        # provide a list of all variables in the netCDF dataset
        grid.grid_times = grid_time
        grid_variables = []
        for nc_variable in nc_variables:
            nc_var = nc_variables[nc_variable]
            nc_var_name = nc_var.name
            nc_var_dtype = nc_var.dtype
            nc_var_dims = nc_var.dimensions
            grid_var = (nc_var_name, nc_var_dtype, nc_var_dims)
            grid_variables.append(grid_var)
        grid.variables = grid_variables
        # provide the angles
        try:
            grid_angles = nc_dataset.variables['angle'][:]
            grid.angles = grid_angles
        except KeyError:
            pass
        # dynamically set variable slicing attributes
        for nc_variable in nc_variables:
            # the slicing implied by the padding for each variable
            # each dimension is sliced in the order they appear in ncdump
            var_slicing = determine_variable_slicing(grid, nc_dataset, nc_variable)
            slice_property_name = '{0}_slice'.format(nc_variable)
            grid.add_property(slice_property_name, var_slicing)  # add slice property
            dim_property_name = '{0}_dim'.format(nc_variable)
            var_dims = nc_variables[nc_variable].dimensions
            grid.add_property(dim_property_name, var_dims)  # add dimension property
        return grid
    else:
        raise SGridNonCompliant(nc_dataset)