'''
Created on Mar 19, 2015

@author: ayan
'''
import netCDF4 as nc4
from custom_exceptions import SGridNonCompliant
from utils import ParsePadding, determine_variable_slicing
from lookup import LAT_GRID_CELL_CENTER_LONG_NAME, LON_GRID_CELL_CENTER_LONG_NAME


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

    def find_grid_cell_center_vars(self):
        nc_vars = self.ncd.variables
        grid_cell_center_lon = None
        grid_cell_center_lat = None
        for nc_var in nc_vars.iterkeys():
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
        return grid_cell_center_lat, grid_cell_center_lon
    
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
        for nc_var in nc_vars.iterkeys():
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


def _set_attributes_from_list(target_object, attribute_val_list):
    for attribute_val_element in attribute_val_list:
        attr_name, attr_val = attribute_val_element
        setattr(target_object, attr_name, attr_val)


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
                grid.face_padding = face_dim_padding  # set face padding
                face_slices = pp.define_recommended_slices(face_dim)
                _set_attributes_from_list(grid, face_slices)
            except AttributeError:
                pass
            try:
                edge1_dim = nc_grid_topology_var.edge1_dimensions
                edge1_dim_padding = pp.parse_padding(edge1_dim)
                grid.edge_1_padding = edge1_dim_padding  # set edge 1 padding
                edge_1_slices = pp.define_recommended_slices(edge1_dim)
                _set_attributes_from_list(grid, edge_1_slices)
            except AttributeError:
                pass
            try:
                edge2_dim = nc_grid_topology_var.edge2_dimensions
                edge2_dim_padding = pp.parse_padding(edge2_dim)
                grid.edge_2_padding = edge2_dim_padding  # set edge 2 padding
                edge_2_slices = pp.define_recommended_slices(edge2_dim)
                _set_attributes_from_list(grid, edge_2_slices)
            except AttributeError:
                pass
            try:
                vertical_dim = nc_grid_topology_var.vertical_dimensions
                vertical_dim_padding = pp.parse_padding(vertical_dim)
                grid.vertical_padding = vertical_dim_padding  # set vertical padding
                vertical_slices = pp.define_recommended_slices(vertical_dim)
                _set_attributes_from_list(grid, vertical_slices)
            except AttributeError:
                pass
        grid_cell_center_vars = ncd.find_grid_cell_center_vars()  # get the variable names for the cell center
        grid.grid_cell_center_vars = grid_cell_center_vars  # set the variables for the grid cell centers in the sgrid object
        grid_cell_center_lat_vals = nc_dataset.variables[grid_cell_center_vars[0]][:]
        grid_cell_center_lon_vals = nc_dataset.variables[grid_cell_center_vars[1]][:]
        grid.grid_cell_center_lat = grid_cell_center_lat_vals  # set the grid cell center latitudes
        grid.grid_cell_center_lon = grid_cell_center_lon_vals  # set the grid cell center longitudes
        grid_time = nc_dataset.variables['time'][:]
        grid.grid_times = grid_time
        nc_variables = nc_dataset.variables
        # dynamically set variable slicing attributes
        for nc_variable in nc_variables:
            # the slicing implied by the padding for each variable
            # each dimension is sliced in the order they appear in ncdump
            var_slicing = determine_variable_slicing(grid, nc_dataset, nc_variable)
            property_name = '{0}_slice'.format(nc_variable)
            grid.add_property(property_name, var_slicing)
        return grid
    else:
        raise SGridNonCompliant(nc_dataset)