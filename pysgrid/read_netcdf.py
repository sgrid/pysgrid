'''
Created on Mar 19, 2015

@author: ayan
'''
import netCDF4 as nc4
from custom_exceptions import SGridNonCompliant


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


def find_grid_topology_vars(nc):
    """
    Get the variables from a netCDF dataset
    that have a cf_role attribute of 'grid_topology'.
    
    :params nc: netCDF dataset
    :type nc: netCDF4.Dataset
    :return: list of variables that contain grid topologies
    :rtype: list
    
    """
    nc_vars = nc.variables
    grid_topology_vars = []
    for nc_var in nc_vars.iterkeys():
        nc_var_obj = nc.variables[nc_var]
        try:
            cf_role = nc_var_obj.cf_role.strip()
            topology_dim = nc_var_obj.topology_dimension
        except AttributeError:
            cf_role = None
            topology_dim = None
        if cf_role == 'grid_topology' and topology_dim == 2:
            grid_topology_vars.append(nc_var)
    return grid_topology_vars


def sgrid_compliant_file(nc):
    """
    Determine whether a dataset is
    SGRID compliant.
    
    :param nc: netCDF dataset
    :type nc: netCDF4.Dataset
    :return: True if dataset is compliant, False if it is not
    :rtype: bool
    
    """
    grid_vars = find_grid_topology_vars(nc)
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
    :param grid:an SGRID object
    :type grid: sgrid.SGrid
    :return: an SGrid object
    :rtype: sgrid.SGrid
    
    """
    is_sgrid_compliant = sgrid_compliant_file(nc_dataset)
    if is_sgrid_compliant:
        if grid_topology_vars is None:
            grid_topology_vars_attr = find_grid_topology_vars(nc_dataset)
        else:
            grid_topology_vars_attr = grid_topology_vars
        grid.grid_topology_vars = grid_topology_vars_attr
        return grid
    else:
        raise SGridNonCompliant(nc_dataset)