'''
Created on Mar 19, 2015

@author: ayan
'''
from .lookup import LAT_GRID_CELL_NODE_LONG_NAME, LON_GRID_CELL_NODE_LONG_NAME


class NetCDFDataset(object):
    
    def __init__(self, nc_dataset_obj):
        self.ncd = nc_dataset_obj
    
    def find_grid_cell_node_vars(self):
        """
        Find the variables for the grid
        cell vertices.
        
        """
        nc_vars = self.ncd.variables
        grid_cell_node_lon = None
        grid_cell_node_lat = None
        for nc_var in nc_vars.keys():
            nc_var_obj = nc_vars[nc_var]
            try:
                nc_var_long_name = nc_var_obj.long_name
            except AttributeError:
                continue
            else:
                if nc_var_long_name in LON_GRID_CELL_NODE_LONG_NAME:
                    grid_cell_node_lon = nc_var
                if nc_var_long_name in LAT_GRID_CELL_NODE_LONG_NAME:
                    grid_cell_node_lat = nc_var
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
            except AttributeError:
                cf_role = None
                topology_dim = None
            else:
                topology_dim = nc_var_obj.topology_dimension
            if cf_role == 'grid_topology' and topology_dim >= 2:
                grid_topology_vars.append(nc_var)
        if len(grid_topology_vars) > 0:
            grid_topology_var = grid_topology_vars[0]
        else:
            grid_topology_var = None
        return grid_topology_var
    
    def search_variables_by_location(self, location_str):
        nc_vars = self.ncd.variables
        search_results = []
        for nc_var in nc_vars.keys():
            nc_var_obj = nc_vars[nc_var]
            try:
                nc_var_location = nc_var_obj.location
            except AttributeError:
                continue
            else:
                if nc_var_location == location_str:
                    search_results.append(nc_var)
        return search_results
    
    def find_coordinates_by_location(self, location_str, topology_dim):
        """
        Find a grid coordinates variables with a location attribute equal
        to location_str. This method can be used to infer edge, face, or
        volume coordinates from the location attribute of a variable.
        
        Location is a required attribute per SGRID conventions.
        
        :param str location_str: the location value to search for
        :param int topology_dim: the topology dimension of the grid
        
        """
        nc_vars = self.ncd.variables
        vars_with_location = self.search_variables_by_location(location_str)
        x_coordinate = None
        y_coordinate = None
        z_coordinate = None
        for var_with_location in vars_with_location:
            location_var = nc_vars[var_with_location]
            location_var_dims = location_var.dimensions
            try:
                location_var_coordinates = location_var.coordinates
            except AttributeError:
                # run through this if a location attributed is defined, but not coordinates
                potential_coordinates = []
                for nc_var in nc_vars.keys():
                    nc_var_obj = nc_vars[nc_var]
                    nc_var_dim_set = set(nc_var_obj.dimensions)
                    if (nc_var_dim_set.issubset(location_var_dims) and 
                        nc_var != var_with_location and 
                        len(nc_var_dim_set) > 0
                        ):
                        potential_coordinates.append(nc_var_obj)
                for potential_coordinate in potential_coordinates:
                    pc_name = potential_coordinate.name
                    if 'lon' in pc_name.lower():
                        x_coordinate = pc_name
                    elif 'lat' in pc_name.lower():
                        y_coordinate = pc_name
                    else:
                        z_coordinate = pc_name  # this might not always work...
            else:
                lvc_split = location_var_coordinates.strip().split(' ')
                for lvc in lvc_split:
                    var_coord = nc_vars[lvc]
                    try:
                        var_coord_standard_name = var_coord.standard_name
                    except AttributeError:
                        if 'lon' in var_coord.name.lower():
                            x_coordinate = lvc
                        elif 'lat' in var_coord.name.lower():
                            y_coordinate = lvc
                    else:
                        if var_coord_standard_name == 'longitude':
                            x_coordinate = lvc
                        elif var_coord_standard_name == 'latitude':
                            y_coordinate = lvc
                if len(lvc_split) == 3:
                    z_coordinate = lvc_split[-1]
                break 
        if topology_dim == 2:
            coordinates = (x_coordinate, y_coordinate)
        else:
            coordinates = (x_coordinate, y_coordinate, z_coordinate)
        if all(coordinates):
            coordinate_result = coordinates
        else:
            coordinate_result = None
        return coordinate_result

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
        if grid_vars is not None:
            sgrid_compliant = True
        else:
            sgrid_compliant = False
        return sgrid_compliant