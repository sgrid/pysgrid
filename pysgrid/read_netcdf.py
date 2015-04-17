'''
Created on Mar 19, 2015

@author: ayan
'''
import netCDF4 as nc4
from .custom_exceptions import SGridNonCompliantError, TopologyDimensionError, deprecated
from .utils import ParsePadding, pair_arrays, determine_variable_slicing
from .variables import SGridVariable
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
                if nc_var_location == location_str:
                    search_results.append(nc_var)
            except AttributeError:
                continue
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
                lvc_split = location_var_coordinates.split(' ')
                if len(lvc_split) == 2:
                    x_coordinate, y_coordinate = lvc_split
                elif len(lvc_split) == 3:
                    x_coordinate, y_coordinate, z_coordinate = lvc_split
                else:
                    raise TopologyDimensionError(topology_dim)
                break
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
        if topology_dim == 2:
            coordinates = (x_coordinate, y_coordinate)
        else:
            coordinates = (x_coordinate, y_coordinate, z_coordinate)
        if all(coordinates):
            coordinate_result = coordinates
        else:
            coordinate_result = None
        return coordinate_result
    
    @deprecated
    def find_coordinations_by_location(self, location_str, topology_dim):
        """
        Find a variable with a location attribute equal
        to location_str.
        
        Location is a required attribute per SGRID conventions.
        
        :param str location_str: the location value to search for
        :param int topology_dim: the topology dimension of the grid
        """
        nc_vars = self.ncd.variables
        vars_with_location = self.search_variables_by_location(location_str)
        result = None
        for nc_var in vars_with_location:
            nc_var_obj = nc_vars[nc_var]
            try:
                nc_var_coordinates = nc_var_obj.coordinates
                nc_var_coord_split = nc_var_coordinates.strip().split(' ')
                x_coordinate = None
                y_coordinate = None
                z_coordinate = None
                for nc_var_coord in nc_var_coord_split:
                    var_coord = nc_vars[nc_var_coord]
                    try:
                        var_coord_standard_name = var_coord.standard_name
                        if var_coord_standard_name == 'longitude':
                            x_coordinate = nc_var_coord
                        elif var_coord_standard_name == 'latitude':
                            y_coordinate = nc_var_coord
                    except AttributeError:
                        continue
                if topology_dim == 2:
                    result = (x_coordinate, y_coordinate)
                elif topology_dim == 3:
                    # finding the z_coordinate isn't fully baked yet
                    result = (x_coordinate, y_coordinate, z_coordinate)
                else:
                    raise Exception('I have no idea what to do....')
                break
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
        if grid_vars is not None:
            sgrid_compliant = True
        else:
            sgrid_compliant = False
        return sgrid_compliant
    
    
class SGridND(object):
    
    topology_dim = None
    
    def __init__(self, sgrid, nc_dataset, topology_variable):
        self._sgrid = sgrid
        self.nc_dataset = nc_dataset
        self.ncd = NetCDFDataset(self.nc_dataset)
        self.topology_variable = topology_variable  # the netCDF variable with a cf_role of 'grid_topology'
        self.topology_var = self.nc_dataset.variables[self.topology_variable]
        self.pp = ParsePadding(self.topology_variable)
        
    @property
    def sgrid(self):
        return self._sgrid
    
    def set_dimensions(self):
        ds_dims = self.nc_dataset.dimensions
        grid_dims = [(ds_dim, len(ds_dims[ds_dim])) for ds_dim in ds_dims]
        self._sgrid.dimensions = grid_dims
        
    def set_topology_dimension(self):
        self._sgrid.topology_dimension = self.topology_dim
    
    def set_edge1_dimensions(self):
        try:
            edge1_dim = self.topology_var.edge1_dimensions
            edge1_dim_padding = self.pp.parse_padding(edge1_dim)
            self._sgrid.edge_1_dimension = edge1_dim
            self._sgrid.edge_1_padding = edge1_dim_padding
        except AttributeError:
            pass
        
    def set_edge1_coordinates(self):
        try:
            edge1_coordinates = self.topology_var.edge1_coordinates
            edge1_coordinates_val = edge1_coordinates.split(' ')
            self._sgrid.edge_1_coordinates = tuple(edge1_coordinates_val)
        except AttributeError:
            pass
        
    def set_edge2_dimensions(self):
        try:
            edge2_dim = self.topology_var.edge2_dimensions
            edge2_dim_padding = self.pp.parse_padding(edge2_dim)
            self._sgrid.edge_2_dimension = edge2_dim
            self._sgrid.edge_2_padding = edge2_dim_padding
        except AttributeError:
            pass
        
    def set_edge2_coordinates(self):
        try:
            edge2_coordinates = self.topology_var.edge2_coordinates
            edge2_coordinates_val = edge2_coordinates.split(' ')
            self._sgrid.edge_2_coordinates = tuple(edge2_coordinates_val)
        except AttributeError:
            pass
        
    def set_all_edge_attributes(self):
        self.set_edge1_dimensions()
        self.set_edge1_coordinates()
        self.set_edge2_dimensions()
        self.set_edge2_coordinates()
        
    def set_sgrid_topology(self):
        topology_dim = self.topology_var.topology_dimension
        self._sgrid.topology_dimension = topology_dim
        
    def set_sgrid_vertical_dimensions(self):
        try:
            vertical_dim = self.topology_var.vertical_dimensions
            vertical_dim_padding = self.pp.parse_padding(vertical_dim)
            self._sgrid.vertical_dimensions = vertical_dim
            self._sgrid.vertical_padding = vertical_dim_padding
        except AttributeError:
            pass
        
    def set_sgrid_node_coordinates(self):
        try:
            node_coordinates = self.topology_var.node_coordinates
            node_coordinate_val = node_coordinates.split(' ')
            self._sgrid.node_coordinates = tuple(node_coordinate_val)
        except AttributeError:
            grid_cell_node_vars = self.ncd.find_grid_cell_node_vars()
            self._sgrid.node_coordinates = grid_cell_node_vars
            
    def set_sgrid_variable_attributes(self):
        dataset_variables = []
        grid_variables = []
        nc_variables = self.nc_dataset.variables
        for nc_variable in nc_variables:
            nc_var = nc_variables[nc_variable]
            nc_var_name = nc_var.name
            dataset_variables.append(nc_var_name)
            sgrid_var = SGridVariable.create_variable(nc_var)
            var_center_slicing = determine_variable_slicing(self._sgrid,
                                                            self.nc_dataset,
                                                            nc_variable,
                                                            method='center')
            sgrid_var.center_slicing = var_center_slicing
            self._sgrid.add_property(sgrid_var.variable, sgrid_var)
            try:
                if nc_var.grid:
                    grid_variables.append(nc_var_name)
            except AttributeError:
                continue
        self._sgrid.variables = dataset_variables
        self._sgrid.grid_variables = grid_variables
        
    def set_sgrid_angles(self):
        try:
            grid_angles = self.nc_dataset.variables['angle'][:]
            self._sgrid.angles = grid_angles
        except KeyError:
            pass
        
    def set_sgrid_time(self):
        try:
            grid_time = self.nc_dataset.variables['time'][:]
        except KeyError:
            grid_time = self.nc_dataset.variables['Times'][:]
        self._sgrid.grid_times = grid_time
    
    def set_sgrid_nd_attributes(self):
        self.set_sgrid_topology()
        # set vertical dimensions
        self.set_sgrid_vertical_dimensions()
        # set node coordinates
        self.set_sgrid_node_coordinates()
        # set variables
        self.set_sgrid_variable_attributes()
        # set the angles
        self.set_sgrid_angles()
        # set time
        self.set_sgrid_time()
        
    def set_all_face_attributes(self):
        raise NotImplementedError
    
    def set_cell_center_lat_lon(self):
        raise NotImplementedError
    
    def set_cell_node_lat_lon(self):
        raise NotImplementedError
        
        
class SGrid2D(SGridND):
    
    topology_dim = 2
    
    def set_face_dimensions(self):
        try:
            face_dim = self.topology_var.face_dimensions
            face_dim_padding = self.pp.parse_padding(face_dim)
            self._sgrid.face_dimensions = face_dim
            self._sgrid.face_padding = face_dim_padding
        except AttributeError:
            pass
        
    def set_face_coordindates(self):
        try:
            face_coordinates = self.topology_var.face_coordinates
            face_coordinate_val = face_coordinates.split(' ')
            self._sgrid.face_coordinates = tuple(face_coordinate_val)
        except AttributeError:
            grid_cell_center_vars = self.ncd.find_coordinates_by_location('face', self.topology_dim)
            self._sgrid.face_coordinates = grid_cell_center_vars
            
    def set_all_face_attributes(self):
        self.set_face_dimensions()
        self.set_face_coordindates()
        
    def set_cell_center_lat_lon(self):
        grid_cell_center_lon_var, grid_cell_center_lat_var = self._sgrid.face_coordinates
        grid_cell_center_lat = self.nc_dataset.variables[grid_cell_center_lat_var][:]
        grid_cell_center_lon = self.nc_dataset.variables[grid_cell_center_lon_var][:]
        self._sgrid.centers = pair_arrays(grid_cell_center_lon, grid_cell_center_lat)
        
    def set_cell_node_lat_lon(self):
        grid_cell_nodes_lat_var, grid_cell_nodes_lon_var = self._sgrid.node_coordinates
        grid_cell_nodes_lat = self.nc_dataset.variables[grid_cell_nodes_lat_var][:]
        grid_cell_nodes_lon = self.nc_dataset.variables[grid_cell_nodes_lon_var][:]
        self._sgrid.nodes = pair_arrays(grid_cell_nodes_lon, grid_cell_nodes_lat)
            
            
class SGrid3D(SGridND):
    
    topology_dim = 3
    
    def set_volume_dimensions(self):
        try:
            vol_dim = self.topology_var.volume_dimensions
            vol_dim_padding = self.pp.parse_padding(vol_dim)
            self._sgrid.volume_dimensions = vol_dim
            self._sgrid.volume_padding = vol_dim_padding
        except AttributeError:
            pass
        
    def set_volume_coordinates(self):
        try:
            volume_coordinates = self.topology_var.volume_coordinates
            volume_coordinates_val = volume_coordinates.split(' ')
            self._sgrid.volume_coordinates = tuple(volume_coordinates_val)
        except AttributeError:
            grid_cell_center_vars = self.ncd.find_coordinates_by_location('volume', self.topology_dim)
            self._sgrid.volume_coordinates = grid_cell_center_vars
            
    def set_edge3_dimensions(self):
        try:
            edge3_dim = self.topology_var.edge3_dimensions
            edge3_dim_padding = self.pp.parse_padding(edge3_dim)
            self._sgrid.edge_3_dimension = edge3_dim
            self._sgrid.edge_3_padding = edge3_dim_padding
        except AttributeError:
            pass
        
    def set_all_edge_attributes(self):
        self.set_edge1_dimensions()
        self.set_edge1_coordinates()
        self.set_edge2_dimensions()
        self.set_edge2_coordinates()
        self.set_edge3_dimensions()
        
    def set_all_face_attributes(self):
        pass
        
    def set_cell_center_lat_lon(self):
        grid_cell_center_lon_var = self._sgrid.volume_coordinates[0]
        grid_cell_center_lat_var = self._sgrid.volume_coordinates[1]
        grid_cell_center_lon = self.nc_dataset.variables[grid_cell_center_lon_var][:]
        grid_cell_center_lat = self.nc_dataset.variables[grid_cell_center_lat_var][:]
        self._sgrid.centers = pair_arrays(grid_cell_center_lon, grid_cell_center_lat)
        
    def set_cell_node_lat_lon(self):
        pass


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
        topology_var = grid_topology_vars_attr
        nc_grid_topology_var = nc_dataset.variables[topology_var]
        if nc_grid_topology_var.topology_dimension == 2:
            sg2 = SGrid2D(grid, nc_dataset, topology_var)
            sg2.set_dimensions()
            sg2.set_topology_dimension()
            sg2.set_sgrid_topology()
            sg2.set_sgrid_vertical_dimensions()
            sg2.set_sgrid_node_coordinates()
            sg2.set_all_edge_attributes()
            sg2.set_all_face_attributes()
            sg2.set_cell_center_lat_lon()
            sg2.set_cell_node_lat_lon()
            sg2.set_sgrid_angles()
            sg2.set_sgrid_time()
            sg2.set_sgrid_variable_attributes()
            result_sgrid = sg2.sgrid    
        elif nc_grid_topology_var.topology_dimension == 3:
            sg3 = SGrid3D(grid, nc_dataset, topology_var)
            sg3.set_dimensions()
            sg3.set_topology_dimension()
            sg3.set_sgrid_topology()
            sg3.set_sgrid_vertical_dimensions()
            sg3.set_sgrid_node_coordinates()
            sg3.set_all_edge_attributes()
            sg3.set_all_face_attributes()
            sg3.set_volume_dimensions()
            sg3.set_volume_coordinates()
            sg3.set_cell_center_lat_lon()
            sg3.set_cell_node_lat_lon()
            sg3.set_sgrid_angles()
            sg3.set_sgrid_time()
            sg3.set_sgrid_variable_attributes()
            result_sgrid = sg3.sgrid
        else:
            raise TopologyDimensionError(nc_grid_topology_var.topology_dimension)
        return result_sgrid
    else:
        raise SGridNonCompliantError(nc_dataset)