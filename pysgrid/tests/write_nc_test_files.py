'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import netCDF4 as nc4
import numpy as np
from pysgrid.lookup import (LON_GRID_CELL_CENTER_LONG_NAME, LAT_GRID_CELL_CENTER_LONG_NAME,
                            LON_GRID_CELL_NODE_LONG_NAME, LAT_GRID_CELL_NODE_LONG_NAME)

test_files = os.path.join(os.path.split(__file__)[0], 'files')
file_name = os.path.join(test_files, 'test_sgrid.nc')

if __name__ == '__main__':
    
    with nc4.Dataset(file_name, 'w') as rg:
        # set dimensions
        z_center = rg.createDimension('z_center', 2)
        z_node = rg.createDimension('z_node', 3)
        time = rg.createDimension('time', 2)
        x_center = rg.createDimension('x_center', 4)
        y_center = rg.createDimension('y_center', 4)
        x_node = rg.createDimension('x_node', 3)
        y_node = rg.createDimension('y_node', 3)
        x_u = rg.createDimension('x_u', 3)
        y_u = rg.createDimension('y_u', 4)
        x_v = rg.createDimension('x_v', 4)
        y_v = rg.createDimension('y_v', 3)
        # create coordinate variables
        z_centers = rg.createVariable('z_center', 'i4', ('z_center',))
        z_nodes = rg.createVariable('z_node', 'i4', ('z_node',))
        times = rg.createVariable('time', 'f8', ('time',))
        x_centers = rg.createVariable('x_center', 'f4', ('x_center',))
        y_centers = rg.createVariable('y_center', 'f4', ('y_center',))
        x_nodes = rg.createVariable('x_node', 'f4', ('x_node',))
        y_nodes = rg.createVariable('y_node', 'f4', ('y_node',))
        x_us = rg.createVariable('x_u', 'f4', ('x_u',))
        y_us = rg.createVariable('y_u', 'f4', ('y_u',))
        x_vs = rg.createVariable('x_v', 'f4', ('x_v',))
        y_vs = rg.createVariable('y_v', 'f4', ('y_v',))
        # create other variables
        grid = rg.createVariable('grid', 'i2')
        u = rg.createVariable('u', 'f4', ('time', 'z_center', 'y_u', 'x_u'))
        v = rg.createVariable('v', 'f4', ('time', 'z_center', 'y_v', 'x_v'))
        lat_u = rg.createVariable('lat_u', 'f4', ('y_u', 'x_u'))
        lon_u = rg.createVariable('lon_u', 'f4', ('y_u', 'x_u'))
        lat_v = rg.createVariable('lat_v', 'f4', ('y_v', 'x_v'))
        lon_v = rg.createVariable('lon_v', 'f4', ('y_v', 'x_v'))
        # create variable attributes
        x_centers.long_name = LON_GRID_CELL_CENTER_LONG_NAME[0]
        y_centers.long_name = LAT_GRID_CELL_CENTER_LONG_NAME[0]
        x_nodes.long_name = LON_GRID_CELL_NODE_LONG_NAME[0]
        y_nodes.long_name = LAT_GRID_CELL_NODE_LONG_NAME[0]
        grid.cf_role = 'grid_topology'
        grid.topology_dimension = 2
        grid.node_dimensions = 'x_node y_node'
        grid.face_dimensions = 'x_center: x_node (padding:both) y_center: y_node (padding:both)'
        grid.edge1_dimensions = 'x_u: x_node (padding:both) y_u: y_node'
        grid.edge2_dimensions = 'x_v: x_node y_v: y_node (padding:both)'
        grid.node_coordinates = 'lon_node lat_node'
        grid.face_coordinates = 'lon_center lat_center'
        grid.edge1_coordinates = 'lon_u lat_u'
        grid.edge2_coordinates = 'lon_v lat_v'
        grid.vertical_dimensions = 'z_center: z_node (padding:none)'
        # create coordinate data
        z_centers[:] = np.ndarray(shape=(2,))
        times[:] = np.ndarray(shape=(2,))
        x_centers[:] = np.ndarray(shape=(4,))
        y_centers[:] = np.ndarray(shape=(4,))
        x_nodes[:] = np.ndarray(shape=(3,))
        y_nodes[:] = np.ndarray(shape=(3,))
        x_us[:] = np.ndarray(shape=(3,))
        y_us[:] = np.ndarray(shape=(4,))
        x_vs[:] = np.ndarray(shape=(4,))
        y_vs[:] = np.ndarray(shape=(3,))
        u[:] = np.ndarray(shape=(2, 2, 4, 3))  # x-directed velocities
        v[:] = np.ndarray(shape=(2, 2, 3, 4))  # y-directed velocities 
        lat_u[:] = np.ndarray(shape=(4, 3))
        lon_u[:] = np.ndarray(shape=(4, 3))
        lat_v[:] = np.ndarray(shape=(3, 4))
        lon_v[:] = np.ndarray(shape=(3, 4))
        
        