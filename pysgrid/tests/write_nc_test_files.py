'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import netCDF4 as nc4
import numpy as np
from pysgrid.lookup import (LON_GRID_CELL_CENTER_LONG_NAME, LAT_GRID_CELL_CENTER_LONG_NAME,
                            LON_GRID_CELL_NODE_LONG_NAME, LAT_GRID_CELL_NODE_LONG_NAME)


TEST_FILES = os.path.join(os.path.split(__file__)[0], 'files')


def fake_dgrid(nc_filename='fake_dgrid.nc'):
    file_name = os.path.join(TEST_FILES, nc_filename)
    with nc4.Dataset(file_name, 'w') as rg:
        # define dims
        rg.createDimension('MMAXZ', 4)
        rg.createDimension('NMAXZ', 4)
        rg.createDimension('MMAX', 4)
        rg.createDimension('NMAX', 4)
        rg.createDimension('KMAX', 2)
        rg.createDimension('KMAX1', 3)
        rg.createDimension('time', 2)     
        # define vars
        xcor = rg.createVariable('XCOR', 'f4', ('MMAX', 'NMAX'))
        ycor = rg.createVariable('YCOR', 'f4', ('MMAX', 'NMAX'))
        xz = rg.createVariable('XZ', 'f4', ('MMAXZ', 'NMAXZ'))
        yz = rg.createVariable('YZ', 'f4', ('MMAXZ', 'NMAXZ'))
        u1 = rg.createVariable('U1', 'f4', ('time', 'KMAX', 'MMAXZ', 'NMAX'))  # dims T, Z, X, Y
        v1 = rg.createVariable('V1', 'f4', ('time', 'KMAX', 'MMAX', 'NMAXZ'))
        times = rg.createVariable('time', 'f8', ('time',))
        grid = rg.createVariable('grid', 'i4')
        # define attributes
        grid.cf_role = 'grid_topology'
        grid.topology_dimension = 2
        grid.node_dimensions = 'MMAX NMAX'
        grid.face_dimensions = 'MMAXZ: MMAX (padding: low) NMAXZ: NMAX (padding: low)'
        grid.node_coordinates = 'XCOR YCOR'
        grid.face_coordinates = 'XZ YZ'
        grid.vertical_dimensions = 'KMAX: KMAX1 (padding: none)'
        u1.grid = 'some grid'
        v1.grid = 'some grid'
        # create variable data
        xcor[:] = np.random.random((4, 4))
        ycor[:] = np.random.random((4, 4))
        xz[:] = np.random.random((4, 4))
        yz[:] = np.random.random((4, 4))
        u1[:] = np.random.random((2, 2, 4, 4))
        v1[:] = np.random.random((2, 2, 4, 4))
        times[:] = np.random.random((2,))
        

def deltares_like_sgrid(nc_filename='test_sgrid_deltares_like.nc'):
    """
    Create a netCDF file that is structurally similar to
    deltares output. Dimension and variable names may differ
    from an actual file.
    
    """
    file_name = os.path.join(TEST_FILES, nc_filename)
    with nc4.Dataset(file_name, 'w') as rg:
        # define dimensions
        y_center = rg.createDimension('MMAXZ', 4)
        x_center = rg.createDimension('NMAXZ', 4)
        y_node = rg.createDimension('MMAX', 4)
        x_node = rg.createDimension('NMAX', 4)
        z_center = rg.createDimension('KMAX', 2)
        z_node = rg.createDimension('KMAX1', 3)
        time = rg.createDimension('time', 2)
        # define variables
        xcor = rg.createVariable('XCOR', 'f4', ('MMAX', 'NMAX'))  # nodes
        ycor = rg.createVariable('YCOR', 'f4', ('MMAX', 'NMAX'))  # nodes
        xz = rg.createVariable('XZ', 'f4', ('MMAXZ', 'NMAXZ'))  # centers
        yz = rg.createVariable('YZ', 'f4', ('MMAXZ', 'NMAXZ'))  # centers
        u1 = rg.createVariable('U1', 'f4', ('time', 'KMAX', 'MMAX', 'NMAXZ'))
        v1 = rg.createVariable('V1', 'f4', ('time', 'KMAX', 'MMAXZ', 'NMAX'))
        times = rg.createVariable('time', 'f8', ('time',))
        grid = rg.createVariable('grid', 'i4')
        latitude = rg.createVariable('latitude', 'f4', ('MMAXZ', 'NMAXZ'))
        longitude = rg.createVariable('longitude', 'f4', ('MMAXZ', 'NMAXZ'))
        grid_latitude = rg.createVariable('grid_latitude', 'f4', ('MMAX', 'NMAX'))
        grid_longitude = rg.createVariable('grid_longitude', 'f4', ('MMAX', 'NMAX')) 
        # define variable attributes
        grid.cf_role = 'grid_topology'
        grid.topology_dimension = 2
        grid.node_dimensions = 'MMAX NMAX'
        grid.face_dimensions = 'MMAXZ: MMAX (padding: low) NMAXZ: NMAX (padding: low)'
        grid.node_coordinates = 'XCOR YCOR'
        grid.face_coordinates = 'XZ YZ'
        grid.vertical_dimensions = 'KMAX: KMAX1 (padding: none)'
        latitude.long_name = LAT_GRID_CELL_CENTER_LONG_NAME[1]
        longitude.long_name = LON_GRID_CELL_CENTER_LONG_NAME[1]
        grid_latitude.long_name = LAT_GRID_CELL_NODE_LONG_NAME[1]
        grid_longitude.long_name = LON_GRID_CELL_NODE_LONG_NAME[1]
        u1.grid = 'some grid'
        v1.grid = 'some grid'
        # create variable data
        xcor[:] = np.random.random((4, 4))
        ycor[:] = np.random.random((4, 4))
        xz[:] = np.random.random((4, 4))
        yz[:] = np.random.random((4, 4))
        u1[:] = np.random.random((2, 2, 4, 4))
        v1[:] = np.random.random((2, 2, 4, 4))
        times[:] = np.random.random((2,))
        latitude[:] = np.random.random((4, 4))
        longitude[:] = np.random.random((4, 4))
        grid_latitude[:] = np.random.random((4, 4))
        grid_longitude[:] = np.random.random((4, 4))
        
        
def roms_like_sgrid(nc_filename='test_sgrid_roms_like.nc'):
    """
    Create a netCDF file that is structurally similar to
    ROMS output. Dimension and variable names may differ
    from an actual file.
    
    """
    file_name = os.path.join(TEST_FILES, nc_filename)
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
        lon_centers = rg.createVariable('lon_center', 'f4', ('y_center', 'x_center'))
        lat_centers = rg.createVariable('lat_center', 'f4', ('y_center', 'x_center'))
        lon_nodes = rg.createVariable('lon_node', 'f4', ('y_node', 'x_node'))
        lat_nodes = rg.createVariable('lat_node', 'f4', ('y_node', 'x_node'))
        lat_u = rg.createVariable('lat_u', 'f4', ('y_u', 'x_u'))
        lon_u = rg.createVariable('lon_u', 'f4', ('y_u', 'x_u'))
        lat_v = rg.createVariable('lat_v', 'f4', ('y_v', 'x_v'))
        lon_v = rg.createVariable('lon_v', 'f4', ('y_v', 'x_v'))
        zeta = rg.createVariable('zeta', 'f4', ('time', 'y_center', 'x_center'))
        # create variable attributes
        lon_centers.long_name = LON_GRID_CELL_CENTER_LONG_NAME[0]
        lon_centers.standard_name = 'longitude'
        lat_centers.long_name = LAT_GRID_CELL_CENTER_LONG_NAME[0]
        lat_centers.standard_name = 'latitude'
        lon_nodes.long_name = LON_GRID_CELL_NODE_LONG_NAME[0]
        lat_nodes.long_name = LAT_GRID_CELL_NODE_LONG_NAME[0]
        grid.cf_role = 'grid_topology'
        grid.topology_dimension = 2
        grid.node_dimensions = 'x_node y_node'
        grid.face_dimensions = 'x_center: x_node (padding: both) y_center: y_node (padding: both)'
        grid.edge1_dimensions = 'x_u: x_node y_u: y_node (padding: both)'
        grid.edge2_dimensions = 'x_v: x_node (padding: both) y_v: y_node'
        grid.node_coordinates = 'lon_node lat_node'
        grid.face_coordinates = 'lon_center lat_center'
        grid.edge1_coordinates = 'lon_u lat_u'
        grid.edge2_coordinates = 'lon_v lat_v'
        grid.vertical_dimensions = 'z_center: z_node (padding: none)'
        zeta.location = 'faces'
        zeta.coordinates = 'time lat_center lon_center'
        u.grid = 'some grid'
        v.grid = 'some grid'
        # create coordinate data
        z_centers[:] = np.random.random(size=(2,))
        times[:] = np.random.random(size=(2,))
        lon_centers[:, :] = np.random.random(size=(4, 4))
        lat_centers[:, :] = np.random.random(size=(4, 4))
        lon_nodes[:] = np.random.random(size=(3, 3))
        lat_nodes[:] = np.random.random(size=(3, 3))
        x_us[:] = np.random.random(size=(3,))
        y_us[:] = np.random.random(size=(4,))
        x_vs[:] = np.random.random(size=(4,))
        y_vs[:] = np.random.random(size=(3,))
        u[:, :, :, :] = np.random.random(size=(2, 2, 4, 3))  # x-directed velocities
        v[:] = np.random.random(size=(2, 2, 3, 4))  # y-directed velocities 
        lat_u[:] = np.random.random(size=(4, 3))
        lon_u[:] = np.random.random(size=(4, 3))
        lat_v[:] = np.random.random(size=(3, 4))
        lon_v[:] = np.random.random(size=(3, 4))
        
        
def wrf_like_sgrid(nc_filename='test_sgrid_wrf_like.nc'):
    file_name = os.path.join(TEST_FILES, nc_filename)
    with  nc4.Dataset(file_name, 'w') as fg:
        # create dimensions
        fg.createDimension('Time', 13)
        fg.createDimension('DateStrLen', 19)
        fg.createDimension('west_east', 73)
        fg.createDimension('south_north', 60)
        fg.createDimension('west_east_stag', 74)
        fg.createDimension('bottom_top', 27)
        fg.createDimension('south_north_stag', 61)
        fg.createDimension('bottom_top_stag', 28)
        # create variables
        times = fg.createVariable('Times', np.dtype(str), ('Time', 'DateStrLen'))
        us = fg.createVariable('U', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east_stag'))
        us.grid = 'grid'
        us.location = 'face1'
        vs = fg.createVariable('V', 'f4', ('Time', 'bottom_top', 'south_north_stag', 'west_east'))
        vs.grid = 'grid'
        vs.location = 'face2'
        ws = fg.createVariable('W', 'f4', ('Time', 'bottom_top_stag', 'south_north', 'west_east'))
        ws.grid = 'grid'
        ws.location = 'face3'
        temps = fg.createVariable('T', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))
        temps.grid = 'grid'
        temps.location = 'volume'
        xlats = fg.createVariable('XLAT', 'f4', ('Time', 'south_north', 'west_east'))
        xlongs = fg.createVariable('XLONG', 'f4', ('Time', 'south_north', 'west_east'))
        znus = fg.createVariable('ZNU', 'f4', ('Time', 'bottom_top'))
        znws = fg.createVariable('ZNW', 'f4', ('Time', 'bottom_top_stag'))
        grid = fg.createVariable('grid', 'i2')
        grid.cf_role = 'grid_topology'
        grid.topology_dimension = 3
        grid.node_dimensions = 'west_east_stag south_north_stag bottom_top_stag'
        grid.volume_dimensions = 'west_east: west_east_stag (padding: none) south_north: south_north_stag (padding: none) bottom_top: bottom_top_stag (padding: none)'
        grid.volume_coordinates = 'XLONG XLAT ZNU'
        # create fake data
        times[:] = np.random.random(size=(13, 19)).astype(str)
        us[:, :, :, :] = np.random.random(size=(13, 27, 60, 74))
        vs[:, :, :, :] = np.random.random(size=(13, 27, 61, 73))
        ws[:, :, :, :] = np.random.random(size=(13, 28, 60, 73))
        temps[:, :, :, :] = np.random.random(size=(13, 27, 60, 73))
        xlats[:, :, :] = np.random.random(size=(13, 60, 73))
        xlongs[:, :, :] = np.random.random(size=(13, 60, 73))
        znus[:, :] = np.random.random(size=(13, 27))
        znws[:, :] = np.random.random(size=(13, 28))
        
        
        
def roms_like_non_compliant_sgrid(nc_filename='test_noncompliant_sgrid_roms_like.nc'):
    """
    Create a netCDF file that is structurally similar to
    ROMS output. Dimension and variable names may differ
    from an actual file.
    
    """
    file_name = os.path.join(TEST_FILES, nc_filename)
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
        lon_centers = rg.createVariable('lon_center', 'f4', ('y_center', 'x_center'))
        lat_centers = rg.createVariable('lat_center', 'f4', ('y_center', 'x_center'))
        lon_nodes = rg.createVariable('lon_node', 'f4', ('y_node', 'x_node'))
        lat_nodes = rg.createVariable('lat_node', 'f4', ('y_node', 'x_node'))
        lat_u = rg.createVariable('lat_u', 'f4', ('y_u', 'x_u'))
        lon_u = rg.createVariable('lon_u', 'f4', ('y_u', 'x_u'))
        lat_v = rg.createVariable('lat_v', 'f4', ('y_v', 'x_v'))
        lon_v = rg.createVariable('lon_v', 'f4', ('y_v', 'x_v'))
        # create variable attributes
        lon_centers.long_name = LON_GRID_CELL_CENTER_LONG_NAME[0]
        lat_centers.long_name = LAT_GRID_CELL_CENTER_LONG_NAME[0]
        lon_nodes.long_name = LON_GRID_CELL_NODE_LONG_NAME[0]
        lat_nodes.long_name = LAT_GRID_CELL_NODE_LONG_NAME[0]
        grid.topology_dimension = 2
        grid.node_dimensions = 'x_node y_node'
        grid.face_dimensions = 'x_center: x_node (padding: both) y_center: y_node (padding: both)'
        grid.edge1_dimensions = 'x_u: x_node y_u: y_node (padding: both)'
        grid.edge2_dimensions = 'x_v: x_node (padding: both) y_v: y_node'
        grid.node_coordinates = 'lon_node lat_node'
        grid.face_coordinates = 'lon_center lat_center'
        grid.edge1_coordinates = 'lon_u lat_u'
        grid.edge2_coordinates = 'lon_v lat_v'
        grid.vertical_dimensions = 'z_center: z_node (padding: none)'
        # create coordinate data
        z_centers[:] = np.random.random(size=(2,))
        times[:] = np.random.random(size=(2,))
        lon_centers[:, :] = np.random.random(size=(4, 4))
        lat_centers[:, :] = np.random.random(size=(4, 4))
        lon_nodes[:] = np.random.random(size=(3, 3))
        lat_nodes[:] = np.random.random(size=(3, 3))
        x_us[:] = np.random.random(size=(3,))
        y_us[:] = np.random.random(size=(4,))
        x_vs[:] = np.random.random(size=(4,))
        y_vs[:] = np.random.random(size=(3,))
        u[:, :, :, :] = np.random.random(size=(2, 2, 4, 3))  # x-directed velocities
        v[:] = np.random.random(size=(2, 2, 3, 4))  # y-directed velocities 
        lat_u[:] = np.random.random(size=(4, 3))
        lon_u[:] = np.random.random(size=(4, 3))
        lat_v[:] = np.random.random(size=(3, 4))
        lon_v[:] = np.random.random(size=(3, 4))

if __name__ == '__main__':
    
    # deltares_like_sgrid()
    # roms_like_sgrid()
    # roms_like_non_compliant_sgrid()
    # fake_dgrid()
    wrf_like_sgrid()
        