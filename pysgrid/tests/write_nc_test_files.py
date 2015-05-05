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


def simulated_dgrid(target_dir=TEST_FILES, nc_filename='fake_dgrid.nc'):
    file_name = os.path.join(target_dir, nc_filename)
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
        # grid.face_coordinates = 'XZ YZ'
        grid.vertical_dimensions = 'KMAX: KMAX1 (padding: none)'
        u1.grid = 'some grid'
        u1.standard_name = 'sea_water_x_velocity'
        u1.location = 'edge2'
        v1.grid = 'some grid'
        v1.standard_name = 'sea_water_y_velocity'
        v1.location = 'edge1'
        # create variable data
        xcor[:] = np.random.random((4, 4))
        ycor[:] = np.random.random((4, 4))
        xz[:] = np.random.random((4, 4))
        yz[:] = np.random.random((4, 4))
        u1[:] = np.random.random((2, 2, 4, 4))
        v1[:] = np.random.random((2, 2, 4, 4))
        times[:] = np.random.random((2,))
        

def deltares_sgrid(target_dir=TEST_FILES, nc_filename='test_sgrid_deltares.nc'):
    """
    Create a netCDF file that is structurally similar to
    deltares output. Dimension and variable names may differ
    from an actual file.
    
    """
    file_name = os.path.join(target_dir, nc_filename)
    with nc4.Dataset(file_name, 'w') as rg:
        # define dimensions
        rg.createDimension('MMAXZ', 4)
        rg.createDimension('NMAXZ', 4)
        rg.createDimension('MMAX', 4)
        rg.createDimension('NMAX', 4)
        rg.createDimension('KMAX', 2)
        rg.createDimension('KMAX1', 3)
        rg.createDimension('time', 2)
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
        latitude.axes = 'X: NMAXZ Y: MMAXZ'
        longitude.long_name = LON_GRID_CELL_CENTER_LONG_NAME[1]
        longitude.axes = 'X: NMAXZ Y: MMAXZ'
        grid_latitude.long_name = LAT_GRID_CELL_NODE_LONG_NAME[1]
        grid_latitude.axes = 'X: NMAX Y: MMAX'
        grid_longitude.long_name = LON_GRID_CELL_NODE_LONG_NAME[1]
        grid_longitude.axes = 'X: NMAX Y: MMAX'
        u1.grid = 'some grid'
        u1.axes = 'X: NMAXZ Y: MMAX Z: KMAX'
        u1.standard_name = 'sea_water_x_velocity'
        v1.grid = 'some grid'
        v1.axes = 'X: NMAX Y: MMAXZ Z: KMAX'
        v1.standard_name = 'sea_water_y_velocity'
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
    return file_name
        
        
def roms_sgrid(target_dir=TEST_FILES, nc_filename='test_sgrid_roms.nc'):
    """
    Create a netCDF file that is structurally similar to
    ROMS output. Dimension and variable names may differ
    from an actual file.
    
    """
    file_name = os.path.join(target_dir, nc_filename)
    with nc4.Dataset(file_name, 'w') as rg:
        # set dimensions
        rg.createDimension('s_rho', 2)
        rg.createDimension('s_w', 3)
        rg.createDimension('time', 2)
        rg.createDimension('xi_rho', 4)
        rg.createDimension('eta_rho', 4)
        rg.createDimension('xi_psi', 3)
        rg.createDimension('eta_psi', 3)
        rg.createDimension('xi_u', 3)
        rg.createDimension('eta_u', 4)
        rg.createDimension('xi_v', 4)
        rg.createDimension('eta_v', 3)
        # create coordinate variables
        z_centers = rg.createVariable('s_rho', 'i4', ('s_rho',))
        rg.createVariable('s_w', 'i4', ('s_w',))
        times = rg.createVariable('time', 'f8', ('time',))
        rg.createVariable('xi_rho', 'f4', ('xi_rho',))
        rg.createVariable('eta_rho', 'f4', ('eta_rho',))
        rg.createVariable('xi_psi', 'f4', ('xi_psi',))
        rg.createVariable('eta_psi', 'f4', ('eta_psi',))
        x_us = rg.createVariable('xi_u', 'f4', ('xi_u',))
        y_us = rg.createVariable('eta_u', 'f4', ('eta_u',))
        x_vs = rg.createVariable('xi_v', 'f4', ('xi_v',))
        y_vs = rg.createVariable('eta_v', 'f4', ('eta_v',))
        # create other variables
        grid = rg.createVariable('grid', 'i2')
        u = rg.createVariable('u', 'f4', ('time', 's_rho', 'eta_u', 'xi_u'))
        v = rg.createVariable('v', 'f4', ('time', 's_rho', 'eta_v', 'xi_v'))
        lon_centers = rg.createVariable('lon_rho', 'f4', ('eta_rho', 'xi_rho'))
        lat_centers = rg.createVariable('lat_rho', 'f4', ('eta_rho', 'xi_rho'))
        lon_nodes = rg.createVariable('lon_psi', 'f4', ('eta_psi', 'xi_psi'))
        lat_nodes = rg.createVariable('lat_psi', 'f4', ('eta_psi', 'xi_psi'))
        lat_u = rg.createVariable('lat_u', 'f4', ('eta_u', 'xi_u'))
        lon_u = rg.createVariable('lon_u', 'f4', ('eta_u', 'xi_u'))
        lat_v = rg.createVariable('lat_v', 'f4', ('eta_v', 'xi_v'))
        lon_v = rg.createVariable('lon_v', 'f4', ('eta_v', 'xi_v'))
        zeta = rg.createVariable('zeta', 'f4', ('time', 'eta_rho', 'xi_rho'))
        # create variable attributes
        lon_centers.long_name = LON_GRID_CELL_CENTER_LONG_NAME[0]
        lon_centers.standard_name = 'longitude'
        lon_centers.axes = 'X: xi_rho Y: eta_rho'
        lat_centers.long_name = LAT_GRID_CELL_CENTER_LONG_NAME[0]
        lat_centers.standard_name = 'latitude'
        lat_centers.axes = 'X: xi_rho Y: eta_rho'
        lon_nodes.long_name = LON_GRID_CELL_NODE_LONG_NAME[0]
        lon_nodes.axes = 'X: xi_psi Y: eta_psi'
        lat_nodes.long_name = LAT_GRID_CELL_NODE_LONG_NAME[0]
        lat_nodes.axes = 'X: xi_psi Y: eta_psi'
        grid.cf_role = 'grid_topology'
        grid.topology_dimension = 2
        grid.node_dimensions = 'xi_psi eta_psi'
        grid.face_dimensions = 'xi_rho: xi_psi (padding: both) eta_rho: eta_psi (padding: both)'
        grid.edge1_dimensions = 'xi_u: xi_psi eta_u: eta_psi (padding: both)'
        grid.edge2_dimensions = 'xi_v: xi_psi (padding: both) eta_v: eta_psi'
        grid.node_coordinates = 'lon_psi lat_psi'
        grid.face_coordinates = 'lon_rho lat_rho'
        grid.edge1_coordinates = 'lon_u lat_u'
        grid.edge2_coordinates = 'lon_v lat_v'
        grid.vertical_dimensions = 's_rho: s_w (padding: none)'
        zeta.location = 'faces'
        zeta.coordinates = 'time lat_rho lon_rho'
        u.grid = 'some grid'
        u.axes = 'X: xi_u Y: eta_u'
        u.location = 'edge1'
        u.standard_name = 'sea_water_x_velocity'
        v.grid = 'some grid'
        v.axes = 'X: xi_v Y: eta_v'
        v.location = 'edge2'
        v.standard_name = 'sea_water_y_velocity'
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
    return file_name


def wrf_sgrid_2d(target_dir=TEST_FILES, nc_filename='test_sgrid_wrf_2.nc'):
    file_name = os.path.join(target_dir, nc_filename)
    with nc4.Dataset(file_name, 'w') as nc:
        nc.createDimension('Time', 2)
        nc.createDimension('DateStrLen', 3)
        nc.createDimension('west_east', 4)
        nc.createDimension('south_north', 5)
        nc.createDimension('west_east_stag', 5)
        nc.createDimension('bottom_top', 3)
        nc.createDimension('south_north_stag', 6)
        nc.createDimension('bottom_top_stag', 4)
        times = nc.createVariable('Times', np.dtype(str), ('Time', 'DateStrLen'))
        us = nc.createVariable('U', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east_stag'))
        us.grid = 'grid'
        us.location = 'edge1'
        vs = nc.createVariable('V', 'f4', ('Time', 'bottom_top', 'south_north_stag', 'west_east'))
        vs.grid = 'grid'
        vs.location = 'edge2'
        ws = nc.createVariable('W', 'f4', ('Time', 'bottom_top_stag', 'south_north', 'west_east'))
        ws.grid = 'grid'
        ws.location = 'face'
        temps = nc.createVariable('T', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))
        temps.grid = 'grid'
        temps.location = 'face'
        xlats = nc.createVariable('XLAT', 'f4', ('Time', 'south_north', 'west_east'))
        xlongs = nc.createVariable('XLONG', 'f4', ('Time', 'south_north', 'west_east'))
        znus = nc.createVariable('ZNU', 'f4', ('Time', 'bottom_top'))
        znws = nc.createVariable('ZNW', 'f4', ('Time', 'bottom_top_stag'))
        grid = nc.createVariable('grid', 'i2')
        grid.cf_role = 'grid_topology'
        grid.topology_dimension = 2
        grid.node_dimensions = 'west_east_stag south_north_stag'
        grid.face_dimensions = ('west_east: west_east_stag (padding: none) '
                                'south_north: south_north_stag (padding: none)'
                                )
        grid.face_coordinates = 'XLONG XLAT'
        grid.vertical_dimensions = 'bottom_top: bottom_top_stag (padding: none)'
        grid.edge1_dimensions = 'west_east_stag south_north: south_north_stag (padding: none)'
        grid.edge2_dimensions = 'west_east: west_east_stag (padding: none) south_north_stag'
        times[:] = np.random.random(size=(2, 3)).astype(str)
        us[:, :, :, :] = np.random.random(size=(2, 3, 5, 5))
        vs[:, :, :, :] = np.random.random(size=(2, 3, 6, 4))
        ws[:, :, :, :] = np.random.random(size=(2, 4, 5, 4))
        temps[:, :, :, :] = np.random.random(size=(2, 3, 5, 4))
        xlats[:, :, :] = np.random.random(size=(2, 5, 4))
        xlongs[:, :, :] = np.random.random(size=(2, 5, 4))
        znus[:, :] = np.random.random(size=(2, 3))
        znws[:, :] = np.random.random(size=(2, 4))
    return file_name
        
        
def wrf_sgrid(target_dir=TEST_FILES, nc_filename='test_sgrid_wrf.nc'):
    """
    Write an SGrid file using 3D conventions.
    
    """
    file_name = os.path.join(target_dir, nc_filename)
    with  nc4.Dataset(file_name, 'w') as fg:
        # create dimensions
        fg.createDimension('Time', 2)
        fg.createDimension('DateStrLen', 3)
        fg.createDimension('west_east', 4)
        fg.createDimension('south_north', 5)
        fg.createDimension('west_east_stag', 5)
        fg.createDimension('bottom_top', 3)
        fg.createDimension('south_north_stag', 6)
        fg.createDimension('bottom_top_stag', 4)
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
        grid.volume_dimensions = ('west_east: west_east_stag (padding: none) '
                                  'south_north: south_north_stag (padding: none) '
                                  'bottom_top: bottom_top_stag (padding: none)')
        grid.volume_coordinates = 'XLONG XLAT ZNU'
        # create fake data
        times[:] = np.random.random(size=(2, 3)).astype(str)
        us[:, :, :, :] = np.random.random(size=(2, 3, 5, 5))
        vs[:, :, :, :] = np.random.random(size=(2, 3, 6, 4))
        ws[:, :, :, :] = np.random.random(size=(2, 4, 5, 4))
        temps[:, :, :, :] = np.random.random(size=(2, 3, 5, 4))
        xlats[:, :, :] = np.random.random(size=(2, 5, 4))
        xlongs[:, :, :] = np.random.random(size=(2, 5, 4))
        znus[:, :] = np.random.random(size=(2, 3))
        znws[:, :] = np.random.random(size=(2, 4))
    return file_name
        
               
def non_compliant_sgrid(target_dir=TEST_FILES, nc_filename='test_noncompliant_sgrid.nc'):
    """
    Create a netCDF file that is structurally similar to
    ROMS output. Dimension and variable names may differ
    from an actual file.
    
    """
    file_name = os.path.join(target_dir, nc_filename)
    with nc4.Dataset(file_name, 'w') as rg:
        # set dimensions
        rg.createDimension('z_center', 2)
        rg.createDimension('z_node', 3)
        rg.createDimension('time', 2)
        rg.createDimension('x_center', 4)
        rg.createDimension('y_center', 4)
        rg.createDimension('x_node', 3)
        rg.createDimension('y_node', 3)
        rg.createDimension('x_u', 3)
        rg.createDimension('y_u', 4)
        rg.createDimension('x_v', 4)
        rg.createDimension('y_v', 3)
        # create coordinate variables
        z_centers = rg.createVariable('z_center', 'i4', ('z_center',))
        rg.createVariable('z_node', 'i4', ('z_node',))
        times = rg.createVariable('time', 'f8', ('time',))
        rg.createVariable('x_center', 'f4', ('x_center',))
        rg.createVariable('y_center', 'f4', ('y_center',))
        rg.createVariable('x_node', 'f4', ('x_node',))
        rg.createVariable('y_node', 'f4', ('y_node',))
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
    return file_name