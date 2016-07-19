'''
Created on Apr 7, 2015

@author: ayan
'''

from __future__ import (absolute_import, division, print_function)

import os
import unittest

from netCDF4 import Dataset

from ..read_netcdf import (NetCDFDataset, parse_axes, parse_padding,
                           parse_vector_axis, find_grid_topology_var)
from .write_nc_test_files import roms_sgrid, wrf_sgrid_2d


class TestParseAxes(unittest.TestCase):

    def setUp(self):
        self.xy = 'X: xi_psi Y: eta_psi'
        self.xyz = 'X: NMAX Y: MMAXZ Z: KMAX'

    def test_xyz_axis_parse(self):
        result = parse_axes(self.xyz)
        expected = ('NMAX', 'MMAXZ', 'KMAX')
        self.assertEqual(result, expected)

    def test_xy_axis_parse(self):
        result = parse_axes(self.xy)
        expected = ('xi_psi', 'eta_psi', None)
        self.assertEqual(result, expected)


class TestParsePadding(unittest.TestCase):

    def setUp(self):
        self.grid_topology = 'some_grid'
        self.with_two_padding = 'xi_rho: xi_psi (padding: both) eta_rho: eta_psi (padding: low)'
        self.with_one_padding = 'xi_v: xi_psi (padding: high) eta_v: eta_psi'
        self.with_no_padding = 'MMAXZ: MMAX NMAXZ: NMAX'

    def test_mesh_name(self):
        result = parse_padding(self.with_one_padding, self.grid_topology)
        mesh_topology = result[0].mesh_topology_var
        expected = 'some_grid'
        self.assertEqual(mesh_topology, expected)

    def test_two_padding_types(self):
        result = parse_padding(self.with_two_padding, self.grid_topology)
        expected_len = 2
        padding_datum_0 = result[0]
        padding_type = padding_datum_0.padding
        sub_dim = padding_datum_0.node_dim
        dim = padding_datum_0.face_dim
        expected_sub_dim = 'xi_psi'
        expected_padding_type = 'both'
        expected_dim = 'xi_rho'
        self.assertEqual(len(result), expected_len)
        self.assertEqual(padding_type, expected_padding_type)
        self.assertEqual(sub_dim, expected_sub_dim)
        self.assertEqual(dim, expected_dim)

    def test_one_padding_type(self):
        result = parse_padding(self.with_one_padding, self.grid_topology)
        expected_len = 1
        padding_datum_0 = result[0]
        padding_type = padding_datum_0.padding
        sub_dim = padding_datum_0.node_dim
        dim = padding_datum_0.face_dim
        expected_padding_type = 'high'
        expected_sub_dim = 'xi_psi'
        expected_dim = 'xi_v'
        self.assertEqual(len(result), expected_len)
        self.assertEqual(padding_type, expected_padding_type)
        self.assertEqual(sub_dim, expected_sub_dim)
        self.assertEqual(dim, expected_dim)

    def test_no_padding(self):
        self.assertRaises(ValueError,
                          parse_padding,
                          padding_str=self.with_no_padding,
                          mesh_topology_var=self.grid_topology
                          )


class TestParseVectorAxis(unittest.TestCase):

    def setUp(self):
        self.standard_name_1 = 'sea_water_y_velocity'
        self.standard_name_2 = 'atmosphere_optical_thickness_due_to_cloud'
        self.standard_name_3 = 'ocean_heat_x_transport_due_to_diffusion'

    def test_std_name_with_velocity_direction(self):
        direction = parse_vector_axis(self.standard_name_1)
        expected_direction = 'Y'
        self.assertEqual(direction, expected_direction)

    def test_std_name_without_direction(self):
        direction = parse_vector_axis(self.standard_name_2)
        self.assertIsNone(direction)

    def test_std_name_with_transport_direction(self):
        direction = parse_vector_axis(self.standard_name_3)
        expected_direction = 'X'
        self.assertEqual(direction, expected_direction)


import pytest
import contextlib
import numpy as np
from pysgrid.lookup import (LON_GRID_CELL_CENTER_LONG_NAME,
                            LAT_GRID_CELL_CENTER_LONG_NAME,
                            LON_GRID_CELL_NODE_LONG_NAME,
                            LAT_GRID_CELL_NODE_LONG_NAME)

@pytest.yield_fixture
def roms(fname='tmp_sgrid_roms.nc'):
    """
    Create a netCDF file that is structurally similar to
    ROMS output. Dimension and variable names may differ
    from an actual file.

    """
    nc = Dataset(fname, 'w')
    # Set dimensions.
    nc.createDimension('s_rho', 2)
    nc.createDimension('s_w', 3)
    nc.createDimension('time', 2)
    nc.createDimension('xi_rho', 4)
    nc.createDimension('eta_rho', 4)
    nc.createDimension('xi_psi', 3)
    nc.createDimension('eta_psi', 3)
    nc.createDimension('xi_u', 3)
    nc.createDimension('eta_u', 4)
    nc.createDimension('xi_v', 4)
    nc.createDimension('eta_v', 3)
    # Create coordinate variables.
    z_centers = nc.createVariable('s_rho', 'i4', ('s_rho',))
    nc.createVariable('s_w', 'i4', ('s_w',))
    times = nc.createVariable('time', 'f8', ('time',))
    nc.createVariable('xi_rho', 'f4', ('xi_rho',))
    nc.createVariable('eta_rho', 'f4', ('eta_rho',))
    nc.createVariable('xi_psi', 'f4', ('xi_psi',))
    nc.createVariable('eta_psi', 'f4', ('eta_psi',))
    x_us = nc.createVariable('xi_u', 'f4', ('xi_u',))
    y_us = nc.createVariable('eta_u', 'f4', ('eta_u',))
    x_vs = nc.createVariable('xi_v', 'f4', ('xi_v',))
    y_vs = nc.createVariable('eta_v', 'f4', ('eta_v',))
    # Create other variables.
    grid = nc.createVariable('grid', 'i2')
    u = nc.createVariable('u', 'f4', ('time', 's_rho', 'eta_u', 'xi_u'))
    v = nc.createVariable('v', 'f4', ('time', 's_rho', 'eta_v', 'xi_v'))
    fake_u = nc.createVariable('fake_u', 'f4', ('time', 's_rho', 'eta_u', 'xi_u'))  # noqa
    lon_centers = nc.createVariable('lon_rho', 'f4', ('eta_rho', 'xi_rho'))
    lat_centers = nc.createVariable('lat_rho', 'f4', ('eta_rho', 'xi_rho'))
    lon_nodes = nc.createVariable('lon_psi', 'f4', ('eta_psi', 'xi_psi'))
    lat_nodes = nc.createVariable('lat_psi', 'f4', ('eta_psi', 'xi_psi'))
    lat_u = nc.createVariable('lat_u', 'f4', ('eta_u', 'xi_u'))
    lon_u = nc.createVariable('lon_u', 'f4', ('eta_u', 'xi_u'))
    lat_v = nc.createVariable('lat_v', 'f4', ('eta_v', 'xi_v'))
    lon_v = nc.createVariable('lon_v', 'f4', ('eta_v', 'xi_v'))
    salt = nc.createVariable('salt', 'f4', ('time', 's_rho', 'eta_rho', 'xi_rho'))  # noqa
    zeta = nc.createVariable('zeta', 'f4', ('time', 'eta_rho', 'xi_rho'))
    # Create variable attributes.
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
    times.standard_name = 'time'
    grid.cf_role = 'grid_topology'
    grid.topology_dimension = 2
    grid.node_dimensions = 'xi_psi eta_psi'
    grid.face_dimensions = 'xi_rho: xi_psi (padding: both) eta_rho: eta_psi (padding: both)'  # noqa
    grid.edge1_dimensions = 'xi_u: xi_psi eta_u: eta_psi (padding: both)'
    grid.edge2_dimensions = 'xi_v: xi_psi (padding: both) eta_v: eta_psi'
    grid.node_coordinates = 'lon_psi lat_psi'
    grid.face_coordinates = 'lon_rho lat_rho'
    grid.edge1_coordinates = 'lon_u lat_u'
    grid.edge2_coordinates = 'lon_v lat_v'
    grid.vertical_dimensions = 's_rho: s_w (padding: none)'
    salt.grid = 'grid'
    zeta.location = 'face'
    zeta.coordinates = 'time lat_rho lon_rho'
    u.grid = 'some grid'
    u.axes = 'X: xi_u Y: eta_u'
    u.coordinates = 'time s_rho lat_u lon_u '
    u.location = 'edge1'
    u.standard_name = 'sea_water_x_velocity'
    v.grid = 'some grid'
    v.axes = 'X: xi_v Y: eta_v'
    v.location = 'edge2'
    v.standard_name = 'sea_water_y_velocity'
    fake_u.grid = 'some grid'
    # Create coordinate data.
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
    u[:] = np.random.random(size=(2, 2, 4, 3))  # x-directed velocities
    v[:] = np.random.random(size=(2, 2, 3, 4))  # y-directed velocities
    fake_u[:] = np.random.random(size=(2, 2, 4, 3))
    lat_u[:] = np.random.random(size=(4, 3))
    lon_u[:] = np.random.random(size=(4, 3))
    lat_v[:] = np.random.random(size=(3, 4))
    lon_v[:] = np.random.random(size=(3, 4))
    salt[:] = np.random.random(size=(2, 2, 4, 4))
    nc.sync()
    yield nc
    nc.close()
    os.remove(fname)


def test_finding_node_variables(roms):
    nc_ds = NetCDFDataset(roms)
    result = nc_ds.find_node_coordinates('xi_psi eta_psi')
    expected = ('lon_psi', 'lat_psi')
    assert result == expected

def test_find_face_coordinates_by_location(roms):
    nc_ds = NetCDFDataset(roms)
    result = nc_ds.find_coordinates_by_location('face', 2)
    expected = ('lon_rho', 'lat_rho')
    assert result == expected

def test_find_edge_coordinates_by_location(roms):
    nc_ds = NetCDFDataset(roms)
    result = nc_ds.find_coordinates_by_location('edge1', 2)
    expected = ('lon_u', 'lat_u')
    assert result == expected

def test_find_grid_topology(roms):
    result = find_grid_topology_var(roms)
    expected = 'grid'
    assert result == expected

def test_find_variables_by_standard_name(roms):
    nc_ds = NetCDFDataset(roms)
    result = nc_ds.find_variables_by_attr(standard_name='time')
    expected = ['time']
    assert result == expected

def test_find_variables_by_standard_name_none(roms):
    nc_ds = NetCDFDataset(roms)
    result = nc_ds.find_variables_by_attr(standard_name='some standard_name')
    assert result == []

def test_sgrid_compliant_check(roms):
    nc_ds = NetCDFDataset(roms)
    result = nc_ds.sgrid_compliant_file()
    assert result == True


class TestNetCDFDatasetWithoutNodes(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = wrf_sgrid_2d()

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)

    def setUp(self):
        self.ds = Dataset(self.sgrid_test_file)
        self.nc_ds = NetCDFDataset(self.ds)

    def tearDown(self):
        self.ds.close()

    def test_node_coordinates(self):
        node_coordinates = self.nc_ds.find_node_coordinates('west_east_stag south_north_stag')
        self.assertIsNone(node_coordinates)

    def test_find_variable_by_attr(self):
        result = self.nc_ds.find_variables_by_attr(cf_role='grid_topology', topology_dimension=2)
        expected = ['grid']
        assert result == expected

    def test_find_variable_by_nonexistant_attr(self):
        result = self.nc_ds.find_variables_by_attr(bird='tufted titmouse')
        self.assertEqual(result, [])
