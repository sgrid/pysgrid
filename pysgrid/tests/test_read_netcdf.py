'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest

import netCDF4 as nc4

from ..custom_exceptions import CannotFindPaddingError
from ..read_netcdf import NetCDFDataset, parse_axes, parse_padding, parse_vector_axis
from .write_nc_test_files import roms_sgrid


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
        self.assertRaises(CannotFindPaddingError, 
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

        
class TestNetCDFDataset(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = roms_sgrid()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
    
    def setUp(self):
        self.ds = nc4.Dataset(self.sgrid_test_file)
        self.nc_ds = NetCDFDataset(self.ds)
        
    def tearDown(self):
        self.ds.close()
        
    def test_finding_node_variables(self):
        result = self.nc_ds.find_grid_cell_node_vars()
        expected = ('lon_psi', 'lat_psi')
        self.assertEqual(result, expected)
        
    def test_find_coordinatates_by_location(self):
        result = self.nc_ds.find_coordinates_by_location('faces', 2)
        expected = ('lon_rho', 'lat_rho')
        self.assertEqual(result, expected)
        
    def test_find_grid_topology(self):
        result = self.nc_ds.find_grid_topology_var()
        expected = 'grid'
        self.assertEqual(result, expected)
        
    def test_sgrid_compliant_check(self):
        result = self.nc_ds.sgrid_compliant_file()
        self.assertTrue(result)