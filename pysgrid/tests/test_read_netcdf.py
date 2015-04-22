'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest
import netCDF4 as nc4
from ..read_netcdf import NetCDFDataset, parse_padding
from ..custom_exceptions import CannotFindPaddingError



CURRENT_DIR = os.path.dirname(__file__)
TEST_FILES = os.path.join(CURRENT_DIR, 'files')


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
        sub_dim = padding_datum_0.sub_dim
        dim = padding_datum_0.dim
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
        sub_dim = padding_datum_0.sub_dim
        dim = padding_datum_0.dim
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


class TestNetCDFDataset(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_roms.nc')
        self.ds = nc4.Dataset(self.sgrid_test_file)
        self.nc_ds = NetCDFDataset(self.ds)
        
    def test_finding_node_variables(self):
        result = self.nc_ds.find_grid_cell_node_vars()
        expected = ('lon_node', 'lat_node')
        self.assertEqual(result, expected)
        
    def test_find_coordinatates_by_location(self):
        result = self.nc_ds.find_coordinates_by_location('faces', 2)
        expected = ('lon_center', 'lat_center')
        self.assertEqual(result, expected)
        
    def test_find_grid_topology(self):
        result = self.nc_ds.find_grid_topology_vars()
        expected = 'grid'
        self.assertEqual(result, expected)
        
    def test_sgrid_compliant_check(self):
        result = self.nc_ds.sgrid_compliant_file()
        self.assertTrue(result)