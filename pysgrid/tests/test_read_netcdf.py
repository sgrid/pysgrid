'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest
import netCDF4 as nc4
from ..read_netcdf import NetCDFDataset


CURRENT_DIR = os.path.dirname(__file__)
TEST_FILES = os.path.join(CURRENT_DIR, 'files')


class TestNetCDFDataset(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_roms_like.nc')
        self.ds = nc4.Dataset(self.sgrid_test_file)
        self.nc_ds = NetCDFDataset(self.ds)
        
    def test_finding_center_variables(self):
        result = self.nc_ds.find_grid_cell_center_vars()
        expected = ('lon_center', 'lat_center')
        self.assertEqual(result, expected)
        
    def test_finding_node_variables(self):
        result = self.nc_ds.find_grid_cell_node_vars()
        expected = ('lon_node', 'lat_node')
        self.assertEqual(result, expected)
        
    def test_find_coordinations_by_location(self):
        result = self.nc_ds.find_coordinations_by_location('faces', 2)
        expected = ('lon_center', 'lat_center')
        self.assertEqual(result, expected)
        
    def test_find_grid_topology(self):
        result = self.nc_ds.find_grid_topology_vars()
        expected = ['grid']
        self.assertEqual(result, expected)
        
    def test_sgrid_compliant_check(self):
        result = self.nc_ds.sgrid_compliant_file()
        self.assertTrue(result)