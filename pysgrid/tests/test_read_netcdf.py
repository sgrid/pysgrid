'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest
import netCDF4 as nc4
from ..read_netcdf import NetCDFDataset


test_files = os.path.join(os.path.split(__file__)[0], 'files')


class TestNetCDFDataset(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(test_files, 'test_sgrid.nc')
        self.ds = nc4.Dataset(self.sgrid_test_file)
        self.nc_ds = NetCDFDataset(self.ds)
        
    def test_finding_center_variables(self):
        result = self.nc_ds.find_grid_cell_center_vars()
        expected = ('y_center', 'x_center')
        self.assertEqual(result, expected)
        
    def test_finding_node_variables(self):
        result = self.nc_ds.find_grid_cell_node_vars()
        expected = ('y_node', 'x_node')
        self.assertEqual(result, expected)
        
    def test_find_grid_topology(self):
        result = self.nc_ds.find_grid_topology_vars()
        expected = ['grid']
        self.assertEqual(result, expected)
        
    def test_sgrid_compliant_check(self):
        result = self.nc_ds.sgrid_compliant_file()
        self.assertTrue(result)