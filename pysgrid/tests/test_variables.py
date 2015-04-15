'''
Created on Apr 15, 2015

@author: ayan
'''
import os
import unittest
import netCDF4 as nc4
from ..variables import SGridVariable


CURRENT_DIR = os.path.dirname(__file__)
TEST_FILES = os.path.join(CURRENT_DIR, 'files')


class TestSGridVariable(unittest.TestCase):
    
    def setUp(self):
        self.test_file = os.path.join(TEST_FILES, 'test_sgrid_roms_like.nc')
        self.dataset = nc4.Dataset(self.test_file)
        self.test_var_1 = self.dataset.variables['u']
        self.test_var_2 = self.dataset.variables['zeta']
        
    def test_create_sgrid_variable_object(self):
        sgrid_var = SGridVariable.create_variable(self.test_var_1)
        self.assertIsInstance(sgrid_var, SGridVariable)
        
    def test_attributes_with_grid(self):
        sgrid_var = SGridVariable.create_variable(self.test_var_1)
        sgrid_var_name = sgrid_var.variable
        sgrid_var_name_expected = 'u'
        sgrid_var_dim = sgrid_var.dimensions
        sgrid_var_dim_expected = ('time', 'z_center', 'y_u', 'x_u')
        sgrid_var_grid = sgrid_var.grid
        sgrid_var_grid_expected = 'some grid'
        sgrid_var_location = sgrid_var.location
        self.assertEqual(sgrid_var_name, sgrid_var_name_expected)
        self.assertEqual(sgrid_var_dim, sgrid_var_dim_expected)
        self.assertEqual(sgrid_var_grid, sgrid_var_grid_expected)
        self.assertIsNone(sgrid_var_location)
        
    def test_attributes_with_location(self):
        sgrid_var = SGridVariable.create_variable(self.test_var_2)
        sgrid_var_name = sgrid_var.variable
        sgrid_var_name_expected = 'zeta'
        sgrid_var_dim = sgrid_var.dimensions
        sgrid_var_grid = sgrid_var.grid
        sgrid_var_location = sgrid_var.location
        sgrid_var_location_expected = 'faces'
        sgrid_var_dim_expected = ('time', 'y_center', 'x_center')
        self.assertEqual(sgrid_var_name, sgrid_var_name_expected)
        self.assertEqual(sgrid_var_dim, sgrid_var_dim_expected)
        self.assertIsNone(sgrid_var_grid)
        self.assertEqual(sgrid_var_location, sgrid_var_location_expected)
        