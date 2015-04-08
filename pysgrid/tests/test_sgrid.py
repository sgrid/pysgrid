'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest
import netCDF4 as nc4
from ..sgrid import SGrid


test_files = os.path.join(os.path.split(__file__)[0], 'files')


class TestSGridCreate(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(test_files, 'test_sgrid.nc')
        self.sg = SGrid()
        
    def test_load_from_file(self):
        sg_obj = self.sg.from_nc_file(self.sgrid_test_file)
        self.assertIsInstance(sg_obj, SGrid)
        
    def test_load_from_dataset(self):
        ds = nc4.Dataset(self.sgrid_test_file)
        sg_obj = self.sg.from_nc_dataset(ds)
        self.assertIsInstance(sg_obj, SGrid)
        

class TestSGridAttributes(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(test_files, 'test_sgrid.nc')
        self.sg_obj = SGrid().from_nc_file(self.sgrid_test_file)
        
    def test_centers(self):
        centers = self.sg_obj.centers
        centers_shape = centers.shape
        expected_shape = (4, 4, 2)
        self.assertEqual(centers_shape, expected_shape)