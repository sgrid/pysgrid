'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest
import netCDF4 as nc4
import numpy as np
from ..sgrid import SGrid


test_files = os.path.join(os.path.split(__file__)[0], 'files')


class TestSGridCreate(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(test_files, 'test_sgrid_roms_like.nc')
        self.sg = SGrid()
  
    def test_load_from_file(self):
        sg_obj = self.sg.from_nc_file(self.sgrid_test_file)
        self.assertIsInstance(sg_obj, SGrid)

    def test_load_from_dataset(self):
        ds = nc4.Dataset(self.sgrid_test_file)
        sg_obj = self.sg.from_nc_dataset(ds)
        self.assertIsInstance(sg_obj, SGrid)
        

class TestSGridWithEdgesAttributes(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(test_files, 'test_sgrid_roms_like.nc')
        self.sg_obj = SGrid().from_nc_file(self.sgrid_test_file)
  
    def test_centers(self):
        centers = self.sg_obj.centers
        centers_shape = centers.shape
        expected_shape = (4, 4, 2)
        self.assertEqual(centers_shape, expected_shape)
        
    def test_variable_slicing(self):
        u_slices = self.sg_obj.u_slice
        v_slices = self.sg_obj.v_slice
        u_expected = (np.s_[:], np.s_[:], np.s_[1:-1], np.s_[:])
        v_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:-1])
        self.assertEqual(u_slices, u_expected)
        self.assertEqual(v_slices, v_expected)
        

class TestSGridWithOutEdgesAttributes(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(test_files, 'test_sgrid_deltares_like.nc')
        self.sg_obj = SGrid().from_nc_file(self.sgrid_test_file)
        
    def test_centers(self):
        centers = self.sg_obj.centers
        centers_shape = centers.shape
        expected_shape = (4, 4, 2)
        self.assertEqual(centers_shape, expected_shape)
        
    def test_variable_slice(self):
        u_slices = self.sg_obj.U1_slice
        v_slices = self.sg_obj.V1_slice
        u_expected = (np.s_[:], np.s_[:], np.s_[1:], np.s_[:])
        v_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:])
        self.assertEqual(u_slices, u_expected)
        self.assertEqual(v_slices, v_expected)