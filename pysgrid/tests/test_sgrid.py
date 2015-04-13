'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest
import netCDF4 as nc4
import numpy as np
import mock
from ..sgrid import SGrid
from ..custom_exceptions import SGridNonCompliant


CURRENT_DIR = os.path.dirname(__file__)
TEST_FILES = os.path.join(CURRENT_DIR, 'files')


class TestSGridCompliant(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_noncompliant_sgrid_roms_like.nc')
        self.sg = SGrid
        
    def test_exception_raised(self):
        self.assertRaises(SGridNonCompliant, 
                          self.sg.from_nc_file, 
                          self.sgrid_test_file
                          )


class TestSGridCreate(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_roms_like.nc')
        self.sg = SGrid
  
    def test_load_from_file(self):
        sg_obj = self.sg.from_nc_file(self.sgrid_test_file)
        self.assertIsInstance(sg_obj, SGrid)

    def test_load_from_dataset(self):
        ds = nc4.Dataset(self.sgrid_test_file)
        sg_obj = self.sg.from_nc_dataset(ds)
        self.assertIsInstance(sg_obj, SGrid)
        

class TestSGridWithOptionalAttributes(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_roms_like.nc')
        self.sg_obj = SGrid.from_nc_file(self.sgrid_test_file)
        self.write_path = os.path.join(CURRENT_DIR, 'test_sgrid_write.nc')
  
    def test_centers(self):
        centers = self.sg_obj.centers
        centers_shape = centers.shape
        expected_shape = (4, 4, 2)
        self.assertEqual(centers_shape, expected_shape)
    
    def test_variables(self):
        dataset_vars = self.sg_obj.variables
        expected_vars = [(u'z_center', np.dtype('int32'), (u'z_center',)), 
                         (u'z_node', np.dtype('int32'), (u'z_node',)), 
                         (u'time', np.dtype('float64'), (u'time',)), 
                         (u'x_center', np.dtype('float32'), (u'x_center',)), 
                         (u'y_center', np.dtype('float32'), (u'y_center',)), 
                         (u'x_node', np.dtype('float32'), (u'x_node',)), 
                         (u'y_node', np.dtype('float32'), (u'y_node',)), 
                         (u'x_u', np.dtype('float32'), (u'x_u',)), 
                         (u'y_u', np.dtype('float32'), (u'y_u',)), 
                         (u'x_v', np.dtype('float32'), (u'x_v',)), 
                         (u'y_v', np.dtype('float32'), (u'y_v',)), 
                         (u'grid', np.dtype('int16'), ()), 
                         (u'u', np.dtype('float32'), (u'time', u'z_center', u'y_u', u'x_u')), 
                         (u'v', np.dtype('float32'), (u'time', u'z_center', u'y_v', u'x_v')), 
                         (u'lon_center', np.dtype('float32'), (u'y_center', u'x_center')), 
                         (u'lat_center', np.dtype('float32'), (u'y_center', u'x_center')), 
                         (u'lon_node', np.dtype('float32'), (u'y_node', u'x_node')), 
                         (u'lat_node', np.dtype('float32'), (u'y_node', u'x_node')), 
                         (u'lat_u', np.dtype('float32'), (u'y_u', u'x_u')), 
                         (u'lon_u', np.dtype('float32'), (u'y_u', u'x_u')), 
                         (u'lat_v', np.dtype('float32'), (u'y_v', u'x_v')), 
                         (u'lon_v', np.dtype('float32'), (u'y_v', u'x_v')),
                         (u'zeta', np.dtype('float32'), (u'time', u'y_center', u'x_center')),
                         ]
        self.assertEqual(dataset_vars, expected_vars)

    def test_variable_slicing(self):
        u_slices = self.sg_obj.u_slice
        v_slices = self.sg_obj.v_slice
        u_expected = (np.s_[:], np.s_[:], np.s_[1:-1], np.s_[:])
        v_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:-1])
        self.assertEqual(u_slices, u_expected)
        self.assertEqual(v_slices, v_expected)
        
    def test_optional_grid_attrs(self):
        face_coordinates = self.sg_obj.face_coordinates
        node_coordinates = self.sg_obj.node_coordinates
        edge_1_coordinates = self.sg_obj.edge_1_coordinates
        edge_2_coordinates = self.sg_obj.edge_2_coordinates
        fc_expected = ('lon_center', 'lat_center')
        nc_expected = ('lon_node', 'lat_node')
        e1c_expected = ('lon_u', 'lat_u')
        e2c_expected = ('lon_v', 'lat_v')
        self.assertEqual(face_coordinates, fc_expected)
        self.assertEqual(node_coordinates, nc_expected)
        self.assertEqual(edge_1_coordinates, e1c_expected)
        self.assertEqual(edge_2_coordinates, e2c_expected)
        
    def test_grid_variables(self):
        grid_variables = self.sg_obj.grid_variables
        expected_grid_variables = ['u', 'v']
        self.assertEqual(grid_variables, expected_grid_variables)
    
    @mock.patch('pysgrid.sgrid.nc4')
    def test_write_sgrid_to_netcdf(self, mock_nc):
        self.sg_obj.save_as_netcdf(self.write_path)
        mock_nc.Dataset.assert_called_with(self.write_path, 'w')
        

class TestSGridWithoutEdgesAttributes(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_deltares_like.nc')
        self.sg_obj = SGrid.from_nc_file(self.sgrid_test_file)
        
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
        xz_slices = self.sg_obj.XZ_slice
        xcor_slices = self.sg_obj.XCOR_slice
        xz_expected = (np.s_[1:], np.s_[1:])
        xcor_expected  = (np.s_[:], np.s_[:])
        self.assertEqual(u_slices, u_expected)
        self.assertEqual(v_slices, v_expected)
        self.assertEqual(xz_slices, xz_expected)
        self.assertEqual(xcor_slices, xcor_expected)
        
    def test_grid_optional_attrs(self):
        face_coordinates = self.sg_obj.face_coordinates
        node_coordinates = self.sg_obj.node_coordinates
        edge_1_coordinates = self.sg_obj.edge_1_coordinates
        edge_2_coordinates = self.sg_obj.edge_2_coordinates
        fc_expected = ('XZ', 'YZ')
        nc_expected = ('XCOR', 'YCOR')
        self.assertEqual(face_coordinates, fc_expected)
        self.assertEqual(node_coordinates, nc_expected)
        self.assertIsNone(edge_1_coordinates)
        self.assertIsNone(edge_2_coordinates)
        
    def test_grid_variables(self):
        grid_variables = self.sg_obj.grid_variables
        expected_grid_variables = ['U1', 'V1']
        self.assertEqual(grid_variables, expected_grid_variables)