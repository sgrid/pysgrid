'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest
import netCDF4 as nc4
import numpy as np
import mock
from pysgrid.sgrid import SGridND, SGrid
from ..utils import GridPadding
from ..custom_exceptions import SGridNonCompliantError


CURRENT_DIR = os.path.dirname(__file__)
TEST_FILES = os.path.join(CURRENT_DIR, 'files')


class TestSGridCompliant(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_noncompliant_sgrid_roms.nc')
        self.sg = SGrid
        
    def test_exception_raised(self):
        self.assertRaises(SGridNonCompliantError, 
                          self.sg.from_nc_file, 
                          self.sgrid_test_file
                          )


class TestSGridCreate(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_roms.nc')
        self.sg = SGrid
  
    def test_load_from_file(self):
        sg_obj = self.sg.from_nc_file(self.sgrid_test_file)
        self.assertIsInstance(sg_obj, SGridND)

    def test_load_from_dataset(self):
        ds = nc4.Dataset(self.sgrid_test_file)
        sg_obj = self.sg.from_nc_dataset(ds)
        self.assertIsInstance(sg_obj, SGridND)
    

class TestSGridWithOptionalAttributes(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_roms.nc')
        self.sg_obj = SGrid.from_nc_file(self.sgrid_test_file)
        self.write_path = os.path.join(CURRENT_DIR, 'test_sgrid_write.nc')
  
    def test_centers(self):
        centers = self.sg_obj.centers
        centers_shape = centers.shape
        expected_shape = (4, 4, 2)
        self.assertEqual(centers_shape, expected_shape)
    
    def test_variables(self):
        dataset_vars = self.sg_obj.variables
        expected_vars = [u'z_center', 
                         u'z_node', 
                         u'time', 
                         u'x_center', 
                         u'y_center', 
                         u'x_node', 
                         u'y_node', 
                         u'x_u', 
                         u'y_u', 
                         u'x_v', 
                         u'y_v', 
                         u'grid', 
                         u'u', 
                         u'v', 
                         u'lon_center', 
                         u'lat_center', 
                         u'lon_node', 
                         u'lat_node', 
                         u'lat_u', 
                         u'lon_u', 
                         u'lat_v', 
                         u'lon_v', 
                         u'zeta'
                         ]
        self.assertEqual(dataset_vars, expected_vars)

    def test_variable_slicing(self):
        u_center_slices = self.sg_obj.u.center_slicing
        v_center_slices = self.sg_obj.v.center_slicing
        u_center_expected = (np.s_[:], np.s_[:], np.s_[1:-1], np.s_[:])
        v_center_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:-1])
        self.assertEqual(u_center_slices, u_center_expected)
        self.assertEqual(v_center_slices, v_center_expected)
        
    def test_optional_grid_attrs(self):
        face_coordinates = self.sg_obj.face_coordinates
        node_coordinates = self.sg_obj.node_coordinates
        edge1_coordinates = self.sg_obj.edge1_coordinates
        edge2_coordinates = self.sg_obj.edge2_coordinates
        fc_expected = ('lon_center', 'lat_center')
        nc_expected = ('lon_node', 'lat_node')
        e1c_expected = ('lon_u', 'lat_u')
        e2c_expected = ('lon_v', 'lat_v')
        self.assertEqual(face_coordinates, fc_expected)
        self.assertEqual(node_coordinates, nc_expected)
        self.assertEqual(edge1_coordinates, e1c_expected)
        self.assertEqual(edge2_coordinates, e2c_expected)
        
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
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_deltares.nc')
        self.sg_obj = SGrid.from_nc_file(self.sgrid_test_file)
        
    def test_centers(self):
        centers = self.sg_obj.centers
        centers_shape = centers.shape
        expected_shape = (4, 4, 2)
        self.assertEqual(centers_shape, expected_shape)
        
    def test_topology_dimension(self):
        topology_dim = self.sg_obj.topology_dimension
        expected_dim = 2
        self.assertEqual(topology_dim, expected_dim)
        
    def test_variable_slice(self):
        u_center_slices = self.sg_obj.U1.center_slicing
        v_center_slices = self.sg_obj.V1.center_slicing
        u_center_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:])
        v_center_expected = (np.s_[:], np.s_[:], np.s_[1:], np.s_[:])
        xz_center_slices = self.sg_obj.XZ.center_slicing
        xcor_center_slices = self.sg_obj.XCOR.center_slicing
        xz_center_expected = (np.s_[1:], np.s_[1:])
        xcor_center_expected  = (np.s_[:], np.s_[:])
        self.assertEqual(u_center_slices, u_center_expected)
        self.assertEqual(v_center_slices, v_center_expected)
        self.assertEqual(xz_center_slices, xz_center_expected)
        self.assertEqual(xcor_center_slices, xcor_center_expected)
        
    def test_grid_optional_attrs(self):
        face_coordinates = self.sg_obj.face_coordinates
        node_coordinates = self.sg_obj.node_coordinates
        edge1_coordinates = self.sg_obj.edge1_coordinates
        edge2_coordinates = self.sg_obj.edge2_coordinates
        fc_expected = ('XZ', 'YZ')
        nc_expected = ('XCOR', 'YCOR')
        self.assertEqual(face_coordinates, fc_expected)
        self.assertEqual(node_coordinates, nc_expected)
        self.assertIsNone(edge1_coordinates)
        self.assertIsNone(edge2_coordinates)
        
    def test_grid_variables(self):
        grid_variables = self.sg_obj.grid_variables
        expected_grid_variables = ['U1', 'V1']
        self.assertEqual(grid_variables, expected_grid_variables)
        
    def test_no_3d_attributes(self):
        self.assertFalse(hasattr(self.sg_obj, 'volume_padding'))
        self.assertFalse(hasattr(self.sg_obj, 'volume_dimensions'))
        self.assertFalse(hasattr(self.sg_obj, 'volume_coordinates'))
        self.assertFalse(hasattr(self.sg_obj, 'face1_padding'))
        self.assertFalse(hasattr(self.sg_obj, 'face1_coordinates'))
        self.assertFalse(hasattr(self.sg_obj, 'face1_dimensions'))
        self.assertFalse(hasattr(self.sg_obj, 'face2_padding'))
        self.assertFalse(hasattr(self.sg_obj, 'face2_coordinates'))
        self.assertFalse(hasattr(self.sg_obj, 'face2_dimensions'))
        self.assertFalse(hasattr(self.sg_obj, 'face3_padding'))
        self.assertFalse(hasattr(self.sg_obj, 'face3_coordinates'))
        self.assertFalse(hasattr(self.sg_obj, 'edge3_padding'))
        self.assertFalse(hasattr(self.sg_obj, 'edge3_coordinates'))
        self.assertFalse(hasattr(self.sg_obj, 'edge3_dimensions'))
        
    def test_2d_attributes(self):
        self.assertTrue(hasattr(self.sg_obj, 'face_padding'))
        self.assertTrue(hasattr(self.sg_obj, 'face_coordinates'))
        self.assertTrue(hasattr(self.sg_obj, 'face_dimensions'))
        self.assertTrue(hasattr(self.sg_obj, 'vertical_padding'))
        self.assertTrue(hasattr(self.sg_obj, 'vertical_dimensions'))


class TestSGridSave(unittest.TestCase):
    """
    Test that SGrid.save_as_netcdf is saving
    content correctly.
    
    There maybe a better way to do this using
    mocks, but this will do for now.
    """
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES,'test_sgrid_deltares.nc' )
        self.sgrid_target = os.path.join(TEST_FILES, 'tmp_sgrid.nc')
        self.sg_obj = SGrid.from_nc_file(self.sgrid_test_file)
        
    def test_save_as_netcdf(self):
        self.sg_obj.save_as_netcdf(self.sgrid_target)
        target = SGrid.from_nc_file(self.sgrid_target)
        target_dims = target.dimensions
        expected_target_dims = [(u'MMAXZ', 4), 
                                (u'NMAXZ', 4), 
                                (u'MMAX', 4), 
                                (u'NMAX', 4), 
                                (u'KMAX', 2), 
                                (u'KMAX1', 3), 
                                (u'time', 2)
                                ]
        target_vars = target.variables
        expected_target_vars = [u'XZ', 
                                u'YZ', 
                                u'XCOR', 
                                u'YCOR', 
                                u'grid', 
                                u'time', 
                                u'U1', 
                                u'V1'
                                ]
        target_grid_vars = target.grid_variables
        expected_target_grid_vars = [u'U1', 
                                     u'V1'
                                     ]
        target_face_coordinates = target.face_coordinates
        expected_target_face_coordinates = (u'XZ', u'YZ')
        self.assertIsInstance(target, SGridND)
        self.assertEqual(target_dims, expected_target_dims)
        self.assertEqual(target_vars, expected_target_vars)
        self.assertEqual(target_grid_vars, expected_target_grid_vars)
        self.assertEqual(target_face_coordinates, expected_target_face_coordinates)      
      
    def tearDown(self):
        os.remove(self.sgrid_target)
        
        
class Test3DimensionalSGrid(unittest.TestCase):
    
    def setUp(self):
        self.sgrid_test_file = os.path.join(TEST_FILES, 'test_sgrid_wrf.nc')
        self.sg_obj = SGrid.from_nc_file(self.sgrid_test_file)
        
    def test_sgrid_instance(self):
        self.assertIsInstance(self.sg_obj, SGridND)
        
    def test_variables(self):
        sg_vars = self.sg_obj.variables
        sg_vars_expected = [u'Times', u'U', u'V', u'W', 
                            u'T', u'XLAT', u'XLONG', 
                            u'ZNU', u'ZNW', u'grid'
                            ]
        self.assertEqual(sg_vars, sg_vars_expected)
        
    def test_volume_padding(self):
        volume_padding = self.sg_obj.volume_padding
        volume_padding_expected = [GridPadding(mesh_topology_var=u'grid', dim=u'west_east', sub_dim=u'west_east_stag', padding=u'none'), 
                                   GridPadding(mesh_topology_var=u'grid', dim=u'south_north', sub_dim=u'south_north_stag', padding=u'none'), 
                                   GridPadding(mesh_topology_var=u'grid', dim=u'bottom_top', sub_dim=u'bottom_top_stag', padding=u'none')
                                   ]
        self.assertEqual(volume_padding, volume_padding_expected)
        
    def test_volume_coordinates(self):
        volume_coordinates = self.sg_obj.volume_coordinates
        volume_coordinates_expected = (u'XLONG', u'XLAT', u'ZNU')
        self.assertEqual(volume_coordinates, volume_coordinates_expected)
        
    def test_slicing_assignment(self):
        u_center_slice = self.sg_obj.U.center_slicing
        u_center_slice_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[:])
        self.assertEqual(u_center_slice, u_center_slice_expected)
        
    def test_sgrid_centers(self):
        centers_shape = self.sg_obj.centers.shape
        expected_shape = (2, 5, 4, 2)
        self.assertEqual(centers_shape, expected_shape)
        
    def test_topology_dimension(self):
        topology_dim = self.sg_obj.topology_dimension
        expected_topology_dim = 3
        self.assertEqual(topology_dim, expected_topology_dim)
        
    def test_no_2d_attributes(self):
        self.assertFalse(hasattr(self.sg_obj, 'face_padding'))
        self.assertFalse(hasattr(self.sg_obj, 'face_coordinates'))
        self.assertFalse(hasattr(self.sg_obj, 'face_dimensions'))
        self.assertFalse(hasattr(self.sg_obj, 'vertical_padding'))
        self.assertFalse(hasattr(self.sg_obj, 'vertical_dimensions'))
        
    def test_3d_attributes(self):
        self.assertTrue(hasattr(self.sg_obj, 'volume_padding'))
        self.assertTrue(hasattr(self.sg_obj, 'volume_dimensions'))
        self.assertTrue(hasattr(self.sg_obj, 'volume_coordinates'))
        self.assertTrue(hasattr(self.sg_obj, 'face1_padding'))
        self.assertTrue(hasattr(self.sg_obj, 'face1_coordinates'))
        self.assertTrue(hasattr(self.sg_obj, 'face1_dimensions'))
        self.assertTrue(hasattr(self.sg_obj, 'face2_padding'))
        self.assertTrue(hasattr(self.sg_obj, 'face2_coordinates'))
        self.assertTrue(hasattr(self.sg_obj, 'face2_dimensions'))
        self.assertTrue(hasattr(self.sg_obj, 'face3_padding'))
        self.assertTrue(hasattr(self.sg_obj, 'face3_coordinates'))
        self.assertTrue(hasattr(self.sg_obj, 'edge3_padding'))
        self.assertTrue(hasattr(self.sg_obj, 'edge3_coordinates'))
        self.assertTrue(hasattr(self.sg_obj, 'edge3_dimensions'))