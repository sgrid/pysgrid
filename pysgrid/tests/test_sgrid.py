'''
Created on Apr 7, 2015

@author: ayan
'''
import os
import unittest

import mock
import netCDF4 as nc4
import numpy as np

from ..custom_exceptions import SGridNonCompliantError
from ..sgrid import SGrid2D, SGrid3D, from_ncfile, from_nc_dataset
from ..utils import GridPadding
from .write_nc_test_files import (deltares_sgrid, deltares_sgrid_no_optional_attr, 
                                  non_compliant_sgrid, roms_sgrid, wrf_sgrid, 
                                  wrf_sgrid_2d)


CURRENT_DIR = os.path.dirname(__file__)
TEST_FILES = os.path.join(CURRENT_DIR, 'files')


class TestSGridCompliant(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = non_compliant_sgrid()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
        
    def test_exception_raised(self):
        self.assertRaises(SGridNonCompliantError, 
                          from_ncfile,
                          self.sgrid_test_file
                          )


class TestSGridCreate(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = roms_sgrid()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
        
    def setUp(self):
        self.ds = nc4.Dataset(self.sgrid_test_file)
        
    def tearDown(self):
        self.ds.close()
  
    def test_load_from_file(self):
        sg_obj = from_ncfile(self.sgrid_test_file)
        self.assertIsInstance(sg_obj, SGrid2D)

    def test_load_from_dataset(self):
        sg_obj = from_nc_dataset(self.ds)
        self.assertIsInstance(sg_obj, SGrid2D)


class TestSGridRomsDataset(unittest.TestCase):
    """
    Test using a representative ROMS file.
    
    """
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = roms_sgrid()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
    
    def setUp(self):
        self.sg_obj = from_ncfile(self.sgrid_test_file)
        self.write_path = os.path.join(CURRENT_DIR, 'test_sgrid_write.nc')
  
    def test_centers(self):
        centers = self.sg_obj.centers
        centers_shape = centers.shape
        expected_shape = (4, 4, 2)
        self.assertEqual(centers_shape, expected_shape)
    
    def test_variables(self):
        dataset_vars = self.sg_obj.variables
        expected_vars = [u's_rho', 
                         u's_w', 
                         u'time', 
                         u'xi_rho', 
                         u'eta_rho', 
                         u'xi_psi', 
                         u'eta_psi', 
                         u'xi_u', 
                         u'eta_u', 
                         u'xi_v', 
                         u'eta_v', 
                         u'grid', 
                         u'u', 
                         u'v', 
                         u'lon_rho', 
                         u'lat_rho', 
                         u'lon_psi', 
                         u'lat_psi', 
                         u'lat_u', 
                         u'lon_u', 
                         u'lat_v', 
                         u'lon_v', 
                         u'zeta'
                         ]
        self.assertEqual(len(dataset_vars), len(expected_vars))
        self.assertEqual(dataset_vars, expected_vars)
        
    def test_grid_variables(self):
        dataset_grid_variables = self.sg_obj.grid_variables
        expected_grid_variables = [u'u', u'v']
        self.assertEqual(len(dataset_grid_variables), len(expected_grid_variables))
        self.assertEqual(set(dataset_grid_variables), set(expected_grid_variables))
        
    def test_non_grid_variables(self):
        dataset_non_grid_variables = self.sg_obj.non_grid_variables
        expected_non_grid_variables = [u's_rho', 
                                       u's_w', 
                                       u'time', 
                                       u'xi_rho', 
                                       u'eta_rho', 
                                       u'xi_psi', 
                                       u'eta_psi', 
                                       u'xi_u', 
                                       u'eta_u', 
                                       u'xi_v', 
                                       u'eta_v', 
                                       u'grid', 
                                       u'lon_rho', 
                                       u'lat_rho', 
                                       u'lon_psi', 
                                       u'lat_psi', 
                                       u'lat_u', 
                                       u'lon_u', 
                                       u'lat_v', 
                                       u'lon_v', 
                                       u'zeta'
                                       ]
        self.assertEqual(len(dataset_non_grid_variables), len(expected_non_grid_variables))
        self.assertEqual(set(dataset_non_grid_variables), set(expected_non_grid_variables))

    def test_variable_slicing(self):
        u_center_slices = self.sg_obj.u.center_slicing
        v_center_slices = self.sg_obj.v.center_slicing
        u_center_expected = (np.s_[:], np.s_[:], np.s_[1:-1], np.s_[:])
        v_center_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:-1])
        self.assertEqual(u_center_slices, u_center_expected)
        self.assertEqual(v_center_slices, v_center_expected)
        
    def test_grid_variable_average_axes(self):
        uc_axis = self.sg_obj.u.center_axis
        uc_axis_expected = 1
        un_axis = self.sg_obj.u.node_axis
        un_axis_expected = 0
        lon_rho_c_axis = self.sg_obj.lon_rho.center_axis
        lon_rho_n_axis = self.sg_obj.lon_rho.node_axis
        self.assertEqual(uc_axis, uc_axis_expected)
        self.assertEqual(un_axis, un_axis_expected)
        self.assertIsNone(lon_rho_c_axis)
        self.assertIsNone(lon_rho_n_axis)
        
    def test_optional_grid_attrs(self):
        face_coordinates = self.sg_obj.face_coordinates
        node_coordinates = self.sg_obj.node_coordinates
        edge1_coordinates = self.sg_obj.edge1_coordinates
        edge2_coordinates = self.sg_obj.edge2_coordinates
        fc_expected = ('lon_rho', 'lat_rho')
        nc_expected = ('lon_psi', 'lat_psi')
        e1c_expected = ('lon_u', 'lat_u')
        e2c_expected = ('lon_v', 'lat_v')
        self.assertEqual(face_coordinates, fc_expected)
        self.assertEqual(node_coordinates, nc_expected)
        self.assertEqual(edge1_coordinates, e1c_expected)
        self.assertEqual(edge2_coordinates, e2c_expected)
    
    @mock.patch('pysgrid.sgrid.nc4')
    def test_write_sgrid_to_netcdf(self, mock_nc):
        self.sg_obj.save_as_netcdf(self.write_path)
        mock_nc.Dataset.assert_called_with(self.write_path, 'w')
        

class TestSGridNoCoordinates(unittest.TestCase):
    """
    Test to make sure that if no coordinates (e.g. face, edge1, etc)
    are specified, those coordinates can be inferred from the dataset.
    
    A file is representing a delft3d dataset is used for this test.
    
    """
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = deltares_sgrid_no_optional_attr()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
        
    def setUp(self):
        self.sgrid_obj = from_ncfile(self.sgrid_test_file)
        
    def test_face_coordinate_inference(self):
        face_coordinates = self.sgrid_obj.face_coordinates
        expected_face_coordinates = (u'XZ', u'YZ')
        self.assertEqual(face_coordinates, expected_face_coordinates)
        
    def test_grid_cell_centers(self):
        centers = self.sgrid_obj.centers
        centers_shape = (4, 4, 2)
        self.assertEqual(centers.shape, centers_shape)
        
    def test_grid_cell_nodes(self):
        nodes = self.sgrid_obj.nodes
        nodes_shape = (4, 4, 2)
        self.assertEqual(nodes.shape, nodes_shape)
        
    def test_grid_angles(self):
        angles = self.sgrid_obj.angles
        angles_shape = (4, 4)
        self.assertEqual(angles.shape, angles_shape)

        
class TestSGridWRFDataset(unittest.TestCase):
    """
    Test a representative WRF file.
    
    """
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = wrf_sgrid_2d()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
        
    def setUp(self):
        self.sg_obj = from_ncfile(self.sgrid_test_file)
        
    def test_topology_dimension(self):
        topology_dim = self.sg_obj.topology_dimension
        expected_dim = 2
        self.assertEqual(topology_dim, expected_dim)
        
    def test_variable_slicing(self):
        u_slice = self.sg_obj.U.center_slicing
        u_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[:])
        v_slice = self.sg_obj.V.center_slicing
        v_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[:])
        self.assertEqual(u_slice, u_expected)
        self.assertEqual(v_slice, v_expected)
        
    def test_variable_average_axes(self):
        u_avg_axis = self.sg_obj.U.center_axis
        u_axis_expected = 1
        v_avg_axis = self.sg_obj.V.center_axis
        v_axis_expected = 0
        self.assertEqual(u_avg_axis, u_axis_expected)
        self.assertEqual(v_avg_axis, v_axis_expected)
        

class TestSGridDelft3dDataset(unittest.TestCase):
    """
    Test using a representative delft3d file.
    
    """
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = deltares_sgrid()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
    
    def setUp(self):
        self.sg_obj = from_ncfile(self.sgrid_test_file)
        
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
        
    def test_averaging_axes(self):
        u1c_axis = self.sg_obj.U1.center_axis
        u1c_expected = 0
        v1n_axis = self.sg_obj.V1.node_axis
        v1n_expected = 0
        latitude_c_axis = self.sg_obj.latitude.center_axis
        latitude_n_axis = self.sg_obj.latitude.node_axis
        self.assertEqual(u1c_axis, u1c_expected)
        self.assertEqual(v1n_axis, v1n_expected)
        self.assertIsNone(latitude_c_axis)
        self.assertIsNone(latitude_n_axis)
        
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
        self.assertEqual(set(grid_variables), set(expected_grid_variables))
        
    def test_angles(self):
        angles = self.sg_obj.angles
        expected_shape = (4, 4)
        self.assertEqual(angles.shape, expected_shape)
        
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
        
        
class TestSGridSaveNoNodeCoordinates(unittest.TestCase):
    """
    Test that SGrid.save_as_netcdf is saving content
    when there are no nodes or node coordinates specified.
    
    This scenario will typically occur with WRF datasets.
    
    """
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = wrf_sgrid_2d()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
        
    def setUp(self):
        self.sgrid_target = os.path.join(TEST_FILES, 'tmp_sgrid.nc')
        self.sg_obj = from_ncfile(self.sgrid_test_file)
        self.sg_obj.save_as_netcdf(self.sgrid_target)
        self.target = from_ncfile(self.sgrid_target)
        
    def tearDown(self):
        os.remove(self.sgrid_target)
        
    def test_sgrid(self):
        self.assertIsInstance(self.target, SGrid2D)
        
    def test_nodes(self):
        nodes = self.target.nodes
        self.assertIsNone(nodes)
        
    def test_node_coordinates(self):
        node_coordinates = self.target.node_coordinates
        self.assertIsNone(node_coordinates)
        
    def test_node_dimesnions(self):
        node_dims = self.target.node_dimensions
        expected = 'west_east_stag south_north_stag'
        self.assertEqual(node_dims, expected)


class TestSGridSaveNodeCoordinates(unittest.TestCase):
    """
    Test that SGrid.save_as_netcdf is saving
    content correctly.
    
    There maybe a better way to do this using
    mocks, but this will do for now.
    
    """
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = deltares_sgrid()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
        
    def setUp(self):
        self.sgrid_target = os.path.join(TEST_FILES, 'tmp_sgrid.nc')
        self.sg_obj = from_ncfile(self.sgrid_test_file)
        self.sg_obj.save_as_netcdf(self.sgrid_target)
        self.target = from_ncfile(self.sgrid_target)
        
    def tearDown(self):
        os.remove(self.sgrid_target)
        
    def test_save_as_netcdf(self):
        """
        Test that the attributes in the
        saved netCDF file are as expected.
        
        """
        target_dims = self.target.dimensions
        expected_target_dims = [(u'MMAXZ', 4), 
                                (u'NMAXZ', 4), 
                                (u'MMAX', 4), 
                                (u'NMAX', 4), 
                                (u'KMAX', 2), 
                                (u'KMAX1', 3), 
                                (u'time', 2)
                                ]
        target_vars = self.target.variables
        expected_target_vars = [u'XZ', 
                                u'YZ', 
                                u'XCOR', 
                                u'YCOR', 
                                u'grid', 
                                u'time', 
                                u'U1', 
                                u'V1', 
                                u'latitude', 
                                u'longitude', 
                                u'grid_latitude', 
                                u'grid_longitude'
                                ]
        target_grid_vars = self.target.grid_variables
        expected_target_grid_vars = [u'U1', 
                                     u'V1'
                                     ]
        target_face_coordinates = self.target.face_coordinates
        expected_target_face_coordinates = (u'XZ', u'YZ')
        self.assertIsInstance(self.target, SGrid2D)
        self.assertEqual(len(target_dims), len(expected_target_dims))
        self.assertEqual(set(target_dims), set(expected_target_dims))
        self.assertEqual(len(target_vars), len(expected_target_vars))
        self.assertEqual(set(target_vars), set(expected_target_vars))
        self.assertEqual(len(target_grid_vars), len(expected_target_grid_vars))
        self.assertEqual(set(target_grid_vars), set(expected_target_grid_vars))
        self.assertEqual(target_face_coordinates, expected_target_face_coordinates)
        
    def test_saved_sgrid_attributes(self):
        """
        Test that calculated/inferred attributes
        are as expected from the saved filed.
        
        """
        u1_var = self.target.U1
        u1_var_center_avg_axis = u1_var.center_axis
        expected_u1_center_axis = 0
        u1_vector_axis = u1_var.vector_axis
        expected_u1_vector_axis = 'X'
        original_angles = self.sg_obj.angles
        saved_angles = self.target.angles
        self.assertEqual(u1_var_center_avg_axis, expected_u1_center_axis)
        self.assertEqual(u1_vector_axis, expected_u1_vector_axis)
        np.testing.assert_almost_equal(original_angles, saved_angles, decimal=3)
        
        
class Test3DimensionalSGrid(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.sgrid_test_file = wrf_sgrid()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.sgrid_test_file)
    
    def setUp(self):
        self.sg_obj = from_ncfile(self.sgrid_test_file)
        
    def test_sgrid_instance(self):
        self.assertIsInstance(self.sg_obj, SGrid3D)
        
    def test_variables(self):
        sg_vars = self.sg_obj.variables
        sg_vars_expected = [u'Times',
                            u'XTIME',
                            u'U',
                            u'V',
                            u'W', 
                            u'T',
                            u'XLAT',
                            u'XLONG', 
                            u'ZNU',
                            u'ZNW',
                            u'grid'
                            ]
        self.assertEqual(len(sg_vars), len(sg_vars_expected))
        self.assertEqual(set(sg_vars), set(sg_vars_expected))
        
    def test_volume_padding(self):
        volume_padding = self.sg_obj.volume_padding
        volume_padding_expected = [GridPadding(mesh_topology_var=u'grid', face_dim=u'west_east', node_dim=u'west_east_stag', padding=u'none'), 
                                   GridPadding(mesh_topology_var=u'grid', face_dim=u'south_north', node_dim=u'south_north_stag', padding=u'none'), 
                                   GridPadding(mesh_topology_var=u'grid', face_dim=u'bottom_top', node_dim=u'bottom_top_stag', padding=u'none')
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