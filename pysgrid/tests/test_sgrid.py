'''
Created on Apr 7, 2015

@author: ayan
'''

from __future__ import (absolute_import, division, print_function)

import os
import unittest

try:
    from unittest import mock
except ImportError:
    import mock

from netCDF4 import Dataset
import numpy as np

from ..sgrid import SGrid, load_sgrid
from .write_nc_test_files import (deltares_sgrid,
                                  deltares_sgrid_no_optional_attr,
                                  non_compliant_sgrid,
                                  roms_sgrid,
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
        self.assertRaises(ValueError,
                          load_sgrid,
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
        self.ds = Dataset(self.sgrid_test_file)

    def tearDown(self):
        self.ds.close()

    def test_load_from_file(self):
        sg_obj = load_sgrid(self.sgrid_test_file)
        self.assertIsInstance(sg_obj, SGrid)

    def test_load_from_dataset(self):
        sg_obj = load_sgrid(self.ds)
        self.assertIsInstance(sg_obj, SGrid)


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
        self.sg_obj = load_sgrid(self.sgrid_test_file)
        self.write_path = os.path.join(CURRENT_DIR, 'test_sgrid_write.nc')

    def test_center_lon(self):
        center_lon = self.sg_obj.center_lon
        self.assertEqual(center_lon.shape, (4, 4))

    def test_center_lat(self):
        center_lat = self.sg_obj.center_lat
        self.assertEqual(center_lat.shape, (4, 4))

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
                         u'fake_u',
                         u'lon_rho',
                         u'lat_rho',
                         u'lon_psi',
                         u'lat_psi',
                         u'lat_u',
                         u'lon_u',
                         u'lat_v',
                         u'lon_v',
                         u'salt',
                         u'zeta'
                         ]
        self.assertEqual(len(dataset_vars), len(expected_vars))
        self.assertEqual(dataset_vars, expected_vars)

    def test_grid_variables(self):
        dataset_grid_variables = self.sg_obj.grid_variables
        expected_grid_variables = [u'u', u'v', u'fake_u', u'salt']
        self.assertEqual(len(dataset_grid_variables),
                         len(expected_grid_variables))
        self.assertEqual(set(dataset_grid_variables),
                         set(expected_grid_variables))

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
        self.assertEqual(len(dataset_non_grid_variables),
                         len(expected_non_grid_variables))
        self.assertEqual(set(dataset_non_grid_variables),
                         set(expected_non_grid_variables))

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

    @mock.patch('pysgrid.sgrid.Dataset')
    def test_write_sgrid_to_netcdf(self, mock_nc):
        self.sg_obj.save_as_netcdf(self.write_path)
        mock_nc.assert_called_with(self.write_path, 'w')


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
        self.sgrid_obj = load_sgrid(self.sgrid_test_file)

    def test_face_coordinate_inference(self):
        face_coordinates = self.sgrid_obj.face_coordinates
        expected_face_coordinates = (u'XZ', u'YZ')
        self.assertEqual(face_coordinates, expected_face_coordinates)

    def test_center_lon(self):
        center_lon = self.sgrid_obj.center_lon
        self.assertEqual(center_lon.shape, (4, 4))

    def test_center_lat(self):
        center_lat = self.sgrid_obj.center_lat
        self.assertEqual(center_lat.shape, (4, 4))

    def test_node_lon(self):
        node_lon = self.sgrid_obj.node_lon
        self.assertEqual(node_lon.shape, (4, 4))

    def test_node_lat(self):
        node_lat = self.sgrid_obj.node_lat
        self.assertEqual(node_lat.shape, (4, 4))

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
        self.sg_obj = load_sgrid(self.sgrid_test_file)

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
        self.sg_obj = load_sgrid(self.sgrid_test_file)

    def test_center_lon(self):
        center_lon = self.sg_obj.center_lon
        self.assertEqual(center_lon.shape, (4, 4))

    def test_center_lat(self):
        center_lat = self.sg_obj.center_lat
        self.assertEqual(center_lat.shape, (4, 4))

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
        xcor_center_expected = (np.s_[:], np.s_[:])
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
        expected_grid_variables = [u'U1', u'V1', u'FAKE_U1', u'W', u'FAKE_W']
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
        self.sg_obj = load_sgrid(self.sgrid_test_file)
        self.sg_obj.save_as_netcdf(self.sgrid_target)
        self.target = load_sgrid(self.sgrid_target)

    def tearDown(self):
        os.remove(self.sgrid_target)

    def test_sgrid(self):
        self.assertIsInstance(self.target, SGrid)

    def test_nodes(self):
        node_lon = self.target.node_lon
        node_lat = self.target.node_lat
        self.assertIsNone(node_lon) and self.assertIsNone(node_lat)

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
        self.sg_obj = load_sgrid(self.sgrid_test_file)
        self.sg_obj.save_as_netcdf(self.sgrid_target)
        self.target = load_sgrid(self.sgrid_target)

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
                                u'U1',
                                u'FAKE_U1',
                                u'V1',
                                u'W',
                                u'FAKE_W',
                                u'time',
                                u'latitude',
                                u'longitude',
                                u'grid_latitude',
                                u'grid_longitude'
                                ]
        target_grid_vars = self.target.grid_variables
        expected_target_grid_vars = [u'U1',
                                     u'FAKE_U1',
                                     u'V1',
                                     u'W',
                                     u'FAKE_W'
                                     ]
        target_face_coordinates = self.target.face_coordinates
        expected_target_face_coordinates = (u'XZ', u'YZ')
        self.assertIsInstance(self.target, SGrid)
        self.assertEqual(len(target_dims), len(expected_target_dims))
        self.assertEqual(set(target_dims), set(expected_target_dims))
        self.assertEqual(len(target_vars), len(expected_target_vars))
        self.assertEqual(set(target_vars), set(expected_target_vars))
        self.assertEqual(len(target_grid_vars), len(expected_target_grid_vars))
        self.assertEqual(set(target_grid_vars), set(expected_target_grid_vars))
        self.assertEqual(target_face_coordinates,
                         expected_target_face_coordinates)

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
        np.testing.assert_almost_equal(original_angles, saved_angles, decimal=3)  # noqa
