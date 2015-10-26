'''
Created on Apr 15, 2015

@author: ayan
'''
import os
import unittest

from netCDF4 import Dataset
import numpy as np

from ..sgrid import SGrid2D
from ..utils import GridPadding
from ..variables import SGridVariable
from .write_nc_test_files import deltares_sgrid, roms_sgrid, wrf_sgrid_2d


class TestSGridVariableROMS(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_file = roms_sgrid()

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.test_file)

    def setUp(self):
        self.face_padding = [GridPadding(mesh_topology_var=u'grid', face_dim=u'xi_rho', node_dim=u'xi_psi', padding=u'both'),
                             GridPadding(mesh_topology_var=u'grid', face_dim=u'eta_rho', node_dim=u'eta_psi', padding=u'both')
                             ]
        self.sgrid = SGrid2D(face_padding=self.face_padding,
                             node_dimensions='xi_psi eta_psi'
                             )
        self.dataset = Dataset(self.test_file)
        self.test_var_1 = self.dataset.variables['u']
        self.test_var_2 = self.dataset.variables['zeta']
        self.test_var_3 = self.dataset.variables['salt']
        self.test_var_4 = self.dataset.variables['fake_u']

    def tearDown(self):
        self.dataset.close()

    def test_create_sgrid_variable_object(self):
        sgrid_var = SGridVariable.create_variable(self.test_var_1, self.sgrid)
        self.assertIsInstance(sgrid_var, SGridVariable)

    def test_attributes_with_grid(self):
        sgrid_var = SGridVariable.create_variable(self.test_var_1, self.sgrid)
        sgrid_var_name = sgrid_var.variable
        sgrid_var_name_expected = 'u'
        sgrid_var_dim = sgrid_var.dimensions
        sgrid_var_dim_expected = ('time', 's_rho', 'eta_u', 'xi_u')
        sgrid_var_grid = sgrid_var.grid
        sgrid_var_grid_expected = 'some grid'
        sgrid_var_location = sgrid_var.location
        sgrid_var_location_expected = 'edge1'
        sgrid_var_dtype = sgrid_var.dtype
        x_axis = sgrid_var.x_axis
        x_axis_expected = 'xi_u'
        y_axis = sgrid_var.y_axis
        y_axis_expected = 'eta_u'
        z_axis = sgrid_var.z_axis
        standard_name = sgrid_var.standard_name
        expected_standard_name = 'sea_water_x_velocity'
        coordinates = sgrid_var.coordinates
        expected_coordinates = ('time', 's_rho', 'lat_u', 'lon_u')
        self.assertEqual(sgrid_var_name, sgrid_var_name_expected)
        self.assertEqual(sgrid_var_dim, sgrid_var_dim_expected)
        self.assertEqual(sgrid_var_grid, sgrid_var_grid_expected)
        self.assertEqual(sgrid_var_location, sgrid_var_location_expected)
        self.assertEqual(sgrid_var_dtype, np.dtype('float32'))
        self.assertEqual(x_axis, x_axis_expected)
        self.assertEqual(y_axis, y_axis_expected)
        self.assertIsNone(z_axis)
        self.assertEqual(standard_name, expected_standard_name)
        self.assertEqual(coordinates, expected_coordinates)

    def test_face_location_inference(self):
        sgrid_var = SGridVariable.create_variable(self.test_var_3, self.sgrid)
        sgrid_var_location = sgrid_var.location
        expected_location = 'face'
        self.assertEqual(sgrid_var_location, expected_location)

    def test_edge_location_inference(self):
        sgrid_var = SGridVariable.create_variable(self.test_var_4, self.sgrid)
        sgrid_var_location = sgrid_var.location
        # representative ROMS sgrid
        # None is expected since edge1 and edge2 attributes are not defined
        self.assertIsNone(sgrid_var_location)

    def test_edge_location_inference_with_defined_edges(self):
        self.sgrid.edge1_padding = [GridPadding(mesh_topology_var=u'grid', face_dim=u'eta_u', node_dim=u'eta_psi', padding=u'both')]
        self.sgrid.edge2_padding = [GridPadding(mesh_topology_var=u'grid', face_dim=u'xi_v', node_dim=u'xi_psi', padding=u'both')]
        sgrid_var = SGridVariable.create_variable(self.test_var_4, self.sgrid)
        sgrid_var_location = sgrid_var.location
        expected_location = 'edge1'
        self.assertEqual(sgrid_var_location, expected_location)

    def test_attributes_with_location(self):
        sgrid_var = SGridVariable.create_variable(self.test_var_2, self.sgrid)
        sgrid_var_name = sgrid_var.variable
        sgrid_var_name_expected = 'zeta'
        sgrid_var_dim = sgrid_var.dimensions
        sgrid_var_grid = sgrid_var.grid
        sgrid_var_location = sgrid_var.location
        sgrid_var_location_expected = 'face'
        sgrid_var_dim_expected = ('time', 'eta_rho', 'xi_rho')
        sgrid_var_dtype = sgrid_var.dtype
        x_axis = sgrid_var.x_axis
        y_axis = sgrid_var.y_axis
        z_axis = sgrid_var.z_axis
        self.assertEqual(sgrid_var_name, sgrid_var_name_expected)
        self.assertEqual(sgrid_var_dim, sgrid_var_dim_expected)
        self.assertIsNone(sgrid_var_grid)
        self.assertEqual(sgrid_var_location, sgrid_var_location_expected)
        self.assertEqual(sgrid_var_dtype, np.dtype('float32'))
        self.assertIsNone(x_axis)
        self.assertIsNone(y_axis)
        self.assertIsNone(z_axis)

    def test_vector_directions(self):
        u_var = SGridVariable.create_variable(self.test_var_1, self.sgrid)
        u_vector_axis = u_var.vector_axis
        expected_u_axis = 'X'
        zeta_var = SGridVariable.create_variable(self.test_var_2, self.sgrid)
        zeta_axis = zeta_var.vector_axis
        self.assertEqual(u_vector_axis, expected_u_axis)
        self.assertIsNone(zeta_axis)


class TestSGridVariablesWRF(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_file = wrf_sgrid_2d()

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.test_file)

    def setUp(self):
        self.face_padding = [GridPadding(mesh_topology_var=u'grid', face_dim=u'west_east', node_dim=u'west_east_stag', padding=u'none'),
                             GridPadding(mesh_topology_var=u'grid', face_dim=u'south_north', node_dim=u'south_north_stag', padding=u'none')
                             ]
        self.node_dimensions = 'west_east_stag south_north_stag'
        self.sgrid = SGrid2D(face_padding=self.face_padding,
                             node_dimensions=self.node_dimensions
                             )
        self.dataset = Dataset(self.test_file)
        self.test_var_1 = self.dataset.variables['SNOW']
        self.test_var_2 = self.dataset.variables['FAKE_U']

    def tearDown(self):
        self.dataset.close()

    def test_face_location_inference(self):
        sg_var = SGridVariable.create_variable(self.test_var_1, self.sgrid)
        sg_var_location = sg_var.location
        expected_location = 'face'
        self.assertEqual(sg_var_location, expected_location)

    def test_edge_location_inference(self):
        sg_var = SGridVariable.create_variable(self.test_var_2, self.sgrid)
        sg_var_location = sg_var.location
        expected_location = 'edge1'
        self.assertEqual(sg_var_location, expected_location)


class TestSGridVariablesDeltares(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_file = deltares_sgrid()

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.test_file)

    def setUp(self):
        self.face_padding = [GridPadding(mesh_topology_var=u'grid', face_dim=u'MMAXZ', node_dim=u'MMAX', padding=u'low'),
                             GridPadding(mesh_topology_var=u'grid', face_dim=u'NMAXZ', node_dim=u'NMAX', padding=u'low')
                             ]
        self.node_dimensions = 'MMAX NMAX'
        self.sgrid = SGrid2D(face_padding=self.face_padding,
                             node_dimensions=self.node_dimensions
                             )
        self.dataset = Dataset(self.test_file)
        self.test_var_1 = self.dataset.variables['FAKE_W']
        self.test_var_2 = self.dataset.variables['FAKE_U1']

    def tearDown(self):
        self.dataset.close()

    def test_face_location_inference(self):
        sg_var = SGridVariable.create_variable(self.test_var_1, self.sgrid)
        sg_var_location = sg_var.location
        expected_location = 'face'
        self.assertEqual(sg_var_location, expected_location)

    def test_edge_location_inference(self):
        sg_var = SGridVariable.create_variable(self.test_var_2, self.sgrid)
        sg_var_location = sg_var.location
        expected_location = 'edge2'
        self.assertEqual(sg_var_location, expected_location)
