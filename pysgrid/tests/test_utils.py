'''
Created on Mar 23, 2015

@author: ayan
'''
import unittest
import numpy as np
from ..utils import ParsePadding, pair_arrays, check_array_dims
from ..custom_exceptions import CannotFindPadding, DimensionMismatch


class TestCheckArrayDims(unittest.TestCase):
    
    def setUp(self):
        self.a = np.array([1, 2, 3])
        self.b = np.array([-1, -2])
        
    def test_exception_raised(self):
        self.assertRaises(DimensionMismatch,
                          check_array_dims,
                          self.a,
                          self.b
                          )


class TestParsePadding(unittest.TestCase):
    
    def setUp(self):
        self.grid_topology = 'some_grid'
        self.pp = ParsePadding(self.grid_topology)
        self.with_two_padding = 'xi_rho: xi_psi (padding: both) eta_rho: eta_psi (padding: low)'
        self.with_one_padding = 'xi_v: xi_psi (padding: high) eta_v: eta_psi'
        self.with_no_padding = 'MMAXZ: MMAX NMAXZ: NMAX'
        
    def test_mesh_name(self):
        result = self.pp.parse_padding(padding_str=self.with_one_padding)
        mesh_topology = result[0].mesh_topology_var
        expected = 'some_grid'
        self.assertEqual(mesh_topology, expected)
    
    def test_two_padding_types(self):
        result = self.pp.parse_padding(padding_str=self.with_two_padding)
        expected_len = 2
        padding_datum_0 = result[0]
        padding_type = padding_datum_0.padding
        node_dim = padding_datum_0.node_dim
        face_dim = padding_datum_0.face_dim
        expected_node_dim = 'xi_psi'
        expected_padding_type = 'both'
        expected_face_dim = 'xi_rho'
        self.assertEqual(len(result), expected_len)
        self.assertEqual(padding_type, expected_padding_type)
        self.assertEqual(node_dim, expected_node_dim)
        self.assertEqual(face_dim, expected_face_dim)
        
    def test_one_padding_type(self):
        result = self.pp.parse_padding(padding_str=self.with_one_padding)
        expected_len = 1
        padding_datum_0 = result[0]
        padding_type = padding_datum_0.padding
        node_dim = padding_datum_0.node_dim
        face_dim = padding_datum_0.face_dim
        expected_padding_type = 'high'
        expected_node_dim = 'xi_psi'
        expected_face_dim = 'xi_v'
        self.assertEqual(len(result), expected_len)
        self.assertEqual(padding_type, expected_padding_type)
        self.assertEqual(node_dim, expected_node_dim)
        self.assertEqual(face_dim, expected_face_dim)
        
    def test_no_padding(self):
        self.assertRaises(CannotFindPadding, 
                          self.pp.parse_padding, 
                          padding_str=self.with_no_padding
                          )
        

class TestPairArrays(unittest.TestCase):
    
    def setUp(self):
        self.a1 = (1, 2)
        self.a2 = (3, 4)
        self.a3 = (5, 6)
        self.b1 = (10, 20)
        self.b2 = (30, 40)
        self.b3 = (50, 60)
        self.a = np.array([self.a1, self.a2, self.a3])
        self.b = np.array([self.b1, self.b2, self.b3])
        
    def test_pair_arrays(self):
        result = pair_arrays(self.a, self.b)
        x = [[(1, 10), (2, 20)],
             [(3, 30), (4, 40)],
             [(5, 50), (6, 60)]]
        expected = np.array(x)
        test_equal = (result == expected).all()
        self.assertTrue(test_equal)
        