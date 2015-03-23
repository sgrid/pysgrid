'''
Created on Mar 23, 2015

@author: ayan
'''
import unittest

from ..utils import ParsePadding


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
        padding_var = padding_datum_0.dim_var
        padding_dim = padding_datum_0.dim_name
        expected_padding_var = 'xi_psi'
        expected_padding_type = 'both'
        expected_padding_dim = 'xi_rho'
        self.assertEqual(len(result), expected_len)
        self.assertEqual(padding_type, expected_padding_type)
        self.assertEqual(padding_var, expected_padding_var)
        self.assertEqual(padding_dim, expected_padding_dim)
        
    def test_one_padding_type(self):
        result = self.pp.parse_padding(padding_str=self.with_one_padding)
        expected_len = 1
        padding_datum_0 = result[0]
        padding_type = padding_datum_0.padding
        padding_var = padding_datum_0.dim_var
        padding_dim = padding_datum_0.dim_name
        expected_padding_type = 'high'
        expected_padding_var = 'xi_psi'
        expected_padding_dim = 'xi_v'
        self.assertEqual(len(result), expected_len)
        self.assertEqual(padding_type, expected_padding_type)
        self.assertEqual(padding_var, expected_padding_var)
        self.assertEqual(padding_dim, expected_padding_dim)
        
    def test_no_padding(self):
        result = self.pp.parse_padding(padding_str=self.with_no_padding)
        self.assertIsNone(result)