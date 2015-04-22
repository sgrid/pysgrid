'''
Created on Mar 23, 2015

@author: ayan
'''
import unittest
import numpy as np
from ..utils import pair_arrays, check_element_equal


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
        np.testing.assert_almost_equal(result, expected, decimal=3)
        

class TestCheckElementEqual(unittest.TestCase):
    
    def setUp(self):
        self.a = [7, 7, 7, 7]
        self.b = [7, 8, 9, 10]
        
    def test_list_with_identical_elements(self):
        result = check_element_equal(self.a)
        self.assertTrue(result)
        
    def test_list_with_different_elements(self):
        result = check_element_equal(self.b)
        self.assertFalse(result)