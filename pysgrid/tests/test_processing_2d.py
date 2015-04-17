'''
Created on Apr 3, 2015

@author: ayan
'''
import unittest
import numpy as np
from ..processing_2d import vector_sum, rotate_vectors, avg_to_cell_center


class TestVectorSum(unittest.TestCase):
    
    def setUp(self):
        self.x_vector = np.array([3, 5, 9, 11])
        self.y_vector = np.array([4, 12, 40, 60])
        self.mismatch_vector = np.array([17, 45])
        
    def test_vector_sum(self):
        sum_result = vector_sum(self.x_vector, self.y_vector)
        expected = np.array([5, 13, 41, 61])
        np.testing.assert_almost_equal(sum_result, expected)
        
    def test_mismatch(self):
        self.assertRaises(ValueError, vector_sum, 
                          self.x_vector, self.mismatch_vector
                          )
        

class TestRotateVectors(unittest.TestCase):
    
    def setUp(self):
        self.x_vector = np.array([3, 5, 9, 11])
        self.y_vector = np.array([4, 12, 40, 60])
        self.angles_simple = np.array([0, np.pi/2, 0, np.pi/2])
        self.angles_complex = np.array([np.pi/6, np.pi/5, np.pi/4, np.pi/3])
        self.angles_mismatch = np.array([np.pi/7])
        
    def test_vector_rotation_simple(self):
        rotated_x, rotated_y = rotate_vectors(self.x_vector, self.y_vector, self.angles_simple)
        expected_x = np.array([3, -12, 9, -60])
        expected_y = np.array([4, 5, 40, 11])
        np.testing.assert_almost_equal(rotated_x, expected_x, decimal=3)
        np.testing.assert_almost_equal(rotated_y, expected_y, decimal=3)
        
    def test_vector_rotation_complex(self):
        rotated_x, rotated_y = rotate_vectors(self.x_vector, self.y_vector, self.angles_complex)
        expected_x = np.array([0.5981, -3.0083, -21.9203, -46.4615])
        expected_y = np.array([4.9641, 12.6471, 34.6482, 39.5263])
        np.testing.assert_almost_equal(rotated_x, expected_x, decimal=3)
        np.testing.assert_almost_equal(rotated_y, expected_y, decimal=3)
        
    def test_mismatch(self):
        self.assertRaises(ValueError, rotate_vectors,
                          self.x_vector, self.y_vector, self.angles_mismatch)
        

class TestAvgToCellCenter(unittest.TestCase):
    
    def setUp(self):
        self.data = np.array([[4, 5, 9, 10], [8, 39, 41, 20]])
        self.avg_dim_0 = 0
        self.avg_dim_1 = 1
        
    def test_no_transpose(self):
        avg_result = avg_to_cell_center(self.data, self.avg_dim_1)
        expected = np.array([[4.5, 7, 9.5], [23.5, 40, 30.5]])
        np.testing.assert_almost_equal(avg_result, expected, decimal=3)
        
    def test_with_transpose(self):
        avg_result = avg_to_cell_center(self.data, self.avg_dim_0)
        expected = np.array([[6, 22, 25, 15]])
        np.testing.assert_almost_equal(avg_result, expected, decimal=3)