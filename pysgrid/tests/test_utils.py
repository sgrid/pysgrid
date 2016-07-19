'''
Created on Mar 23, 2015

@author: ayan
'''

from __future__ import (absolute_import, division, print_function)

import unittest

import numpy as np

from pysgrid.utils import (calculate_bearing,
                           calculate_angle_from_true_east,
                           check_element_equal,
                           does_intersection_exist,
                           pair_arrays,
                           )


class TestDoesIntersectionExist(unittest.TestCase):

    def setUp(self):
        self.tuple_a = (718, 903, 1029, 1701)
        self.tuple_b = (718, 828)
        self.tuple_c = (15, 20)

    def test_intersect_exists(self):
        result = does_intersection_exist(self.tuple_a, self.tuple_b)
        self.assertTrue(result)

    def test_intersect_does_not_exist(self):
        result = does_intersection_exist(self.tuple_a, self.tuple_c)
        self.assertFalse(result)


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


class TestCalculateBearing(unittest.TestCase):

    def setUp(self):
        self.points = np.array([(-93.51105439, 11.88846735),
                                (-93.46607342, 11.90917952)])
        self.point_1 = self.points[:-1, :]
        self.point_2 = self.points[1:, :]

    def test_bearing_calculation(self):
        result = calculate_bearing(self.point_1, self.point_2)
        expected = 64.7947
        np.testing.assert_almost_equal(result, expected, decimal=3)


class TestCalculateAngleFromTrueEast(unittest.TestCase):

    def setUp(self):
        self.vertical_1 = np.array([[-122.41, 37.78], [-122.33, 37.84], [-122.22, 37.95]])  # noqa
        self.vertical_2 = np.array([[-90.07, 29.95], [-89.97, 29.84], [-89.91, 29.76]])  # noqa
        self.vertical_3 = np.array([[-89.40, 43.07], [-89.49, 42.93], [-89.35, 42.84]])  # noqa
        self.vertical_4 = np.array([[-122.41, 37.78], [-122.53, 37.84], [-122.67, 37.95]])  # noqa
        self.centers = np.array((self.vertical_1,
                                 self.vertical_2,
                                 self.vertical_3,
                                 self.vertical_4))

    def test_angle_from_true_east_calculation(self):
        bearing_start_points = self.centers[:, :-1, :]
        bearing_end_points = self.centers[:, 1:, :]
        angle_from_true_east = calculate_angle_from_true_east(bearing_start_points, bearing_end_points)  # noqa
        expected_values = np.array([[0.7598, 0.9033, 0.9033],
                                    [-0.903, -0.994, -0.994],
                                    [-2.011, -0.719, -0.719],
                                    [-3.706, -3.926, -3.926]
                                    ])
        expected_shape = (4, 3)
        np.testing.assert_almost_equal(angle_from_true_east, expected_values, decimal=3)  # noqa
        self.assertEqual(angle_from_true_east.shape, expected_shape)
