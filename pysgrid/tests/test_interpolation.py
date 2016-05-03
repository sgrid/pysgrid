'''
Created on Feb 17, 2016

@author: jay.hennen
'''

import numpy as np
from netCDF4 import Dataset
from pysgrid import SGrid, load_grid

node_lon = np.array(([1, 3, 5], [1, 3, 5], [1, 3, 5]))
node_lat = np.array(([1, 1, 1], [3, 3, 3], [5, 5, 5]))
edge2_lon = np.array(([0, 2, 4, 6], [0, 2, 4, 6], [0, 2, 4, 6]))
edge2_lat = np.array(([1, 1, 1, 1], [3, 3, 3, 3], [5, 5, 5, 5]))
edge1_lon = np.array(([1, 3, 5], [1, 3, 5], [1, 3, 5], [1, 3, 5]))
edge1_lat = np.array(([0, 0, 0], [2, 2, 2], [4, 4, 4], [6, 6, 6]))
center_lon = np.array(([0, 2, 4, 6], [0, 2, 4, 6], [0, 2, 4, 6], [0, 2, 4, 6]))
center_lat = np.array(([0, 0, 0, 0], [2, 2, 2, 2], [4, 4, 4, 4], [6, 6, 6, 6]))

sgrid = SGrid(node_lon=node_lon,
              node_lat=node_lat,
              edge1_lon=edge1_lon,
              edge1_lat=edge1_lat,
              edge2_lon=edge2_lon,
              edge2_lat=edge2_lat,
              center_lon=center_lon,
              center_lat=center_lat)

c_var = np.array(([0, 0, 0, 0], [0, 1, 2, 0], [0, 2, 1, 0], [0, 0, 0, 0]))
e2_var = np.array(([1, 0, 0, 1], [0, 1, 2, 0], [0, 0, 0, 0]))
e1_var = np.array(([1, 1, 0], [0, 1, 0], [0, 2, 0], [1, 1, 0]))
n_var = np.array(([0, 1, 0], [1, 0, 1], [0, 1, 0]))

ptsx, ptsy = np.mgrid[0:6:600j, 0:6:600j]
pts = np.stack((ptsx, ptsy), axis=-1)


def test_locate_faces():
    diagonal = [[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [5, 5], [6, 6]]
    ind_ans = np.ma.masked_equal(
        [[-1, -1], [0, 0], [0, 0], [1, 1], [1, 1], [-1, -1], [-1, -1]], -1)
    indices = sgrid.locate_faces(diagonal)
    res = (indices.data == ind_ans.data).all() and (
        indices.mask == ind_ans.mask).all()
    assert(res)


def test_points_in_polys():
    from pysgrid.utils import points_in_polys
    points = np.array(
        [[0, 0],
         [1, 0],
         [2, 0],
         [0, 1],
         [0, 2],
         [1, 2],
         [2, 2],
         [2, 1]])
    polygon = np.array(([0, 0], [2, 0], [2, 2], [0, 2])).reshape(1, 4, 2)
    pinp = np.array([points_in_polys(point.reshape(1, 2), polygon)
                     for point in points]).reshape(-1)
    answer = sgrid.locate_faces(points + 1) == [0, 0]
    answer = np.logical_and(answer[:, 0], answer[:, 1])
    res = (answer == pinp).all()
    print pinp
    print sgrid.locate_faces(points + 1).data
    print sgrid.locate_faces(points + 1).mask
    print answer
    print res
    assert 0
    assert(res)


def test_index_translation():
    ptsx, ptsy = np.mgrid[0:6:7j, 0:6:7j]
    points = np.stack((ptsx, ptsy), axis=-1).reshape(-1, 2)
    indices = sgrid.locate_faces(points)
    # to rho
    c_ind = sgrid.translate_index(points, indices, 'center')
    cgrid = SGrid(node_lon=center_lon,
                  node_lat=center_lat)
    answer = cgrid.locate_faces(points)
    res = (c_ind == answer).all()
    assert(res)

    e1_ind = sgrid.translate_index(points, indices, 'edge1')
    e1_grid = SGrid(node_lon=edge1_lon,
                    node_lat=edge1_lat)
    answer = e1_grid.locate_faces(points)
    res = (e1_ind == answer).all()
    assert(res)

    e2_ind = sgrid.translate_index(points, indices, 'edge2')
    e2_grid = SGrid(node_lon=edge2_lon,
                    node_lat=edge2_lat)
    answer = e2_grid.locate_faces(points)
    res = (e2_ind == answer).all()
    assert(res)

    # Test case for places where edge1 increases slope clockwise
    te1_node_lon = np.array(([1, 3, 5], [1, 3, 5], [2, 4, 6]))
    te1_node_lat = np.array(([1, 1, 1], [3, 3, 3], [5, 5, 5]))
    te1_edge1_lon = np.array(
        ([1, 3, 5], [1, 3, 5], [1.5, 3.5, 5.5], [2.5, 5.5, 7.5]))
    te1_edge1_lat = np.array(([0, 0, 0], [2, 2, 2], [4, 4, 4], [6, 6, 6]))
    test_points = np.array(([1.01, 3], [3.01, 3.01], [3.01, 2.99], [5.01, 3]))

    g = SGrid(node_lon=te1_node_lon,
              node_lat=te1_node_lat,
              edge1_lon=te1_edge1_lon,
              edge1_lat=te1_edge1_lat)

    indices = g.locate_faces(test_points)
    translated_indices = g.translate_index(test_points, indices, 'edge1')
    answer = np.ma.masked_array(([-1, -1], [1, 0], [1, 0], [-1, -1]),
                                mask=[[True, True], [False, False], [False, False], [True, True]])
    res = (translated_indices == answer).all()
    assert(res)

    # Test case for places where edge1 increases slope counterclockwise
    te1_node_lon = np.array(([2, 4, 6], [2, 4, 6], [1, 3, 5]))
    te1_node_lat = np.array(([1, 1, 1], [3, 3, 3], [5, 5, 5]))
    te1_edge1_lon = np.array(
        ([2, 4, 6], [2, 4, 6], [1.5, 3.5, 5.5], [0.5, 2.5, 4.5]))
    te1_edge1_lat = np.array(([0, 0, 0], [2, 2, 2], [4, 4, 4], [6, 6, 6]))
    test_points = np.array(([1.99, 3], [3.99, 3.01], [3.99, 2.99], [5.99, 3]))

    g = SGrid(node_lon=te1_node_lon,
              node_lat=te1_node_lat,
              edge1_lon=te1_edge1_lon,
              edge1_lat=te1_edge1_lat)

    indices = g.locate_faces(test_points)
    translated_indices = g.translate_index(test_points, indices, 'edge1')
    answer = np.ma.masked_array(([-1, -1], [1, 1], [1, 1], [-1, -1]),
                                mask=[[True, True], [False, False], [False, False], [True, True]])
    res = (translated_indices == answer).all()
    assert(res)

    # Test case for places where edge2 increases slope counterclockwise
    te2_node_lon = np.array(([1, 3, 5], [1, 3, 5], [1, 3, 5]))
    te2_node_lat = np.array(([1, 1, 2], [3, 3, 4], [5, 5, 6]))
    te2_edge2_lon = np.array(
        ([0, 2, 4, 6], [0, 2, 4, 6], [0, 2, 4, 6]))
    te2_edge2_lat = np.array(
        ([1, 1, 1.5, 2.5], [3, 3, 3.5, 4.5], [5, 5, 5.5, 6.5]))
    test_points = np.array(([3, 1.01], [3.01, 3.01], [2.99, 3.01], [3, 5.01]))

    g = SGrid(node_lon=te2_node_lon,
              node_lat=te2_node_lat,
              edge2_lon=te2_edge2_lon,
              edge2_lat=te2_edge2_lat)

    indices = g.locate_faces(test_points)
    translated_indices = g.translate_index(test_points, indices, 'edge2')
    answer = np.ma.masked_array(([-1, -1], [0, 1], [0, 1], [-1, -1]),
                                mask=[[True, True], [False, False], [False, False], [True, True]])
    res = (translated_indices == answer).all()
    assert(res)

    # Test case for places where edge2 increases slope clockwise
    te2_node_lon = np.array(([1, 3, 5], [1, 3, 5], [1, 3, 5]))
    te2_node_lat = np.array(([2, 2, 1], [4, 4, 3], [6, 6, 5]))
    te2_edge2_lon = np.array(
        ([0, 2, 4, 6], [0, 2, 4, 6], [0, 2, 4, 6]))
    te2_edge2_lat = np.array(
        ([2, 2, 1.5, 0.5], [4, 4, 3.5, 2.5], [6, 6, 5.5, 4.5]))
    test_points = np.array(([3, 1.99], [3.01, 3.99], [2.99, 3.99], [3, 5.99]))

    g = SGrid(node_lon=te2_node_lon,
              node_lat=te2_node_lat,
              edge2_lon=te2_edge2_lon,
              edge2_lat=te2_edge2_lat)

    indices = g.locate_faces(test_points)
    translated_indices = g.translate_index(test_points, indices, 'edge2')
    answer = np.ma.masked_array(([-1, -1], [1, 1], [1, 1], [-1, -1]),
                                mask=[[True, True], [False, False], [False, False], [True, True]])
    res = (translated_indices == answer).all()
    assert(res)


def test_interpolation_alphas():
    points = np.array(([2, 2], [2, 4], [4, 2], [4, 4]))
    alphas_c = sgrid.interpolation_alphas(points, location='center')
    alphas_e1 = sgrid.interpolation_alphas(points, location='edge1')
    alphas_e2 = sgrid.interpolation_alphas(points, location='edge2')
    alphas_n = sgrid.interpolation_alphas(points, location='node')

    answer_c = np.array([[1., 0., 0., 0.],
                         [1., 0., 0., 0.],
                         [1., 0., 0., 0.],
                         [1., 0., 0., 0.]])
    answer_e1 = np.array([[0.5, 0., 0., 0.5],
                          [0.5, 0., 0., 0.5],
                          [0.5, 0., 0., 0.5],
                          [0.5, 0., 0., 0.5]])
    answer_e2 = np.array([[0.5, 0.5, 0., 0.],
                          [0.5, 0.5, 0., 0.],
                          [0.5, 0.5, 0., 0.],
                          [0.5, 0.5, 0., 0.]])
    answer_n = np.array([[0.25, 0.25, 0.25, 0.25],
                         [0.25, 0.25, 0.25, 0.25],
                         [0.25, 0.25, 0.25, 0.25],
                         [0.25, 0.25, 0.25, 0.25]])

    assert(np.all(alphas_c == answer_c))
    assert(np.all(alphas_n == answer_n))
    assert(np.all(alphas_e1 == answer_e1))
    assert(np.all(alphas_e2 == answer_e2))
