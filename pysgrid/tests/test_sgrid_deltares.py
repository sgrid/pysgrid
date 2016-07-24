"""
Created on Apr 7, 2015

@author: ayan

"""

from __future__ import (absolute_import, division, print_function)

import pytest
import numpy as np

from ..sgrid import SGrid, load_grid

from .write_nc_test_files import (deltares_sgrid,
                                  deltares_sgrid_no_optional_attr)


"""
Test SGrid No Coordinates.

Test to make sure that if no coordinates (e.g. face, edge1, etc)
are specified, those coordinates can be inferred from the dataset.

A file is representing a delft3d dataset is used for this test.

"""


@pytest.fixture
def sgrid_deltares(deltares_sgrid_no_optional_attr):
    return load_grid(deltares_sgrid_no_optional_attr)


def test_face_coordinate_inference(sgrid_deltares):
    face_coordinates = sgrid_deltares.face_coordinates
    expected_face_coordinates = (u'XZ', u'YZ')
    assert face_coordinates == expected_face_coordinates


def test_center_lon_deltares_no_coord(sgrid_deltares):
    center_lon = sgrid_deltares.center_lon
    assert center_lon.shape == (4, 4)


def test_center_lat_deltares_no_coord(sgrid_deltares):
    center_lat = sgrid_deltares.center_lat
    assert center_lat.shape == (4, 4)


def test_node_lon(sgrid_deltares):
    node_lon = sgrid_deltares.node_lon
    assert node_lon.shape == (4, 4)


def test_node_lat(sgrid_deltares):
    node_lat = sgrid_deltares.node_lat
    assert node_lat.shape == (4, 4)


def test_grid_angles(sgrid_deltares):
    angles = sgrid_deltares.angles
    angles_shape = (4, 4)
    assert angles.shape == angles_shape


"""
Test SGrid Delft3d Dataset.

"""


@pytest.fixture
def sgrid_obj_deltares(deltares_sgrid):
    return load_grid(deltares_sgrid)


def test_center_lon_deltares(sgrid_obj_deltares):
    center_lon = sgrid_obj_deltares.center_lon
    assert center_lon.shape == (4, 4)


def test_center_lat_deltares(sgrid_obj_deltares):
    center_lat = sgrid_obj_deltares.center_lat
    assert center_lat.shape == (4, 4)


def test_topology_dimension_deltares(sgrid_obj_deltares):
    topology_dim = sgrid_obj_deltares.topology_dimension
    expected_dim = 2
    assert topology_dim == expected_dim


def test_variable_slice(sgrid_obj_deltares):
    u_center_slices = sgrid_obj_deltares.U1.center_slicing
    v_center_slices = sgrid_obj_deltares.V1.center_slicing
    u_center_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:])
    v_center_expected = (np.s_[:], np.s_[:], np.s_[1:], np.s_[:])
    xz_center_slices = sgrid_obj_deltares.XZ.center_slicing
    xcor_center_slices = sgrid_obj_deltares.XCOR.center_slicing
    xz_center_expected = (np.s_[1:], np.s_[1:])
    xcor_center_expected = (np.s_[:], np.s_[:])
    assert u_center_slices == u_center_expected
    assert v_center_slices == v_center_expected
    assert xz_center_slices == xz_center_expected
    assert xcor_center_slices == xcor_center_expected


def test_averaging_axes(sgrid_obj_deltares):
    u1c_axis = sgrid_obj_deltares.U1.center_axis
    u1c_expected = 0
    v1n_axis = sgrid_obj_deltares.V1.node_axis
    v1n_expected = 0
    latitude_c_axis = sgrid_obj_deltares.latitude.center_axis
    latitude_n_axis = sgrid_obj_deltares.latitude.node_axis
    assert u1c_axis == u1c_expected
    assert v1n_axis == v1n_expected
    assert latitude_c_axis is None
    assert latitude_n_axis is None


def test_grid_optional_attrs(sgrid_obj_deltares):
    face_coordinates = sgrid_obj_deltares.face_coordinates
    node_coordinates = sgrid_obj_deltares.node_coordinates
    edge1_coordinates = sgrid_obj_deltares.edge1_coordinates
    edge2_coordinates = sgrid_obj_deltares.edge2_coordinates
    fc_expected = ('XZ', 'YZ')
    nc_expected = ('XCOR', 'YCOR')
    assert face_coordinates == fc_expected
    assert node_coordinates == nc_expected
    assert edge1_coordinates is None
    assert edge2_coordinates is None


def test_grid_variables_deltares(sgrid_obj_deltares):
    grid_variables = sgrid_obj_deltares.grid_variables
    expected_grid_variables = [u'U1', u'V1', u'FAKE_U1', u'W', u'FAKE_W']
    assert set(grid_variables) == set(expected_grid_variables)


def test_angles(sgrid_obj_deltares):
    angles = sgrid_obj_deltares.angles
    expected_shape = (4, 4)
    assert angles.shape == expected_shape


def test_no_3d_attributes(sgrid_obj_deltares):
    assert not hasattr(sgrid_obj_deltares, 'volume_padding')
    assert not hasattr(sgrid_obj_deltares, 'volume_dimensions')
    assert not hasattr(sgrid_obj_deltares, 'volume_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'face1_padding')
    assert not hasattr(sgrid_obj_deltares, 'face1_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'face1_dimensions')
    assert not hasattr(sgrid_obj_deltares, 'face2_padding')
    assert not hasattr(sgrid_obj_deltares, 'face2_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'face2_dimensions')
    assert not hasattr(sgrid_obj_deltares, 'face3_padding')
    assert not hasattr(sgrid_obj_deltares, 'face3_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'edge3_padding')
    assert not hasattr(sgrid_obj_deltares, 'edge3_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'edge3_dimensions')


def test_2d_attributes(sgrid_obj_deltares):
    assert hasattr(sgrid_obj_deltares, 'face_padding')
    assert hasattr(sgrid_obj_deltares, 'face_coordinates')
    assert hasattr(sgrid_obj_deltares, 'face_dimensions')
    assert hasattr(sgrid_obj_deltares, 'vertical_padding')
    assert hasattr(sgrid_obj_deltares, 'vertical_dimensions')


"""
Test SGrid Save Node Coordinates.

Test that `SGrid.save_as_netcdf `is saving the content correctly.

"""


def test_round_trip(deltares_sgrid, tmpdir):
    """
    TODO: add more "round-trip" tests.

    """
    fname = tmpdir.mkdir('files').join('deltares_roundtrip.nc')
    sg_obj = load_grid(deltares_sgrid)
    sg_obj.save_as_netcdf(fname)


@pytest.fixture
def sgrid_deltares_sgrid(deltares_sgrid):
    return load_grid(deltares_sgrid)


def test_save_as_netcdf(sgrid_deltares_sgrid):
    """
    Test that the attributes in the
    saved netCDF file are as expected.

    """
    target_dims = sgrid_deltares_sgrid.dimensions
    expected_target_dims = [(u'MMAXZ', 4),
                            (u'NMAXZ', 4),
                            (u'MMAX', 4),
                            (u'NMAX', 4),
                            (u'KMAX', 2),
                            (u'KMAX1', 3),
                            (u'time', 2)
                            ]
    target_vars = sgrid_deltares_sgrid.variables
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
    target_grid_vars = sgrid_deltares_sgrid.grid_variables
    expected_target_grid_vars = [u'U1',
                                 u'FAKE_U1',
                                 u'V1',
                                 u'W',
                                 u'FAKE_W']
    target_face_coordinates = sgrid_deltares_sgrid.face_coordinates
    expected_target_face_coordinates = (u'XZ', u'YZ')
    assert isinstance(sgrid_deltares_sgrid, SGrid)
    assert len(target_dims) == len(expected_target_dims)
    assert set(target_dims) == set(expected_target_dims)
    assert len(target_vars) == len(expected_target_vars)
    assert set(target_vars) == set(expected_target_vars)
    assert len(target_grid_vars) == len(expected_target_grid_vars)
    assert set(target_grid_vars) == set(expected_target_grid_vars)
    assert target_face_coordinates == expected_target_face_coordinates


def test_saved_sgrid_attributes(sgrid_deltares_sgrid):
    """
    Test that calculated/inferred attributes
    are as expected from the saved filed.

    """
    u1_var = sgrid_deltares_sgrid.U1
    u1_var_center_avg_axis = u1_var.center_axis
    expected_u1_center_axis = 0
    u1_vector_axis = u1_var.vector_axis
    expected_u1_vector_axis = 'X'
    original_angles = sgrid_deltares_sgrid.angles
    saved_angles = sgrid_deltares_sgrid.angles
    assert u1_var_center_avg_axis == expected_u1_center_axis
    assert u1_vector_axis == expected_u1_vector_axis
    np.testing.assert_almost_equal(original_angles, saved_angles, decimal=3)  # noqa
